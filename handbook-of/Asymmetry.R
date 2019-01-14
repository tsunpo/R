# =============================================================================
# Library      : Transcritional Strand Asymmetry
# Name         : handbook-of/Asymmetry.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 29/01/18
# =============================================================================

# -----------------------------------------------------------------------------
# Method: Finding mutations locate within Ensembl genes (Step 1)
# Last Modified: 29/01/18
# -----------------------------------------------------------------------------
getEnsGeneSNV <- function(pos, ensGene.chr) {
   ensGene.chr.start <- subset(ensGene.chr, pos >= start_position)
   ensGene.chr.start.end <- subset(ensGene.chr.start, pos <= end_position)
 
   return(ensGene.chr.start.end)
}

getEnsGeneSNVs <- function(vcf.chr, ensGene.chr.gene) {
   vcf.chr.start <- subset(vcf.chr, POS >= ensGene.chr.gene$start_position)
   vcf.chr.start.end <- subset(vcf.chr.start, POS <= ensGene.chr.gene$end_position)
 
   return(vcf.chr.start.end)
}

getEnsGeneTSS <- function(vcf.chr, ensGene.chr.gene, distance) {
   if (ensGene.chr.gene$strand == 1) {
      vcf.chr.start <- subset(vcf.chr, POS >= ensGene.chr.gene$start_position - distance)
      vcf.chr.start.end <- subset(vcf.chr.start, POS <= ensGene.chr.gene$start_position)
   } else {
      vcf.chr.start <- subset(vcf.chr, POS >= ensGene.chr.gene$end_position)
      vcf.chr.start.end <- subset(vcf.chr.start, POS <= ensGene.chr.gene$end_position + distance)
   }
 
   return(vcf.chr.start.end)
}

getMergedTable <- function(ensGene.chr.gene, vcf.gene.chr0) {
   ensGene.chr.genes <- ensGene.chr.gene
   if (nrow(vcf.gene.chr0) >= 2)
      for (v in 2:nrow(vcf.gene.chr0))
      ensGene.chr.genes <- rbind(ensGene.chr.genes, ensGene.chr.gene)
  
   return(cbind(ensGene.chr.genes, vcf.gene.chr0))
}

# -----------------------------------------------------------------------------
# Method: Read in all the SNVs (Step 3.1)
# Last Modified: 25/01/18
# -----------------------------------------------------------------------------
initMergedReport <- function(isExon) {
   colnames <- c("SAMPLE", "CHROM", "POS", "REF", "ALT", "ensembl_gene_id", "strand")
   if (isExon)
      colnames <- c(colnames, "exon")
   
   return(toTable("", length(colnames), 0, colnames))
}

getMergedReport <- function(sample, vcf) {
   vcf <- cbind(toTable(sample, 1, nrow(vcf), "SAMPLE"), vcf)
   
   return(vcf[,c("SAMPLE", "CHROM", "POS", "REF", "ALT", "ensembl_gene_id", "strand")])
}

## Define transcriptional asymmetry
## Sig.    s1        s2        s3        s4        s5        s6
## idx:    1         3         5         7         9         11
## w/c:    w   c     w   c     w   c     w   c     w   c     w   c   ## Reference or anti-reference strand
## s6:     C>A/G>T   C>G/G>C   C>T/G>A   T>A/A>T   T>C/A>G   T>G/A>C
REFS <- c("C","G",  "C","G",  "C","G",  "T","A",  "T","A",  "T","A")
ALTS <- c("A","T",  "G","C",  "T","A",  "A","T",  "C","G",  "G","C")
idxs <- seq(1, 12, 2)

###                               s   w/c  tx(strand)
## Data structure: E.g tx.snv.s6[[1]][[1]][[1]] is C>A Tx(+)
##                     tx.snv.s6[[1]][[1]][[2]] is C>A Tx(-)
##                     tx.snv.s6[[1]][[2]][[1]] is G>T Tx(+)
##                     tx.snv.s6[[1]][[2]][[2]] is G>T Tx(-)
##
initTableS6 <- function(isExon) {
   tx.snv.s6 <- list(list(list(), list()), list(list(), list()), list(list(), list()), list(list(), list()), list(list(), list()), list(list(), list()))
   for (i in 1:6)
      for (j in 1:2)
         for (k in 1:2)
            tx.snv.s6[[i]][[j]][[k]] <- initMergedReport(isExon)
 
   return(tx.snv.s6)
}

getTableS6 <- function(tx.snv, isExon) {
   tx.snv.s6 <- initTableS6(isExon)

   for (i in 1:length(idxs)) {
      idx <- idxs[i]
 
      watson <- subset(subset(tx.snv, REF == REFS[idx]), ALT == ALTS[idx])
      tx.snv.s6[[i]][[1]][[1]] <- rbind(tx.snv.s6[[i]][[1]][[1]], subset(watson, strand == 1))
      tx.snv.s6[[i]][[1]][[2]] <- rbind(tx.snv.s6[[i]][[1]][[2]], subset(watson, strand == -1))
 
      crick  <- subset(subset(tx.snv, REF == REFS[idx+1]), ALT == ALTS[idx+1])
      tx.snv.s6[[i]][[2]][[1]] <- rbind(tx.snv.s6[[i]][[2]][[1]], subset(crick, strand == 1))
      tx.snv.s6[[i]][[2]][[2]] <- rbind(tx.snv.s6[[i]][[2]][[2]], subset(crick, strand == -1))
   }
 
   return(tx.snv.s6)
}

# -----------------------------------------------------------------------------
# Method: Determine SNVs are on exons or introns (Step 3.2)
# Last Modified: 15/03/18
# -----------------------------------------------------------------------------
isSNVonExon <- function(chr, pos, ensGene.transcript.exon) {
   ensGene.transcript.exon.chr <- subset(ensGene.transcript.exon, chrom == chr)
   ensGene.transcript.exon.start <- subset(ensGene.transcript.exon.chr, pos >= exonStart)
   ensGene.transcript.exon.start.end <- subset(ensGene.transcript.exon.start, pos <= exonEnd)
 
   if (nrow(ensGene.transcript.exon.start.end) > 0)
      return(T)
   return(F)
}

# -----------------------------------------------------------------------------
# Method: Keep only SNVs on expressed genes (Step 4.1)
# Last Modified: 25/01/18
# -----------------------------------------------------------------------------
getTxQ4 <- function(tpm.gene.log2, gene.list) {
   if (!is.na(gene.list)[1])
      tpm.gene.log2 <- tpm.gene.log2[gene.list,]
   tpm.gene.log2 <- tpm.gene.log2[order(tpm.gene.log2$MEDIAN),]
   
   q <- quantile(tpm.gene.log2$MEDIAN)
   tx.q4 <- list()
   tx.q4[[4]] <- rownames(subset(tpm.gene.log2, MEDIAN > as.numeric(q[4])))
   tx.q4[[3]] <- rownames(subset(subset(tpm.gene.log2, MEDIAN > as.numeric(q[3])), MEDIAN <= as.numeric(q[4])))
   tx.q4[[2]] <- rownames(subset(subset(tpm.gene.log2, MEDIAN > as.numeric(q[2])), MEDIAN <= as.numeric(q[3])))
   tx.q4[[1]] <- rownames(subset(tpm.gene.log2, MEDIAN <= as.numeric(q[2])))
 
   return(tx.q4)
}

getLengthQ4 <- function(ensGene, gene.list) {
   ensGene$LENGTH <- abs(ensGene$start_position - ensGene$end_position)
 
   if (!is.na(gene.list)[1])
      ensGene <- ensGene[gene.list,]
   ensGene <- ensGene[order(ensGene$LENGTH),]
 
   q <- quantile(ensGene$LENGTH)
   tx.q4 <- list()
   tx.q4[[4]] <- rownames(subset(ensGene, LENGTH > as.numeric(q[4])))
   tx.q4[[3]] <- rownames(subset(subset(ensGene, LENGTH > as.numeric(q[3])), LENGTH <= as.numeric(q[4])))
   tx.q4[[2]] <- rownames(subset(subset(ensGene, LENGTH > as.numeric(q[2])), LENGTH <= as.numeric(q[3])))
   tx.q4[[1]] <- rownames(subset(ensGene, LENGTH <= as.numeric(q[2])))
 
   return(tx.q4)
}


getTxQ4Fixed <- function(tpm.gene.log2, tx.q4.fix) {
   tx.q4 <- list()
   for (q in 1:4)
      tx.q4[[q]] <- intersect(tx.q4.fix[[q]], rownames(tpm.gene.log2))

   return(tx.q4)
}

# -----------------------------------------------------------------------------
# Method: SNV asymetrey (Step 5.1)
# Last Modified: 25/01/18
# -----------------------------------------------------------------------------
getMutGenes <- function(s1, tpm.gene) {
   return(intersect(unique(s1$ensembl_gene_id), rownames(tpm.gene)))
} 

getMutPerMb <- function(s1, tpm.gene) {
   txs <- getMutGenes(s1, tpm.gene)
   s1 <- subset(s1, ensembl_gene_id %in% txs)
   
   mb <- 0
   for (t in 1:length(txs))
      mb <- mb + abs(ensGene[txs[t],]$end_position - ensGene[txs[t],]$start_position)
 
   return(nrow(s1) / mb * 1000000)
}

getMain <- function(rownames) {
   return(paste0(rownames[1], " vs ", rownames[2]))
}

getLog2Colours <- function(q4) {
   cols <- c()
   for (c in 1:ncol(q4))
      if (log2(q4[1,c]/q4[2,c]) < 0)
         cols <- c(cols, "lightskyblue")
      else
         cols <- c(cols, "sandybrown")  
  
   return(cols)
}

# -----------------------------------------------------------------------------
# Method: Transcription-associated SNV asymmetry (Step 5.2 and 8.1)
# Last Modified: 01/02/18
# -----------------------------------------------------------------------------
getLengthTx <- function(tx) {
   return(abs(ensGene[tx,]$end_position - ensGene[tx,]$start_position))
}

getMutPerMbTx <- function(s1, tx) {
   s1 <- subset(s1, ensembl_gene_id == tx)
   mb <- getLengthTx(tx)
 
   return(nrow(s1) / mb * 1000000)
}

getMutPerMbTxs <- function(s1, txs) {
   #s1 <- subset(s1, ensembl_gene_id %in% txs)          ## BUG?
   #txs <- intersect(unique(s1$ensembl_gene_id), txs)   ## BUG?
   txs <- intersect(unique(s1$ensembl_gene_id), txs)    ## SWAP 23/06/18
   s1 <- subset(s1, ensembl_gene_id %in% txs)           ## SWAP 23/06/18
   
   mb <- 0
   for (t in 1:length(txs))
      mb <- mb + abs(ensGene[txs[t],]$end_position - ensGene[txs[t],]$start_position)
 
   return(nrow(s1) / mb * 1000000)
}

# =============================================================================
# Inner Class  : Replicative Strand Asymmetry
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 29/01/18
# =============================================================================

# -----------------------------------------------------------------------------
# Method: Divide genes into two groups (hand-on and co-directional; Step 6.2)
# Last Modified: 20/02/18
# -----------------------------------------------------------------------------
getHeadOnCollision <- function(ensGene.tx.rt, headon) {
   if (headon == "ho")
      return(rbind(subset(subset(ensGene.tx.rt, strand == -1), RT == 1), subset(subset(ensGene.tx.rt, strand == 1), RT == -1)))
   else
      return(rbind(subset(subset(ensGene.tx.rt, strand == -1), RT == -1), subset(subset(ensGene.tx.rt, strand == 1), RT == 1)))
}

getTxQ4RT <- function(ensGene.tx.rt, headon, tpm.gene.log2) {
   gene.list <- rownames(getHeadOnCollision(ensGene.tx.rt, headon))
   #tpm.gene.log2 <- tpm.gene.log2[gene.list,]   ## REMOVED 07/10/18
   
   return(getTxQ4(tpm.gene.log2, gene.list))
}

getStrand <- function(genes.list, strand) {
   if (strand == "fw")
      return(rownames(subset(ensGene[genes.list,], strand == 1)))
   else
      return(rownames(subset(ensGene[genes.list,], strand == -1)))
}

getTxQ4RTStrand <- function(ensGene.tx.rt.st, headon, strand, tpm.gene.log2) {
   gene.list <- rownames(getHeadOnCollision(ensGene.tx.rt.st, headon))
   gene.list <- getStrand(gene.list, strand)
   #tpm.gene.log2 <- tpm.gene.log2[gene.list,]
 
   return(getTxQ4(tpm.gene.log2, gene.list))
}

getTxRTStrand <- function(ensGene.tx.rt.st, headon, strand) {
   gene.list <- rownames(getHeadOnCollision(ensGene.tx.rt.st, headon))
   gene.list <- getStrand(gene.list, strand)
 
   return(gene.list)
}

# -----------------------------------------------------------------------------
# Method: Divide genes into two groups (hand-on and co-directional; Step 7)
# Last Modified: 05/04/18
# -----------------------------------------------------------------------------
getRTTxQ4 <- function(j, i) {
   q4 <- q4s.rt[[j]][[i]]
   if (j == 1)
      colnames(q4) <- c("", "            Tx", "", "")
   if (j == 2)
      colnames(q4) <- c("                                     RT-Tx (Co-directional)", "", "", "")
   if (j == 3)
      colnames(q4) <- c("                                     RT-Tx (Head-on)", "", "", "")
   
   return(q4)
}

getRTTxMain <- function(i, asyms) {
   main <- getMain(rownames(asyms[[i]]))
   if (j == 2)
      main <- "Leading strand"
   else if (j == 3)
      main <- "Lagging strand"
 
   return(main)
}

getRTTxQ4_2 <- function(q4s.rt.ho, j, i, headon) {
   q4 <- q4s.rt.ho[[j]][[i]]
   if (headon == "ho") {
      if (j == 1)
         colnames(q4) <- c("                                     RT-Tx (Head-on)", "", "", "")
      if (j == 2)
         colnames(q4) <- c("                                      RT(L)-Tx(+)", "", "", "")
      if (j == 3)
         colnames(q4) <- c("                                      RT(R)-Tx(-)", "", "", "")
   } else {
      if (j == 1)
         colnames(q4) <- c("                                     RT-Tx (Co-directional)", "", "", "")
      if (j == 2)
         colnames(q4) <- c("                                      RT(R)-Tx(+)", "", "", "")
      if (j == 3)
         colnames(q4) <- c("                                      RT(L)-Tx(-)", "", "", "")
   }
   
   return(q4)
}

getRTTxMain_2 <- function(i, asyms, headon) {
   main <- getMain(rownames(asyms[[i]]))
   if (j == 1)
      if (headon == "ho")
         main <- "Lagging strand"
      else
         main <- "Leading strand"
   if (j == 2)
      main <- "Forward transcription"
   if (j == 3)
      main <- "Reverse transcription"
 
   return(main)
}

# -----------------------------------------------------------------------------
# Method: Final reports (s6.rt.st; Step 7)
# Last Modified: 05/04/18
# -----------------------------------------------------------------------------
initTableS6RTSt <- function() {
   return(list(list(list(list(list(), list())), list(list(list(), list()))), list(list(list(list(), list())), list(list(list(), list()))), list(list(list(list(), list())), list(list(list(), list()))), list(list(list(list(), list())), list(list(list(), list()))), list(list(list(list(), list())), list(list(list(), list()))), list(list(list(list(), list())), list(list(list(), list())))))
}








# =============================================================================
# Inner Class: Collections of test/obsolete/deprecated methods
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 01/02/18
# =============================================================================

# -----------------------------------------------------------------------------
# Method: Step 4.2
# Last Modified: 16/03/18
# -----------------------------------------------------------------------------
# getTxQ4Exon <- function(tx.q4, ensGene.tx.exon) {
#    tx.q4.exon <- list()
#    tx.q4.exon[[4]] <- intersect(tx.q4[[4]], rownames(ensGene.tx.exon))
#    tx.q4.exon[[3]] <- intersect(tx.q4[[3]], rownames(ensGene.tx.exon))
#    tx.q4.exon[[2]] <- intersect(tx.q4[[2]], rownames(ensGene.tx.exon))
#    tx.q4.exon[[1]] <- intersect(tx.q4[[1]], rownames(ensGene.tx.exon))
#  
#   return(tx.q4.exon)
# }

#getMutPerMb0 <- function(s1) {
#   txs <- unique(s1$ensembl_gene_id)
#   mpm <- 0
#   for (t in 1:length(txs))
#      mpm <- mpm + ((nrow(subset(s1, ensembl_gene_id == txs[t])) / abs(ensGene[txs[t],]$end_position - ensGene[txs[t],]$start_position)) * 1000000)
# 
#   return(mpm / length(txs))
#}

#getMutPerMbTx0 <- function(s1, txs) {
#   mpm <- 0
#   for (t in 1:length(txs))
#      mpm <- mpm + ((nrow(subset(s1, ensembl_gene_id == txs[t])) / abs(ensGene[txs[t],]$end_position - ensGene[txs[t],]$start_position)) * 1000000)
# 
#   return(mpm / length(txs))
#}
