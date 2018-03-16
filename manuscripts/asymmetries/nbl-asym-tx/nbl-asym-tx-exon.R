# =============================================================================
# Manuscript   : The dangeous case of transcription in NBL
# Chapter II   : 
# Name         : manuscripts/asymmetries/nbl-asym-tx-exon.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 16/03/18
# =============================================================================
#wd.source <- "/projects/cangen/tyang2/dev/R"              ## tyang2@cheops
#wd.source <- "/ngs/cangen/tyang2/dev/R"                   ## tyang2@gauss
wd.source <- "/Users/tpyang/Work/dev/R"                    ## tpyang@localhost

handbooks <- c("Common.R", "Asymmetry.R", "ReplicationTiming.R", "DifferentialExpression.R")   ## Required handbooks/libraries for the manuscript
invisible(sapply(handbooks, function(b) source(file.path(wd.source, "handbook-of", b))))

load(file.path(wd.source, "guide-to-the", "hg19.RData"))   ## The bioinformatician's guide to the reference genome

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 29/01/18
# -----------------------------------------------------------------------------
BASE <- "NBL"
#wd     <- paste0("/ngs/cangen/tyang2/", BASE, "/analysis/")                   ## tyang2@gauss
#wd.ngs <- paste0("/ngs/cangen/tyang2/", BASE, "/ngs/WGS/")
wd     <- paste0("/Users/tpyang/Work/uni-koeln/tyang2/", BASE, "/analysis/")   ## tpyang@localhost
wd.ngs <- paste0("/Users/tpyang/Work/uni-koeln/tyang2/", BASE, "/ngs/WGS/")

wd.asym       <- paste0(wd, "asymmetries/", tolower(BASE), "-wgs-asym/")
wd.asym.data  <- paste0(wd.asym, "data/")
wd.asym.plots <- paste0(wd.asym, "plots/")
setwd(wd.asym)

samples <- readTable(paste0(wd.ngs, "nbl_wgs_n56-1.list"), header=F, rownames=F, sep="")

# -----------------------------------------------------------------------------
# Step 8: Read in all the SNVs
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

###
##
load(file=paste0(wd.asym.data, "nbl_asym_tx_snv.RData"))   ## All SNVs on "expressed" genes, regardless protein coding or not
ens.tx.pcg.snv <- intersect(unique(tx.snv$ensembl_gene_id), rownames(tpm.gene.nbl.pcg))   ## From Step 4.0; ADD 16/03/18
tx.pcg.snv <- subset(tx.snv, ensembl_gene_id %in% ens.tx.pcg.snv)
# > nrow(tx.snv)
# [1] 47207
# > nrow(tx.pcg.snv)
# [1] 45909

ensGene.transcript.exon <- subset(ensGene.transcript.exon, exonFrame != -1)   ## exonFrames value of -1 means that the exon is UTR
# > nrow(ensGene.transcript.exon)
# [1] 667505

tx.pcg.snv$exon <- F
tx.pcg.snv$exon <- mapply(x = 1:nrow(tx.pcg.snv), function(x) isSNVonExon(tx.pcg.snv$CHROM[x], tx.pcg.snv$POS[x], ensGene.transcript.exon))
save(tx.pcg.snv, file=paste0(wd.asym.data, "nbl_asym_tx_pcg_snv_exon.RData"))   ## All SNVs on expressed, protein coding genes
#for (x in 1:nrow(tx.pcg.snv))
#   tx.pcg.snv$exon[x] <- isSNVonExon(tx.pcg.snv$CHROM[x], tx.pcg.snv$POS[x], ensGene.transcript.exon)
# > nrow(subset(tx.pcg.snv, exon == T))
# [1] 2418
# > nrow(subset(tx.pcg.snv, exon == F))
# [1] 43491
# > 2418/43491
# [1] 0.05559771

# > length(unique(tx.pcg.snv$ensembl_gene_id))
# [1] 10083
# > length(unique(subset(tx.pcg.snv, exon == T)$ensembl_gene_id))
# [1] 2076
# > 2076/10083
# [1] 0.2058911
# > length(unique(subset(tx.pcg.snv, exon == T)$SAMPLE))
# [1] 55

# -----------------------------------------------------------------------------
# Step 2: Divide SNVs into two groups (exons and introns)
# Last Modified: 16/03/18
# -----------------------------------------------------------------------------
getTxQ4Exon <- function(tx.q4, ensGene.tx.exon) {
   tx.q4.exon <- list()
   tx.q4.exon[[4]] <- intersect(tx.q4[[4]], rownames(ensGene.tx.exon))
   tx.q4.exon[[3]] <- intersect(tx.q4[[3]], rownames(ensGene.tx.exon))
   tx.q4.exon[[2]] <- intersect(tx.q4[[2]], rownames(ensGene.tx.exon))
   tx.q4.exon[[1]] <- intersect(tx.q4[[1]], rownames(ensGene.tx.exon))
 
   return(tx.q4.exon)
}

exon <- F
prefix <- "intron"
tx.input <- subset(tx.pcg.snv, exon == F)

ens.tx.input <- unique(tx.input$ensembl_gene_id)
# > length(ens.tx.input)   ## Exon
# [1] 2076
# > length(ens.tx.input)   ## Intron
# [1] 9482
ensGene.protein_coding <- subset(ensGene, gene_biotype == "protein_coding")
length(intersect(ens.tx.input, ensGene.protein_coding$ensembl_gene_id))
# [1] 2076

ensGene.tx.input <- ensGene[ens.tx.input,]
tx.q4.input <- getTxQ4Exon(tx.q4, ensGene.tx.input)
### 
##
for (q in 1:4)   ## Exons
   print(length(intersect(ens.tx.input, tx.q4[[q]])))
# [1] 565
# [1] 561
# [1] 479
# [1] 471

for (q in 1:4)   ## Introns
   print(length(intersect(ens.tx.input, tx.q4[[q]])))
# [1] 2311
# [1] 2522
# [1] 2449
# [1] 2200

### 
## (Figure S4)
q4s <- list()
pdf(paste0(wd.asym.plots, "nbl_asym_tx_pcg_snv_q4s_", prefix, ".pdf"), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   idx <- idxs[i]
 
   q4 <- as.matrix(toTable(0, 4, 2, c("n=", "n=", "n=", "n=")))
   rownames(q4) <- c(paste0(REFS[idx], "ntx:", REFS[idx+1], "tx"), paste0(REFS[idx+1], "ntx:", REFS[idx], "tx"))
   for (q in 1:4) {
      txs <- tx.q4.input[[q]]
  
      ## Cntx:Gtx  =   tx.snv.s6[[1]][[1]][[1]]) C>A Tx(+)  +  tx.snv.s6[[1]][[2]][[2]] G>T Tx(-)
      CntxGtx <- rbind(tx.snv.s6[[i]][[1]][[1]],               tx.snv.s6[[i]][[2]][[2]])
      #CntxGtx <- subset(CntxGtx, ensembl_gene_id %in% txs)
      #q4[1, q] <- getMutPerMbTxs(bp, unique(CntxGtx$ensembl_gene_id))
      txs.cg <- intersect(unique(CntxGtx$ensembl_gene_id), txs)
      #q4[1, q] <- getMutPerMbTxs(tx.input, txs.cg)
      q4[1, q] <- getMutPerMbTxs(CntxGtx, txs.cg)
      
      ## Gntx:Ctx  =   tx.snv.s6[[1]][[2]][[1]]) G>T Tx(+)  +  tx.snv.s6[[1]][[1]][[2]] C>A Tx(-)       
      GntxCtx <- rbind(tx.snv.s6[[i]][[2]][[1]],               tx.snv.s6[[i]][[1]][[2]])
      #GntxCtx <- subset(GntxCtx, ensembl_gene_id %in% txs) 
      #q4[2, q] <- getMutPerMbTxs(bp, unique(GntxCtx$ensembl_gene_id))
      txs.gc <- intersect(unique(GntxCtx$ensembl_gene_id), txs)
      #q4[2, q] <- getMutPerMbTxs(tx.input, txs.gc)
      q4[2, q] <- getMutPerMbTxs(GntxCtx, txs.gc)
        
      colnames(q4)[q] <- paste0(colnames(q4)[q], length(unique(c(txs.cg, txs.gc))), " (", length(txs.cg), "+", length(txs.gc), ")")
   }
 
   q4s[[i]] <- q4
   barplot(q4, ylab="SNVs/Mb", main=paste(rownames(q4), collapse="/"), beside=TRUE, width=.3, col=c("#E6E6E6", "#4D4D4D"))
}
dev.off()
save(q4s, file=paste0(wd.asym.data, "nbl_asym_tx_pcg_snv_q4s_", prefix, ".RData"))

## ADD 16/03/18
#q4s.exon[[2]] <- q4s
q4s.exon[[3]] <- q4s
save(q4s.exon, file=paste0(wd.asym.data, "nbl_asym_tx_pcg_snv_q4s_exon.RData"))

###
## Refining the plot
pdf(paste0(wd.asym.plots, "nbl_asym_tx_pcg_snv_q4s_", prefix, ".pdf"), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   q4 <- q4s[[i]]
   if (exon)
      colnames(q4) <- c("                                     Tx (Exon)", "", "", "")
   else
      colnames(q4) <- c("                                     Tx (Intron)", "", "", "")
   #for (c in 1:ncol(q4))
   #   colnames(q4)[c] <- gsub("n=", "", unlist(strsplit(colnames(q4)[c], " "))[1])
 
   barplot(q4, ylab="SNVs/Mb", main=getMain(rownames(asyms[[i]])), beside=TRUE, width=.3, col=c("#E6E6E6", "#4D4D4D"))   ##, xlab="Number of genes")
   mtext(paste(rownames(q4), collapse=" vs "), cex=0.55, font=3, line=0.5)
}
dev.off()

# -----------------------------------------------------------------------------
# Step 6.2: Log2 transcription-associated SNV asymmetry (Figure S5)
# Last Modified: 26/01/18
# -----------------------------------------------------------------------------
pdf(paste0(wd.asym.plots, "nbl_asym_tx_pcg_snv_q4s_", prefix, "_log2_ylim0.5.pdf"), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   q4 <- q4s[[i]]
   if (exon)
      colnames(q4) <- c("                                    Tx (Exon)", "", "", "")
   else
      colnames(q4) <- c("                                    Tx (Intron)", "", "", "")
 
   cols <- c()   
   for (c in 1:ncol(q4))
      if (log2(q4[1,c]/q4[2,c]) > 0)
         cols <- c(cols, "#4D4D4D")
      else
         cols <- c(cols, "#E6E6E6")    
 
   barplot(log2(q4[1,]/q4[2,]), ylab="Asymmetry", ylim=c(-0.5, 0.5), main=getMain(rownames(asyms[[i]])), beside=TRUE, width=.3, col=cols)
   mtext(paste0("log2(", paste(rownames(q4), collapse="/"), ")"), cex=0.55, font=3, line=0.5)
}
dev.off()

# -----------------------------------------------------------------------------
# Step 7: Replication–transcription collision-induced SNV asymmetry (Figure 1)
# Last Modified: 22/02/18
# -----------------------------------------------------------------------------
getQ4 <- function(j, i) {
   q4 <- q4s.exon[[j]][[i]]
   if (j == 1)
      colnames(q4) <- c("", "            Tx", "", "")
   else if (j == 2)
      colnames(q4) <- c("                                     Tx (Exon)", "", "", "")
   else if (j ==3)
      colnames(q4) <- c("                                     Tx (Intron)", "", "", "")
   
   return(q4)
}

###
##
idx <- 0
for (i in 1:6) {
   pdf(paste0(wd.asym.plots, "nbl_asym_tx_pcg_snv_q4s_exonintron_", REFS[i+idx], ">", ALTS[i+idx], ".pdf"), height=4.5, width=7)
   par(mfrow=c(2, 3)) 
   
   for (j in 1:3) {
      q4 <- getQ4(j, i)
      
      #barplot(q4, ylab="SNVs/Mb", main=getMain(rownames(asyms[[i]])), beside=TRUE, width=.3, col=c("#E6E6E6", "#4D4D4D"), ylim=c(0, q4s.exon[[2]][[i]][2, 1]))   ##, xlab="Number of genes")   ## For C>A only
      #barplot(q4, ylab="SNVs/Mb", main=getMain(rownames(asyms[[i]])), beside=TRUE, width=.3, col=c("#E6E6E6", "#4D4D4D"), ylim=c(0, q4s.exon[[1]][[i]][1, 4]))   ##, xlab="Number of genes")   ## For T>A
      barplot(q4, ylab="SNVs/Mb", main=getMain(rownames(asyms[[i]])), beside=TRUE, width=.3, col=c("#E6E6E6", "#4D4D4D"), ylim=c(0, q4s.exon[[1]][[i]][2, 4]))   ##, xlab="Number of genes")   ## For T>C
      #barplot(q4, ylab="SNVs/Mb", main=getMain(rownames(asyms[[i]])), beside=TRUE, width=.3, col=c("#E6E6E6", "#4D4D4D"), ylim=c(0, q4s.exon[[2]][[i]][1, 4]))   ##, xlab="Number of genes")
      mtext(paste(rownames(q4), collapse=" vs "), cex=0.55, font=3, line=0.5)      
   }
 
   for (j in 1:3) {
      q4 <- getQ4(j, i)
  
      cols <- c()   
      for (c in 1:ncol(q4))
         if (log2(q4[1,c]/q4[2,c]) > 0)
            cols <- c(cols, "#4D4D4D")
         else
            cols <- c(cols, "#E6E6E6")    
      
      barplot(log2(q4[1,]/q4[2,]), ylab="Asymmetry", ylim=c(-0.5, 0.5), main=getMain(rownames(asyms[[i]])), beside=TRUE, width=.3, col=cols)
      mtext(paste0("log2(", paste(rownames(q4), collapse="/"), ")"), cex=0.55, font=3, line=0.5)
   }
   
   dev.off()
   idx <- idx + 1
}














# -----------------------------------------------------------------------------
# Step 6: Read in all the SVs
# Last Modified: 31/01/18
# -----------------------------------------------------------------------------
sv.bp.tpm.cna <- NA
for (s in 1:length(samples)) {
   sample <- samples[s]
 
   ## INPUT: INV
   #sv.file <- paste0(wd.ngs, sample, "/", sample, "_ANALYSIS/", sample, "_unsnarl.txt")
   sv.file <- paste0(wd.ngs, "rearrangements/", sample, "_unsnarl.txt")
   sv.bp <- initBreakPoint(read.peiflyne.unsnarl.txt(sv.file, cutoff=10))   ## Filter 1: TO-DO
   #sv.bp <- subset(sv.bp, Length_BP >= 10000)                              ## Filter 2: TO-DO
 
   ## INPUT: Expression
   sv.bp.tpm <- setBreakPointTranscriptionEnsGene(sample, sv.bp,     tpm.gene.sclc.log2, ensGene.tx.rt, pair=1)
   sv.bp.tpm <- setBreakPointTranscriptionEnsGene(sample, sv.bp.tpm, tpm.gene.sclc.log2, ensGene.tx.rt, pair=2)
 
   ## INPUT: CNA
   cna.file <- paste0(wd.ngs, sample, "/", sample, "_ANALYSIS/", sample, "_iCN.seg")
   sv.bp.tpm <- setBreakPointCNA(sample, sv.bp.tpm, read.peiflyne.iCN.seg(cna.file), pair=1)
   sv.bp.tpm <- setBreakPointCNA(sample, sv.bp.tpm, read.peiflyne.iCN.seg(cna.file), pair=2)
 
   ##
   if (is.na(sv.bp.tpm.cna))
      sv.bp.tpm.cna <- sv.bp.tpm
   else
      sv.bp.tpm.cna <- rbind(sv.bp.tpm.cna, sv.bp.tpm)
} 

break1 <- sv.bp.tpm.cna[,c("BP_ENS_ID_1", "BP_iCN_1")]   ## TO-DO 06/02/18: add sample IDs!!!
break1 <- subset(break1, !is.na(BP_ENS_ID_1))
colnames(break1) <- c("ensembl_gene_id", "iCN")

break2 <- sv.bp.tpm.cna[,c("BP_ENS_ID_2", "BP_iCN_2")]
break2 <- subset(break2, !is.na(BP_ENS_ID_2))
colnames(break2) <- c("ensembl_gene_id", "iCN")

tx.bp <- rbind(break1, break2)

###
##
ens.tx.pcg.snv.bp <- intersect(unique(tx.bp$ensembl_gene_id), unique(tx.pcg.snv$ensembl_gene_id))
tx.pcg.snv.bp <- subset(tx.bp, ensembl_gene_id %in% ens.tx.pcg.snv.bp)
# > nrow(tx.pcg.snv.bp)
# [1] 10639

ens.tx.pcg.snv.indel.bp <- intersect(unique(tx.bp$ensembl_gene_id), unique(tx.pcg.snv.indel$ensembl_gene_id))
tx.pcg.snv.indel.bp <- subset(tx.bp, ensembl_gene_id %in% ens.tx.pcg.snv.indel.bp)
# > nrow(tx.pcg.snv.indel.bp)
# [1] 10303









###
##
ens.tx.pcg.snv.indel.bp <- intersect(unique(tx.bp$ensembl_gene_id), unique(tx.pcg.snv.indel$ensembl_gene_id))
tx.pcg.snv.indel.bp <- subset(tx.bp, ensembl_gene_id %in% ens.tx.pcg.snv.indel.bp)
# > nrow(tx.pcg.snv.indel.bp)
# [1] 10303

nhej <- setdiff(tx.pcg.snv.bp[,1], tx.pcg.snv.indel.bp[,1])
for (q in 1:4)
   print(length(intersect(nhej, tx.q4[[q]])))
# [1] 41
# [1] 47
# [1] 59
# [1] 60
# > length(nhej)
# [1] 207

### 
##
test <- setdiff(ens.tx.pcg.snv, ens.tx.pcg.snv.indel)
for (q in 1:4)
   print(length(intersect(test, tx.q4[[q]])))
# [1] 959
# [1] 706
# [1] 666
# [1] 946

for (q in 1:4)
   print(length(intersect(ens.tx.pcg.snv.indel.bp, tx.q4[[q]])))






# -----------------------------------------------------------------------------
# Step 7:
# Last Modified: 30/01/18
# -----------------------------------------------------------------------------
getHeadOnCollision <- function(ensGene.tx.rt, headon) {
   if (headon)
      return(rbind(subset(subset(ensGene.tx.rt, strand == -1), SLOPE_START > 0), subset(subset(ensGene.tx.rt, strand == 1), SLOPE_END < 0)))
   else
      return(rbind(subset(subset(ensGene.tx.rt, strand == -1), SLOPE_END < 0), subset(subset(ensGene.tx.rt, strand == 1), SLOPE_START > 0)))
}

getHeadOnCollision0 <- function(ensGene.tx.rt, headon) {
   if (headon)
      return(rbind(subset(subset(ensGene.tx.rt, SLOPE_START > 0), strand == -1), subset(subset(ensGene.tx.rt, SLOPE_START < 0), strand == 1)))
   else
      return(rbind(subset(subset(ensGene.tx.rt, SLOPE_START > 0), strand == 1), subset(subset(ensGene.tx.rt, SLOPE_START < 0), strand == -1)))
}

getTxQ4RT <- function(tx.q4, ensGene.tx.rt, headon) {
   tx.q4.rt <- list()
   tx.q4.rt[[4]] <- intersect(tx.q4[[4]], rownames(getHeadOnCollision(ensGene.tx.rt, headon)))
   tx.q4.rt[[3]] <- intersect(tx.q4[[3]], rownames(getHeadOnCollision(ensGene.tx.rt, headon)))
   tx.q4.rt[[2]] <- intersect(tx.q4[[2]], rownames(getHeadOnCollision(ensGene.tx.rt, headon)))
   tx.q4.rt[[1]] <- intersect(tx.q4[[1]], rownames(getHeadOnCollision(ensGene.tx.rt, headon)))
 
   return(tx.q4.rt)
}

#getLaggingStrandColour <- function(headon) {
#   if (headon)
#      return(c("salmon", "#E6E6E6"))
#   else
#      return(c("#E6E6E6", "salmon"))
#}

###
##
ensGene.tx.rt.nona <- subset(subset(ensGene.tx.rt, !is.na(SLOPE_START)), !is.na(SLOPE_END))
ensGene.tx.rt.nona$SIGN <- sign(ensGene.tx.rt.nona$SLOPE_START / ensGene.tx.rt.nona$SLOPE_END)
ensGene.tx.rt.nona.sign <- subset(ensGene.tx.rt.nona, SIGN == 1)[,-10]
# > nrow(ensGene.tx.rt.nona)
# [1] 18093
# > nrow(ensGene.tx.rt.nona.sign)
# [1] 16896

tx.input <- tx.pcg.snv.indel.bp
headon <- F   ## Replication–transcription Head-on (T) / Co-directional (F) collisions
prefix <- "cd"
ens.tx.pcg.snv.indel.bp.rt <- intersect(unique(tx.input$ensembl_gene_id), unique(ensGene.tx.rt.nona.sign$ensembl_gene_id))
ensGene.tx.rt.input <- ensGene.tx.rt.nona.sign[ens.tx.pcg.snv.indel.bp.rt,]
tx.q4.rt.input <- getTxQ4RT(tx.q4, ensGene.tx.rt.input, headon=headon)
# > length(unique(tx.input$ensembl_gene_id))
# [1] 2885
# > nrow(ensGene.tx.rt.input)
# [1] 2483

q4s <- list()
pdf(paste0(wd.asym.plots, "sclc_asym_tx_rt_snv_indel_bp_q4s_", prefix, ".pdf"), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   idx <- idxs[i]
 
   q4 <- as.matrix(toTable(0, 4, 2, c("n=", "n=", "n=", "n=")))
   rownames(q4) <- c(paste0(REFS[idx], "ntx:", REFS[idx+1], "tx"), paste0(REFS[idx+1], "ntx:", REFS[idx], "tx"))
   for (q in 1:4) {
      txs <- tx.q4.rt.input[[q]]
      
      ## Cntx:Gtx  =   tx.snv.s6[[1]][[1]][[1]]) C>A Tx(+)  +  tx.snv.s6[[1]][[2]][[2]] G>T Tx(-)
      CntxGtx <- rbind(tx.snv.s6[[i]][[1]][[1]],               tx.snv.s6[[i]][[2]][[2]])
      #CntxGtx <- subset(CntxGtx, ensembl_gene_id %in% txs)
      #q4[1, q] <- getMutPerMbTxs(bp, unique(CntxGtx$ensembl_gene_id))
      txs.cg <- intersect(unique(CntxGtx$ensembl_gene_id), txs)
      q4[1, q] <- getMutPerMbTxs(tx.input, txs.cg)

      ## Gntx:Ctx  =   tx.snv.s6[[1]][[2]][[1]]) G>T Tx(+)  +  tx.snv.s6[[1]][[1]][[2]] C>A Tx(-)       
      GntxCtx <- rbind(tx.snv.s6[[i]][[2]][[1]],               tx.snv.s6[[i]][[1]][[2]])
      #GntxCtx <- subset(GntxCtx, ensembl_gene_id %in% txs) 
      #q4[2, q] <- getMutPerMbTxs(bp, unique(GntxCtx$ensembl_gene_id))
      txs.gc <- intersect(unique(GntxCtx$ensembl_gene_id), txs)
      q4[2, q] <- getMutPerMbTxs(tx.input, txs.gc)
      
      colnames(q4)[q] <- paste0(colnames(q4)[q], length(unique(c(txs.cg, txs.gc))), " (", length(txs.cg), "+", length(txs.gc), ")")
   }
   
   q4s[[i]] <- q4
   barplot(q4, ylab="Breaks/Mb", main=paste(rownames(q4), collapse="/"), beside=TRUE, width=.3, col=c("#E6E6E6", "#4D4D4D"))
}
dev.off()
save(q4s, file=paste0(wd.asym.data, "sclc_asym_tx_rt_snv_indel_bp_q4s_", prefix, ".RData"))

###
## Refining the plot
pdf(paste0(wd.asym.plots, "sclc_asym_tx_rt_snv_indel_bp_q4s_", prefix, ".pdf"), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   q4 <- q4s[[i]]
   if (headon)
      colnames(q4) <- c("                                     RT-Tx (Head-on)", "", "", "")
   else
      colnames(q4) <- c("                                     RT-Tx (Co-directional)", "", "", "")
   #for (c in 1:ncol(q4))
   #   colnames(q4)[c] <- gsub("n=", "", unlist(strsplit(colnames(q4)[c], " "))[1])
 
   barplot(q4, ylab="Breaks/Mb", main=getMain(rownames(asyms[[i]])), beside=TRUE, width=.3, col=c("#E6E6E6", "#4D4D4D"))   ##, xlab="Number of genes")
   mtext(paste(rownames(q4), collapse=" vs "), cex=0.55, font=3, line=0.5)
}
dev.off()

# -----------------------------------------------------------------------------
# Step 8: Log2 transcription-associated SNV asymmetry (Figure 3)
# Last Modified: 26/01/18
# -----------------------------------------------------------------------------
pdf(paste0(wd.asym.plots, "sclc_asym_tx_rt_snv_indel_bp_q4s_", prefix, "_log2_ylim0.2.pdf"), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   q4 <- q4s[[i]]
   if (headon)
      colnames(q4) <- c("                                     RT-Tx (Head-on)", "", "", "")
   else
      colnames(q4) <- c("                                     RT-Tx (Co-directional)", "", "", "")
   
   cols <- c()   
   for (c in 1:ncol(q4))
      if (log2(q4[1,c]/q4[2,c]) > 0)
         cols <- c(cols, "#4D4D4D")
      else
         cols <- c(cols, "#E6E6E6")    
   
   barplot(log2(q4[1,]/q4[2,]), ylab="Asymmetry", ylim=c(-0.2, 0.2), main=getMain(rownames(asyms[[i]])), beside=TRUE, width=.3, col=cols)
   mtext(paste0("log2(", paste(rownames(q4), collapse="/"), ")"), cex=0.55, font=3, line=0.5)
}
dev.off()







