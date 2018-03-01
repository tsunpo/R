# =============================================================================
# Manuscript   : The dangeous case of DNA replication in SCLC
# Chapter I    : Transcriptional strand asymmetry
# Name         : manuscripts/asymmetries/sclc-asym-tx.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 29/01/18
# =============================================================================
#wd.source <- "/projects/cangen/tyang2/dev/R"              ## tyang2@cheops
wd.source <- "/ngs/cangen/tyang2/dev/R"                    ## tyang2@gauss
#wd.source <- "/Users/tpyang/Work/dev/R"                   ## tpyang@localhost

handbooks <- c("Common.R", "Asymmetry.R", "ReplicationTiming.R")   ## Required handbooks/libraries for the manuscript
invisible(sapply(handbooks, function(b) source(file.path(wd.source, "handbook-of", b))))

load(file.path(wd.source, "guide-to-the", "hg19.RData"))   ## The bioinformatician's guide to the reference genome
load(file.path(wd.source, "guide-to-the", "hg19.1kb.gc.RData"))

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 29/01/18
# -----------------------------------------------------------------------------
BASE <- "SCLC"
wd     <- paste0("/ngs/cangen/tyang2/", BASE, "/analysis/")                     ## tyang2@gauss
wd.ngs <- paste0("/ngs/cangen/tyang2/", BASE, "/ngs/WGS/")
#wd     <- paste0("/Users/tpyang/Work/uni-koeln/tyang2/", BASE, "/analysis/")   ## tpyang@localhost
#wd.ngs <- paste0("/Users/tpyang/Work/uni-koeln/tyang2/", BASE, "/ngs/WGS/")

wd.asym       <- paste0(wd, "asymmetries/", tolower(BASE), "-wgs-asym/")
wd.asym.data  <- paste0(wd.asym, "data/")
wd.asym.plots <- paste0(wd.asym, "plots/")
setwd(wd.asym)

samples <- readTable(paste0(wd.ngs, "sclc_wgs_n101.list"), header=F, rownames=F, sep="")

# -----------------------------------------------------------------------------
# Step 1: Finding mutations locate within Ensembl genes
# Last Modified: 23/01/18
# -----------------------------------------------------------------------------
for (s in 1:length(samples)) {
   #sample <- SAMPLE 
   sample <- samples[s]
   vcf <- read.peiflyne.mutcall.filtered.vcf(paste0(wd.ngs, sample, "/", sample, "_ANALYSIS/", sample, "_mutcall_filtered.vcf"), pass=T, rs=F)

   colnames <- c(colnames(ensGene), colnames(vcf))
   vcf.gene <- toTable("", length(colnames), 0, colnames)
   #vcf.gene.10kb <- NA
   for (c in 1:24) {
      chr <- chrs[c]
      vcf.chr <- subset(vcf, CHROM == chr)
      ensGene.chr <- subset(ensGene, chromosome_name == chr)  
      
      vcf.gene.chr <- toTable("", length(colnames), 0, colnames)
      for (g in 1:nrow(ensGene.chr)) {
         ensGene.chr.gene <- ensGene.chr[g,]
         
         vcf.gene.chr0 <- getEnsGeneSNVs(vcf.chr, ensGene.chr.gene)
         if (nrow(vcf.gene.chr0) != 0) {
            if (nrow(vcf.gene.chr) == 0)
               vcf.gene.chr <- getMergedTable(ensGene.chr.gene, vcf.gene.chr0)
            else
               vcf.gene.chr <- rbind(vcf.gene.chr, getMergedTable(ensGene.chr.gene, vcf.gene.chr0))
            
            #vcf.chr.gene.10kb <- getEnsGeneTSS(vcf.chr, ensGene.chr.gene, 10000)
         }
      }
      
      if (nrow(vcf.gene.chr) != 0) {
         vcf.gene.chr <- vcf.gene.chr[order(vcf.gene.chr$POS),]
      
         if (nrow(vcf.gene) == 0)
            vcf.gene <- vcf.gene.chr
         else
            vcf.gene <- rbind(vcf.gene, vcf.gene.chr)
      }
   }
   
   writeTable(vcf.gene, gzfile(paste0(wd.asym.data, sample, "_mut.ens.vcf.gz")), colnames=T, rownames=F, sep="\t")
}

# -----------------------------------------------------------------------------
# Step 2: Filtering out genes that is transcribed in SCLC 
# Last Modified: 24/01/18
# -----------------------------------------------------------------------------   
load("/ngs/cangen/tyang2/dev/R/manuscripts/luad-lcnec-sclc-rnaseq-de/de-sclc-tpm/data/sclc_kallisto_0.43.1_tpm.gene_r5_p47.RData")              ## tyang2@gauss
#load("/Users/tpyang/Work/uni-koeln/tyang2/SCLC/analysis/expression/kallisto/sclc-rnaseq-de/data/sclc_kallisto_0.43.1_tpm.gene_r5_p47.RData")   ## tpyang@localhost
tpm.gene.sclc <- tpm.gene
tpm.gene.sclc.pcg <- getEnsGeneFiltered(tpm.gene.sclc, ensGene, autosomeOnly=T, proteinCodingOnly=T, withHLA=F)
tpm.gene.sclc.hla <- getEnsGeneHLA(tpm.gene.sclc, ensGene)
tpm.gene.sclc.tcr <- getEnsGeneTCR(tpm.gene.sclc, ensGene)
tpm.gene.sclc.ig  <- getEnsGeneIG(tpm.gene.sclc, ensGene)

for (s in 1:length(samples)) {
   #sample <- SAMPLE 
   sample <- samples[s]
   
   vcf <- readTable(paste0(wd.asym.data, sample, "_mut.ens.vcf.gz"), header=T, rownames=F, sep="")
   ens.tx <- intersect(unique(vcf$ensembl_gene_id), rownames(tpm.gene.sclc))   ## Get "expressed" genes
   vcf.tx <- subset(vcf, ensembl_gene_id %in% ens.tx)
   writeTable(vcf.tx, gzfile(paste0(wd.asym.data, sample, "_mut.ens.tx.vcf.gz")), colnames=T, rownames=F, sep="\t")
   
   ## Seperate SNVs
   vcf.tx.snv <- subset(vcf.tx,     REF %in% c("A", "T", "C", "G"))
   vcf.tx.snv <- subset(vcf.tx.snv, ALT %in% c("A", "T", "C", "G"))   
   writeTable(vcf.tx.snv, gzfile(paste0(wd.asym.data, sample, "_mut.ens.tx.snv.vcf.gz")), colnames=T, rownames=F, sep="\t")
   
   ## Seperate Indels
   vcf.tx$TMP <- paste0(vcf.tx$REF, vcf.tx$ALT)
   vcf.tx.indel <- vcf.tx[which(nchar(vcf.tx$TMP) != 2),]
   writeTable(vcf.tx.indel[,-16], gzfile(paste0(wd.asym.data, sample, "_mut.ens.tx.indel.vcf.gz")), colnames=T, rownames=F, sep="\t")
}

# -----------------------------------------------------------------------------
# Step 3.0: Read in all the SNVs
# Last Modified: 01/02/18
# -----------------------------------------------------------------------------
getLog2andMedian <- function(tpm.gene) {
   tpm.gene.log2 <- log2(tpm.gene + 0.01)
   tpm.gene.log2$MEDIAN <- mapply(x = 1:nrow(tpm.gene.log2), function(x) median(as.numeric(tpm.gene.log2[x, ])))
 
   return(tpm.gene.log2)
}

getTxQ4 <- function(tpm.gene.sclc.log2) {
   q <- quantile(tpm.gene.sclc.log2$MEDIAN)
 
   tx.q4 <- list()
   tx.q4[[4]] <- rownames(subset(tpm.gene.sclc.log2, MEDIAN > as.numeric(q[4])))
   tx.q4[[3]] <- rownames(subset(subset(tpm.gene.sclc.log2, MEDIAN > as.numeric(q[3])), MEDIAN <= as.numeric(q[4])))
   tx.q4[[2]] <- rownames(subset(subset(tpm.gene.sclc.log2, MEDIAN > as.numeric(q[2])), MEDIAN <= as.numeric(q[3])))
   tx.q4[[1]] <- rownames(subset(tpm.gene.sclc.log2, MEDIAN <= as.numeric(q[2])))
 
   return(tx.q4)
}

###
##
tx.snv <- getMergedReport(samples[1], readTable(paste0(wd.asym.data, samples[1], "_mut.ens.tx.snv.vcf.gz"), header=T, rownames=F, sep=""))[1,]
tx.snv <- tx.snv[-1,]
for (s in 1:length(samples)) {
   sample <- samples[s]
 
   vcf <- readTable(paste0(wd.asym.data, sample, "_mut.ens.tx.snv.vcf.gz"), header=T, rownames=F, sep="")
   vcf <- subset(vcf, CHROM %in% paste0("chr", c(1:22)))   ## Keep only autosomes
   snv <- getMergedReport(sample, vcf)
   tx.snv <- rbind(tx.snv, snv)
}
save(tx.snv, file=paste0(wd.asym.data, "sclc_asym_tx_snv.RData"))

## Keep only SNVs on protein coding genes
ens.tx.snv.pcg <- intersect(unique(tx.snv$ensembl_gene_id), rownames(tpm.gene.sclc.pcg))
tx.snv.pcg <- subset(tx.snv, ensembl_gene_id %in% ens.tx.snv.pcg)
# > nrow(tx.snv)
# [1] 1441258
# > nrow(tx.snv.pcg)
# [1] 1414362
# > nrow(tpm.gene.sclc.pcg)
# [1] 16395
# > length(ens.tx.snv.pcg)
# [1] 15928

###
## Use only protein coding genes
tpm.gene.sclc.pcg.log2 <- getLog2andMedian(tpm.gene.sclc.pcg)
tpm.gene.sclc.pcg.log2 <- tpm.gene.sclc.pcg.log2[ens.tx.snv.pcg,]

tx.q4 <- getTxQ4(tpm.gene.sclc.pcg.log2)



# -----------------------------------------------------------------------------
# Step 3.1: 
# Last Modified: 01/02/18
# -----------------------------------------------------------------------------
getTableS6 <- function(vcf) {
   colnames <- c("SAMPLE", colnames(vcf))
   tx.snv.s6 <- list(list(list(), list()), list(list(), list()), list(list(), list()), list(list(), list()), list(list(), list()), list(list(), list()))
   
   for (i in 1:6)
      for (j in 1:2)
         for (k in 1:2)
            tx.snv.s6[[i]][[j]][[k]] <- toTable("", length(colnames), 0, colnames)
   
   return(tx.snv.s6)
}

###
## Define transcriptional asymmetry
## Sig.    s1(S4)    s2        s3(S2/13) s4        s5(S5)    s6
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
tx.snv.s6 <- getTableS6(getMergedReport(samples[1], readTable(paste0(wd.asym.data, samples[1], "_mut.ens.tx.snv.vcf.gz"), header=T, rownames=F, sep="")))
for (s in 1:length(samples)) {
   sample <- samples[s]
   vcf <- readTable(paste0(wd.asym.data, sample, "_mut.ens.tx.snv.vcf.gz"), header=T, rownames=F, sep="")
   vcf <- subset(vcf, CHROM %in% paste0("chr", c(1:22)))   ## Keep only autosomes
   vcf <- getMergedReport(sample, vcf)
   
   for (i in 1:length(idxs)) {
      idx <- idxs[i]
      
      watson <- subset(subset(vcf, REF == REFS[idx]), ALT == ALTS[idx])
      tx.snv.s6[[i]][[1]][[1]] <- rbind(tx.snv.s6[[i]][[1]][[1]], subset(watson, strand == 1))
      tx.snv.s6[[i]][[1]][[2]] <- rbind(tx.snv.s6[[i]][[1]][[2]], subset(watson, strand == -1))
      
      crick  <- subset(subset(vcf, REF == REFS[idx+1]), ALT == ALTS[idx+1])
      tx.snv.s6[[i]][[2]][[1]] <- rbind(tx.snv.s6[[i]][[2]][[1]], subset(crick, strand == 1))
      tx.snv.s6[[i]][[2]][[2]] <- rbind(tx.snv.s6[[i]][[2]][[2]], subset(crick, strand == -1))
   }
}
save(tx.snv.s6, file=paste0(wd.asym.data, "sclc_asym_tx_snv_s6.RData"))

# -----------------------------------------------------------------------------
# Step 3.2: SNV asymetrey (Figure 1)
# Last Modified: 25/01/18
# -----------------------------------------------------------------------------
asyms <- list()
pdf(paste0(wd.asym.plots, "sclc_asym_tx_pcg_snv_s6_asyms.pdf"), height=4, width=6)
par(mfrow=c(2, 3))
#par(mar=c(3.5, 3.5, 1,1), mgp=c(2, 0.65, 0), cex=0.9)
for (i in 1:length(idxs)) {
   idx <- idxs[i]

   asym <- as.matrix(toTable(0, 2, 2, c("Tx(+)", "Tx(-)")))
   rownames(asym) <- c(paste0(REFS[idx], ">", ALTS[idx]), paste0(REFS[idx+1], ">", ALTS[idx+1]))
   for (j in 1:2)
      for (k in 1:2)
         asym[j, k] <- getMutPerMb(tx.snv.s6[[i]][[j]][[k]], tpm.gene.sclc.pcg)
               ## E.g. getMutPerMb(tx.snv.s6[[1]][[1]][[1]])   ## E.g. C>A Tx(+)
               ##      getMutPerMb(tx.snv.s6[[1]][[1]][[2]])   ##      C>A Tx(-)
               ##      getMutPerMb(tx.snv.s6[[1]][[2]][[1]])   ##      G>T Tx(+)
               ##      getMutPerMb(tx.snv.s6[[1]][[2]][[2]])   ##      G>T Tx(-)
   
   asyms[[i]] <- asym
   barplot(asym, ylab="SNVs/Mb", main=getMain(rownames(asym)), beside=TRUE, width=.3, col=c("#E6E6E6", "#4D4D4D", "#4D4D4D", "#E6E6E6"))
   mtext(paste(c(paste0(REFS[idx], ":", REFS[idx+1]), paste0(REFS[idx+1], ":", REFS[idx])), collapse=" vs "), cex=0.55, font=3, line=0.5)
}
dev.off()
save(tx.snv.s6, asyms, file=paste0(wd.asym.data, "sclc_asym_tx_pcg_snv_s6_asyms.RData"))

###
## Refining the plot
pdf(paste0(wd.asym.plots, "sclc_asym_tx_pcg_snv_s6_asyms_ylim300.pdf"), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   asym <- asyms[[i]]
   idx <- idxs[i]
   
   barplot(asym, ylab="SNVs/Mb", main=getMain(rownames(asym)), beside=TRUE, width=.3, col=c("#E6E6E6", "#4D4D4D", "#4D4D4D", "#E6E6E6"), ylim=c(0, 300))
   mtext(paste(c(paste0(REFS[idx], ":", REFS[idx+1]), paste0(REFS[idx+1], ":", REFS[idx])), collapse=" vs "), cex=0.55, font=3, line=0.5)
}
dev.off()

# -----------------------------------------------------------------------------
# Step 4.0: Reports
# Last Modified: 01/02/18
# -----------------------------------------------------------------------------
tpm.gene.sclc.pcg.log2 <- getLog2Median(tpm.gene.sclc.pcg)
tpm.gene.sclc.hla.log2 <- getLog2Median(tpm.gene.sclc.hla)
tpm.gene.sclc.tcr.log2 <- getLog2Median(tpm.gene.sclc.tcr)
tpm.gene.sclc.ig.log2  <- getLog2Median(tpm.gene.sclc.ig)

report <- toTable(0, 4, 4, c("TPM_MEDIAN", "SNV_MB", "INV_MB", "SV_MB"))
rownames(report) <- c("Protein_Coding", "HLA", "TCR", "IG")

median(tpm.gene.sclc.pcg.log2$MEDIAN)

getMutPerMbTxs(tx.snv, rownames(tpm.gene.sclc.pcg.log2))
getMutPerMbTxs(tx.snv, rownames(tpm.gene.sclc.hla.log2))
getMutPerMbTxs(tx.snv, rownames(tpm.gene.sclc.tcr.log2))
getMutPerMbTxs(tx.snv, rownames(tpm.gene.sclc.ig.log2))

getMutPerMbTxs(tx.indel, rownames(tpm.gene.sclc.pcg.log2))
getMutPerMbTxs(tx.indel, rownames(tpm.gene.sclc.hla.log2))
getMutPerMbTxs(tx.indel, rownames(tpm.gene.sclc.tcr.log2))
getMutPerMbTxs(tx.indel, rownames(tpm.gene.sclc.ig.log2))

getMutPerMbTxs(tx.bp, rownames(tpm.gene.sclc.pcg.log2))
getMutPerMbTxs(tx.bp, rownames(tpm.gene.sclc.hla.log2))
getMutPerMbTxs(tx.bp, rownames(tpm.gene.sclc.tcr.log2))
getMutPerMbTxs(tx.bp, rownames(tpm.gene.sclc.ig.log2))

# -----------------------------------------------------------------------------
# Step 4.1: Transcription-associated SNV asymmetry (Figure 2)
# Last Modified: 26/01/18
# -----------------------------------------------------------------------------
q4s <- list()
pdf(paste0(wd.asym.plots, "sclc_asym_tx_pcg_snv_s6_asyms_q4s.pdf"), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   idx <- idxs[i]
 
   q4 <- as.matrix(toTable(0, 4, 2, c("n=", "n=", "n=", "n=")))
   rownames(q4) <- c(paste0(REFS[idx], "ntx:", REFS[idx+1], "tx"), paste0(REFS[idx+1], "ntx:", REFS[idx], "tx"))
   for (q in 1:4) {
      txs <- tx.q4[[q]]
      ## Cntx:Gtx  =   tx.snv.s6[[1]][[1]][[1]]) C>A Tx(+)  +  tx.snv.s6[[1]][[2]][[2]] G>T Tx(-)
      CntxGtx <- rbind(tx.snv.s6[[i]][[1]][[1]],               tx.snv.s6[[i]][[2]][[2]])
      txs.cg <- intersect(unique(CntxGtx$ensembl_gene_id), txs)
      q4[1, q] <- getMutPerMbTxs(CntxGtx, txs.cg)

      ## Gntx:Ctx  =   tx.snv.s6[[1]][[2]][[1]]) G>T Tx(+)  +  tx.snv.s6[[1]][[1]][[2]] C>A Tx(-)       
      GntxCtx <- rbind(tx.snv.s6[[i]][[2]][[1]],               tx.snv.s6[[i]][[1]][[2]])
      txs.gc <- intersect(unique(GntxCtx$ensembl_gene_id), txs)
      q4[2, q] <- getMutPerMbTxs(GntxCtx, txs.gc)
      
      colnames(q4)[q] <- paste0(colnames(q4)[q], length(unique(c(txs.cg, txs.gc))), " (", length(txs.cg), "+", length(txs.gc), ")")
   }
   
   q4s[[i]] <- q4
   barplot(q4, ylab="SNVs/Mb", main=getMain(rownames(q4)), beside=TRUE, width=.3, col=c("#E6E6E6", "#4D4D4D"))
}
dev.off()
save(tx.snv.s6, asyms, q4s, file=paste0(wd.asym.data, "sclc_asym_tx_pcg_snv_s6_asyms_q4s.RData"))

###
## Refining the plot
pdf(paste0(wd.asym.plots, "sclc_asym_tx_pcg_snv_s6_asyms_q4s.pdf"), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   q4 <- q4s[[i]]
   colnames(q4) <- c("", "            Tx", "", "")
   #for (c in 1:ncol(q4))
   #   colnames(q4)[c] <- gsub("n=", "", unlist(strsplit(colnames(q4)[c], " "))[1])
   
   barplot(q4, ylab="SNVs/Mb", main=getMain(rownames(asyms[[i]])), beside=TRUE, width=.3, col=c("#E6E6E6", "#4D4D4D"))   ##, ylim=c(0, 410))
   mtext(paste(rownames(q4), collapse=" vs "), cex=0.55, font=3, line=0.5)
}
dev.off()

# -----------------------------------------------------------------------------
# Step 4.2: Log2 transcription-associated SNV asymmetry (Figure 3)
# Last Modified: 26/01/18
# -----------------------------------------------------------------------------
pdf(paste0(wd.asym.plots, "sclc_asym_tx_pcg_snv_s6_asyms_q4s_log2_ylim1.6.pdf"), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   q4 <- q4s[[i]]
   colnames(q4) <- c("", "           Tx", "", "")
   
   barplot(log2(q4[1,]/q4[2,]), ylab="Asymmetry", ylim=c(-1.6, 1.6), main=getMain(rownames(asyms[[i]])), beside=TRUE, width=.3, col=c("#E6E6E6", "#E6E6E6", "#E6E6E6", "#E6E6E6"))
   mtext(paste0("log2(", paste(rownames(q4), collapse="/"), ")"), cex=0.55, font=3, line=0.5)
}
dev.off()










# =============================================================================
# Inner Class  : Collections of test/obsolete/deprecated methods
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified:
# =============================================================================
for (s in 1:length(samples)) {
   sample <- samples[s]
   
   pdf(paste0(wd.asym.plots, "sclc_snv_s6_", sample, ".pdf"))
   par(mfrow=c(3,1))
   #par(mar=c(3.5, 3.5, 1,1), mgp=c(2, 0.65, 0), cex=0.9)
   for (i in 1:length(idxs)) {
      idx <- idxs[i]
  
      hist(wt)
      hist(mpg)
      hist(disp)
      
      tx.snv.s6[[i]][[1]][[1]]   ## E.g. C>A Tx(+)
      tx.snv.s6[[i]][[1]][[2]]   ##      C>A Tx(-)
  
      tx.snv.s6[[i]][[2]][[1]]   ##      G>T Tx(+)
      tx.snv.s6[[i]][[2]][[2]]   ##      G>T Tx(-)
   }
   dev.off()
}










# > getMutPerMb(tx.snv.s6[[i]][[1]][[1]])   ## E.g. C>A Tx(+)
# [1] 162.091
# >    getMutPerMb(tx.snv.s6[[i]][[1]][[2]])   ##      C>A Tx(-)
# [1] 293.3805
# >    getMutPerMb(tx.snv.s6[[i]][[2]][[1]])   ##      G>T Tx(+)
# [1] 281.7484
# >    getMutPerMb(tx.snv.s6[[i]][[2]][[2]])   ##      G>T Tx(-)
# [1] 175.6149
# >    
#  >    getMutPerMb0(tx.snv.s6[[i]][[1]][[1]])   ## E.g. C>A Tx(+)
# [1] 176.2969
# >    getMutPerMb0(tx.snv.s6[[i]][[1]][[2]])   ##      C>A Tx(-)
# [1] 264.0553
# >    getMutPerMb0(tx.snv.s6[[i]][[2]][[1]])   ##      G>T Tx(+)
# [1] 262.0042
# >    getMutPerMb0(tx.snv.s6[[i]][[2]][[2]])   ##      G>T Tx(-)
# [1] 194.1061








 
   snv.s6[[1]][[1]][[1]] <- vcf.CA.p   ## C>A/G>T
   snv.s6[[1]][[1]][[2]] <- vcf.CA.m
   snv.s6[[1]][[2]][[1]] <- vcf.GT.p
   snv.s6[[1]][[2]][[2]] <- vcf.GT.m
   
   snv.s6[[2]][[1]][[1]] <- vcf.CG.p   ## C>G/G>C
   snv.s6[[2]][[1]][[2]] <- vcf.CG.m
   snv.s6[[2]][[2]][[1]] <- vcf.GC.p
   snv.s6[[2]][[2]][[2]] <- vcf.GC.m
   
   snv.s6[[3]][[1]][[1]] <- vcf.CT.p   ## C>T/G>A
   snv.s6[[3]][[1]][[2]] <- vcf.CT.m
   snv.s6[[3]][[2]][[1]] <- vcf.GA.p
   snv.s6[[3]][[2]][[2]] <- vcf.GA.m
   
   snv.s6[[4]][[1]][[1]] <- vcf.CG.p
   snv.s6[[4]][[1]][[2]] <- vcf.CG.m
   snv.s6[[4]][[2]][[1]] <- vcf.GC.p
   snv.s6[[4]][[2]][[2]] <- vcf.GC.m
   
   
   snv.s6[[1]][[1]] <- vcf.CA
   snv.s6[[1]][[2]] <- vcf.GT
   snv.s6[[2]][[1]] <- vcf.CG
   snv.s6[[2]][[2]] <- vcf.GC
   
   
   ## C>A/G>T (Smoking)
   vcf.CA <- subset(subset(vcf, REF == "C"), ALT == "A")
   vcf.CA.p <- subset(vcf.CA, strand == 1)
   vcf.CA.m <- subset(vcf.CA, strand == -1)
 
   vcf.GT <- subset(subset(vcf, REF == "G"), ALT == "T")
   vcf.GT.p <- subset(vcf.GT, strand == 1)
   vcf.GT.m <- subset(vcf.GT, strand == -1)



##for (s in 1:length(samples)) {
for (s in 1:10) {
   #sample <- SAMPLE
   sample <- samples[s]
   vcf <- readTable(paste0(wd.asym.data, sample, "_mut.ens.tx.snv.vcf.gz"), header=T, rownames=F, sep="")
  
   ## C>A/G>T (Smoking)
   vcf.CA <- subset(subset(vcf, REF == "C"), ALT == "A")
   vcf.CA.p <- subset(vcf.CA, strand == 1)
   vcf.CA.m <- subset(vcf.CA, strand == -1)
      
   vcf.GT <- subset(subset(vcf, REF == "G"), ALT == "T")
   vcf.GT.p <- subset(vcf.GT, strand == 1)
   vcf.GT.m <- subset(vcf.GT, strand == -1)

   ## C>G/G>C
   vcf.CG <- subset(subset(vcf, REF == "C"), ALT == "G")
   vcf.CG.p <- subset(vcf.CG, strand == 1)
   vcf.CG.m <- subset(vcf.CG, strand == -1)
   
   vcf.GC <- subset(subset(vcf, REF == "G"), ALT == "C")
   vcf.GC.p <- subset(vcf.GC, strand == 1)
   vcf.GC.m <- subset(vcf.GC, strand == -1)
   
   ## C>T/G>A (APOBEC/AID)
   vcf.CT <- subset(subset(vcf, REF == "C"), ALT == "T")
   vcf.CT.p <- subset(vcf.CT, strand == 1)
   vcf.CT.m <- subset(vcf.CT, strand == -1)
   
   vcf.GA <- subset(subset(vcf, REF == "G"), ALT == "A")
   vcf.GA.p <- subset(vcf.GA, strand == 1)
   vcf.GA.m <- subset(vcf.GA, strand == -1)

   ## T>A/A>T
   vcf.TA <- subset(subset(vcf, REF == "T"), ALT == "A")
   vcf.TA.p <- subset(vcf.TA, strand == 1)
   vcf.TA.m <- subset(vcf.TA, strand == -1)
   
   vcf.AT <- subset(subset(vcf, REF == "A"), ALT == "T")
   vcf.AT.p <- subset(vcf.AT, strand == 1)
   vcf.AT.m <- subset(vcf.AT, strand == -1)
   
   ## T>C/A>G (Sig. 5)
   vcf.TC <- subset(subset(vcf, REF == "T"), ALT == "C")
   vcf.TC.p <- subset(vcf.TC, strand == 1)
   vcf.TC.m <- subset(vcf.TC, strand == -1)
   
   vcf.AG <- subset(subset(vcf, REF == "A"), ALT == "G")
   vcf.AG.p <- subset(vcf.AG, strand == 1)
   vcf.AG.m <- subset(vcf.AG, strand == -1)
   
   ## T>G/A>C
   vcf.TG <- subset(subset(vcf, REF == "T"), ALT == "G")
   vcf.TG.p <- subset(vcf.TG, strand == 1)
   vcf.TG.m <- subset(vcf.TG, strand == -1)
   
   vcf.AC <- subset(subset(vcf, REF == "A"), ALT == "C")
   vcf.AC.p <- subset(vcf.AC, strand == 1)
   vcf.AC.m <- subset(vcf.AC, strand == -1)
}
