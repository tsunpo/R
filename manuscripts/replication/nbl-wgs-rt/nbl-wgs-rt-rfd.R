# =============================================================================
# Manuscript   : 
# Chapter II   : 
# Name         : manuscripts/replicaiton/nbl-wgs-rt-rfd.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 11/09/19; 12/11/18
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/re/home/tyang2/dev/R"                ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "ReplicationForkDirectionality.R", "ReplicationTiming.R", "DifferentialExpression.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.bed.gc.1kb.RData"))
load(file.path(wd.src.ref, "hg19.ensGene.bed.1kb.RData"))

# -----------------------------------------------------------------------------
# 
# Last Modified: 31/10/18
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
#wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE <- "NBL"
base <- tolower(BASE)
method <- "rpkm"

wd.ngs   <- file.path(wd, BASE, "ngs/WGS")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.rt    <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data/bstrps")
#wd.rt.plots <- file.path(wd.rt, "plots/bstrps")
wd.rt.plots <- file.path(wd.rt, "plots/nrfd")

# -----------------------------------------------------------------------------
# Bootstrap distribution
# Last Modified: 02/11/18
# -----------------------------------------------------------------------------
nrds.RT.BSTRPS <- getBootstrap(base, "SLOPE")
save(nrds.RT.BSTRPS, file=file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.d.rt.RT.SLOPE.RData")))
# > nrow(nrds.RT.BSTRPS)
# [1] 2652467

file.name <- file.path(wd.rt.plots, paste0("hist_", base, "_rpkm_SLOPE.pdf"))
main.text <- c(paste0(BASE, " bootstrap distribution"), paste0("Chr1-22 (1-kbs)"))
xlab.text <- "Number of right-leading resamplings"
plotBootstrapHist(nrds.RT.BSTRPS, file.name, main.text, xlab.text, 100, boundary.break)

# -----------------------------------------------------------------------------
# Create RT + RFD data
# Last Modified: 22/10/1
# -----------------------------------------------------------------------------
#load(file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.d.rt.RT.SLOPE.RData")))
load(file.path(wd.anlys, "replication", paste0(base, "-wgs-rt-m2"), "data", paste0(base, "_", method, ".gc.cn.d.rt.log2s_", "m2-m1", ".RData")))
# nrow(nrds)
# [1] 2659570

#install.packages("zoo", method="wget")
library("zoo")
kb <- 20
nrds.RT.NRFD <- getRTNRFD(nrds, nrds.RT.BSTRPS, bed.gc, kb)

save(nrds.RT.NRFD, file=file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.d.rt.log2s.nrfd.", kb, "kb_", "m2-m1", ".RData")))
writeTable(nrds.RT.NRFD, gzfile(file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.d.rt.log2s.nrfd.", kb, "kb_", "m2-m1", ".txt.gz"))), colnames=T, rownames=T, sep="\t")
nrds.RT.NRFD.nbl <- nrds.RT.NRFD
# > nrow(nrds.RT.NRFD.nbl)
# [1] 2652467

snr$S[3] <- sd(nrds.RT.NRFD.nbl$SPLINE)
snr$N[3] <- sd(nrds.RT.NRFD.nbl$RT - nrds.RT.NRFD.nbl$SPLINE)

# -----------------------------------------------------------------------------
# Plot bootstrap RFD data
# Last Modified: 04/11/18
# -----------------------------------------------------------------------------
boundary.upper <- 950   ## 500-520 breaks
boundary.lower <-  50   ## 480-500 breaks
boundary.break <-  45   ## 45 breaks each centering 500

###
##
kb <- 20
load(file=file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.d.rt.log2s.nrfd.", kb, "kb_", "m2-m1", ".RData")))

genes <- c("MRPL47", "ANAPC11", "C20orf24", "MED10", "RPL38", "TRIAP1", "ICT1", "SDHAF2")
genes <- c("TERT", "MDM2", "ALK", "ATRX")
for (g in 1:length(genes)) {
   chr <- subset(ensGene, external_gene_name == genes[g])$chromosome_name
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   start_position <- subset(ensGene, external_gene_name == genes[g])$start_position
   end_position   <- subset(ensGene, external_gene_name == genes[g])$end_position
   strand <- subset(ensGene, external_gene_name == genes[g])$strand
 
   size  <- 2000000
   start <- start_position + size
   end   <- start_position - size
   if (strand < 0) {
      start <- end_position + size
      end   <- end_position - size
   }
 
   file.name <- file.path(wd.rt.plots, paste0("NRFD_", BASE, "_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, "_", genes[g]))
   plotBootstrapRFD(file.name, BASE, chr, end, start, nrds.RT.NRFD, bed.gc.chr, boundary.upper, boundary.lower, "png", width=5, kb, F, genes[g])
}





## Chr17 (BIRC5)
c <- 17
chr <- chrs[c]
bed.gc.chr <- subset(bed.gc, CHR == chrs[c])

file.name <- file.path(wd.rt.plots, paste0("NRFD_", base, "_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, "_TTR_BIRC5"))
plotBootstrapRFD(file.name, BASE, chr,  74210267, 78210267, nrds.RT.NRFD, bed.gc.chr, boundary.upper, boundary.lower, "png", width=5, kb, gene="BIRC5")

## Chr17 (RDM1)
c <- 17
chr <- chrs[c]
bed.gc.chr <- subset(bed.gc, CHR == chrs[c])

file.name <- file.path(wd.rt.plots, paste0("NRFD_", base, "_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, "_TTR_RDM1"))
plotBootstrapRFD(file.name, BASE, chr,  32257777, 36257777, nrds.RT.NRFD, bed.gc.chr, boundary.upper, boundary.lower, "png", width=5, kb, gene="RDM1")

## Chr17 (TP53)
c <- 17
chr <- chrs[c]
bed.gc.chr <- subset(bed.gc, CHR == chrs[c])

file.name <- file.path(wd.rt.plots, paste0("NRFD_", base, "_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, "_TTR_TP53"))
plotBootstrapRFD(file.name, BASE, chr,  5590856, 9590856, nrds.RT.NRFD, bed.gc.chr, boundary.upper, boundary.lower, "png", width=5, kb, gene="TP53")

## Chr4
c <- 4
chr <- chrs[c]
bed.gc.chr <- subset(bed.gc, CHR == chr)

#file.name <- file.path(wd.rt.plots, paste0("NRFD_", base, "_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, "_TTR"))
#plotBootstrapRFD(file.name, BASE, chr, 142575001, 172575001, nrds.RT.NRFD, bed.gc.chr, boundary.upper, boundary.lower, "png", width=5, kb, withUnclassified=T)
file.name <- file.path(wd.rt.plots, paste0("NRFD_", base, "_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, "_TTR_IRF2"))
plotBootstrapRFD(file.name, BASE, chr, 183308867, 187395734, nrds.RT.NRFD, bed.gc.chr, boundary.upper, boundary.lower, "png", width=5, kb, gene="IRF2")

## Chr5
c <- 5
chr <- chrs[c]
bed.gc.chr <- subset(bed.gc, CHR == chr)

#file.name <- file.path(wd.rt.plots, paste0("NRFD_", base, "_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, "_TTR"))
#plotBootstrapRFD(file.name, BASE, chr, 142575001, 172575001, nrds.RT.NRFD, bed.gc.chr, boundary.upper, boundary.lower, "png", width=5, kb, withUnclassified=T)
file.name <- file.path(wd.rt.plots, paste0("NRFD_", base, "_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, "_TTR_TERT"))
plotBootstrapRFD(file.name, BASE, chr, 0, 3295184, nrds.RT.NRFD, bed.gc.chr, boundary.upper, boundary.lower, "png", width=5, kb, gene="TERT")

## Chr1
c <- 1
chr <- chrs[c]
bed.gc.chr <- subset(bed.gc, CHR == chr)

#file.name <- file.path(wd.rt.plots, paste0("NRFD_", base, "_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, "_TTR"))
#plotBootstrapRFD(file.name, BASE, chr, 142575001, 172575001, nrds.RT.NRFD, bed.gc.chr, boundary.upper, boundary.lower, "png", width=5, kb, withUnclassified=T)
file.name <- file.path(wd.rt.plots, paste0("NRFD_", base, "_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, "_TTR_TOR1AIP1"))
plotBootstrapRFD(file.name, BASE, chr, 174851177, 184851177, nrds.RT.NRFD, bed.gc.chr, boundary.upper, boundary.lower, "png", width=5, kb)

## Chr13 (BRCA2)
c <- 13
chr <- chrs[c]
bed.gc.chr <- subset(bed.gc, CHR == chrs[c])

file.name <- file.path(wd.rt.plots, paste0("NRFD_", base, "_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, "_TTR_BRCA2"))
plotBootstrapRFD(file.name, BASE, chr,  30889611, 34889611, nrds.RT.NRFD, bed.gc.chr, boundary.upper, boundary.lower, "png", width=5, kb)

## Chr15 (PIF1)
c <- 15
chr <- chrs[c]
bed.gc.chr <- subset(bed.gc, CHR == chr)

file.name <- file.path(wd.rt.plots, paste0("NRFD_", base, "_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, "_TTR_PIF1"))
plotBootstrapRFD(file.name, BASE, chr, 63117867, 67117867, nrds.RT.NRFD, bed.gc.chr, boundary.upper, boundary.lower, "png", width=5, kb)

# -----------------------------------------------------------------------------
# 
# Last Modified: 12/08/20
# -----------------------------------------------------------------------------
nrds.RT.NRFD.nbl.ttr <- getBootstrapTTR(nrds.RT.NRFD.nbl, 0.9)

nrds.RT.NRFD.nbl.ctr <- getBootstrapCTR(nrds.RT.NRFD.nbl, 0.9)
nrds.RT.NRFD.nbl.ctr.iz <- subset(nrds.RT.NRFD.nbl.ctr, NRFD > 0)
nrds.RT.NRFD.nbl.ctr.tz <- subset(nrds.RT.NRFD.nbl.ctr, NRFD < 0)

nrds.RT.NRFD.nbl.ctr.iz.e <- subset(nrds.RT.NRFD.nbl.ctr.iz, RT > 0)
nrds.RT.NRFD.nbl.ctr.iz.l <- subset(nrds.RT.NRFD.nbl.ctr.iz, RT < 0)
nrds.RT.NRFD.nbl.ctr.tz.e <- subset(nrds.RT.NRFD.nbl.ctr.tz, RT > 0)
nrds.RT.NRFD.nbl.ctr.tz.l <- subset(nrds.RT.NRFD.nbl.ctr.tz, RT < 0)

nrow(nrds.RT.NRFD.nbl.ctr.iz.e)/nrow(nrds.RT.NRFD.nbl)
# [1] 0.07562017
nrow(nrds.RT.NRFD.nbl.ctr.iz.l)/nrow(nrds.RT.NRFD.nbl)
# [1] 0.04072021
nrow(nrds.RT.NRFD.nbl.ctr.tz.e)/nrow(nrds.RT.NRFD.nbl)
# [1] 0.0571664
nrow(nrds.RT.NRFD.nbl.ctr.tz.l)/nrow(nrds.RT.NRFD.nbl)
# [1] 0.05736848









## Chr2
c <- 2
chr <- chrs[c]
bed.gc.chr <- subset(bed.gc, CHR == chr)

file.name <- file.path(wd.rt.plots, paste0("NRFD_", base, "_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, "_TTR"))
plotBootstrapRFD(file.name, BASE, chr, 110000000, 130000000, nrds.RT.NRFD, bed.gc.chr, boundary.upper, boundary.lower, "png", width=5, kb)

## Chr12
c <- 12
chr <- chrs[c]
bed.gc.chr <- subset(bed.gc, CHR == chrs[c])

file.name <- file.path(wd.rt.plots, paste0("NRFD_", base, "_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, "_TTR"))
plotBootstrapRFD(file.name, BASE, chr,  97500000, 105000000, nrds.RT.NRFD, bed.gc.chr, boundary.upper, boundary.lower, "png", width=5, kb)

## Chr1
c <- 1
chr <- chrs[c]
bed.gc.chr <- subset(bed.gc, CHR == chr)

file.name <- file.path(wd.rt.plots, paste0("NRFD_", base, "_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, "_TTR"))
plotBootstrapRFD(file.name, BASE, chr, 142575001, 172575001, nrds.RT.NRFD, bed.gc.chr, boundary.upper, boundary.lower, "png", width=5, kb, withUnclassified=T)










load(file.path(wd.anlys, "replication", paste0(base, "-wgs-rt-m2"), "data", paste0(base, "_", method, ".gc.cn.d.rt.log2s_", "m2-m1", ".RData")))
nrds.RT <- getRT(nrds, bed.gc)
# > nrow(nrds.RT)
# [1] 2652467

nrds.RT.RFD <- getRTRFD(nrds.RT, nrds.RT.BSTRPS)
save(nrds.RT.RFD, file=file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.d.rt.log2s.rfd_", "m2-m1", ".RData")))
writeTable(nrds.RT.RFD, gzfile(file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.d.rt.log2s.rfd_", "m2-m1", ".txt.gz"))), colnames=T, rownames=T, sep="\t")
nrds.RT.RFD.nbl <- nrds.RT.RFD
# > nrow(nrds.RT.RFD.nbl)
# [1] 2652467











# -----------------------------------------------------------------------------
# Wilcoxon rank sum test (non-parametric; T vs N)
# Last Modified: 22/10/19
# -----------------------------------------------------------------------------
nrds.RT.RFD.nbl.c.e <- nrds.RT.RFD.nbl[overlaps.c.e,]
# > nrow(nrds.RT.RFD.nbl.c.e)
# [1] 40451

nrds.RT.RFD.nbl.c.e.q4 <- subset(nrds.RT.RFD.nbl.c.e, RT >= quantile(nrds.RT.RFD.nbl.c.e$RT)[4])











# -----------------------------------------------------------------------------
# Plot bootstrap RFD data
# Last Modified: 04/11/18
# -----------------------------------------------------------------------------
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]
 
   ## RFD   
   load(file=file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.d.rt.RT.SLOPE_", chr, ".RData")))
   nrds.RT.BSTRPS.chr$RFD <- getRFD(nrds.RT.BSTRPS.chr)

   ## Chr12
   file.name <- file.path(wd.rt.plots, paste0("RFD_", base, "_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ""))
   plotBootstrapRFD(file.name, BASE, chr, NA, NA, nrds.chr, bed.gc.chr, nrds.RT.BSTRPS.chr, boundary.upper, boundary.lower, "png", width=10)
   plotBootstrapRFD(file.name, BASE, chr,  97500000, 105000000, nrds.chr, bed.gc.chr, nrds.RT.BSTRPS.chr, boundary.upper, boundary.lower, "png", width=5)
 
   ## Chr2
   file.name <- file.path(wd.rt.plots, paste0("RFD_", base, "_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ""))
   plotBootstrapRFD(file.name, BASE, chr, NA, NA, nrds.chr, bed.gc.chr, nrds.RT.BSTRPS.chr, boundary.upper, boundary.lower, "png", width=10)
   plotBootstrapRFD(file.name, BASE, chr,  37000000,  40000000, nrds.chr, bed.gc.chr, nrds.RT.BSTRPS.chr, boundary.upper, boundary.lower, "png", width=5)
   plotBootstrapRFD(file.name, BASE, chr, 215000000, 220000000, nrds.chr, bed.gc.chr, nrds.RT.BSTRPS.chr, boundary.upper, boundary.lower, "png", width=5)
   plotBootstrapRFD(file.name, BASE, chr, 110000000, 130000000, nrds.chr, bed.gc.chr, nrds.RT.BSTRPS.chr, boundary.upper, boundary.lower, "png", width=10)
   
   ## Chr13
   file.name <- file.path(wd.rt.plots, paste0("RFD_", base, "_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ""))
   plotBootstrapRFD(file.name, BASE, chr, NA, NA, nrds.chr, bed.gc.chr, nrds.RT.BSTRPS.chr, boundary.upper, boundary.lower, "png", width=10)
}









# -----------------------------------------------------------------------------
# Distributions in bootstrapped data
# Last Modified: 02/11/18
# -----------------------------------------------------------------------------
nrds.RT.BSTRPS <- NULL
for (c in 1:22) {
   chr <- chrs[c]

   load(file=file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.d.rt.RT.SPLINE_", chr, ".RData")))
   if (is.null(nrds.RT.BSTRPS))
      nrds.RT.BSTRPS <- nrds.RT.BSTRPS.chr
   else
      nrds.RT.BSTRPS <- rbind(nrds.RT.BSTRPS, nrds.RT.BSTRPS.chr)
}
save(nrds.RT.BSTRPS, file=file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.d.rt.RT.SPLINE.RData")))
file.name <- file.path(wd.rt.plots, paste0("hist_", base, "_rpkm_SPLINE_BSTRPS.pdf"))
main.text <- paste0("Chr1-22 (", BASE, ")")
xlab.text <- c("Number of early replication counts", "(out of 1,000 bootstrappings)")
plotBootstrapHist(nrds.RT.BSTRPS, file.name, main.text, xlab.text, 200, boundary.break)
# > nrow(nrds.RT.BSTRPS)
# [1] 2652467












# -----------------------------------------------------------------------------
# Plot RO and RT from bootstrapped data
# Last Modified: 04/11/18
# -----------------------------------------------------------------------------
for (c in 1:22) {
   #c <- 8
   chr <- chrs[c]

   ## Replication fork directionality (RFD)  
   load(file=file.path(wd.rt.data, paste0("bed.gc.rt_", base, "_bstrps", bstrps, "_", chr, ".RData")))
   bed.gc.rt.chr$RFD <- mapply(x = 1:nrow(bed.gc.rt.chr), function(x) as.numeric(getRFD(bed.gc.rt.chr[x, ])))   ## ADD 05/12/18
   right.idx  <- which(bed.gc.rt.chr$RIGHT_LEADING > boundary.upper)
   left.idx   <- which(bed.gc.rt.chr$RIGHT_LEADING < boundary.lower)
   boundary.idx <- setdiff(c(1:nrow(bed.gc.rt.chr)), c(right.idx, left.idx))
   bed.gc.chr <- bed.gc[rownames(bed.gc.rt.chr),]
   
   bed.gc.chr.sub <- subset(bed.gc.chr[boundary.idx,], START > 85000000)
   bed.gc.chr.sub <- subset(bed.gc.chr.sub, START < 87500000)
   rpkms.chr.rt[rownames(bed.gc.rt.chr[rownames(bed.gc.chr.sub),]),]
   
   start <- 128747680 - 500000
   end   <- 128753674 + 500000
   file.name <- file.path(wd.rt.plots, paste0(base, "_RFD_bstrps1000_", chr))
   plotBootstrapsRFD(file.name, BASE, chr, start, end, bed.gc.chr, bed.gc.rt.chr, right.idx, left.idx, boundary.idx, "png")      ## see ReplicationTiming.R
   
   ## Replication timing
   rt.chr <- readTable(file.path(wd.rt, "data", paste0(base, "_rpkm.corr.gc.d.rt_", chr, "_", BASE, "-", BASE, "_n101-92.txt.gz")), header=T, rownames=T, sep="\t") 
   rt.chr <- rt.chr[rownames(bed.gc.rt.chr),]
   
   file.name  <- file.path(wd.rt.plots, paste0(base, "_RT_bstrps1000_", chr))
   plotBootstrapsRT(file.name, BASE, chr, start, end, rt.chr, bed.gc.chr, right.idx, left.idx, boundary.idx, ymax=0.15, "png")   ## see ReplicationTiming.R
}

# -----------------------------------------------------------------------------
# Transcription vs. replication time
# Last Modified: 08/11/18
# -----------------------------------------------------------------------------
load(file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
tpm.gene      <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=T, proteinCodingOnly=F)
tpm.gene.log2 <- getLog2andMedian(tpm.gene, pseudocount=0.01)
# > nrow(tpm.gene.log2)
# [1] 18440

load(file=file.path(wd.rt.data, paste0("ensGene.rt_", base, "_bstrps", bstrps, ".RData")))   ## Load objects ensGene.rt.start and ensGene.rt.end
ensGene.rt.tx <- getEnsGeneTxRFD(ensGene, ensGene.rt.start, ensGene.rt.end, tpm.gene.log2)

##
file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_RFD_TSS+TTS_pc1e2"))
plotEnsGeneTxRFD(file.name, BASE, ensGene.rt.tx, boundary.upper, "png")

#leading.ratio <- log2(boundary.upper/500)
leading.ratio <- (boundary.upper - (1000 - boundary.upper))/1000
sclc.tx.right.0 <- rownames(subset(ensGene.rt.tx, RFD_TSS >  leading.ratio))
sclc.tx.left.0  <- rownames(subset(ensGene.rt.tx, RFD_TSS < -leading.ratio))
sclc.tx.5050    <- setdiff(rownames(ensGene.rt.tx), c(sclc.tx.right.0, sclc.tx.left.0))
sclc.tx.inconsist.0 <- rownames(subset(ensGene.rt.tx, CONSIST < 0))   ## This may include "50/50" if cutoff for leading-count ratio change
# length(intersect(sclc.tx.5050, sclc.tx.inconsist.0))
# [1] 42

## ADD 20/11/18: 50/50
sclc.tx.5050.right <- rownames(subset(ensGene.rt.tx[sclc.tx.5050,], RFD_TTS >  leading.ratio))
sclc.tx.5050.left  <- rownames(subset(ensGene.rt.tx[sclc.tx.5050,], RFD_TTS < -leading.ratio))
sclc.tx.5050.5050  <- setdiff(sclc.tx.5050, c(sclc.tx.5050.right, sclc.tx.5050.left))
sclc.tx.5050.5050.right <- rownames(subset(ensGene.rt.tx[sclc.tx.5050.5050,], RFD_TTS > 0))
sclc.tx.5050.5050.left  <- rownames(subset(ensGene.rt.tx[sclc.tx.5050.5050,], RFD_TTS < 0))

## ADD 20/11/18: Remove inconsistent genes
sclc.tx.consist.right <- setdiff(sclc.tx.right.0, sclc.tx.inconsist.0)
sclc.tx.consist.left  <- setdiff(sclc.tx.left.0,  sclc.tx.inconsist.0)
sclc.tx.consist <- c(sclc.tx.consist.right, sclc.tx.consist.left)

sclc.tx.inconsist <- setdiff(sclc.tx.inconsist.0, sclc.tx.5050)
sclc.tx.inconsist.right <- rownames(subset(ensGene.rt.tx[sclc.tx.inconsist,], RFD_TSS >  leading.ratio))
sclc.tx.inconsist.left  <- rownames(subset(ensGene.rt.tx[sclc.tx.inconsist,], RFD_TSS < -leading.ratio))
testWbyEnsGeneRTTx(ensGene.rt.tx, sclc.tx.inconsist, sclc.tx.consist)
# [1] 0.01322831
## ADD 20/11/18: How many boundarys are inconsistent (which will be visulised in the boxplot later)
# > 8494+8432+1056+96   ## Right + Left + Inconsistent + 50/50
# [1] 18078

###
##
file.name  <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_RFD_TSS+TTS_consist+5050.pdf"))
main.text <- "Replication fork directionality (RFD) in SCLC"
xlab.text <- c("RFD = (R-L)/(R+L)", "TSS")
ylab.text <- "TTS"
sclc.tx <- c(sclc.tx.consist, sclc.tx.inconsist)
cols <- rep("mediumpurple1", length(sclc.tx))

pdf(file.name, height=6, width=6)#, units="in", res=300)
plot(ensGene.rt.tx[sclc.tx,]$RFD_TTS ~ ensGene.rt.tx[sclc.tx,]$RFD_TSS, xlab=xlab.text, ylab=ylab.text, main=main.text, col=cols)

points(ensGene.rt.tx[sclc.tx.consist.right,]$RFD_TTS ~ ensGene.rt.tx[sclc.tx.consist.right,]$RFD_TSS, col="sandybrown")
points(ensGene.rt.tx[sclc.tx.consist.left,]$RFD_TTS  ~ ensGene.rt.tx[sclc.tx.consist.left,]$RFD_TSS,  col="steelblue1")
abline(v=0, lty=5, lwd=0.85, col="black")
text( 0.5,  0.5, paste0("94.04% (", separator(length(sclc.tx.consist.right)), ")"), cex=1, col="black")
text(-0.5, -0.5, paste0("94.21% (", separator(length(sclc.tx.consist.left)), ")"), cex=1, col="black") 

#legend("bottom", c("Left-leading in SCLC", paste0("L/R (n=", length(sclc.tx.boundary.o.c), ") in SCLC"), "Right-leading in SCLC"), col=c("steelblue1", "red", "sandybrown"), pch=1, cex=0.8, horiz=T)
mtext(paste0("Consistent and inconsistent genes (n=", separator(length(sclc.tx)), ")"), cex=1, line=0.4)
dev.off()

# -----------------------------------------------------------------------------
# LM for sclc.tx.consist.right and sclc.tx.consist.left
# Last Modified: 07/12/18
# -----------------------------------------------------------------------------
plotRFD <- function(file.name, main.text, ensGene.rt.tx, tx.list1, col1, tx.list2=NULL, col2=NULL, ext) {
   xlab.text <- c("RFD = (R-L)/(R+L)", "TSS")
   ylab.text <- "TTS"
   tx.list1 <- intersect(tx.list1, rownames(ensGene.rt.tx))
   ensGene.rt.tx.plot <- ensGene.rt.tx[tx.list1,]
   if (!is.null(tx.list2) && length(tx.list2) != 0) {
      tx.list2 <- intersect(tx.list2, rownames(ensGene.rt.tx))
      ensGene.rt.tx.plot <- ensGene.rt.tx[c(tx.list1, tx.list2),]
   }
   
   if (max(ensGene.rt.tx.plot[tx.list1,]$RFD_TSS) > 0) {
      xlim <- c(0, 1)
      ylim <- c(0, 1)
   } else {
      xlim <- c(-1, 0)
      ylim <- c(-1, 0) 
   }
   xmin <- xlim[1]
   xmax <- xlim[2]
   ymin <- ylim[1]
   ymax <- ylim[2]
   
   if (ext == "pdf") {
      pdf(paste0(file.name, ".pdf"), height=6, width=6)
   } else if (ext == "png")
      png(paste0(file.name, ".png"), height=6, width=6, units="in", res=300)   ## ADD 16/05/17: res=300
   plot(RFD_TTS ~ RFD_TSS, data=ensGene.rt.tx.plot, ylim=ylim, xlim=xlim, ylab=ylab.text, xlab=xlab.text, main=main.text[1], col=col1)
   if (!is.null(tx.list2))
      points(RFD_TTS ~ RFD_TSS, data=ensGene.rt.tx.plot[tx.list2,], col=col2)
   
   ## Linear regression
   ## https://www.statlect.com/fundamentals-of-statistics/R-squared-of-a-linear-regression
   ensGene.rt.tx.plot$SLOPE <- ensGene.rt.tx.plot$AS / ensGene.rt.tx.plot$LENGTH * 2   ## BUG *2 23/01/19
   lm.fit <- lm(ensGene.rt.tx.plot$RFD_TTS ~ ensGene.rt.tx.plot$RFD_TSS + ensGene.rt.tx.plot$LENGTH + ensGene.rt.tx.plot$TPM + ensGene.rt.tx.plot$SLOPE)
   intercept <- coef(lm.fit)[1]
   r2 <- summary(lm.fit)$r.squared
   
   if (max(ensGene.rt.tx.plot[tx.list1,]$RFD_TTS) > 0) {   ## 0.98, 0.92, 0.80, 0.74
      text(0.18, 0.98, paste0("LM's R^2 = ", round0(r2*100, digits=2), "%"), cex=1.2)
      if (intercept > 0)
         text(0.142, 0.92, paste0("Intercept = ", round0(intercept, digits=2)), cex=1.2)
      else
         text(0.155, 0.92, paste0("Intercept = ", round0(intercept, digits=2)), cex=1.2)         
   } else {                                     ## -0.74, -0.80, -0.92, -0.98
      text(-0.18, -0.92, paste0("LM's R^2 = ", round0(r2*100, digits=2), "%"), cex=1.2)
      if (intercept > 0)
         text(-0.218, -0.98, paste0("Intercept = ", round0(intercept, digits=2)), cex=1.2)
      else
         text(-0.205, -0.98, paste0("Intercept = ", round0(intercept, digits=2)), cex=1.2)
   }
   
   abline(lm.fit)
   mtext(main.text[2], cex=1.2, line=0.3)
   dev.off()
   
   return(as.numeric(intercept))
}

## Consistent
file.name  <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_RFD_TSS+TTS_consist_right"))
main.text <- c("Right-leading consistent genes in SCLC and LCL", paste0("(n=", separator(length(sclc.tx.consist.right)), ")"))
plotRFD(file.name, main.text, ensGene.rt.tx, sclc.tx.consist.right, "sandybrown", ext="pdf")

file.name  <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_RFD_TSS+TTS_consist_left"))
main.text <- c("Left-leading consistent genes in SCLC", paste0("(n=", separator(length(sclc.tx.consist.left)), ")"))
plotRFD(file.name, main.text, ensGene.rt.tx, sclc.tx.consist.left, "steelblue1", ext="pdf")

## Right (CD and HO)
sclc.tx.consist.right.cd <- rownames(subset(ensGene.rt.tx[sclc.tx.consist.right,], CD > 0))
sclc.tx.consist.right.ho <- rownames(subset(ensGene.rt.tx[sclc.tx.consist.right,], CD < 0))

file.name  <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_RFD_TSS+TTS_consist_right_CD"))
main.text <- c("Right-leading consistent genes in SCLC", paste0("(Co-directional; n=", separator(length(sclc.tx.consist.right.cd)), ")"))
plotRFD(file.name, main.text, ensGene.rt.tx, sclc.tx.consist.right.cd, "sandybrown", ext="pdf")

file.name  <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_RFD_TSS+TTS_consist_right_HO"))
main.text <- c("Right-leading consistent genes in SCLC", paste0("(Head-on; n=", separator(length(sclc.tx.consist.right.ho)), ")"))
plotRFD(file.name, main.text, ensGene.rt.tx, sclc.tx.consist.right.ho, "sandybrown", ext="pdf")

## Left (CD and HO)
sclc.tx.consist.left.cd <- rownames(subset(ensGene.rt.tx[sclc.tx.consist.left,], CD > 0))
sclc.tx.consist.left.ho <- rownames(subset(ensGene.rt.tx[sclc.tx.consist.left,], CD < 0))

file.name  <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_RFD_TSS+TTS_consist_left_CD"))
main.text <- c("Left-leading consistent genes in SCLC", paste0("(Co-directional; n=", separator(length(sclc.tx.consist.left.cd)), ")"))
plotRFD(file.name, main.text, ensGene.rt.tx, sclc.tx.consist.left.cd, "steelblue1", ext="pdf")

file.name  <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_RFD_TSS+TTS_consist_left_HO"))
main.text <- c("Left-leading consistent genes in SCLC", paste0("(Head-on; n=", separator(length(sclc.tx.consist.left.ho)), ")"))
plotRFD(file.name, main.text, ensGene.rt.tx, sclc.tx.consist.left.ho, "steelblue1", ext="pdf")

###
## Consistent genes in (Q1 to Q4); use ALL expressed genes
load(file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
#tpm.gene      <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=F, proteinCodingOnly=F)
tpm.gene.log2 <- getLog2andMedian(tpm.gene, pseudocount=0.01)

tx.q4 <- getTxQ4(tpm.gene.log2, NA)
for (q in 1:4)
   print(length(tx.q4[[q]]))
# [1] 4783
# [1] 4783
# [1] 4782
# [1] 4783
# for (q in 1:4)
#    print(length(intersect(tx.q4.fix[[q]], tx.q4.fix.all[[q]])))
# [1] 4561
# [1] 4554
# [1] 4580
# [1] 4610

tx.length.q4 <- getLengthQ4(ensGene[rownames(tpm.gene.log2),], NA)
for (q in 1:4)
   print(length(tx.length.q4[[q]]))

genes.rb1.up   <- rownames(subset(de.tpm.gene, LOG2_FC > 0))   ## ADD 13/01/19
genes.rb1.down <- rownames(subset(de.tpm.gene, LOG2_FC < 0))   ## ADD 13/01/19

## AS: (TSS - TTS)/2
textP <- function(p, at, sth) {
   text <- paste0("P-value=", sprintf("%4.2e", p))
   text <- paste0(text, " at ", at, "=", sth)
   if (p <= 1e-4)
      text <- paste0("***", text, "    ")
   else if (p <= 1e-3)
      text <- paste0("**", text, "   ")
   else if (p <= 1e-2)
      text <- paste0("*", text, "  ")
   
   return(text)
}

#rho <- 0.62333
#textP(0.001, "rho", round0(rho, digits=2))

plotLengthAS <- function(shifts, lengths, file.name, main.text, xlab.text) {
   ylab.text <- "RFD shift (initiation efficiency)"
   lm.fit <- lm(shifts ~ lengths)
   r2 <- summary(lm.fit)$r.squared
   
   pdf(paste0(file.name, ".pdf"), height=6, width=6)
   plot(shifts ~ lengths, ylab=ylab.text, xlab=xlab.text, main=main.text[1])
   abline(lm.fit)
   mtext(paste0("R^2=", paste0(round0(r2*100, digits=2), "%")), cex=1.2, line=0.3)
   dev.off()
}

tpm.gene$MEDIAN <- mapply(x = 1:nrow(tpm.gene), function(x) median(as.numeric(tpm.gene[x,])))
ensGene.rt.tx$TPM <- tpm.gene[rownames(ensGene.rt.tx),]$MEDIAN
ensGene.rt.tx$RFD_START <- mapply(x = 1:nrow(ensGene.rt.tx), function(x) as.numeric(getRFD(ensGene.rt.tx[x,  7:10])))   ## START = start_position
ensGene.rt.tx$RFD_END   <- mapply(x = 1:nrow(ensGene.rt.tx), function(x) as.numeric(getRFD(ensGene.rt.tx[x, 11:14])))   ## END   = end_position
ensGene.rt.tx$AS <- (ensGene.rt.tx$RFD_END - ensGene.rt.tx$RFD_START) / 2
#ensGene.rt.tx$AS[which(ensGene.rt.tx$strand == -1)] <- ensGene.rt.tx$AS[which(ensGene.rt.tx$strand == -1)] * -1
ensGene.rt.tx.con <- subset(ensGene.rt.tx, CONSIST == 1)

ensGene.rt.tx.input <- subset(ensGene.rt.tx.con, AS < 0)
ensGene.rt.tx.input <- subset(ensGene.rt.tx.input, CD > 0)
file.name <- file.path(wd.rt.plots, paste0("plot_", base, "_ensGene.rt.tx_bstrps1000_RFD_shift-length_DS_CD.pdf"))
main.text <- paste0("Decending genes (CD)")
xlab.text <- "Gene length [log10]"
plotLengthAS(ensGene.rt.tx.input$AS, log10(ensGene.rt.tx.input$LENGTH), file.name, main.text, xlab.text)

#left.idx <- which(ensGene.rt.tx.con$RT < 0)
#ensGene.rt.tx.con$SLOPE[left.idx] <- (ensGene.rt.tx.con$RFD_TSS[left.idx] - ensGene.rt.tx.con$RFD_TTS[left.idx]) / (ensGene.rt.tx.con$LENGTH[left.idx] / 1000)
#ensGene.rt.tx.con.no0 <- subset(ensGene.rt.tx.con, SLOPE != 0)

##
#ensGene.rt.tx.con.as <- subset(ensGene.rt.tx.con, AS > 0)
intercepts <- list(list(list(list(), list(), list(), list()), list(list(), list(), list(), list())), list(list(list(), list(), list(), list()), list(list(), list(), list(), list())))

tx.list1 <- sclc.tx.consist.right.ho
tx.list2 <- sclc.tx.5050.5050.right.ho
for (q in 1:4) {
   tx.list1.input <- intersect(tx.list1, tx.q4[[q]])
   tx.list2.input <- intersect(tx.list2, tx.q4[[q]])
   #tx.list1.input <- intersect(tx.list1.input, genes.rb1.down)   ## ADD 13/01/19
   #tx.list2.input <- intersect(tx.list2.input, genes.rb1.down)   ## ADD 13/01/19
   #tx.list1.input <- intersect(tx.list1.input, rownames(ensGene.rt.tx.con.as))
   
   file.name  <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_RFD_TTS~TSS+LENGTH+TPM+SLOPE_consist+5050_right_HO_Q", q))
   main.text <- c(paste0("Q", q, " right-moving consistent genes in SCLC"), paste0("(Head-on; n=", separator(length(tx.list1.input)), "+", length(tx.list2.input), ")"))
   intercepts[[1]][[2]][[q]] <- plotRFD(file.name, main.text, ensGene.rt.tx.con, tx.list1.input, "sandybrown", tx.list2=tx.list2.input, col2="red", ext="pdf")
}

tx.list1 <- sclc.tx.consist.left.ho
tx.list2 <- sclc.tx.5050.5050.left.ho
for (q in 1:4) {
   tx.list1.input <- intersect(tx.list1, tx.q4[[q]])
   tx.list2.input <- intersect(tx.list2, tx.q4[[q]])
   #tx.list1.input <- intersect(tx.list1.input, genes.rb1.down)   ## ADD 13/01/19
   #tx.list2.input <- intersect(tx.list2.input, genes.rb1.down)   ## ADD 13/01/19
   
   file.name  <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_RFD_TTS~TSS+LENGTH+TPM+SLOPE_consist+5050_left_HO_Q", q))
   main.text <- c(paste0("Q", q, " left-moving consistent genes in SCLC"), paste0("(Head-on; n=", separator(length(tx.list1.input)), "+", length(tx.list2.input), ")"))
   intercepts[[2]][[2]][[q]] <- plotRFD(file.name, main.text, ensGene.rt.tx.con, tx.list1.input, "steelblue1", tx.list2=tx.list2.input, col2="red", ext="pdf")
}

###
## Plot RFD_END constant
# http://blog.minitab.com/blog/adventures-in-statistics-2/regression-analysis-how-to-interpret-the-constant-y-intercept
#TTS.TSS.intercepts <- intercepts

qs <- c(1, 2, 3, 4)
is <- c()
for (rt in 1:2)
   for (st in 1:2)
      for (q in 1:4)
         is <- c(is, intercepts[[rt]][[st]][[q]])

file.name  <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_RFD_TTS~TSS+LENGTH+TPM+SLOPE_consist+5050_Constants3"))
main.text <- c("lm(TTS~TSS+LENGTH+TPM+SLOPE+Constant)", "among consistent genes in SCLC")
xlab.text <- "Expression"
ylab.text <- "TTS constant (when TSS=0)"
xlim <- c(1, 4)
max <- max(abs(min(is)), max(is))
max <- 0.2367002   ## L+TPM
#max <- 0.3701957   ## AS
#max <- 0.4
ylim <- c(-max, max)
pdf(paste0(file.name, ".pdf"), height=6, width=6)
plot(as.numeric(intercepts[[1]][[1]]) ~ qs, ylim=ylim, xlim=xlim, ylab=ylab.text, xlab=xlab.text, main=main.text[1], type="l", col="sandybrown", xaxt="n", lwd=2)
points(as.numeric(intercepts[[1]][[1]]) ~ qs, col="blue", pch=19)
lines(as.numeric(intercepts[[1]][[2]]) ~ qs, col="sandybrown", lwd=2)
points(as.numeric(intercepts[[1]][[2]]) ~ qs, col="red", pch=19)

lines(as.numeric(intercepts[[2]][[1]]) ~ qs, col="steelblue1", lwd=2)
points(as.numeric(intercepts[[2]][[1]]) ~ qs, col="blue", pch=19)
lines(as.numeric(intercepts[[2]][[2]]) ~ qs, col="steelblue1", lwd=2)
points(as.numeric(intercepts[[2]][[2]]) ~ qs, col="red", pch=19)

abline(h=0, lty=5, lwd=0.85)
legend("topleft", legend=c("Co-directional", "Head-on"), col=c("blue", "red"), pch=19)
legend("bottomleft", legend=c("Right-moving", "Left-moving"), col=c("sandybrown", "steelblue1"), lty=1, lwd=2, bty="n")
axis(side=1, at=seq(1, 4, by=1), labels=paste0("Q", qs))
mtext(main.text[2], cex=1.2, line=0.3)
dev.off()

###
## 15/01/19
tx.list1 <- sclc.tx.consist.right.cd
tx.list2 <- sclc.tx.5050.5050.right.cd
tx.list3 <- sclc.tx.consist.right.ho
tx.list4 <- sclc.tx.5050.5050.right.ho
for (q in 1:4) {
   tx.list1.input <- c(tx.list1 , tx.list2)
   tx.list2.input <- c(tx.list3 , tx.list4)
   
   tx.list1.input <- intersect(tx.list1.input, tx.q4.fix.all[[q]])
   tx.list2.input <- intersect(tx.list2.input, tx.q4.fix.all[[q]])
   
   print(paste0("Q", q, ", Length, CD vs. HO: ", testW(ensGene.rt.tx[tx.list1.input,]$LENGTH, ensGene.rt.tx[tx.list2.input,]$LENGTH)))
}

tx.list <- sclc.tx.5050.5050.right.ho
for (q in 1:4) {
   genes <- intersect(tx.list, tx.q4.fix.all[[q]])
 
   for (g in 1:length(genes))
     print(paste0("Q", q, ": ", ensGene[genes[g],]$external_gene_name))
}

save(ensGene.rt.tx, sclc.tx.consist.right.cd, sclc.tx.consist.right.ho, sclc.tx.5050.5050.right.cd, sclc.tx.5050.5050.right.ho, sclc.tx.consist.left.cd, sclc.tx.consist.left.ho, sclc.tx.5050.5050.left.cd, sclc.tx.5050.5050.left.ho, file=file.path(wd.rt.plot, "test.RData"))


###
## B2M; 14/01/19
genes <- c("SPG11", "B2M", "PATL2", "TRIM69", "BAX", "DHDH", "NUCB1")
for (g in 1:length(genes))
   plotTxDensityHistogram(genes[g], BASE, wd.rt.plots, tpm.gene.log2, ensGene.rt.tx, tx.q4.fix.all, 0.01)  ## See DifferentialExpression.R

## Plot only samples with CN=2
plotTxDensityHistogram("SPG11", BASE, wd.rt.plots, tpm.gene.log2[,c(data[cn2.idx,]$Sample, "MEDIAN")], ensGene.rt.tx, tx.q4, 0.01)  ## See DifferentialExpression.R

###
## CNs in replication-stress region near B2M
samples.rna <- colnames(tpm.gene.log2)[-ncol(tpm.gene.log2)]
samples.wgs <- readTable(file.path(wd.ngs, "sclc_wgs_n101.list"), header=F, rownames=F, sep="")
samples <- intersect(samples.rna, samples.wgs)
# length(samples)
# [1] 70

getCNFromSegment <- function(ensGene.gene, segs) {
   segs.chr <- segs[which(segs$CHR == ensGene.gene$chromosome_name),]
   segs.chr.end <- segs.chr[which(segs.chr$END >= ensGene.gene$start_position),]
   segs.chr.end.start <- segs.chr.end[which(segs.chr.end$START <= ensGene.gene$end_position),]
 
   return(median(segs.chr.end.start$Ratio))
}

plotTxCN <- function(cns, txs, file.name, main.text, xlab.text, ylab.text) {
   pdf(paste0(file.name, ".pdf"), height=4, width=6)
   plot(cns ~ txs, ylab=ylab.text, xlab=xlab.text, main=main.text[1])
 
   ##
   lm.fit <- lm(cns ~ txs)
   intercept <- coef(lm.fit)[1]
   r2 <- summary(lm.fit)$r.squared

   #text(0.18, 0.98, paste0("LM's R^2 = ", round0(r2*100, digits=2), "%"), cex=1.2)
   #text(0.142, 0.92, paste0("Intercept = ", round0(intercept, digits=2)), cex=1.2)
   abline(lm.fit)
   #mtext(main.text[2], cex=1.2, line=0.3)
   dev.off()
}

#samples.loh <- samples
genes <- c("SPG11", "B2M", "PATL2", "TRIM69", "BAX", "DHDH", "NUCB1")
genes <- c("CTDSPL2", "EIF3J")
for (g in 1:length(genes)) {
   gene <- genes[g]
   ensGene.gene <- subset(ensGene, external_gene_name == gene)
   ensembl_gene_id <- rownames(ensGene.gene)
   txs <- as.numeric(tpm.gene.log2[ensembl_gene_id, samples])
   xlab.text <- "log2(TPM+0.01)"
 
   cns <- c()   
   icns <- c()
   for (s in 1:length(samples)) {
      sample <- samples[s]

      ## INPUT: RPKM (from 1a) or SEQ
      segs <- read.peiflyne.cn.seg(file.path(wd.ngs, "peiflyne", sample, paste0(sample, "_ANALYSIS/", sample, "_cn.seg")))   ## See ReplicationTiming.R / Inner Class / PeifLyne File Reader
      cns <- c(cns, getCNFromSegment(ensGene.gene, segs))
      
      segs <- read.peiflyne.cn.seg(file.path(wd.ngs, "peiflyne", sample, paste0(sample, "_ANALYSIS/", sample, "_iCN.seg")))   ## See ReplicationTiming.R / Inner Class / PeifLyne File Reader
      icns <- c(icns, getCNFromSegment(ensGene.gene, segs))
   }
   main.text <- paste0(gene, " (", ensembl_gene_id, ") in ", BASE, " (n=", length(samples), ")")
   file.name <- file.path(wd.rt.plots, paste0("plot_", base, "_cn-tx_", gene, ".pdf"))
   plotTxCN(cns, txs, file.name, main.text, xlab.text, "Copy number")
   file.name <- file.path(wd.rt.plots, paste0("plot_", base, "_icn-tx_", gene, ".pdf"))
   plotTxCN(icns, txs, file.name, main.text, xlab.text, "iCN")
   
   data <- toTable(2, 5, length(icns), c("Sample", "Tx", "CN", "iCN", "iCN.ylab"))
   data$Sample <- samples
   data$Tx  <- txs
   data$CN  <- cns
   data$iCN <- icns
   data <- data[order(data$iCN),]
  
   loh.idx <- which(data$iCN < 2)
   cn2.idx <- which(data$iCN == 2)
   amp.idx <- which(data$iCN > 2)
   data$iCN.ylab[loh.idx] <- 1
   data$iCN.ylab[amp.idx] <- 3
   names <- c(paste0("iCN<2\nn=", length(loh.idx)), paste0("iCN=2\nn=", length(cn2.idx)), paste0("iCN>2\nn=", length(amp.idx)))
   #samples.loh <- intersect(samples.loh, samples[loh.idx])
   
   ## SRC
   #p <- testU(data$Tx[loh.idx], data$Tx[-loh.idx])
   rho <- cor.test(data$CN, data$Tx, method="spearman", exact=F)[[4]]
   p   <- cor.test(data$CN, data$Tx, method="spearman", exact=F)[[3]]
   
   file.name <- file.path(wd.rt.plots, paste0("boxplot_", base, "_icn-tx_", gene, ".pdf"))
   pdf(file.name, height=4, width=6)   #, units="in", res=300)
   boxplot(Tx ~ iCN.ylab, data=data, horizontal=T, las=1, xlab=xlab.text, names=names, main=main.text, outline=T)
   mtext.text <- paste0("P-value=", sprintf("%4.2e", p))
   mtext.text <- paste0(mtext.text, " at rho=", round0(rho, digits=2))
   if (p <= 1e-4) {
      mtext.text <- paste0("***", mtext.text, "    ")
   } else if (p <= 1e-3) {
      mtext.text <- paste0("**", mtext.text, "   ")
   } else if (p <= 1e-2) {
      mtext.text <- paste0("*", mtext.text, "  ")
   }
   mtext(mtext.text, cex=1.2, line=0.3)
   dev.off()
}
# > samples.loh 
# [1] "S00022" "S00035" "S00050" "S00356" "S00472" "S00501" "S00825" "S00827" "S00829" "S00830" "S00831" "S00832" "S00837" "S00838" "S01297" "S01366" "S01453" "S01494" "S01512" "S01524" "S01542"
# [22] "S01563"








plot(RFD_TTS ~ RFD_TSS, data=ensGene.rt.tx.plot, ylim=ylim, xlim=xlim, ylab=ylab.text, xlab=xlab.text, main=main.text[1], col=col1)
lines(x2,y2,col="green")




 
pdf(paste0(file.name, ".pdf"), height=6, width=6)
hist()
dev.off()




###
## Consistent + 50/50 (consistent)
file.name  <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_RFD_TTS~TSS+LENGTH_consist+5050_right"))
main.text <- c("Right-leading consistent genes in SCLC", paste0("(n=", separator(length(sclc.tx.consist.right)), "+", length(sclc.tx.5050.5050.right), ")"))
plotRFD(file.name, main.text, ensGene.rt.tx, sclc.tx.consist.right, "sandybrown", tx.list2=sclc.tx.5050.5050.right, col2="red", ext="pdf")

file.name  <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_RFD_TSS+TTS_consist+5050_left"))
main.text <- c("Left-leading consistent genes in SCLC", paste0("(n=", separator(length(sclc.tx.consist.left)), "+", length(sclc.tx.5050.5050.left), ")"))
plotRFD(file.name, main.text, ensGene.rt.tx, sclc.tx.consist.left, "steelblue1", tx.list2=sclc.tx.5050.5050.left, col2="red", ext="pdf")

## Right (CD and HO)
sclc.tx.5050.5050.right.cd <- rownames(subset(ensGene.rt.tx[sclc.tx.5050.5050.right,], CD > 0))
sclc.tx.5050.5050.right.ho <- rownames(subset(ensGene.rt.tx[sclc.tx.5050.5050.right,], CD < 0))

file.name  <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_RFD_TSS+TTS_consist+5050_right_CD"))
main.text <- c("Right-leading consistent genes in SCLC", paste0("(Co-directional; n=", separator(length(sclc.tx.consist.right.cd)), "+", length(sclc.tx.5050.5050.right.cd), ")"))
plotRFD(file.name, main.text, ensGene.rt.tx, sclc.tx.consist.right.cd, "sandybrown", tx.list2=sclc.tx.5050.5050.right.cd, col2="red", ext="pdf")

file.name  <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_RFD_TSS+TTS_consist+5050_right_HO"))
main.text <- c("Right-leading consistent genes in SCLC", paste0("(Head-on; n=", separator(length(sclc.tx.consist.right.ho)), "+", length(sclc.tx.5050.5050.right.ho), ")"))
plotRFD(file.name, main.text, ensGene.rt.tx, sclc.tx.consist.right.ho, "sandybrown", tx.list2=sclc.tx.5050.5050.right.ho, col2="red", ext="pdf")

## Left (CD and HO)
sclc.tx.5050.5050.left.cd <- rownames(subset(ensGene.rt.tx[sclc.tx.5050.5050.left,], CD > 0))
sclc.tx.5050.5050.left.ho <- rownames(subset(ensGene.rt.tx[sclc.tx.5050.5050.left,], CD < 0))

file.name  <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_RFD_TSS+TTS_consist+5050_left_CD"))
main.text <- c("Left-leading consistent genes in SCLC", paste0("(Co-directional; n=", separator(length(sclc.tx.consist.left.cd)), "+", length(sclc.tx.5050.5050.left.cd), ")"))
plotRFD(file.name, main.text, ensGene.rt.tx, sclc.tx.consist.left.cd, "steelblue1", tx.list2=sclc.tx.5050.5050.left.cd, col2="red", ext="pdf")

file.name  <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_RFD_TSS+TTS_consist+5050_left_HO"))
main.text <- c("Left-leading consistent genes in SCLC", paste0("(Head-on; n=", separator(length(sclc.tx.consist.left.ho)), "+", length(sclc.tx.5050.5050.left.ho), ")"))
plotRFD(file.name, main.text, ensGene.rt.tx, sclc.tx.consist.left.ho, "steelblue1", tx.list2=sclc.tx.5050.5050.left.ho, col2="red", ext="pdf")

###
## Consistent in SCLC and LCL
## Right (CD and HO)
getTSSTTSbetweenTissues <- function(tx.list, ensGene.rt.tx, ensGene.rt.tx.lcl) {
   overlaps <- intersect(rownames(ensGene.rt.tx), rownames(ensGene.rt.tx.lcl))
   overlaps <- intersect(tx.list, overlaps)
    
   tss <- which((ensGene.rt.tx[overlaps,]$RFD_TSS * ensGene.rt.tx.lcl[overlaps,]$RFD_TSS) > 0)
   tts <- which((ensGene.rt.tx[overlaps,]$RFD_TTS * ensGene.rt.tx.lcl[overlaps,]$RFD_TTS) > 0)
   
   return(overlaps[intersect(tss, tts)])
}

###
## Consistent + 50/50 (consistent) in SCLC and LCL
sclc.tx.consist.right.lcl <- getTSSTTSbetweenTissues(sclc.tx.consist.right, ensGene.rt.tx, ensGene.rt.tx.lcl)
sclc.tx.consist.left.lcl  <- getTSSTTSbetweenTissues(sclc.tx.consist.left, ensGene.rt.tx, ensGene.rt.tx.lcl)
sclc.tx.5050.5050.right.lcl <- getTSSTTSbetweenTissues(sclc.tx.5050.5050.right, ensGene.rt.tx, ensGene.rt.tx.lcl)
sclc.tx.5050.5050.left.lcl  <- getTSSTTSbetweenTissues(sclc.tx.5050.5050.left, ensGene.rt.tx, ensGene.rt.tx.lcl)

file.name  <- file.path(wd.rt.plots, paste0(base, "+lcl_ensGene.rt.tx_bstrps1000_RFD_TSS+TTS_consist+5050_right"))
main.text <- c("Right-leading consistent genes in SCLC and LCL", paste0("(n=", separator(length(sclc.tx.consist.right.lcl)), "+", length(sclc.tx.5050.5050.right.lcl), ")"))
plotRFD(file.name, main.text, ensGene.rt.tx, sclc.tx.consist.right.lcl, "sandybrown", tx.list2=sclc.tx.5050.5050.right.lcl, col2="red", ext="pdf")

file.name  <- file.path(wd.rt.plots, paste0(base, "+lcl_ensGene.rt.tx_bstrps1000_RFD_TSS+TTS_consist+5050_left"))
main.text <- c("Left-leading consistent genes in SCLC and LCL", paste0("(n=", separator(length(sclc.tx.consist.left.lcl)), "+", length(sclc.tx.5050.5050.left.lcl), ")"))
plotRFD(file.name, main.text, ensGene.rt.tx, sclc.tx.consist.left.lcl, "steelblue1", tx.list2=sclc.tx.5050.5050.left.lcl, col2="red", ext="pdf")

## Right (CD and HO)
sclc.tx.consist.right.cd.lcl <- getTSSTTSbetweenTissues(sclc.tx.consist.right.cd, ensGene.rt.tx, ensGene.rt.tx.lcl)
sclc.tx.consist.right.ho.lcl <- getTSSTTSbetweenTissues(sclc.tx.consist.right.ho, ensGene.rt.tx, ensGene.rt.tx.lcl)
sclc.tx.5050.5050.right.cd.lcl <- getTSSTTSbetweenTissues(sclc.tx.5050.5050.right.cd, ensGene.rt.tx, ensGene.rt.tx.lcl)
sclc.tx.5050.5050.right.ho.lcl <- getTSSTTSbetweenTissues(sclc.tx.5050.5050.right.ho, ensGene.rt.tx, ensGene.rt.tx.lcl)

file.name  <- file.path(wd.rt.plots, paste0(base, "+lcl_ensGene.rt.tx_bstrps1000_RFD_TSS+TTS_consist+5050_right_CD"))
main.text <- c("Right-leading consistent genes in SCLC", paste0("(Co-directional; n=", separator(length(sclc.tx.consist.right.cd.lcl)), "+", length(sclc.tx.5050.5050.right.cd.lcl), ")"))
plotRFD(file.name, main.text, ensGene.rt.tx, sclc.tx.consist.right.cd.lcl, "sandybrown", tx.list2=sclc.tx.5050.5050.right.cd.lcl, col2="red", ext="pdf")

file.name  <- file.path(wd.rt.plots, paste0(base, "+lcl_ensGene.rt.tx_bstrps1000_RFD_TSS+TTS_consist+5050_right_HO"))
main.text <- c("Right-leading consistent genes in SCLC", paste0("(Head-on; n=", separator(length(sclc.tx.consist.right.ho.lcl)), "+", length(sclc.tx.5050.5050.right.ho.lcl), ")"))
plotRFD(file.name, main.text, ensGene.rt.tx, sclc.tx.consist.right.ho.lcl, "sandybrown", tx.list2=sclc.tx.5050.5050.right.ho.lcl, col2="red", ext="pdf")

## Left (CD and HO)
sclc.tx.consist.left.cd.lcl <- getTSSTTSbetweenTissues(sclc.tx.consist.left.cd, ensGene.rt.tx, ensGene.rt.tx.lcl)
sclc.tx.consist.left.ho.lcl <- getTSSTTSbetweenTissues(sclc.tx.consist.left.ho, ensGene.rt.tx, ensGene.rt.tx.lcl)
sclc.tx.5050.5050.left.cd.lcl <- getTSSTTSbetweenTissues(sclc.tx.5050.5050.left.cd, ensGene.rt.tx, ensGene.rt.tx.lcl)
sclc.tx.5050.5050.left.ho.lcl <- getTSSTTSbetweenTissues(sclc.tx.5050.5050.left.ho, ensGene.rt.tx, ensGene.rt.tx.lcl)

file.name  <- file.path(wd.rt.plots, paste0(base, "+lcl_ensGene.rt.tx_bstrps1000_RFD_TSS+TTS_consist+5050_left_CD"))
main.text <- c("Left-leading consistent genes in SCLC", paste0("(Co-directional; n=", separator(length(sclc.tx.consist.left.cd.lcl)), "+", length(sclc.tx.5050.5050.left.cd.lcl), ")"))
plotRFD(file.name, main.text, ensGene.rt.tx, sclc.tx.consist.left.cd.lcl, "steelblue1", tx.list2=sclc.tx.5050.5050.left.cd.lcl, col2="red", ext="pdf")

file.name  <- file.path(wd.rt.plots, paste0(base, "+lcl_ensGene.rt.tx_bstrps1000_RFD_TSS+TTS_consist+5050_left_HO"))
main.text <- c("Left-leading consistent genes in SCLC", paste0("(Head-on; n=", separator(length(sclc.tx.consist.left.ho.lcl)), "+", length(sclc.tx.5050.5050.left.ho.lcl), ")"))
plotRFD(file.name, main.text, ensGene.rt.tx, sclc.tx.consist.left.ho.lcl, "steelblue1", tx.list2=sclc.tx.5050.5050.left.ho.lcl, col2="red", ext="pdf")









sclc.tx.5050.5050.lcl <- getTSSTTSbetweenTissues(sclc.tx.5050.5050, ensGene.rt.tx, ensGene.rt.tx.lcl)
sclc.tx.5050.5050.lcl.cd <- rownames(subset(ensGene.rt.tx[sclc.tx.5050.5050.lcl,], CD > 0))
sclc.tx.5050.5050.lcl.ho <- rownames(subset(ensGene.rt.tx[sclc.tx.5050.5050.lcl,], CD < 0))

sclc.tx.consist.right.cd <- rownames(subset(ensGene.rt.tx[sclc.tx.consist.right,], CD > 0))
sclc.tx.consist.right.ho <- rownames(subset(ensGene.rt.tx[sclc.tx.consist.right,], CD < 0))





file.name  <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_RFD_TSS+TTS_consist_right_CD"))
main.text <- c("Right-leading consistent genes in SCLC", paste0("(Co-directional; n=", separator(length(sclc.tx.consist.right.cd)), ")"))
plotRFD(file.name, main.text, ensGene.rt.tx, sclc.tx.consist.right.cd, "sandybrown", ext="pdf")

file.name  <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_RFD_TSS+TTS_consist_right_HO"))
main.text <- c("Right-leading consistent genes in SCLC", paste0("(Head-on; n=", separator(length(sclc.tx.consist.right.ho)), ")"))
plotRFD(file.name, main.text, ensGene.rt.tx, sclc.tx.consist.right.ho, "sandybrown", ext="pdf")

## Left (CD and HO)
sclc.tx.consist.left.cd <- rownames(subset(ensGene.rt.tx[sclc.tx.consist.left,], CD > 0))
sclc.tx.consist.left.ho <- rownames(subset(ensGene.rt.tx[sclc.tx.consist.left,], CD < 0))

file.name  <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_RFD_TSS+TTS_consist_left_CD"))
main.text <- c("Left-leading consistent genes in SCLC", paste0("(Co-directional; n=", separator(length(sclc.tx.consist.left.cd)), ")"))
plotRFD(file.name, main.text, ensGene.rt.tx, sclc.tx.consist.left.cd, "steelblue1", ext="pdf")

file.name  <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_RFD_TSS+TTS_consist_left_HO"))
main.text <- c("Left-leading consistent genes in SCLC", paste0("(Head-on; n=", separator(length(sclc.tx.consist.left.ho)), ")"))
plotRFD(file.name, main.text, ensGene.rt.tx, sclc.tx.consist.left.ho, "steelblue1", ext="pdf")







ensGene.rt.tx[sclc.tx.consist.right,]$RFD_TTS ~ ensGene.rt.tx[sclc.tx.consist.right,]$RFD_TSS

xlim <- 1
ylim <- 1
file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_RFD_TSS+TTS_consist_right.png"))
main.text <- c("Right-leading consistent genes in SCLC", paste0("(n=", separator(length(sclc.tx.inconsist)), ")"))
plotSRCTxLength(file.name, main.text, ensGene.rt.tx, sclc.tx.inconsist, "purple1", xlim, ylim)







###
##
file.name  <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_RFD_TSS+TTS_5050.pdf"))
main.text <- "Replication fork directionality (RFD) in SCLC"
xlab.text <- c("RFD = (R-L)/(R+L)", "TSS")
ylab.text <- "TTS"
sclc.tx <- c(sclc.tx.5050)
cols <- rep("red", length(sclc.tx.5050))

pdf(file.name, height=6, width=6)
plot(ensGene.rt.tx[sclc.tx,]$RFD_TTS ~ ensGene.rt.tx[sclc.tx,]$RFD_TSS, xlab=xlab.text, ylab=ylab.text, main=main.text, col=cols)

points(ensGene.rt.tx[sclc.tx.5050.right,]$RFD_TTS ~ ensGene.rt.tx[sclc.tx.5050.right,]$RFD_TSS, col="sandybrown")
points(ensGene.rt.tx[sclc.tx.5050.left,]$RFD_TTS  ~ ensGene.rt.tx[sclc.tx.5050.left,]$RFD_TSS,  col="steelblue1")
points(ensGene.rt.tx[sclc.tx.5050.5050,]$RFD_TTS  ~ ensGene.rt.tx[sclc.tx.5050.5050,]$RFD_TSS,  col="red")
abline(h=0, lty=5, lwd=0.85, col="black")
text(0,  0.5, paste0("32.29% (", separator(length(sclc.tx.5050.right)), ")"), cex=0.85, col="black")
text(0, -0.5, paste0("36.46% (", separator(length(sclc.tx.5050.left)), ")"), cex=0.85, col="black") 
text(0,    0, paste0("31.25% (", separator(length(sclc.tx.5050.5050)), ")"), cex=0.85, col="black") 

#legend("bottom", c("Left-leading in SCLC", paste0("L/R (n=", length(sclc.tx.boundary.o.c), ") in SCLC"), "Right-leading in SCLC"), col=c("steelblue1", "red", "sandybrown"), pch=1, cex=0.8, horiz=T)
mtext(paste0("R/L-leading genes (n=", separator(length(sclc.tx.5050)), ")"), cex=1, line=0.4)
dev.off()

# -----------------------------------------------------------------------------
# Plot RO and RT for these 30 candidate genes (sclc.tx.5050.5050)
# Last Modified: 05/12/18
# -----------------------------------------------------------------------------
wd.rt.plots.5050 <- file.path(wd.rt.plots, "5050")
for (g in 1:length(sclc.tx.5050.5050)) {
   g <- 19
   gene <- ensGene[sclc.tx.5050.5050[g],]
   chr <- gene$chromosome_name
   length <- gene$end_position - gene$start_position
   start <- gene$start_position - 10000
   end <- gene$end_position + 10000
   
   ## Replication fork directionality (RFD)  
   load(file=file.path(wd.rt.data, paste0("bed.gc.rt_", base, "_bstrps", bstrps, "_", chr, ".RData")))
   bed.gc.rt.chr$RFD <- mapply(x = 1:nrow(bed.gc.rt.chr), function(x) as.numeric(getRFD(bed.gc.rt.chr[x, ])))   ## ADD 05/12/18
   right.idx  <- which(bed.gc.rt.chr$RIGHT_LEADING > boundary.upper)
   left.idx   <- which(bed.gc.rt.chr$RIGHT_LEADING < boundary.lower)
   boundary.idx <- setdiff(c(1:nrow(bed.gc.rt.chr)), c(right.idx, left.idx))
   bed.gc.chr <- bed.gc[rownames(bed.gc.rt.chr),]
 
   start <- gene$start_position - 500000
   end <- gene$end_position + 500000
   file.name <- file.path(wd.rt.plots.5050, paste0(base, "_RFD_bstrps1000_", chr))
   plotBootstrapsRFD(file.name, BASE, chr, start, end, bed.gc.chr, bed.gc.rt.chr, right.idx, left.idx, boundary.idx, "png")       ## see ReplicationTiming.R
 
   ## Replication timing
   rt.chr <- readTable(file.path(wd.rt, "data", paste0(base, "_rpkm.corr.gc.d.rt_", chr, "_", BASE, "-", BASE, "_n101-92.txt.gz")), header=T, rownames=T, sep="\t") 
   rt.chr <- rt.chr[rownames(bed.gc.rt.chr),]
 
   file.name  <- file.path(wd.rt.plots.5050, paste0(base, "_RT_bstrps1000_", chr))
   plotBootstrapsRT(file.name, gene=gene$external_gene_name, BASE, chr, start, end, rt.chr, bed.gc.chr, right.idx, left.idx, boundary.idx, ymax=0.15, "png")   ## see ReplicationTiming.R
}


# -----------------------------------------------------------------------------
# Inconsistent genes have lower expression levels
# Last Modified: 12/11/18
# -----------------------------------------------------------------------------
testWbyEnsGeneRTTx(ensGene.rt.tx, sclc.tx.inconsist, sclc.tx.consist)
# [1] 0.01322831
testWbyEnsGeneRTTx(ensGene.rt.tx, sclc.tx.inconsist.right, sclc.tx.inconsist.left)
# [1] 0.1486169
testWbyEnsGeneRTTx(ensGene.rt.tx, sclc.tx.inconsist, sclc.tx.50.50)
# [1] 0.857977
testWbyEnsGeneRTTx(ensGene.rt.tx, sclc.tx.right, sclc.tx.left)
# [1] 0.66278
testWbyEnsGeneRTTx(ensGene.rt.tx, sclc.tx.consist, sclc.tx.50.50)
# [1] 0.5220981

file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_inconsistent_boxplot_pc1e2.pdf"))
main.text <- c("Inconsistent genes in SCLC", paste0("n=", separator(length(sclc.tx.inconsist))))
xlab.text <- "Replication time"
ylab.text <- "log2(TPM+0.01)"
ymin <- min(ensGene.rt.tx$MEDIAN)
ymax <- max(ensGene.rt.tx$MEDIAN)
overlaps <- intersect(sclc.tx.inconsist, sclc.tx.50.50)
names <- c(paste0("Consistent (n=", separator(length(sclc.tx.consist)),")"), paste0("Inconsistent genes (n=", separator(length(sclc.tx.inconsist)),")"))

pdf(file.name, height=6, width=3.5)   #, units="in", res=300)
boxplot( MEDIAN ~ CONSIST, data=ensGene.rt.tx[c(sclc.tc.consist, sclc.tc.inconsist),], ylim=c(ymin, ymax), ylab=ylab.text, xlab=xlab.text, names=names, main=main.text, border=c("black", "purple1"), outline=T)
#beeswarm(MEDIAN ~ RT, data=ensGene.rt.start.sclc.inconsist[overlaps,], col="red", pch=1, add=T)
legend("topleft", paste0("L/R (n=", separator(length(overlaps)), ")"), col="red", pch=1, cex=0.8)
dev.off()

# -----------------------------------------------------------------------------
# Inconsistent genes have lower expression levels
# Last Modified: 12/11/18
# -----------------------------------------------------------------------------
plotSRCTxLength <- function(file.name, main.text, ensGene.rt.tx, tx.list, col, xlim, ylim) {
   xlab.text <- "Gene length (log10"
   ylab.text <- "log2(TPM+0.01)"
   ensGene.rt.tx <- ensGene.rt.tx[tx.list,]
   ensGene.rt.tx$LENGTH <- log10(ensGene.rt.tx$LENGTH)
   
   if (is.null(xlim)) xlim <- c(min(ensGene.rt.tx$LENGTH), max(ensGene.rt.tx$LENGTH))
   if (is.null(ylim)) ylim <- c(min(ensGene.rt.tx$MEDIAN), max(ensGene.rt.tx$MEDIAN))
   xmin <- xlim[1]
   xmax <- xlim[2]
   ymin <- ylim[1]
   ymax <- ylim[2]
   pdf(file.name, height=6, width=6)
   plot(MEDIAN ~ LENGTH, data=ensGene.rt.tx, ylim=ylim, xlim=xlim, ylab=ylab.text, xlab=xlab.text, main=main.text[1], col=col)

   ## Spearman's rank correlation
   ## https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cor.test.html
   rho    <- cor.test(ensGene.rt.tx$MEDIAN, ensGene.rt.tx$LENGTH, method="spearman", exact=F)[[4]]
   pvalue <- cor.test(ensGene.rt.tx$MEDIAN, ensGene.rt.tx$LENGTH, method="spearman", exact=F)[[3]]
   
   if (rho > 0)
      text(xmin + (xmax - xmin)/6.5, ymax, paste0("Spearman's rho = ", sprintf("%.2f", round(rho, digits=2))), cex=0.95, col=col)
   else
      text(xmin + (xmax - xmin)/6.1, ymax, paste0("Spearman's rho = ", sprintf("%.2f", round(rho, digits=2))), cex=0.95, col=col)
   text(xmin + (xmax - xmin)/7.5, ymax - (ymax-ymin)/17, paste0("p-value = ", formatC(pvalue, format="E", digits=2)), cex=0.95, col=col)

   ## Linear regression
   ## https://www.statlect.com/fundamentals-of-statistics/R-squared-of-a-linear-regression
   lm.fit <- lm(MEDIAN ~ LENGTH, data=ensGene.rt.tx)
   intercept <- coef(lm.fit)[1]
   slope <- summary(lm.fit)$coefficients[2,1]
   p  <- anova(lm.fit)$'Pr(>F)'[1]
   r2 <- summary(lm.fit)$r.squared
   
   if (slope > 0)
      text(xmin + (xmax - xmin)/8.6, intercept + slope*xmin - (ymax-ymin)/18, paste0("LM's slope = ", sprintf("%.2f", round(slope, digits=2))), cex=0.95)   #col="gray55" 
   else
      text(xmin + (xmax - xmin)/8.1, intercept + slope*xmin - (ymax-ymin)/18, paste0("LM's slope = ", sprintf("%.2f", round(slope, digits=2))), cex=0.95) 
   text(xmin + (xmax - xmin)/7.5,  intercept + slope*xmin - (ymax-ymin)/18*2, paste0("p-value = ", formatC(p, format="E", digits=2)), cex=0.95)
   text(xmin + (xmax - xmin)/13.2, intercept + slope*xmin - (ymax-ymin)/18*3, paste0("R^2 = ", sprintf("%.2f", round(r2*100, digits=2)), "%"), cex=0.95)
   
   abline(lm.fit)
   mtext(main.text[2], cex=1.2, line=0.5)
   dev.off()
}

xlim <- c(log10(min(ensGene.rt.tx$LENGTH)), log10(max(ensGene.rt.tx$LENGTH)))
ylim <- c(min(ensGene.rt.tx$MEDIAN), max(ensGene.rt.tx$MEDIAN))
file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_length_inconsist_pc1e2.pdf"))
main.text <- c("Inconsistent genes in SCLC", paste0("(n=", separator(length(sclc.tx.inconsist)), ")"))
plotSRCTxLength(file.name, main.text, ensGene.rt.tx, sclc.tx.inconsist, "purple1", xlim, ylim)

file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_length_right_pc1e2.pdf"))
main.text <- c("Right-leading genes in SCLC", paste0("(n=", separator(length(sclc.tx.right)), ")"))
plotSRCTxLength(file.name, main.text, ensGene.rt.tx, sclc.tx.right, "sandybrown", xlim, ylim)

file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_length_left_pc1e2.pdf"))
main.text <- c("Left-leading genes in SCLC", paste0("(n=", separator(length(sclc.tx.left)), ")"))
plotSRCTxLength(file.name, main.text, ensGene.rt.tx, sclc.tx.left, "steelblue1", xlim, ylim)

file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_length_5050_pc1e2.pdf"))
main.text <- c("Left-/Right-leading genes in SCLC", paste0("(n=", separator(length(sclc.tx.50.50)), ")"))
plotSRCTxLength(file.name, main.text, ensGene.rt.tx, sclc.tx.50.50, "red", xlim, ylim)

# -----------------------------------------------------------------------------
# Inconsistent genes have lower expression levels (CUT3; Final TCR gene list)
# Last Modified: 26/11/18
# -----------------------------------------------------------------------------
ensGene.input <- ensGene.tx.rt.nona.sign.ca
ensGene.input <- subset(ensGene.input, MUT_CG >= 3)
ensGene.input <- subset(ensGene.input, MUT_GC >= 3)

overlaps <- intersect(rownames(ensGene.rt.tx), rownames(ensGene.input))
ensGene.rt.tx.ca <- cbind(ensGene.rt.tx[overlaps,], ensGene.input[overlaps, 13:16])
sclc.tx.right.ca <- intersect(sclc.tx.right, overlaps)
sclc.tx.left.ca  <- intersect(sclc.tx.left, overlaps)
sclc.tx.50.50.ca <- intersect(sclc.tx.50.50, overlaps)
sclc.tx.consist.ca <- intersect(sclc.tx.consist, overlaps)
sclc.tx.inconsist.ca <- intersect(sclc.tx.inconsist, overlaps)

xlim <- c(log10(min(ensGene.rt.tx.ca$LENGTH)), log10(max(ensGene.rt.tx.ca$LENGTH)/1000000))
ylim <- c(min(ensGene.rt.tx.ca$MEDIAN), max(ensGene.rt.tx.ca$MEDIAN))
file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_length_inconsist_pc1e2_CA_Mb.pdf"))
main.text <- c("Inconsistent genes in SCLC", paste0("(n=", separator(length(sclc.tx.inconsist.ca)), "; C>A/G>T)"))
plotSRCTxLength(file.name, main.text, ensGene.rt.tx.ca, sclc.tx.inconsist.ca, "purple1", xlim, ylim)

file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_length_right_pc1e2_CA_Mb.pdf"))
main.text <- c("Right-leading genes in SCLC", paste0("(n=", separator(length(sclc.tx.right.ca)), "; C>A/G>T)"))
plotSRCTxLength(file.name, main.text, ensGene.rt.tx.ca, sclc.tx.right.ca, "sandybrown", xlim, ylim)

file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_length_left_pc1e2_CA_Mb.pdf"))
main.text <- c("Left-leading genes in SCLC", paste0("(n=", separator(length(sclc.tx.left.ca)), "; C>A/G>T)"))
plotSRCTxLength(file.name, main.text, ensGene.rt.tx.ca, sclc.tx.left.ca, "steelblue1", xlim, ylim)

file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_length_5050_pc1e2_CA_Mb.pdf"))
main.text <- c("Left-/Right-leading genes in SCLC", paste0("(n=", separator(length(sclc.tx.50.50.ca)), "; C>A/G>T)"))
plotSRCTxLength(file.name, main.text, ensGene.rt.tx.ca, sclc.tx.50.50.ca, "red", xlim, ylim)

# -----------------------------------------------------------------------------
# Inconsistent genes have lower expression levels (CUT3; Final TCR gene list)
# Last Modified: 12/11/18
# -----------------------------------------------------------------------------
plotSRCMutLength <- function(file.name, main.text, ensGene.rt.tx, tx.list, col, xlim, ylim) {
   xlab.text <- "Gene length (log10)"
   ylab.text <- "Mutation number (log10)"
   ensGene.rt.tx <- ensGene.rt.tx[tx.list,]
   ensGene.rt.tx$LENGTH <- log10(ensGene.rt.tx$LENGTH)
   ensGene.rt.tx$MUT <- log10(ensGene.rt.tx$MUT)
   
   if (is.null(xlim)) xlim <- c(min(ensGene.rt.tx$LENGTH), max(ensGene.rt.tx$LENGTH))
   if (is.null(ylim)) ylim <- c(min(ensGene.rt.tx$MUT), max(ensGene.rt.tx$MUT))
   xmin <- xlim[1]
   xmax <- xlim[2]
   ymin <- ylim[1]
   ymax <- ylim[2]
   pdf(file.name, height=6, width=6)
   plot(MUT ~ LENGTH, data=ensGene.rt.tx, ylim=ylim, xlim=xlim, ylab=ylab.text, xlab=xlab.text, main=main.text[1], col=col)
 
   ## Spearman's rank correlation
   ## https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cor.test.html
   rho    <- cor.test(ensGene.rt.tx$MUT, ensGene.rt.tx$LENGTH, method="spearman", exact=F)[[4]]
   pvalue <- cor.test(ensGene.rt.tx$MUT, ensGene.rt.tx$LENGTH, method="spearman", exact=F)[[3]]
 
   if (rho > 0)
      text(xmin + (xmax - xmin)/6.5, ymax, paste0("Spearman's rho = ", sprintf("%.2f", round(rho, digits=2))), cex=0.95, col=col)
   else
      text(xmin + (xmax - xmin)/6.1, ymax, paste0("Spearman's rho = ", sprintf("%.2f", round(rho, digits=2))), cex=0.95, col=col)
   
   if (pvalue >= 1e-100)
      text(xmin + (xmax - xmin)/7.5, ymax - (ymax-ymin)/17, paste0("p-value = ", formatC(pvalue, format="E", digits=2)), cex=0.95, col=col)
   else if (pvalue == 0)
      text(xmin + (xmax - xmin)/6.7, ymax - (ymax-ymin)/17, "p-value < 2.2E-16 ***", cex=0.95, col=col)
   else
      text(xmin + (xmax - xmin)/6.95, ymax - (ymax-ymin)/17, paste0("p-value = ", formatC(pvalue, format="E", digits=2)), cex=0.95, col=col)
 
   ## Linear regression
   ## https://www.statlect.com/fundamentals-of-statistics/R-squared-of-a-linear-regression
   lm.fit <- lm(MUT ~ LENGTH, data=ensGene.rt.tx)
   intercept <- coef(lm.fit)[1]
   slope <- summary(lm.fit)$coefficients[2,1]
   p  <- anova(lm.fit)$'Pr(>F)'[1]
   r2 <- summary(lm.fit)$r.squared
 
   if (slope > 0)
      text(xmin + (xmax - xmin)/8.6, (ymax-ymin)*3/4 - (ymax-ymin)/18, paste0("LM's slope = ", sprintf("%.2f", round(slope, digits=2))), cex=0.95)   #col="gray55" 
   else
      text(xmin + (xmax - xmin)/8.1, (ymax-ymin)*3/4 - (ymax-ymin)/18, paste0("LM's slope = ", sprintf("%.2f", round(slope, digits=2))), cex=0.95) 
   
   if (p >= 1e-100)
      text(xmin + (xmax - xmin)/7.45, (ymax-ymin)*3/4 - (ymax-ymin)/18*2, paste0("p-value = ", formatC(p, format="E", digits=2)), cex=0.95)
   else if (p == 0)
      text(xmin + (xmax - xmin)/6.7, (ymax-ymin)*3/4 - (ymax-ymin)/18*2, "p-value < 2.2E-16 ***", cex=0.95)
   else
      text(xmin + (xmax - xmin)/6.95, (ymax-ymin)*3/4 - (ymax-ymin)/18*2, paste0("p-value = ", formatC(p, format="E", digits=2)), cex=0.95)
   
   if (r2 < 0.1)  
      text(xmin + (xmax - xmin)/13.2, (ymax-ymin)*3/4 - (ymax-ymin)/18*3, paste0("R^2 = ", sprintf("%.2f", round(r2*100, digits=2)), "%"), cex=0.95)
   else
      text(xmin + (xmax - xmin)/11.6, (ymax-ymin)*3/4 - (ymax-ymin)/18*3, paste0("R^2 = ", sprintf("%.2f", round(r2*100, digits=2)), "%"), cex=0.95)
   
   abline(lm.fit)
   mtext(main.text[2], cex=1.2, line=0.5)
   dev.off()
}

xlim <- log10(c(min(ensGene.rt.tx.ca$LENGTH), max(ensGene.rt.tx.ca$LENGTH)))
ylim <- log10(c(min(ensGene.rt.tx.ca$MUT), max(ensGene.rt.tx.ca$MUT)))
file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_mut_inconsist_CA.pdf"))
main.text <- c("Inconsistent genes in SCLC", paste0("(n=", separator(length(sclc.tx.inconsist.ca)), "; C>A/G>T)"))
plotSRCMutLength(file.name, main.text, ensGene.rt.tx.ca, sclc.tx.inconsist.ca, "purple1", xlim, ylim)

file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_mut_right_CA.pdf"))
main.text <- c("Right-leading genes in SCLC", paste0("(n=", separator(length(sclc.tx.right.ca)), "; C>A/G>T)"))
plotSRCMutLength(file.name, main.text, ensGene.rt.tx.ca, sclc.tx.right.ca, "sandybrown", xlim, ylim)

file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_mut_left_CA.pdf"))
main.text <- c("Left-leading genes in SCLC", paste0("(n=", separator(length(sclc.tx.left.ca)), "; C>A/G>T)"))
plotSRCMutLength(file.name, main.text, ensGene.rt.tx.ca, sclc.tx.left.ca, "steelblue1", xlim, ylim)

file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_mut_5050_CA.pdf"))
main.text <- c("Left-/Right-leading genes in SCLC", paste0("(n=", separator(length(sclc.tx.50.50.ca)), "; C>A/G>T)"))
plotSRCMutLength(file.name, main.text, ensGene.rt.tx.ca, sclc.tx.50.50.ca, "red", xlim, ylim)

##
ensGene.rt.tx.ca.right <- ensGene.rt.tx.ca[sclc.tx.right.ca,]
cor <- cor.test(ensGene.rt.tx.ca.right$MUT, ensGene.rt.tx.ca.right$LENGTH, method="spearman", exact=F)
lm.fit <- lm(MUT ~ LENGTH, data=ensGene.rt.tx.ca.right)

ensGene.rt.tx.ca.inconsist <- ensGene.rt.tx.ca[sclc.tx.inconsist.ca,]
cor <- cor.test(ensGene.rt.tx.ca.inconsist$MUT, ensGene.rt.tx.ca.inconsist$LENGTH, method="spearman", exact=F)
lm.fit <- lm(MUT ~ LENGTH, data=ensGene.rt.tx.ca.inconsist)

# -----------------------------------------------------------------------------
# Inconsistent genes have lower expression levels (CUT3; Final TCR gene list)
# Last Modified: 12/11/18
# -----------------------------------------------------------------------------
plotSRCMutTx <- function(file.name, main.text, ensGene.rt.tx, tx.list, col, xlim, ylim) {
   xlab.text <- "log2(TPM+0.01)"
   ylab.text <- "Mutation number (log10)"
   ensGene.rt.tx <- ensGene.rt.tx[tx.list,]
   ensGene.rt.tx$MUT <- log10(ensGene.rt.tx$MUT)
 
   if (is.null(xlim)) xlim <- c(min(ensGene.rt.tx$MEDIAN), max(ensGene.rt.tx$MEDIAN))
   if (is.null(ylim)) ylim <- c(min(ensGene.rt.tx$MUT), max(ensGene.rt.tx$MUT))
   xmin <- xlim[1]
   xmax <- xlim[2]
   ymin <- ylim[1]
   ymax <- ylim[2]
   pdf(file.name, height=6, width=6)
   plot(MUT ~ MEDIAN, data=ensGene.rt.tx, ylim=ylim, xlim=xlim, ylab=ylab.text, xlab=xlab.text, main=main.text[1], col=col)
 
   ## Spearman's rank correlation
   ## https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cor.test.html
   rho    <- cor.test(ensGene.rt.tx$MUT, ensGene.rt.tx$MEDIAN, method="spearman", exact=F)[[4]]
   pvalue <- cor.test(ensGene.rt.tx$MUT, ensGene.rt.tx$MEDIAN, method="spearman", exact=F)[[3]]
 
   if (rho > 0)
      text(xmin + (xmax - xmin)/6.5, ymax, paste0("Spearman's rho = ", sprintf("%.2f", round(rho, digits=2))), cex=0.95, col=col)
   else
      text(xmin + (xmax - xmin)/6.1, ymax, paste0("Spearman's rho = ", sprintf("%.2f", round(rho, digits=2))), cex=0.95, col=col)
 
   if (pvalue >= 1e-100)
      text(xmin + (xmax - xmin)/7.5, ymax - (ymax-ymin)/17, paste0("p-value = ", formatC(pvalue, format="E", digits=2)), cex=0.95, col=col)
   else if (pvalue == 0)
      text(xmin + (xmax - xmin)/6.75, ymax - (ymax-ymin)/17, "p-value < 2.2E-16 ***", cex=0.95, col=col)
   else
      text(xmin + (xmax - xmin)/7, ymax - (ymax-ymin)/17, paste0("p-value = ", formatC(pvalue, format="E", digits=2)), cex=0.95, col=col)
 
   ## Linear regression
   ## https://www.statlect.com/fundamentals-of-statistics/R-squared-of-a-linear-regression
   lm.fit <- lm(MUT ~ MEDIAN, data=ensGene.rt.tx)
   intercept <- coef(lm.fit)[1]
   slope <- summary(lm.fit)$coefficients[2,1]
   p  <- anova(lm.fit)$'Pr(>F)'[1]
   r2 <- summary(lm.fit)$r.squared
 
   if (slope > 0)
      text(xmin + (xmax - xmin)/8.6, (ymax-ymin)*3/4 - (ymax-ymin)/18, paste0("LM's slope = ", sprintf("%.2f", round(slope, digits=2))), cex=0.95)   #col="gray55" 
   else
      text(xmin + (xmax - xmin)/8.1, (ymax-ymin)*3/4 - (ymax-ymin)/18, paste0("LM's slope = ", sprintf("%.2f", round(slope, digits=2))), cex=0.95) 
 
   if (p >= 1e-100)
      text(xmin + (xmax - xmin)/7.5, (ymax-ymin)*3/4 - (ymax-ymin)/18*2, paste0("p-value = ", formatC(p, format="E", digits=2)), cex=0.95)
   else if (p == 0)
      text(xmin + (xmax - xmin)/6.75, (ymax-ymin)*3/4 - (ymax-ymin)/18*2, "p-value < 2.2E-16 ***", cex=0.95)
   else
      text(xmin + (xmax - xmin)/7, (ymax-ymin)*3/4 - (ymax-ymin)/18*2, paste0("p-value = ", formatC(p, format="E", digits=2)), cex=0.95)
 
   if (r2 < 0.1)  
      text(xmin + (xmax - xmin)/13.2, (ymax-ymin)*3/4 - (ymax-ymin)/18*3, paste0("R^2 = ", sprintf("%.2f", round(r2*100, digits=2)), "%"), cex=0.95)
   else
      text(xmin + (xmax - xmin)/11.5, (ymax-ymin)*3/4 - (ymax-ymin)/18*3, paste0("R^2 = ", sprintf("%.2f", round(r2*100, digits=2)), "%"), cex=0.95)
 
   abline(lm.fit)
   mtext(main.text[2], cex=1.2, line=0.5)
   dev.off()
}

xlim <- c(min(ensGene.rt.tx.ca$MEDIAN), max(ensGene.rt.tx.ca$MEDIAN))
ylim <- log10(c(min(ensGene.rt.tx.ca$MUT), max(ensGene.rt.tx.ca$MUT)))
file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_tx_inconsist_CA.pdf"))
main.text <- c("Inconsistent genes in SCLC", paste0("(n=", separator(length(sclc.tx.inconsist.ca)), "; C>A/G>T)"))
plotSRCMutTx(file.name, main.text, ensGene.rt.tx.ca, sclc.tx.inconsist.ca, "purple1", xlim, ylim)

file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_tx_right_CA.pdf"))
main.text <- c("Right-leading genes in SCLC", paste0("(n=", separator(length(sclc.tx.right.ca)), "; C>A/G>T)"))
plotSRCMutTx(file.name, main.text, ensGene.rt.tx.ca, sclc.tx.right.ca, "sandybrown", xlim, ylim)

file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_tx_left_CA.pdf"))
main.text <- c("Left-leading genes in SCLC", paste0("(n=", separator(length(sclc.tx.left.ca)), "; C>A/G>T)"))
plotSRCMutTx(file.name, main.text, ensGene.rt.tx.ca, sclc.tx.left.ca, "steelblue1", xlim, ylim)

file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_tx_5050_CA.pdf"))
main.text <- c("Left-/Right-leading genes in SCLC", paste0("(n=", separator(length(sclc.tx.50.50.ca)), "; C>A/G>T)"))
plotSRCMutTx(file.name, main.text, ensGene.rt.tx.ca, sclc.tx.50.50.ca, "red", xlim, ylim)







xlab.text <- "Gene length (log10)"
ylab.text <- "Expression (log2)"

pdf(file.name, height=6, width=6)
plot(  MEDIAN ~ log10(LENGTH), data=subset(ensGene.rt.start.sclc.inconsist, RATIO > 0), ylab=ylab.text, xlab=xlab.text, main=main.text, col="sandybrown")
points(MEDIAN ~ log10(LENGTH), data=subset(ensGene.rt.start.sclc.inconsist, RATIO < 0), col="steelblue1", pch=1)
abline(lm(MEDIAN ~ log10(LENGTH), data=subset(ensGene.rt.start.sclc.inconsist, RATIO > 0)), col="sandybrown")
abline(lm(MEDIAN ~ log10(LENGTH), data=subset(ensGene.rt.start.sclc.inconsist, RATIO < 0)), col="steelblue1")
legend("topright", c("Left-leading", "Right-leading"), col=c("steelblue1", "sandybrown"), pch=1, cex=0.8, horiz=F)
dev.off()

##
xlab.text <- "Gene length (log10)"
ylab.text <- "Expression (log2)"
ymin <- min(ensGene.rt.tx$MEDIAN)
ymax <- max(ensGene.rt.tx$MEDIAN)
xmin <- min(ensGene.rt.tx$MEDIAN)
xmax <- max(ensGene.rt.tx$MEDIAN)

file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_length_inconsistent.pdf"))
main.text <- c("Right-leading genes in SCLC", paste0("n=", separator(length(sclc.tx.right))))
pdf(file.name, height=6, width=6)
plot(MEDIAN ~ log10(LENGTH), data=ensGene.rt.tx[sclc.tx.inconsist,], ylab=ylab.text, xlab=xlab.text, main=main.text, col="purple1")
abline(lm(MEDIAN ~ log10(LENGTH), data=ensGene.rt.tx[sclc.tx.inconsist,]), ylim=c(ymin, ymax), col="black")
dev.off()
cor.test(ensGene.rt.tx[sclc.tx.inconsist,]$MEDIAN, log10(ensGene.rt.tx[sclc.tx.inconsist,]$LENGTH), method="spearman", exact=F)[[4]]
cor.test(ensGene.rt.tx[sclc.tx.inconsist,]$MEDIAN, log10(ensGene.rt.tx[sclc.tx.inconsist,]$LENGTH), method="spearman", exact=F)[[3]]


file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_length_right.pdf"))
main.text <- c("Right-leading genes in SCLC", paste0("n=", separator(length(sclc.tx.right))))
pdf(file.name, height=6, width=6)
plot(MEDIAN ~ log10(LENGTH), data=ensGene.rt.tx[sclc.tx.right,], ylab=ylab.text, xlab=xlab.text, main=main.text, col="sandybrown")
abline(lm(MEDIAN ~ log10(LENGTH), data=ensGene.rt.tx[sclc.tx.right,]), col="black")
dev.off()
cor.test(ensGene.rt.tx[sclc.tx.right,]$MEDIAN, log10(ensGene.rt.tx[sclc.tx.right,]$LENGTH), method="spearman", exact=F)[[4]]
cor.test(ensGene.rt.tx[sclc.tx.right,]$MEDIAN, log10(ensGene.rt.tx[sclc.tx.right,]$LENGTH), method="spearman", exact=F)[[3]]


##
file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_length_left.pdf"))
main.text <- c("Left-leading genes in SCLC", paste0("n=", separator(length(sclc.tx.left))))
pdf(file.name, height=6, width=6)
plot(MEDIAN ~ log10(LENGTH), data=ensGene.rt.tx[sclc.tx.left, ], ylab=ylab.text, xlab=xlab.text, main=main.text, col="steelblue1")
abline(lm(MEDIAN ~ log10(LENGTH), data=ensGene.rt.tx[sclc.tx.left,]), col="black")
dev.off()
cor.test(ensGene.rt.tx[sclc.tx.left,]$MEDIAN, log10(ensGene.rt.tx[sclc.tx.left,]$LENGTH), method="spearman", exact=F)[[4]]
cor.test(ensGene.rt.tx[sclc.tx.left,]$MEDIAN, log10(ensGene.rt.tx[sclc.tx.left,]$LENGTH), method="spearman", exact=F)[[3]]


##
file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_length_5050.pdf"))
main.text <- c("Left/Right genes in SCLC", paste0("n=", separator(length(sclc.tx.left))))
pdf(file.name, height=6, width=6)
plot(MEDIAN ~ log10(LENGTH), data=ensGene.rt.tx[sclc.tx.50.50,], ylab=ylab.text, xlab=xlab.text, main=main.text, col="red")
abline(lm(MEDIAN ~ log10(LENGTH), data=ensGene.rt.tx[sclc.tx.50.50,]), col="black")
dev.off()
cor.test(ensGene.rt.tx[sclc.tx.50.50,]$MEDIAN, log10(ensGene.rt.tx[sclc.tx.50.50,]$LENGTH), method="spearman", exact=F)[[4]]
cor.test(ensGene.rt.tx[sclc.tx.50.50,]$MEDIAN, log10(ensGene.rt.tx[sclc.tx.50.50,]$LENGTH), method="spearman", exact=F)[[3]]



##
file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_Start+End.png"))
plotEnsGeneRTTxRO(file.name, ensGene.rt.tx, boundary.upper, "png")

sclc.tx.right  <- rownames(subset(ensGene.rt.start.sclc, RIGHT_LEADING > cutoff.upper))
sclc.tx.left   <- rownames(subset(ensGene.rt.start.sclc, RIGHT_LEADING < cutoff.lower))
sclc.tx.boundary <- setdiff(rownames(ensGene.rt.start.sclc), c(sclc.tx.right, sclc.tx.left))
sclc.tx.inconsist <- rownames(subset(ensGene.rt.start.sclc, SIGN < 0))
## ADD 20/11/18: Remove inconsistent genes
sclc.tx.inconsist <- setdiff(sclc.tx.inconsist, sclc.tx.boundary)
sclc.tx.right     <- setdiff(sclc.tx.right, sclc.tx.inconsist)
sclc.tx.left      <- setdiff(sclc.tx.left,  sclc.tx.inconsist)
## ADD 20/11/18: How many boundarys are inconsistent (which will be visulised in the boxplot later)
length(intersect(sclc.tx.boundary, sclc.tx.inconsist))
# [1] 45
# > 8486+8442+1098+97-45   ## Right+Left+Inconsistent+boundary-Overlaps
# [1] 18078

file.name  <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_Start+End.png"))
main.text  <- "Transcription vs. replication time in SCLC"
mtext.text <- paste0("Expressed genes (n=", separator(nrow(ensGene.rt.start.sclc)), ")")
xlab.text <- "Right-leading ratio (log2)"
ylab.text <- "log2(TPM+1)"

pdf(file.name, height=6, width=6)   #, units="in", res=300)
cols <- rep("steelblue1", nrow(ensGene.rt.start.sclc))
plot(  MEDIAN ~ RATIO, data=ensGene.rt.start.sclc, xlab=xlab.text, ylab=ylab.text, main=main.text, col=cols)
points(MEDIAN ~ RATIO, data=ensGene.rt.start.sclc[sclc.tx.right,],     col="sandybrown")
points(MEDIAN ~ RATIO, data=ensGene.rt.start.sclc[sclc.tx.inconsist,], col="purple1")
points(MEDIAN ~ RATIO, data=ensGene.rt.start.sclc[sclc.tx.boundary,],    col="red")
legend("top", c("Left-leading", paste0("L/R (n=", separator(length(sclc.tx.boundary)), ")"), "Inconsistent", "Right-leading"), col=c("steelblue1", "red", "purple1", "sandybrown"), pch=1, cex=0.75, horiz=T)
mtext(mtext.text, cex=1, line=0.5)
dev.off()

# -----------------------------------------------------------------------------
# Inconsistent genes have lower expression levels
# Last Modified: 12/11/18
# -----------------------------------------------------------------------------
sclc.tx.inconsist <- subset(ensGene.rt.start.sclc, SIGN == -1)$MEDIAN
sclc.tx.consist   <- subset(ensGene.rt.start.sclc, SIGN == 1)$MEDIAN
testW(sclc.tx.inconsist, sclc.tx.consist)
#[1] 0.007379227

file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_inconsistent_beeswarm.pdf"))
main.text <- c("Inconsistent genes in SCLC", paste0("n=", separator(length(sclc.tx.inconsist))))
xlab.text <- "Replication time"
ylab.text <- "log2(TPM+1)"
ymin <- min(ensGene.rt.start.sclc$MEDIAN)
ymax <- max(ensGene.rt.start.sclc$MEDIAN)
overlaps <- intersect(sclc.tx.inconsist, sclc.tx.boundary)
names <- c(paste0("Consistent (n=", separator(length(sclc.tx.consist)),")"), paste0("Inconsistent (n=", separator(length(sclc.tx.inconsist)),")"))

pdf(file.name, height=6, width=3.5)   #, units="in", res=300)
boxplot( MEDIAN ~ RT, data=ensGene.rt.start.sclc.inconsist, ylim=c(ymin, ymax), ylab=ylab.text, xlab=xlab.text, names=names, main=main.text, border="purple1", outline=T)
beeswarm(MEDIAN ~ RT, data=ensGene.rt.start.sclc.inconsist[overlaps,], col="red", pch=1, add=T)
legend("topleft", paste0("L/R (n=", separator(length(overlaps)), ")"), col="red", pch=1, cex=0.8)
dev.off()



# -----------------------------------------------------------------------------
# Inconsistent genes have lower expression levels (TO-DO)
# Last Modified: 12/11/18
# -----------------------------------------------------------------------------
ensGene.rt.start.sclc.inconsist <- ensGene.rt.start.sclc[sclc.tx.inconsist,]
sclc.tx.inconsist.right <- subset(ensGene.rt.start.sclc.inconsist, RATIO > 0)$MEDIAN
sclc.tx.inconsist.left  <- subset(ensGene.rt.start.sclc.inconsist, RATIO < 0)$MEDIAN
testW(sclc.tx.inconsist.right, sclc.tx.inconsist.left)
# [1] 1.235423e-15

file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_inconsistent_beeswarm.pdf"))
main.text <- c("Inconsistent genes in SCLC", paste0("n=", separator(length(sclc.tx.inconsist))))
xlab.text <- "Replication time"
ylab.text <- "log2(TPM+1)"
ymin <- min(ensGene.rt.start.sclc$MEDIAN)
ymax <- max(ensGene.rt.start.sclc$MEDIAN)
overlaps <- intersect(sclc.tx.inconsist, sclc.tx.boundary)
names <- c(paste0("L (n=", separator(length(sclc.tx.inconsist.left)),")"), paste0("R (n=", separator(length(sclc.tx.inconsist.right)),")"))

pdf(file.name, height=6, width=3.5)   #, units="in", res=300)
boxplot( MEDIAN ~ RT, data=ensGene.rt.start.sclc.inconsist, ylim=c(ymin, ymax), ylab=ylab.text, xlab=xlab.text, names=names, main=main.text, border="purple1", outline=T)
beeswarm(MEDIAN ~ RT, data=ensGene.rt.start.sclc.inconsist[overlaps,], col="red", pch=1, add=T)
legend("topleft", paste0("L/R (n=", separator(length(overlaps)), ")"), col="red", pch=1, cex=0.8)
dev.off()

##
sclc.tx.all <- ensGene.rt.start.sclc$MEDIAN
testW(sclc.tx.all, c(sclc.tx.inconsist.right, sclc.tx.inconsist.left))
# [1] 0.01134909
testW(sclc.tx.all, sclc.tx.inconsist.right)
# [1] 8.092788e-05
testW(sclc.tx.all, sclc.tx.inconsist.left)
# [1] 3.372877e-11

##
sclc.tx.inconsist.right.length <- subset(ensGene.rt.start.sclc.inconsist, RATIO > 0)$LENGTH
sclc.tx.inconsist.left.length  <- subset(ensGene.rt.start.sclc.inconsist, RATIO < 0)$LENGTH
testW(sclc.tx.inconsist.right.length, sclc.tx.inconsist.left.length)
# [1] 2.826203e-13

##
file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_length_inconsistent.pdf"))
main.text <- c("Inconsistent genes in SCLC", paste0("n=", separator(length(sclc.tx.inconsist))))
xlab.text <- "Gene length (log10)"
ylab.text <- "log2(TPM+1)"

pdf(file.name, height=6, width=6)
plot(  MEDIAN ~ log10(LENGTH), data=subset(ensGene.rt.start.sclc.inconsist, RATIO > 0), ylab=ylab.text, xlab=xlab.text, main=main.text, col="sandybrown")
points(MEDIAN ~ log10(LENGTH), data=subset(ensGene.rt.start.sclc.inconsist, RATIO < 0), col="steelblue1", pch=1)
abline(lm(MEDIAN ~ log10(LENGTH), data=subset(ensGene.rt.start.sclc.inconsist, RATIO > 0)), col="sandybrown")
abline(lm(MEDIAN ~ log10(LENGTH), data=subset(ensGene.rt.start.sclc.inconsist, RATIO < 0)), col="steelblue1")
legend("topright", c("Left-leading", "Right-leading"), col=c("steelblue1", "sandybrown"), pch=1, cex=0.8, horiz=F)
dev.off()

writeTable(rownames(ensGene.rt.start.sclc.inconsist), file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_inconsistent.list")), colnames=F, rownames=F, sep="")
writeTable(rownames(subset(ensGene.rt.start.sclc.inconsist, RATIO > 0)), file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_inconsistent_right.list")), colnames=F, rownames=F, sep="")
writeTable(rownames(subset(ensGene.rt.start.sclc.inconsist, RATIO < 0)), file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_inconsistent_left.list")), colnames=F, rownames=F, sep="")

# -----------------------------------------------------------------------------
# 
# Last Modified: 12/11/18
# -----------------------------------------------------------------------------
ensGene.rt.start.sclc <- cbind(ensGene[rownames(ensGene.rt.start.sclc), 1:7], ensGene.rt.start.sclc)
ensGene.rt.start.sclc$LENGTH <- abs(ensGene.rt.start.sclc$start_position - ensGene.rt.start.sclc$end_position)

right <- subset(ensGene.rt.start.sclc, RIGHT_LEADING > cutoff.upper)$MEDIAN
left  <- subset(ensGene.rt.start.sclc, RIGHT_LEADING < cutoff.lower)$MEDIAN
testW(right, left)
# [1] 0.1335412

##
xlab.text <- "Gene length (log10)"
ylab.text <- "log2(TPM+1)"

file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_length_right.pdf"))
main.text <- c("Right-leading genes in SCLC", paste0("n=", separator(length(sclc.tx.right))))
pdf(file.name, height=6, width=6)
plot(MEDIAN ~ log10(LENGTH), data=subset(ensGene.rt.start.sclc, RIGHT_LEADING > cutoff.upper), ylab=ylab.text, xlab=xlab.text, main=main.text, col="sandybrown")
abline(lm(MEDIAN ~ log10(LENGTH), data=subset(ensGene.rt.start.sclc, RIGHT_LEADING > cutoff.upper)), col="orange")
dev.off()

##
file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_length_left.pdf"))
main.text <- c("Left-leading genes in SCLC", paste0("n=", separator(length(sclc.tx.left))))
pdf(file.name, height=6, width=6)
plot(MEDIAN ~ log10(LENGTH), data=subset(ensGene.rt.start.sclc, RIGHT_LEADING < cutoff.lower), ylab=ylab.text, xlab=xlab.text, main=main.text, col="steelblue1")
abline(lm(MEDIAN ~ log10(LENGTH), data=subset(ensGene.rt.start.sclc, RIGHT_LEADING < cutoff.lower)), col="steelblue2")
dev.off()

# -----------------------------------------------------------------------------
# Transcription vs. replication time (Cell cycle genes; Dominguez et al 2016)
# Last Modified: 08/11/18
# -----------------------------------------------------------------------------
file.name  <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_periodic.png"))
main.text <- "Transcription vs. replication time in SCLC"
xlab.text <- "Right-leading ratio (log2)"
ylab.text <- "log2(TPM+1)"
cols <- rep("grey", nrow(ensGene.rt.start.sclc))

png(file.name, height=6, width=6, units="in", res=300)
plot(  MEDIAN ~ RATIO, data=ensGene.rt.start.sclc, xlab=xlab.text, ylab=ylab.text, main=main.text, col=cols)
points(MEDIAN ~ RATIO, data=ensGene.rt.start.sclc[periodic.G2M,], col="forestgreen")
points(MEDIAN ~ RATIO, data=ensGene.rt.start.sclc[periodic.G1S,], col="orange")
legend("bottom", c("G1/S", "G2/M"), col=c("orange", "forestgreen"), pch=1, cex=1, horiz=T)
mtext("Periodic gene lists (Dominguez 2016)", cex=1, line=0.5)
dev.off()

