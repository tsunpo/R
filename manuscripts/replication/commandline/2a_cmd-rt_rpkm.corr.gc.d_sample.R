#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
BASE1 <- args[1]   ## Cancer type
PAIR1 <- args[2]   ## T(umour) or S (phase)
LIST1 <- args[3]
BASE0 <- args[4]   ## Normal type
PAIR0 <- args[5]   ## N(ormal) or B(lood) or G1 (phase)
LIST0 <- args[6]
SAMPLE <- args[7]
SAMPLE <- gsub("-", ".", SAMPLE)
base1 <- tolower(BASE1)
base0 <- tolower(BASE0)

# =============================================================================
# Name: 2a_cmd-rt_rpkm.corr.gc.d_sample.R (commandline mode)
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 22/10/18
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"         ## tyang2@cheops
#wd.src <- "/re/home/tyang2/dev/R"                ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"             ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "DifferentialExpression.R", "ReplicationTiming.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.1kb.gc.RData"))

# -----------------------------------------------------------------------------
# Replication timing
# Last Modified: 21/04/19; 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"   ## tyang2@cheops
wd1.ngs    <- file.path(wd, BASE1, "ngs/WGS")
wd1.ngs.data <- file.path(wd1.ngs, "data") 
wd0.ngs    <- file.path(wd, BASE0, "ngs/WGS")
wd0.ngs.data <- file.path(wd0.ngs, "data")

wd1.anlys  <- file.path(wd, BASE1, "analysis")
wd1.rt     <- file.path(wd1.anlys, "replication", paste0(base1, "-wgs-rt"))
wd1.rt.data   <- file.path(wd1.rt, "data")

samples1 <- readTable(file.path(wd1.ngs, LIST1), header=F, rownames=F, sep="")
samples1 <- gsub("-", ".", samples1)   ## ADD 03/07/17 for LCL (e.g. NA19240-2 to NA19240.2)
#samples0 <- readTable(file.path(wd0.ngs, LIST0), header=F, rownames=F, sep="")
#samples0 <- gsub("-", ".", samples0)   ## ADD 03/07/17 for LCL (e.g. NA19240.2 to NA19240.2)
n1 <- length(samples1)
#n0 <- length(samples0)
#if (BSTRP != 0) {
#   samples1 <- samples1[sort(sample(1:n1, n1, replace=T))]
#   samples0 <- samples0[sort(sample(1:n0, n0, replace=T))]
#}

rpkms.T.chr.d.sample.T.all <- NULL
rpkms.chr.rt.lcl.RT.all <- NULL
cors <- toTable(0, 2, 22, c("chr", "cor"))
cors$chr <- 1:22
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   
   ## Read depth
   rpkms.T.chr.d <- pipeGetDetectedRD(wd1.ngs.data, BASE1, chr, PAIR1)[, c("BED", samples1)]   ## ADD BACK 09/04/19; REMOVED 15/02/19; if length(samples) == 1

   ## Individual read depth
   rpkms.T.chr.d.sample <- rpkms.T.chr.d[,c("BED", SAMPLE)]
   colnames(rpkms.T.chr.d.sample) <- c("BED", "T")
   rpkms.T.chr.d.sample$T <- log2(rpkms.T.chr.d.sample$T + 0.01)
   rpkms.T.chr.d.sample.T <- setSpline(rpkms.T.chr.d.sample, bed.gc.chr, "T")

   ## Replication timing
   rpkms.chr.rt.lcl <-readTable(paste0("/projects/cangen/tyang2/LCL/analysis/replication/lcl-wgs-rt/data/lcl_rpkm.corr.gc.d.rt.lcl_", chr, "_LCL-LCL_n7-7.txt.gz"), header=T, rownames=T, sep="\t")
   #rpkms.chr.rt.lcl <-readTable(paste0("/projects/cangen/tyang2/SCLC/analysis/replication/sclc-wgs-rt/data/sclc_rpkm.corr.gc.d.rt_", chr, "_SCLC-SCLC_n101-92.txt.gz"), header=T, rownames=T, sep="\t")
   rpkms.chr.rt.lcl <- setScaledRT(rpkms.chr.rt.lcl, pseudocount=0.01, recaliRT=T, scaledRT=T) 
   rpkms.chr.rt.lcl.RT <- setSpline(rpkms.chr.rt.lcl, bed.gc.chr, "RT")
   
   ## Keep 1kb slopes based on overlapping windows
   overlaps <- intersect(rpkms.T.chr.d.sample.T$BED, rpkms.chr.rt.lcl.RT$BED)
   cors$cor[c] <- getCor(rpkms.chr.rt.lcl.RT[overlaps,]$SPLINE, rpkms.T.chr.d.sample.T[overlaps,]$SPLINE, method="spearman")
   
   ##
   if (is.null(rpkms.T.chr.d.sample.T.all)) {
      rpkms.T.chr.d.sample.T.all <- rpkms.T.chr.d.sample.T[overlaps,]
   } else {
      rpkms.T.chr.d.sample.T.all <- rbind(rpkms.T.chr.d.sample.T.all, rpkms.T.chr.d.sample.T[overlaps,])
   }
   
   if (is.null(rpkms.chr.rt.lcl.RT.all)) {
      rpkms.chr.rt.lcl.RT.all <- rpkms.chr.rt.lcl.RT[overlaps,]
   } else {
      rpkms.chr.rt.lcl.RT.all <- rbind(rpkms.chr.rt.lcl.RT.all, rpkms.chr.rt.lcl.RT[overlaps,])
   }
}
cor <- getCor(rpkms.chr.rt.lcl.RT.all$SPLINE, rpkms.T.chr.d.sample.T.all$SPLINE, method="spearman")
save(cor, cors, file=file.path(wd1.rt.data, "samples", paste0("rd-vs-rt_", SAMPLE, "-vs-lcl_spearman_spline.RData")))
