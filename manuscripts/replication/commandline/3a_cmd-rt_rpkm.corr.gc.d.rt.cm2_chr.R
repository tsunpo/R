#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
BASE <- args[1]   ## Cancer type
PAIR <- args[2]   ## T(umour)
LIST <- args[3]
CHR  <- as.numeric(args[4])
N    <- as.numeric(args[5])
base <- tolower(BASE)

# =============================================================================
# Name: 3a_cmd-rt_rpkm.corr.gc.d.cm2.rt_chr.R (commandline mode)
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 22/05/19
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"         ## tyang2@cheops
#wd.src <- "/re/home/tyang2/dev/R"                ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"             ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "ReplicationTiming.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.1kb.gc.RData"))

# -----------------------------------------------------------------------------
# Replication timing
# Last Modified: 22/05/19
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"   ## tyang2@cheops
wd.ngs    <- file.path(wd, BASE, "ngs/WGS")
wd.ngs.data <- file.path(wd.ngs, "data")

wd.anlys   <- file.path(wd, BASE, "analysis")
wd.rt      <- file.path(wd.anlys, "replication", paste0(base, "-cl-rt"))
wd.rt.data <- file.path(wd.rt, "data")

samples <- readTable(file.path(wd.ngs, LIST), header=T, rownames=T, sep="")

#for (c in 1:22) {
   chr <- chrs[CHR]

   ## Replication timing
   rpkms.chr.d <- readTable(file.path(wd.rt.data, paste0(base, "_rpkm.corr.gc.d_", chr, "_", PAIR, "_n", N, ".txt.gz")), header=T, rownames=T, sep="\t")
   colnames(rpkms.chr.d) <- toupper(gsub("\\.", "", colnames(rpkms.chr.d)))
   rpkms.chr.d <- rpkms.chr.d[, rownames(samples)]
   
   rpkms.T.chr.d <- rpkms.chr.d[, which(samples[,chr] == 1)]
   rpkms.N.chr.d <- rpkms.chr.d[, which(samples[,chr] == 0)]
   n1 <- ncol(rpkms.T.chr.d)
   n0 <- ncol(rpkms.N.chr.d)
   
   rpkms.T.chr.d$MEDIAN <- mapply(x = 1:nrow(rpkms.T.chr.d), function(x) median(as.numeric(rpkms.T.chr.d[x,])))
   rpkms.N.chr.d$MEDIAN <- mapply(x = 1:nrow(rpkms.N.chr.d), function(x) median(as.numeric(rpkms.N.chr.d[x,])))
   
   rpkms.chr <- toTable(NA, 3, nrow(rpkms.T.chr.d), c("T", "N", "RT"))
   rownames(rpkms.chr) <- rownames(rpkms.T.chr.d)
   rpkms.chr$T  <- rpkms.T.chr.d$MEDIAN
   rpkms.chr$N  <- rpkms.N.chr.d$MEDIAN
   rpkms.chr$RT <- mapply(x = 1:nrow(rpkms.chr), function(x) getLog2RDRatio(rpkms.chr[x,]$T, rpkms.chr[x,]$N, pseudocount=0.01))
   
   writeTable(outputRT(rpkms.chr), gzfile(file.path(wd.rt.data, paste0(base, "_rpkm.corr.gc.d.rt.cm2_", chr, "_", PAIR, "-", PAIR, "_n", n1, "-", n0, ".txt.gz"))), colnames=T, rownames=F, sep="\t")
#}
