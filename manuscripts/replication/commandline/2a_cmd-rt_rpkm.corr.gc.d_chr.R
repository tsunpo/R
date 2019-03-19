#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
BASE1 <- args[1]   ## Cancer type
PAIR1 <- args[2]   ## T(umour) or S (phase)
LIST1 <- args[3]
BASE0 <- args[4]   ## Normal type
PAIR0 <- args[5]   ## N(ormal) or B(lood) or G1 (phase)
LIST0 <- args[6]
#BSTRP <- as.numeric(args[7])
CHR  <- as.numeric(args[7])
N <- as.numeric(args[8])
base1 <- tolower(BASE1)
base0 <- tolower(BASE0)

# =============================================================================
# Name: 2a_cmd-rt_rpkm.corr.gc.d_bstrap.R (commandline mode)
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 22/10/18
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"         ## tyang2@cheops
#wd.src <- "/re/home/tyang2/dev/R"                ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"             ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Common.R", "DifferentialExpression.R", "ReplicationTiming.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.1kb.gc.RData"))

# -----------------------------------------------------------------------------
# Replication timing
# Last Modified: 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"   ## tyang2@cheops
wd1.ngs    <- file.path(wd, BASE1, "ngs/WGS")
wd1.ngs.data <- file.path(wd1.ngs, "data") 
wd0.ngs    <- file.path(wd, BASE0, "ngs/WGS")
wd0.ngs.data <- file.path(wd0.ngs, "data")

wd1.anlys  <- file.path(wd, BASE1, "analysis")
wd1.rt     <- file.path(wd1.anlys, "replication", paste0(base1, "-wgs-rt"))
wd1.rt.data   <- file.path(wd1.rt, "data")
wd0.anlys  <- file.path(wd, BASE0, "analysis")
wd0.rt     <- file.path(wd0.anlys, "replication", paste0(base1, "-wgs-rt"))
wd0.rt.data   <- file.path(wd0.rt, "data")
#if (BSTRP != 0)
#   wd1.rt.data <- file.path(wd1.rt, "data/bstrps", BSTRP)

samples1 <- readTable(file.path(wd1.ngs, LIST1), header=F, rownames=F, sep="")
samples1 <- gsub("-", ".", samples1)   ## ADD 03/07/17 for LCL (e.g. NA19240-2 to NA19240.2)
samples0 <- readTable(file.path(wd0.ngs, LIST0), header=F, rownames=F, sep="")
samples0 <- gsub("-", ".", samples0)   ## ADD 03/07/17 for LCL (e.g. NA19240.2 to NA19240.2)
n1 <- length(samples1)
n0 <- length(samples0)

#for (c in 1:22) {
   chr <- chrs[CHR]

   ## Read depth
   #rpkms.T.chr.d <- pipeGetDetectedRD(wd1.ngs.data, BASE1, chr, PAIR1)
   #rpkms.N.chr.d <- pipeGetDetectedRD(wd0.ngs.data, BASE0, chr, PAIR0)   
   #rpkms.T.chr.d$MEDIAN <- mapply(x = 1:nrow(rpkms.T.chr.d), function(x) median(as.numeric(rpkms.T.chr.d[x, -1])))   ## ADD 15/02/19; To skip the first column "BED"
   #rpkms.N.chr.d$MEDIAN <- mapply(x = 1:nrow(rpkms.N.chr.d), function(x) median(as.numeric(rpkms.N.chr.d[x, -1])))   ## ADD 15/02/19; To skip the first column "BED"
   #rpkms.T.chr.d$BED <- rownames(rpkms.T.chr.d)   ## REMOVED 15/02/19; Already keep the first column "BED" in previous step; ## CHANGE 23/10/18: After mapply(median()) otherwise warnings due to $BED
   #rpkms.N.chr.d$BED <- rownames(rpkms.N.chr.d)
   
   #writeTable(rpkms.T.chr.d[, c("BED", "MEDIAN")], gzfile(file.path(wd1.rt.data, paste0(base1, "_rpkm.corr.gc.d_", chr, "_", PAIR1, "_n", n1, ".txt.gz"))), colnames=T, rownames=F, sep="\t")
   #writeTable(rpkms.N.chr.d[, c("BED", "MEDIAN")], gzfile(file.path(wd1.rt.data, paste0(base0, "_rpkm.corr.gc.d_", chr, "_", PAIR0, "_n", n0, ".txt.gz"))), colnames=T, rownames=F, sep="\t")

   ## Replication timing
   if (N != 0) {
      rpkms.T.chr.d <- readTable(file.path(wd1.rt.data, paste0(base1, "_rpkm.corr.gc.d_", chr, "_", PAIR1, "_n", N, ".txt.gz")), header=T, rownames=T, sep="\t")
      rpkms.T.chr.d <- rpkms.T.chr.d[,c("BED", samples1)]
      rpkms.N.chr.d <- readTable(file.path(wd0.rt.data, paste0(base0, "_rpkm.corr.gc.d_", chr, "_", PAIR0, "_n", N, ".txt.gz")), header=T, rownames=T, sep="\t")
      rpkms.N.chr.d <- rpkms.N.chr.d[,c("BED", samples0)]
   } else {
      rpkms.T.chr.d <- readTable(file.path(wd1.rt.data, paste0(base1, "_rpkm.corr.gc.d_", chr, "_", PAIR1, "_n", n1, ".txt.gz")), header=T, rownames=T, sep="\t")
      rpkms.N.chr.d <- readTable(file.path(wd0.rt.data, paste0(base0, "_rpkm.corr.gc.d_", chr, "_", PAIR0, "_n", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   }
   rpkms.T.chr.d$MEDIAN <- mapply(x = 1:nrow(rpkms.T.chr.d), function(x) median(as.numeric(rpkms.T.chr.d[x, -1])))   ## ADD 15/02/19; To skip the first column "BED"
   rpkms.N.chr.d$MEDIAN <- mapply(x = 1:nrow(rpkms.N.chr.d), function(x) median(as.numeric(rpkms.N.chr.d[x, -1])))   ## ADD 15/02/19; To skip the first column "BED"
   
   overlaps <- intersect(rownames(rpkms.T.chr.d), rownames(rpkms.N.chr.d))
   rpkms.T.chr.d.rt <- rpkms.T.chr.d[overlaps,]
   rpkms.N.chr.d.rt <- rpkms.N.chr.d[overlaps,]
   
   rpkms.chr <- toTable(NA, 3, length(overlaps), c("T", "N", "RT"))
   rownames(rpkms.chr) <- overlaps
   rpkms.chr$T  <- rpkms.T.chr.d.rt$MEDIAN
   rpkms.chr$N  <- rpkms.N.chr.d.rt$MEDIAN
   rpkms.chr$RT <- mapply(x = 1:nrow(rpkms.chr), function(x) getLog2RDRatio(rpkms.chr[x,]$T, rpkms.chr[x,]$N, pseudocount=0.01))   ## CHANGED 11/02/19: Was pseudocount=0.01
   
   #plotRD(wd.rt.plots, samples[s], "_rpkm.corr.gc.d_", chr, PAIR, rpkms.chr[,s], bed.gc.chr, "png")
   writeTable(outputRT(rpkms.chr), gzfile(file.path(wd1.rt.data, paste0(base1, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ".txt.gz"))), colnames=T, rownames=F, sep="\t")
#}

