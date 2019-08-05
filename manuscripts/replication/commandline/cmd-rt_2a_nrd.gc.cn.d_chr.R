#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
BASE1 <- args[1]   ## Cancer type
PAIR1 <- args[2]   ## T(umour) or S (phase)
TXT1  <- args[3]
BASE0 <- args[4]   ## Normal type
PAIR0 <- args[5]   ## N(ormal) or B(lood) or G1 (phase)
TXT0  <- args[6]
CHR   <- as.numeric(args[7])
N     <- as.numeric(args[8])
METHOD <- args[9]
base1  <- tolower(BASE1)
base0  <- tolower(BASE0)
method <- tolower(METHOD)

# =============================================================================
# Name: 2a_cmd-rt_rpkm.corr.gc.d_chr.R (commandline mode)
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 23/04/19; 22/10/18
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
# Last Modified: 23/04/19; 31/08/18; 13/06/17
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

samples1 <- readTable(file.path(wd1.ngs, TXT1), header=T, rownames=T, sep="")
samples0 <- readTable(file.path(wd0.ngs, TXT0), header=T, rownames=T, sep="")
if (N != 0) {
   samples1 <- subset(samples1, M2 == 1)[,1]
   samples0 <- subset(samples0, M2 == 0)[,1]
} else {
   samples1 <- samples1[,1]
   samples0 <- samples0[,1]
}
samples1 <- gsub("-", ".", samples1)
samples0 <- gsub("-", ".", samples0)
n1 <- length(samples1)
n0 <- length(samples0)

#for (c in 1:22) {
   chr <- chrs[CHR]

   ## Read depth
   rds.T.chr.d <- pipeGetDetectedRD(wd1.ngs.data, BASE1, chr, PAIR1, method)[, c("BED", samples1)]   ## ADD BACK 09/04/19; REMOVED 15/02/19; if length(samples) == 1
   rds.N.chr.d <- pipeGetDetectedRD(wd0.ngs.data, BASE0, chr, PAIR0, method)[, c("BED", samples0)]
   writeTable(rds.T.chr.d[, c("BED", samples1)], gzfile(file.path(wd1.rt.data, paste0(base1, "_", method, ".gc.cn.d_", chr, "_", PAIR1, "_n", n1, ".txt.gz"))), colnames=T, rownames=F, sep="\t")
   writeTable(rds.N.chr.d[, c("BED", samples0)], gzfile(file.path(wd1.rt.data, paste0(base0, "_", method, ".gc.cn.d_", chr, "_", PAIR0, "_n", n0, ".txt.gz"))), colnames=T, rownames=F, sep="\t")
   
   ## Replication timing
   rds.T.chr.d <- readTable(file.path(wd1.rt.data, paste0(base1, "_", method, ".gc.cn.d_", chr, "_", PAIR1, "_n", n1, ".txt.gz")), header=T, rownames=T, sep="\t")
   rds.N.chr.d <- readTable(file.path(wd0.rt.data, paste0(base0, "_", method, ".gc.cn.d_", chr, "_", PAIR0, "_n", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   if (N != 0) {
      rds.T.chr.d <- rds.T.chr.d[,c("BED", samples1)]
      rds.N.chr.d <- rds.N.chr.d[,c("BED", samples0)]
   }
   rds.T.chr.d$MEDIAN <- mapply(x = 1:nrow(rds.T.chr.d), function(x) median(as.numeric(rds.T.chr.d[x, -1])))   ## ADD 15/02/19; To skip the first column "BED"
   rds.N.chr.d$MEDIAN <- mapply(x = 1:nrow(rds.N.chr.d), function(x) median(as.numeric(rds.N.chr.d[x, -1])))   ## ADD 15/02/19; To skip the first column "BED"

   overlaps <- intersect(rownames(rds.T.chr.d), rownames(rds.N.chr.d))
   rds.T.chr.d.rt <- rds.T.chr.d[overlaps,]
   rds.N.chr.d.rt <- rds.N.chr.d[overlaps,]
   
   rds.chr <- toTable(NA, 3, length(overlaps), c("T", "N", "RT"))
   rownames(rds.chr) <- overlaps
   rds.chr$T  <- rds.T.chr.d.rt$MEDIAN
   rds.chr$N  <- rds.N.chr.d.rt$MEDIAN
   rds.chr$RT <- mapply(x = 1:nrow(rds.chr), function(x) getLog2RDRatio(rds.chr[x,]$T, rds.chr[x,]$N, pseudocount=0.01))   ## CHANGED 11/02/19: Was pseudocount=0.01
   
   #plotRD(wd.rt.plots, samples[s], "_rpkm.corr.gc.d_", chr, PAIR, rpkms.chr[,s], bed.gc.chr, "png")
   writeTable(outputRT(rds.chr), gzfile(file.path(wd1.rt.data, paste0(base1, "_", method, ".gc.cn.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ".txt.gz"))), colnames=T, rownames=F, sep="\t")
#}

