#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
BASE  <- args[1]   ## Cancer type
PAIR1 <- args[2]   ## T(umour) or S (phase)
LIST1 <- args[3]
PAIR0 <- args[4]   ## N(ormal) or B(lood) or G1 (phase)
LIST0 <- args[5]
#CHR   <- as.numeric(args[6])
METHOD <- args[6]
#RT     <- args[7]  ## NA or RT
BSTRP  <- as.numeric(args[7])
base   <- tolower(BASE)
method <- tolower(METHOD)

# =============================================================================
# Name: cmd-rt_3a_nrd.gc.cn.d.rt_bstrp.R (adopted from cmd-rt_2a_nrd.gc.cn.d_chr.R)
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 12/11/19; 23/04/19; 22/10/18
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"         ## tyang2@cheops
#wd.src <- "/re/home/tyang2/dev/R"                ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"             ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "ReplicationTiming.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.bed.gc.1kb.RData"))

# -----------------------------------------------------------------------------
# Replication timing
# Last Modified: 23/04/19; 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"   ## tyang2@cheops
wd.ngs    <- file.path(wd, BASE, "ngs/WGS")
wd.ngs.data <- file.path(wd.ngs, "data")

wd.anlys <- file.path(wd, BASE, "analysis")
wd.rt    <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data <- file.path(wd.rt, "data")
if (BSTRP != 0)
   wd.rt.data <- file.path(wd.rt, "data/bstrps", BSTRP)

samples1 <- readTable(file.path(wd.ngs, LIST1), header=F, rownames=F, sep="")
samples0 <- readTable(file.path(wd.ngs, LIST0), header=F, rownames=F, sep="")
samples1 <- gsub("-", ".", samples1)
samples0 <- gsub("-", ".", samples0)
n1 <- length(samples1)
n0 <- length(samples0)
if (BSTRP != 0) {
   samples1 <- samples1[sort(sample(1:n1, n1, replace=T))]
   samples0 <- samples0[sort(sample(1:n0, n0, replace=T))]
}

for (c in 1:22) {
   #chr <- chrs[CHR]
   chr <- chrs[c]
   
   ## Read depth
   nrds.T.chr.d <- pipeGetDetectedRD(wd.ngs.data, BASE, chr, PAIR1, method)[, c("BED", samples1)]   ## ADD BACK 09/04/19; REMOVED 15/02/19; if length(samples) == 1
   nrds.N.chr.d <- pipeGetDetectedRD(wd.ngs.data, BASE, chr, PAIR0, method)[, c("BED", samples0)]
   #writeTable(nrds.T.chr.d[, c("BED", samples1)], gzfile(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d_", chr, "_", PAIR1, "_n", n1, ".txt.gz"))), colnames=T, rownames=F, sep="\t")
   #writeTable(nrds.N.chr.d[, c("BED", samples0)], gzfile(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d_", chr, "_", PAIR0, "_n", n0, ".txt.gz"))), colnames=T, rownames=F, sep="\t")
   
   ## Replication timing
   nrds.T.chr.d$MEDIAN <- mapply(x = 1:nrow(nrds.T.chr.d), function(x) median(as.numeric(nrds.T.chr.d[x, -1])))   ## ADD 15/02/19; To skip the first column "BED"
   nrds.N.chr.d$MEDIAN <- mapply(x = 1:nrow(nrds.N.chr.d), function(x) median(as.numeric(nrds.N.chr.d[x, -1])))   ## ADD 15/02/19; To skip the first column "BED"

   overlaps <- intersect(rownames(nrds.T.chr.d), rownames(nrds.N.chr.d))
   nrds.T.chr.d.rt <- nrds.T.chr.d[overlaps,]
   nrds.N.chr.d.rt <- nrds.N.chr.d[overlaps,]
   
   nrds.chr <- toTable(NA, 3, length(overlaps), c("T", "N", "RT"))
   rownames(nrds.chr) <- overlaps
   nrds.chr$T  <- nrds.T.chr.d.rt$MEDIAN
   nrds.chr$N  <- nrds.N.chr.d.rt$MEDIAN
   nrds.chr$RT <- mapply(x = 1:nrow(nrds.chr), function(x) (nrds.chr[x,]$T / nrds.chr[x,]$N))

   writeTable(outputRT(nrds.chr), gzfile(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ".txt.gz"))), colnames=T, rownames=F, sep="\t")
}

# -----------------------------------------------------------------------------
# Adapted from sclc-wgs-rt-m2.R
# -----------------------------------------------------------------------------
nrds <- getLog2ScaledRT(wd.rt.data, base, method, PAIR1, PAIR0, n1, n0, chrs, bed.gc)

nrds.RT <- toTable(0, 4, 0, c("BED", "RT", "SPLINE", "SLOPE"))
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
 
   ## Read in replicaiton time
   overlaps <- intersect(rownames(bed.gc.chr), rownames(nrds))   ## Changed 29/11/19: From intersect(rownames(nrds), rownames(bed.gc.chr))
   nrds.chr <- nrds[overlaps,]
   nrds.chr.RT <- setSpline(nrds.chr, bed.gc.chr, "RT", removeCentromere=F, returnAll=F)   ## Changed 29/11/19

   nrds.RT <- rbind(nrds.RT, nrds.chr.RT)
}
save(nrds.RT, file=file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d.rt.RT.RData")))
