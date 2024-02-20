#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
BASE <- args[1]   ## Cancer type
PAIR <- args[2]   ## T(umour) or N(ormal)
TXT  <- args[3]
#CHR <- as.numeric(args[4])
N    <- as.numeric(args[4])
METHOD   <- args[5]
INSILICO <- args[6]
BSTRP  <- as.numeric(args[7])
#CV     <- args[8]
base   <- tolower(BASE)
method <- tolower(METHOD)
#cv     <- tolower(CV)

# =============================================================================
# Name: 3a_cmd-rt_nrd.gc.cn.d_bstrp_insilico.R (adapted from 2c_cmd-rt_nrd.gc.cn.d_chr_insilico.R)
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 18/09/19; 23/04/19; 22/10/18
# =============================================================================
wd.src <- "/nfs/users/nfs_t/ty2/dev/R"         ## tyang2@cheops
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
wd <- "/lustre/scratch127/casm/team294rr/ty2"   ## tyang2@cheops
wd.ngs    <- file.path(wd, BASE, "ngs/WGS")
wd.ngs.data <- file.path(wd.ngs, "data")

wd.anlys  <- file.path(wd, BASE, "analysis")
wd.rt     <- file.path(wd.anlys, "replication", paste0(base, "-nl-wgs-rt"))
wd.rt.data   <- file.path(wd.rt, "data")
if (BSTRP != 0)
   wd.rt.data <- file.path(wd.rt, "data/bstrps", BSTRP)
#if (CV != "NA")
#   wd.rt.data <- file.path(wd.rt, "data", cv, BSTRP)

samples  <- readTable(file.path(wd.ngs, TXT), header=T, rownames=T, sep="")
samples1 <- samples[which(samples[, INSILICO] == 2),][,1]
samples0 <- samples[which(samples[, INSILICO] == 1),][,1]
samples1 <- gsub("-", ".", samples1)
samples0 <- gsub("-", ".", samples0)
n1 <- length(samples1)
n0 <- length(samples0)
if (BSTRP != 0) {
   samples1 <- samples1[sort(sample(1:n1, n1, replace=T))]
   samples0 <- samples0[sort(sample(1:n0, n0, replace=T))]
}

# -----------------------------------------------------------------------------
# Adapted from sclc-wgs-rt-m2.R
# -----------------------------------------------------------------------------
nrds <- getLog2ScaledRT(wd.rt.data, base, method, PAIR, PAIR, n1, n0, chrs, bed.gc)

nrds.RT <- toTable(0, 4, 0, c("BED", "RT", "SPLINE", "SLOPE"))
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
 
   ## Read in replicaiton time
   overlaps <- intersect(rownames(bed.gc.chr), rownames(nrds))   ## ## Changed 29/11/19: From intersect(rownames(nrds), rownames(bed.gc.chr))
   nrds.chr <- nrds[overlaps,]
   nrds.chr.RT <- setSpline2(nrds.chr, bed.gc.chr, "RT")   ## Changed 29/11/19

   nrds.RT <- rbind(nrds.RT, nrds.chr.RT)
}
save(nrds.RT, file=file.path(wd.rt.data, paste0(base, "_", method, ".cn.m.rt.RData")))

