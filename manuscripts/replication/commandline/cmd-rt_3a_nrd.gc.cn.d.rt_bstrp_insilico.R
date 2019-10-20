#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
BASE <- args[1]   ## Cancer type
PAIR <- args[2]   ## T(umour) or N(ormal)
TXT  <- args[3]
#CHR  <- as.numeric(args[4])
N    <- as.numeric(args[4])
METHOD   <- args[5]
INSILICO <- args[6]
BSTRP  <- as.numeric(args[7])
base   <- tolower(BASE)
method <- tolower(METHOD)

# =============================================================================
# Name: 3a_cmd-rt_nrd.gc.cn.d_bstrp_insilico.R (adapted from 2c_cmd-rt_nrd.gc.cn.d_chr_insilico.R)
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 18/09/19; 23/04/19; 22/10/18
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

wd.anlys  <- file.path(wd, BASE, "analysis")
wd.rt     <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data   <- file.path(wd.rt, "data")
if (BSTRP != 0)
   wd.rt.data <- file.path(wd.rt, "data/bstrps", BSTRP)

samples  <- readTable(file.path(wd.ngs, TXT), header=T, rownames=T, sep="")
samples1 <- samples[which(samples[, INSILICO] == 1),][,1]
samples0 <- samples[which(samples[, INSILICO] == 0),][,1]
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

   ## Replication timing
   nrds.chr.d <- readTable(file.path(wd.rt, "data", paste0(base, "_", method, ".gc.cn.d_", chr, "_", PAIR, "_n", N, ".txt.gz")), header=T, rownames=T, sep="\t")
   colnames(nrds.chr.d) <- gsub("\\.", "", colnames(nrds.chr.d))   ## ADD 05/10/19
   colnames(nrds.chr.d) <- toupper(colnames(nrds.chr.d))           ## ADD 05/10/19
   nrds.T.chr.d <- nrds.chr.d[,c("BED", samples1)]
   nrds.N.chr.d <- nrds.chr.d[,c("BED", samples0)]

   nrds.T.chr.d$MEDIAN <- mapply(x = 1:nrow(nrds.T.chr.d), function(x) median(as.numeric(nrds.T.chr.d[x, -1])))   ## ADD 15/02/19; To skip the first column "BED"
   nrds.N.chr.d$MEDIAN <- mapply(x = 1:nrow(nrds.N.chr.d), function(x) median(as.numeric(nrds.N.chr.d[x, -1])))   ## ADD 15/02/19; To skip the first column "BED"

   nrds.chr <- toTable(NA, 3, nrow(nrds.chr.d), c("T", "N", "RT"))
   rownames(nrds.chr) <- rownames(nrds.chr.d)   ## BUG 18/09/19
   nrds.chr$T  <- nrds.T.chr.d$MEDIAN
   nrds.chr$N  <- nrds.N.chr.d$MEDIAN
   nrds.chr$RT <- mapply(x = 1:nrow(nrds.chr), function(x) (nrds.chr[x,]$T / nrds.chr[x,]$N))
   
   writeTable(outputRT(nrds.chr), gzfile(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d.rt_", chr, "_", PAIR, "-", PAIR, "_n", n1, "-", n0, ".txt.gz"))), colnames=T, rownames=F, sep="\t")
}
