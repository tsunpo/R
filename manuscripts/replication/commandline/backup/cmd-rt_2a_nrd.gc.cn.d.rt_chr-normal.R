#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
BASE  <- args[1]   ## Cancer type
PAIR0 <- args[2]   ## N(ormal) or B(lood) or G1 (phase)
LIST0 <- args[3]
CHR   <- as.numeric(args[4])
METHOD <- args[5]
base   <- tolower(BASE)
method <- tolower(METHOD)

# =============================================================================
# Name: cmd-rt_2a_nrd.gc.cn.d_chr.R (commandline mode)
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

samples0 <- readTable(file.path(wd.ngs, LIST0), header=F, rownames=F, sep="")
samples0 <- gsub("-", ".", samples0)
n0 <- length(samples0)

#for (c in 1:22) {
   chr <- chrs[CHR]

   ## Read depth
   nrds.N.chr.d <- pipeGetDetectedRD(wd.ngs.data, BASE, chr, PAIR0, method)[, c("BED", samples0)]
   writeTable(nrds.N.chr.d[, c("BED", samples0)], gzfile(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d_", chr, "_", PAIR0, "_n", n0, ".txt.gz"))), colnames=T, rownames=F, sep="\t")
#}
