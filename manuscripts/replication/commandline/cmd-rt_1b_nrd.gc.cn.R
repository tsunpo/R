#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
BASE   <- args[1]   ## Cancer type
SAMPLE <- args[2]   ## Sample ID
PAIR   <- args[3]   ## T(umour) or N(ormal) pair
CN     <- as.logical(args[4])   ## T(rue) or F(alse) to correct CN
METHOD <- args[5]
base   <- tolower(BASE)
method <- tolower(METHOD)

# =============================================================================
# Name: 1b_cmd-rt_bam.rpkm.corr.gc.R (commandline mode)
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 11/06/17
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"         ## tyang2@cheops
#wd.src <- "/re/home/tyang2/dev/R"                ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"             ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "ReplicationTiming.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.1kb.bed.gc.RData"))

# -----------------------------------------------------------------------------
# Calculate absolute RPKM
# Last Modified: 14/05/17
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"   ## tyang2@cheops
wd.ngs    <- file.path(wd, BASE, "ngs/WGS")
wd.ngs.data <- file.path(wd.ngs, "data")
wd.anlys  <- file.path(wd, BASE, "analysis")
wd.rt     <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data   <- file.path(wd.rt, "data")

## Read copy number information from PeifLyne
#for (s in 1:length(samples)) {
   #sample <- samples[s]
   sample <- SAMPLE 

   ## INPUT: RPKM (from 1a) & SEQ
   nrd <- readTable(file.path(wd.ngs.data, sample, paste0(sample, "_", PAIR, ".", method, ".txt.gz")), header=T, rownames=T, sep="")
   segs <- NA
   if (CN)
      segs <- read.peiflyne.cn.seg(file.path(wd.ngs, "coverage", sample, paste0(sample, "_ANALYSIS/", sample, "_cn.seg")))   ## See ReplicationTiming.R / Inner Class / PeifLyne File Reader
   
   ## Correcte read counts fpr copy number
   nrd.gc.cn <- getNRDGCCN(nrd, segs, bed.gc, PAIR, CN)   ## See ReplicationTiming.R
   writeTable(nrd.gc.cn, gzfile(file.path(wd.ngs.data, sample, paste0(sample, "_", PAIR, ".", method, ".gc.cn.txt.gz"))), colnames=T, rownames=F, sep="\t")
#}
