#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
BASE <- args[1]   ## Cancer type
LIST <- args[2]   ## Sample list
base <- tolower(BASE)

# =============================================================================
# Name: 1a_cmd-rt_bam.rpkm.R (Commandline mode)
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 13/06/17
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
# Calculate normalised read counts (RPKM)
# Last Modified: 01/05/17
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"   ## tyang2@cheops
wd.ngs    <- file.path(wd, BASE, "ngs/WGS")
wd.ngs.data <- file.path(wd.ngs, "data")

samples <- readTable(file.path(wd.ngs, LIST), header=F, rownames=F, sep="")

## Read coverage from PeifLyne
for (s in 1:length(samples)) {
   sample <- samples[s]
 
   ## INPUT: *_cn.txt.gz (read counts)
   bam <- read.peiflyne.cn.txt(file.path(wd.ngs, "coverage", sample, paste0(sample, "_ANALYSIS/", sample, "_cn.txt.gz")))   ## See ReplicationTiming.R / Inner Class / PeifLyne File Reader
   rpkm.T <- initRPKM(bam, bed.gc, pair="T")   ## See ReplicationTiming.R
   rpkm.N <- initRPKM(bam, bed.gc, pair="N")
 
   writeTable(rpkm.T, gzfile(paste0(wd.ngs.data, sample, "/", sample, "_T.bam.rpkm.txt.gz")), colnames=T, rownames=F, sep="\t")
   writeTable(rpkm.N, gzfile(paste0(wd.ngs.data, sample, "/", sample, "_N.bam.rpkm.txt.gz")), colnames=T, rownames=F, sep="\t")
}
