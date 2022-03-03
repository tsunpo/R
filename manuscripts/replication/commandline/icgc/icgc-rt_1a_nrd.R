#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
BASE <- args[1]     ## Cancer type
LIST <- args[2]     ## Sample list
METHOD <- args[3]   ## RPKB, RPKM, MEAN
base   <- tolower(BASE)
method <- tolower(METHOD) 
 
# =============================================================================
# Name: 1a_cmd-rt_rd.nrd.R (commandline mode)
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 02/08/19; 13/06/17
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"         ## tyang2@cheops
#wd.src <- "/re/home/tyang2/dev/R"                ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"             ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "ReplicationTiming.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.bed.gc.icgc.RData"))

# -----------------------------------------------------------------------------
# To calculate normalised read counts (NRD)
# Last Modified: 02/08/19; 01/05/17
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"   ## tyang2@cheops
wd.ngs <- file.path(wd, "ICGC", BASE, "ngs/WGS")
wd.ngs.data <- file.path(wd.ngs, "data")

samples <- readTable(file.path(wd.ngs, LIST), header=T, rownames=2, sep="")

## Read coverage from PeifLyne
for (s in 1:length(samples$wgs_id)) {
   wgs_id      <- samples$wgs_id[s]
   specimen_id <- samples$specimen_id[s]
   
   ## INPUT: *_cn.txt.gz (read counts)
   rd <- read.peiflyne.icgc.cn.txt(file.path("/projects/cangen/PCAWG-Repository/PCAWG.raw.data/copy_number/converted_data", paste0(wgs_id, "_CONVERTED"), paste0(wgs_id, "_cn.txt")))   ## See ReplicationTiming.R / Inner Class / PeifLyne File Reader
   nrd.T <- initNRD(rd, bed.gc, pair="T", method=METHOD)   ## See ReplicationTiming.R
   #nrd.N <- initNRD(rd, bed.gc, pair="N", method=METHOD)
   rm(rd)
   
   writeTable(nrd.T, gzfile(file.path(wd.ngs.data, specimen_id, paste0(specimen_id, "_T.", method, ".txt.gz"))), colnames=T, rownames=F, sep="\t")
   #writeTable(nrd.N, gzfile(file.path(wd.ngs.data, sample, paste0(sample, "_N.", method, ".txt.gz"))), colnames=T, rownames=F, sep="\t")
   rm(nrd.T)
}
