#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
BASE <- args[1]   ## Cancer type
LIST <- args[2]   ## Samples
CHR  <- as.numeric(args[3]) 
PAIR <- args[4]   ## T(umour) or N(ormal) pair
METHOD <- args[5]
base   <- tolower(BASE)
method <- tolower(METHOD)

# =============================================================================
# Name: 1c_cmd-rt_rd_nrd_corr_chr.R (commandline mode)
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
load(file.path(wd.src.ref, "hg19.bed.gc.peiflyne.wes.new.RData"))

# -----------------------------------------------------------------------------
# Calculate absolute RPKM
# Last Modified: 14/05/17
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"   ## tyang2@cheops
wd.ngs    <- file.path(wd, BASE, "ngs/WES_NEW")
wd.ngs.data <- file.path(wd.ngs, "data")

samples <- readTable(file.path(wd.ngs, LIST), header=F, rownames=F, sep="")

##
#for (c in 1:length(chrs)) {
   chr <- chrs[CHR]
   bed.gc.chr <- rownames(subset(bed.gc, CHR == chr))

   nrds.chr <- initReadDepthPerChromosome(samples, bed.gc.chr) 
   for (s in 1:length(samples)) {
      sample <- samples[s]
    
      nrd <- readTable(file.path(wd.ngs.data, sample, paste0(sample, "_", PAIR, ".", method, ".gc.cn.txt.gz")), header=T, rownames=T, sep="")
      nrd.chr <- nrd[bed.gc.chr, c("BED", "NRD_GC_CN")]
      
      nrds.chr[nrd.chr$BED, s] <- nrd.chr$NRD_GC_CN
   }

   nrds.chr$BED <- rownames(nrds.chr)
   writeTable(nrds.chr[,c("BED", samples)], gzfile(file.path(wd.ngs.data, paste0(tolower(BASE), "_", method, ".gc.cn_", chr, "_", PAIR, ".txt.gz"))), colnames=T, rownames=F, sep="\t")
#}
