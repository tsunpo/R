#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
BASE   <- args[1]   ## Cancer type

# =============================================================================
# Manuscript   : 
# Chapter II   : 
# Name         : manuscripts/replicaiton/cll-wgs-rt-rfd.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 11/09/19; 12/11/18
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/re/home/tyang2/dev/R"                ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "ReplicationForkDirectionality.R", "ReplicationTiming.R", "Transcription.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.bed.gc.icgc.RData"))
#load(file.path(wd.src.ref, "hg19.ensGene.bed.1kb.RData"))

# -----------------------------------------------------------------------------
# 
# Last Modified: 31/10/18
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
#wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
#BASE <- ""
base <- tolower(BASE)
method <- "rpkm"

wd.ngs   <- file.path(wd, "ICGC", BASE, "ngs/WGS")
wd.anlys <- file.path(wd, "ICGC",BASE, "analysis")

wd.rt    <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data/resampling")
#wd.rt.plots <- file.path(wd.rt, "plots/resampling")
wd.rt.plots <- file.path(wd.rt, "plots/nrfd")

# -----------------------------------------------------------------------------
# Bootstrap distribution
# Last Modified: 02/11/18
# -----------------------------------------------------------------------------
nrds.RT.BSTRPS <- getBootstrap(base, "SLOPE")
save(nrds.RT.BSTRPS, file=file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.m.rt.RT.SLOPE.RData")))
# > nrow(nrds.RT.BSTRPS)
# [1]

##
#nrds.RT.BSTRPS$RFD <- NA
#nrds.RT.BSTRPS$RFD <- getRFD(nrds.RT.BSTRPS)
#
#file.name <- file.path(wd.rt.plots, paste0("hist_", base, "_rpkm_SLOPE_RFD>0.9.pdf"))
#main.text <- c(paste0(BASE, " bootstrap distribution"), paste0(""))
#xlab.text <- "Number of rightward forks per kb window"
#plotBootstrapHist(nrds.RT.BSTRPS, file.name, main.text, xlab.text, 100, boundary.break)

# -----------------------------------------------------------------------------
# Create RT + RFD data
# Last Modified: 22/10/1
# -----------------------------------------------------------------------------
#load(file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.d.rt.RT.SLOPE.RData")))
load(file.path(wd.anlys, "replication", paste0(base, "-wgs-rt-m2"), "data", paste0(base, "_", method, ".gc.cn.m.rt.log2s_", "m2-m1", ".RData")))
# nrow(nrds)
# [1] 

#install.packages("zoo", method="wget")
library("zoo")
kb <- 20
nrds.RT.NRFD <- getRTNRFD(nrds, nrds.RT.BSTRPS, bed.gc, kb)

save(nrds.RT.NRFD, file=file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.m.rt.log2s.nrfd.", kb, "kb_", "m2-m1", ".RData")))
writeTable(nrds.RT.NRFD, gzfile(file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.m.rt.log2s.nrfd.", kb, "kb_", "m2-m1", ".txt.gz"))), colnames=T, rownames=T, sep="\t")
#nrds.RT.NRFD.clle.es <- nrds.RT.NRFD
# > nrow(nrds.RT.NRFD.clle.es)
# [1] 
