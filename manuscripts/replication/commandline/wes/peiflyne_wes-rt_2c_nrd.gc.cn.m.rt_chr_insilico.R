#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
BASE <- args[1]   ## Cancer type
PAIR <- args[2]   ## T(umour) or N(ormal)
TXT  <- args[3]
CHR  <- as.numeric(args[4])
N    <- as.numeric(args[5])
METHOD   <- args[6]
INSILICO <- args[7]
#BSTRP  <- as.numeric(args[8])
base   <- tolower(BASE)
method <- tolower(METHOD)

# =============================================================================
# Name: 2c_cmd-rt_nrd.gc.cn.d.rt_chr_insilico.R (commandline mode)
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 30/11/20; 23/04/19; 22/10/18
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
# Replication timing
# Last Modified: 23/04/19; 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"   ## tyang2@cheops
wd.ngs    <- file.path(wd, BASE, "ngs/WES_NEW")
wd.ngs.data <- file.path(wd.ngs, "data")

wd.anlys  <- file.path(wd, BASE, "analysis")
wd.rt     <- file.path(wd.anlys, "replication", paste0(base, "-wes-rt"))
wd.rt.data   <- file.path(wd.rt, "data_wes_new")

samples  <- readTable(file.path(wd.ngs, TXT), header=T, rownames=T, sep="")
samples1 <- samples[which(samples[, INSILICO] == 2),][,1]
samples0 <- samples[which(samples[, INSILICO] == 1),][,1]
samples1 <- gsub("\\.", "-", samples1)
samples0 <- gsub("\\.", "-", samples0)
n1 <- length(samples1)
n0 <- length(samples0)

#for (c in 1:22) {
   chr <- chrs[CHR]

   ## Replication timing
   nrds.chr.m <- readTable(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.m_", chr, "_", PAIR, "_n", N, ".txt.gz")), header=T, rownames=T, sep="\t")
   colnames(nrds.chr.m) <- gsub("\\.", "-", colnames(nrds.chr.m))   ## ADD 29/11/20; 05/10/19
   nrds.chr.m <- nrds.chr.m[which(nrds.chr.m$MEDIAN != 0),]     ## ADD 20/02/22; 29/11/20
   
   ## S-like/G1-like ratio
   nrds.T.chr.m <- nrds.chr.m[,c("BED", samples1)]
   nrds.N.chr.m <- nrds.chr.m[,c("BED", samples0)]
   rm(nrds.chr.m)   ## ADD 22/02/22
   nrds.T.chr.m$MEDIAN <- mapply(x = 1:nrow(nrds.T.chr.m), function(x) median(as.numeric(nrds.T.chr.m[x, -1])))   ## ADD 15/02/19; To skip the first column "BED"
   nrds.N.chr.m$MEDIAN <- mapply(x = 1:nrow(nrds.N.chr.m), function(x) median(as.numeric(nrds.N.chr.m[x, -1])))   ## ADD 15/02/19; To skip the first column "BED"
   
   ## Same as cmd-rt_2a_nrd.gc.cn.m.rt_chr.R
   nrds.T.chr.m <- nrds.T.chr.m[which(nrds.T.chr.m$MEDIAN != 0),]   ## ADD 20/02/22; 29/11/20
   nrds.N.chr.m <- nrds.N.chr.m[which(nrds.N.chr.m$MEDIAN != 0),]   ## ADD 20/02/22; 29/11/20
   
   overlaps <- intersect(rownames(nrds.T.chr.m), rownames(nrds.N.chr.m))
   nrds.T.chr.m.o <- nrds.T.chr.m[overlaps,]
   nrds.N.chr.m.o <- nrds.N.chr.m[overlaps,]
   rm(nrds.T.chr.m, nrds.N.chr.m)   ## ADD 22/02/22
   
   nrds.chr <- toTable(NA, 3, nrow(nrds.T.chr.m.o), c("T", "N", "RT"))
   rownames(nrds.chr) <- rownames(nrds.T.chr.m.o)
   nrds.chr$T  <- nrds.T.chr.m.o$MEDIAN
   nrds.chr$N  <- nrds.N.chr.m.o$MEDIAN
   rm(nrds.T.chr.m.o, nrds.N.chr.m.o)   ## ADD 22/02/22
   nrds.chr$RT <- mapply(x = 1:nrow(nrds.chr), function(x) (nrds.chr[x,]$T / nrds.chr[x,]$N))

   writeTable(outputRT(nrds.chr), gzfile(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.m.rt_", chr, "_", PAIR, "-", PAIR, "_n", n1, "-", n0, ".txt.gz"))), colnames=T, rownames=F, sep="\t")
#}
