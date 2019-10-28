#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
BASE  <- "NBL"
PAIR1 <- "T"
#PAIR0 <- "N"
TXT   <- "nbl_wgs_n57-1.txt"
METHOD <- "RPKM"
base   <- tolower(BASE)
method <- tolower(METHOD)
n <- 56

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
overlaps.c.e <- readTable(file.path(wd.src.ref, "overlaps.c.e_n40451.txt.gz"), header=F, rownames=F, sep="\t")

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

samples  <- readTable(file.path(wd.ngs, TXT), header=T, rownames=T, sep="")
samples1 <- samples[which(samples[, "M2"] == 1),][,1]
samples0 <- samples[which(samples[, "M2"] == 0),][,1]
samples1 <- gsub("-", ".", samples1)
samples0 <- gsub("-", ".", samples0)
n1 <- length(samples1)
n0 <- length(samples0)

nrds.T.d <- NULL
nrds.N.d <- NULL
for (c in 1:22) {
   chr <- chrs[c]
 
   nrds.T.chr.d <- readTable(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d_", chr, "_", PAIR1, "_n", n, ".txt.gz")), header=T, rownames=F, sep="\t")
   nrds.N.chr.d <- readTable(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d_", chr, "_", PAIR1, "_n", n, ".txt.gz")), header=T, rownames=F, sep="\t")
   rownames(nrds.T.chr.d) <- nrds.T.chr.d$BED
   rownames(nrds.N.chr.d) <- nrds.N.chr.d$BED
   nrds.T.chr.d <- nrds.T.chr.d[intersect(overlaps.c.e, rownames(nrds.T.chr.d)), samples1]
   nrds.N.chr.d <- nrds.N.chr.d[intersect(overlaps.c.e, rownames(nrds.N.chr.d)), samples0]
   
   if (is.null(nrds.T.d)) {
      nrds.T.d <- nrds.T.chr.d
      nrds.N.d <- nrds.N.chr.d
   } else {
      nrds.T.d <- rbind(nrds.T.d, nrds.T.chr.d)
      nrds.N.d <- rbind(nrds.N.d, nrds.N.chr.d)
   }
}
nrds.T.d <- nrds.T.d[overlaps.c.e,]
nrds.N.d <- nrds.N.d[overlaps.c.e,]

##
test <- toTable(NA, 2, length(overlaps.c.e), c("BED", "P"))
test$BED <- overlaps.c.e
rownames(test) <- overlaps.c.e

test$P <- mapply(x = 1:nrow(test), function(x) testU(log2(as.numeric(nrds.T.d[x,])), log2(as.numeric(nrds.N.d[x,]))))
library(qvalue)
test$FDR <- qvalue(test$P)$qvalue
test <- test[order(test$P),]

save(test, file=file.path(wd.rt.data, paste0("overlaps.c.e.", base, ".RData")))
