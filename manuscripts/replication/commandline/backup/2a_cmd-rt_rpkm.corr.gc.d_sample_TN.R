#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
BASE1 <- args[1]   ## Cancer type
PAIR1 <- args[2]   ## T(umour) or S (phase)
LIST1 <- args[3]
BASE0 <- args[4]   ## Normal type
PAIR0 <- args[5]   ## N(ormal) or B(lood) or G1 (phase)
LIST0 <- args[6]
SAMPLE <- args[7]
base1 <- tolower(BASE1)
base0 <- tolower(BASE0)

# =============================================================================
# Name: 2a_cmd-rt_rpkm.corr.gc.d_sample.R (commandline mode)
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 22/10/18
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"         ## tyang2@cheops
#wd.src <- "/re/home/tyang2/dev/R"                ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"             ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "DifferentialExpression.R", "ReplicationTiming.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.1kb.gc.RData"))

# -----------------------------------------------------------------------------
# Replication timing
# Last Modified: 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"   ## tyang2@cheops
wd1.ngs    <- file.path(wd, BASE1, "ngs/WGS")
wd1.ngs.data <- file.path(wd1.ngs, "data") 
wd0.ngs    <- file.path(wd, BASE0, "ngs/WGS")
wd0.ngs.data <- file.path(wd0.ngs, "data")

wd1.anlys  <- file.path(wd, BASE1, "analysis")
wd1.rt     <- file.path(wd1.anlys, "replication", paste0(base1, "-wgs-rt"))
wd1.rt.data   <- file.path(wd1.rt, "data")

samples1 <- readTable(file.path(wd1.ngs, LIST1), header=F, rownames=F, sep="")
samples1 <- gsub("-", ".", samples1)   ## ADD 03/07/17 for LCL (e.g. NA19240-2 to NA19240.2)
samples0 <- readTable(file.path(wd0.ngs, LIST0), header=F, rownames=F, sep="")
samples0 <- gsub("-", ".", samples0)   ## ADD 03/07/17 for LCL (e.g. NA19240.2 to NA19240.2)
n1 <- length(samples1)
n0 <- length(samples0)
#if (BSTRP != 0) {
#   samples1 <- samples1[sort(sample(1:n1, n1, replace=T))]
#   samples0 <- samples0[sort(sample(1:n0, n0, replace=T))]
#}

cors <- toTable(0, 2, 22, c("chr", "cor"))
cors$chr <- 1:22
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   
   ## Read depth
   #rpkms.T.chr.d <- readTable(file.path(wd1.rt.data, paste0(base1, "_rpkm.corr.gc.d_", chr, "_", PAIR1, "_n", n1, ".txt.gz")), header=T, rownames=T, sep="\t")
   #rpkms.N.chr.d <- readTable(file.path(wd1.rt.data, paste0(base0, "_rpkm.corr.gc.d_", chr, "_", PAIR0, "_n", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.T.chr.d <- pipeGetDetectedRD(wd1.ngs.data, BASE1, chr, PAIR1)[, c("BED", samples1)]   ## ADD BACK 09/04/19; REMOVED 15/02/19; if length(samples) == 1
   rpkms.N.chr.d <- pipeGetDetectedRD(wd0.ngs.data, BASE0, chr, PAIR0)[, c("BED", samples0)]   
   overlaps <- intersect(rpkms.T.chr.d$BED, rpkms.N.chr.d$BED)
   
   ## Individual replication timing
   rpkms.chr.rt.sample <- toTable(0, 3, length(overlaps), c("BED", "T", "N"))
   rownames(rpkms.chr.rt.sample) <- overlaps
   rpkms.chr.rt.sample$BED <- overlaps
   rpkms.chr.rt.sample$T <- rpkms.T.chr.d[overlaps, SAMPLE]
   rpkms.chr.rt.sample$N <- rpkms.N.chr.d[overlaps, SAMPLE]
   rpkms.chr.rt.sample <- setScaledRT(rpkms.chr.rt.sample, pseudocount=0.01, recaliRT=T, flipRT=F, scaledRT=T)   ## Be careful!! flipRT=T in CLL
   rpkms.chr.rt.sample.RT <- setSlope(rpkms.chr.rt.sample, bed.gc.chr, "RT")
   
   ## Replication timing
   rpkms.chr.rt <- readTable(file.path(wd1.rt.data, paste0(base1, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, flipRT=F, scaledRT=T)
   rpkms.chr.rt.RT <- setSlope(rpkms.chr.rt, bed.gc.chr, "RT")
   
   #rpkms.chr.rt.lcl <-readTable(paste0("/projects/cangen/tyang2/LCL/analysis/replication/lcl-wgs-rt-lcl/data/lcl_rpkm.corr.gc.d.rt.lcl_", chr, "_LCL-LCL_n7-7.txt.gz"), header=T, rownames=T, sep="\t")
   rpkms.chr.rt.lcl <-readTable(paste0("/projects/cangen/tyang2/SCLC/analysis/replication/sclc-wgs-rt/data/sclc_rpkm.corr.gc.d.rt_", chr, "_SCLC-SCLC_n101-92.txt.gz"), header=T, rownames=T, sep="\t")
   rpkms.chr.rt.lcl <- setScaledRT(rpkms.chr.rt.lcl, pseudocount=0.01, recaliRT=T, flipRT=F, scaledRT=T) 
   rpkms.chr.rt.lcl.RT <- setSlope(rpkms.chr.rt.lcl, bed.gc.chr, "RT")
   
   ## Keep 1kb slopes based on overlapping windows
   overlaps <- intersect(rpkms.chr.rt.RT$BED, rpkms.chr.rt.lcl.RT$BED)
   overlaps <- intersect(overlaps, rpkms.chr.rt.sample.RT$BED)
   #rpkms.chr.rt.RT     <- rpkms.chr.rt.RT[overlaps,]
   #rpkms.chr.rt.lcl.RT <- rpkms.chr.rt.lcl.RT[overlaps,]
   #rpkms.chr.rt.sample.RT <- rpkms.chr.rt.sample.RT[overlaps,]
   
   cors$cor[c] <- getCor(rpkms.chr.rt.lcl.RT[overlaps,]$SLOPE, rpkms.chr.rt.sample.RT[overlaps,]$SLOPE)
}
save(cors, file=file.path(wd1.rt.data, "samples", paste0("rt-vs-rt_", SAMPLE, "-vs-sclc_cors-pearson.RData")))
