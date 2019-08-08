#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
BASE <- args[1]   ## Cancer type
PAIR <- args[2]   ## T(umour) or S (phase)
SAMPLE <- args[3]
#SAMPLE <- toupper(gsub("-", "", SAMPLE))   ## ADDED for NBL CCL
METHOD <- args[4]
base   <- tolower(BASE)
method <- tolower(METHOD)

# =============================================================================
# Name: 2a_cmd-rt_rpkm.corr.gc.d_sample.R (commandline mode)
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 22/10/18
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"         ## tyang2@cheops
#wd.src <- "/re/home/tyang2/dev/R"                ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"             ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "ReplicationTiming.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.1kb.gc.RData"))

# -----------------------------------------------------------------------------
# Replication timing
# Last Modified: 21/04/19; 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"   ## tyang2@cheops
wd.ngs    <- file.path(wd, BASE, "ngs/WGS")
wd.ngs.data <- file.path(wd.ngs, "data") 

wd.anlys <- file.path(wd, BASE, "analysis")
wd.rt    <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data <- file.path(wd.rt, "data")

#samples <- readTable(file.path(wd1.ngs, LIST), header=F, rownames=F, sep="")
#samples <- toupper(gsub("-", "", samples))   ## ADD 03/07/17 for LCL (e.g. NA19240-2 to NA19240.2)   ## ADDED for NBL CCL
#n <- length(samples)

## nrds
load(file.path(wd, "LCL/analysis/replication/lcl-wgs-rt/data", paste0("nrds_lcl-s-g1_", method, ".RData")))

nrds.T.chr.d.sample.T.all <- NULL
nrds.chr.RT.all <- NULL
cors <- toTable(0, 2, 22, c("chr", "cor"))
cors$chr <- 1:22
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   
   ## Read depth
   nrds.T.chr.d <- pipeGetDetectedRD(wd.ngs.data, BASE, chr, PAIR, method)
   colnames(nrds.T.chr.d) <- toupper(gsub("\\.", "", colnames(nrds.T.chr.d)))
   #nrds.T.chr.d <- nrds.T.chr.d[, c("BED", samples1)]   ## ADD BACK 09/04/19; REMOVED 15/02/19; if length(samples) == 1

   ## Individual read depth
   nrds.T.chr.d.sample <- nrds.T.chr.d[,c("BED", SAMPLE)]
   colnames(nrds.T.chr.d.sample) <- c("BED", "T")
   #nrds.T.chr.d.sample$T <- log2(nrds.T.chr.d.sample$T + 0.01)
   nrds.T.chr.d.sample.T <- setSpline(nrds.T.chr.d.sample, bed.gc.chr, "T")

   ## Replication timing
   nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]  ## Reference LCL S/G1 ratio
   nrds.chr.RT <- setSpline(nrds.chr, bed.gc.chr, "RT")
   
   ## Keep 1kb slopes based on overlapping windows
   overlaps <- intersect(nrds.T.chr.d.sample.T$BED, nrds.chr.RT$BED)
   cors$cor[c] <- getCor(nrds.chr.RT[overlaps,]$SPLINE, nrds.T.chr.d.sample.T[overlaps,]$SPLINE, method="spearman")
   
   ##
   if (is.null(nrds.T.chr.d.sample.T.all)) {
      nrds.T.chr.d.sample.T.all <- nrds.T.chr.d.sample.T[overlaps,]
   } else {
      nrds.T.chr.d.sample.T.all <- rbind(nrds.T.chr.d.sample.T.all, nrds.T.chr.d.sample.T[overlaps,])
   }
   
   if (is.null(nrds.chr.RT.all)) {
      nrds.chr.RT.all <- nrds.chr.RT[overlaps,]
   } else {
      nrds.chr.RT.all <- rbind(nrds.chr.RT.all, nrds.chr.RT[overlaps,])
   }
}
cor <- getCor(nrds.chr.RT.all$SPLINE, nrds.T.chr.d.sample.T.all$SPLINE, method="spearman")
save(cor, cors, file=file.path(wd.rt.data, "samples", paste0("rd-vs-rt_", SAMPLE, "-vs-lcl_spline_spearman.RData")))
