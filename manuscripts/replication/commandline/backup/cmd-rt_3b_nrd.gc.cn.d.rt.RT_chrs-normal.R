#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
BASE  <- args[1]   ## Cancer type
PAIR1 <- args[2]   ## N(ormal) or B(lood) or G1 (phase)
PAIR0 <- args[3]   ## N(ormal) or B(lood) or G1 (phase)
NUM1  <- as.numeric(args[4])
NUM0  <- as.numeric(args[5])
BSTRP <- as.numeric(args[6])
base   <- tolower(BASE)
method <- "rpkm"

# =============================================================================
# Name: cmd-rt_3b_nrd.gc.cs.d.rt.RT_chrs.R (commandline mode)
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 18/09/19; 22/10/18
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
# Last Modified: 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"   ## tyang2@cheops
wd.anlys <- file.path(wd, BASE, "analysis")
wd.rt    <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt-normal"))
wd.rt.data  <- file.path(wd.rt, "data/bstrps", BSTRP)

# -----------------------------------------------------------------------------
# Adapted from sclc-wgs-rt-m2.R
# -----------------------------------------------------------------------------
nrds <- getLog2ScaledRT(wd.rt.data, base, method, PAIR1, PAIR0, NUM1, NUM0, chrs, bed.gc)

nrds.RT <- toTable(0, 4, 0, c("BED", "RT", "SPLINE", "SLOPE"))
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)

   ## Read in replicaiton time
   overlaps <- intersect(rownames(nrds), rownames(bed.gc.chr))
   nrds.chr <- nrds[overlaps,]
   nrds.chr.RT <- setSpline(nrds.chr, bed.gc.chr, "RT")
   
   ## Determine replication direction for each expressed gene
   ## https://stackoverflow.com/questions/43615469/how-to-calculate-the-slope-of-a-smoothed-curve-in-r
   overlaps <- intersect(rownames(nrds.chr.RT), rownames(bed.gc.chr))
   bed.gc.chr <- bed.gc.chr[overlaps,]
   
   slopes <- diff(nrds.chr.RT$SPLINE)/diff(bed.gc.chr$START/1E6)   ## ADD 31/10/18
   nrds.chr.RT$SLOPE <- NA
   nrds.chr.RT <- nrds.chr.RT[1:length(slopes),]
   nrds.chr.RT$SLOPE <- slopes   ## length(slopes) is 1 less than nrow(bed.gc.chr.rt), as no slope for the last 1KB window 
   
   nrds.RT <- rbind(nrds.RT, nrds.chr.RT)
}
save(nrds.RT, file=file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d.rt.RT.RData")))
