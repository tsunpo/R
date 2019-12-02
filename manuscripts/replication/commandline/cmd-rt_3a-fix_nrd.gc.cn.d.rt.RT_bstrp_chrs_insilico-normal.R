#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
BASE <- args[1]   ## Cancer type
PAIR <- args[2]   ## T(umour) or N(ormal)
TXT  <- args[3]
#CHR <- as.numeric(args[4])
N    <- as.numeric(args[4])
METHOD   <- args[5]
INSILICO <- args[6]
BSTRP  <- as.numeric(args[7])
CV     <- args[8]
base   <- tolower(BASE)
method <- tolower(METHOD)
cv     <- tolower(CV)

# =============================================================================
# Name: 3a_cmd-rt_nrd.gc.cn.d_bstrp_insilico.R (adapted from 2c_cmd-rt_nrd.gc.cn.d_chr_insilico.R)
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 18/09/19; 23/04/19; 22/10/18
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
# Adapted from sclc-wgs-rt-m2.R
# -----------------------------------------------------------------------------
nrds <- getLog2ScaledRT(wd.rt.data, base, method, PAIR, PAIR, n1, n0, chrs, bed.gc)

nrds.RT <- toTable(0, 4, 0, c("BED", "RT", "SPLINE"))
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
 
   ## Read in replicaiton time
   overlaps <- intersect(rownames(bed.gc.chr), rownames(nrds))   ## Changed 29/11/19: From intersect(rownames(nrds), rownames(bed.gc.chr))
   nrds.chr <- nrds[overlaps,]
   nrds.chr.RT <- setSpline(nrds.chr, bed.gc.chr, "RT", removeCentromere=F, returnAll=F)   ## Changed 29/11/19

   nrds.RT <- rbind(nrds.RT, nrds.chr.RT)
}
save(nrds.RT, file=file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d.rt.RT.RData")))

