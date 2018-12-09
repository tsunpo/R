#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
BASE   <- args[1]   ## Cancer type
BSTRPS <- as.numeric(args[2])
CHR    <- as.numeric(args[3])
base   <- tolower(BASE)

# =============================================================================
# Name: 2d_cmd-rt_rpkm.corr.gc.d.rt_bstrp1000_bed.R (Commandline mode)
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 01/11/18
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
# Read in bootstrapped data (in 2c and 2d)
# Last Modified: 01/11/18
# -----------------------------------------------------------------------------
wd <- file.path("/projects/cangen/tyang2", BASE)   ## tyang2@cheops
wd.anlys   <- file.path(wd, "analysis")
wd.rt      <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data <- file.path(wd.rt, "data/bstrps")
bed.gc <- bed[which(bed$GC > 0),]   ## Only keep partitions (in the BED file) with valid GC content

bed.gc.rt <- loadBEDGCRT(wd.rt.data, 1, tolower(base))
bed.gc.rt <- getBEDGCRT(bed.gc.rt, "SLOPE", 1)
for (b in 2:BSTRPS) {
   bed.gc.rt.b <- loadBEDGCRT(wd.rt.data, b, tolower(base))
   overlaps <- intersect(rownames(bed.gc.rt), rownames(bed.gc.rt.b))

   bed.gc.rt <- cbind(bed.gc.rt[overlaps,], bed.gc.rt.b[overlaps, "SLOPE"])
   colnames(bed.gc.rt)[1+b] <- paste0("SLOPE_", b)
}
bed.gc.rt <- bed.gc.rt[,-1]
#save(bed.gc.rt, file=file.path(wd.rt.data, paste0("bed.gc.rt_", base, "_bstrps", BSTRPS, ".RData")))

# -----------------------------------------------------------------------------
# Determin replication timing from bootstrapped data (in 2c and 2d)
# Last Modified: 01/11/18
# -----------------------------------------------------------------------------
#load(file.path(wd.rt.data, "bed.gc.rt_sclc_bstrps1000.RData"))
#for (c in 1:22) {
   chr <- chrs[CHR]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   
   overlaps <- intersect(rownames(bed.gc.rt), rownames(bed.gc.chr))
   bed.gc.rt.chr <- bed.gc.rt[overlaps,]
   rm(bed.gc.rt)   ## ADD 01/11/18
   
   bed.gc.rt.chr <- pipeLeading(bed.gc.rt.chr, BSTRPS)
   save(bed.gc.rt.chr, file=file.path(wd.rt.data, paste0("bed.gc.rt_", base, "_bstrps", BSTRPS, "_", chr, "_full.RData")))
   
   bed.gc.rt.chr <- bed.gc.rt.chr[,c("RIGHT_LEADING", "LEFT_LEADING", "RFD")]
   save(bed.gc.rt.chr, file=file.path(wd.rt.data, paste0("bed.gc.rt_", base, "_bstrps", BSTRPS, "_", chr, ".RData")))
#}
