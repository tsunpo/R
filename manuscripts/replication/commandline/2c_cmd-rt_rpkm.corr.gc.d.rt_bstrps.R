#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
BASE   <- args[1]   ## Cancer type
BSTRPS <- as.numeric(args[2])
base   <- tolower(BASE)

# =============================================================================
# Name: 2c_cmd-rt_rpkm.corr.gc.d.rt_bstrp1000.R (Commandline mode)
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 31/10/18
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
# Last Modified: 31/10/18
# -----------------------------------------------------------------------------
wd <- file.path("/projects/cangen/tyang2", BASE)   ## tyang2@cheops
wd.anlys   <- file.path(wd, "analysis")
wd.rt      <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data <- file.path(wd.rt, "data/bstrps")

ensGene.rt       <- loadEnsGeneRT(wd.rt.data, 1, base)
ensGene.rt.start <- getEnsGeneRT(ensGene.rt, "SLOPE_START", 1)
ensGene.rt.end   <- getEnsGeneRT(ensGene.rt, "SLOPE_END", 1)
for (b in 2:BSTRPS) {
   ensGene.rt     <- loadEnsGeneRT(wd.rt.data, b, base)
   overlaps.start <- intersect(rownames(ensGene.rt.start), rownames(ensGene.rt))
   overlaps.end   <- intersect(rownames(ensGene.rt.end), rownames(ensGene.rt))
   
   ensGene.rt.start <- cbind(ensGene.rt.start[overlaps.start,], ensGene.rt[overlaps.start, "SLOPE_START"])
   ensGene.rt.end   <- cbind(ensGene.rt.end[overlaps.end,], ensGene.rt[overlaps.end, "SLOPE_END"])
   colnames(ensGene.rt.start)[1+b] <- paste0("SLOPE_START_", b)
   colnames(ensGene.rt.end)[1+b]   <- paste0("SLOPE_END_", b)
}
ensGene.rt.start <- ensGene.rt.start[,-1]
ensGene.rt.end <- ensGene.rt.end[,-1]
#save(ensGene.rt.start, ensGene.rt.end, file=file.path(wd.rt.data, paste0("ensGene.rt_", base, "_bstrps", BSTRPS, ".RData")))

# -----------------------------------------------------------------------------
# Determine replication timing data from bootstrapped data (in 2c and 2d)
# Last Modified: 31/10/18
# -----------------------------------------------------------------------------
#load(file.path(wd.rt, "ensGene.rt_sclc_bstrp1000.RData"))
ensGene.rt.start <- pipeLeading(ensGene.rt.start, BSTRPS)
length(which(ensGene.rt.start$RT == 1))
length(which(ensGene.rt.start$RT == -1))
length(which(ensGene.rt.start$RT == 0))

ensGene.rt.end   <- pipeLeading(ensGene.rt.end, BSTRPS)
length(which(ensGene.rt.end$RT == 1))
length(which(ensGene.rt.end$RT == -1))
length(which(ensGene.rt.end$RT == 0))

save(ensGene.rt.start, ensGene.rt.end, file=file.path(wd.rt.data, paste0("ensGene.rt_", base, "_bstrps", BSTRPS, "_full.RData")))

# -----------------------------------------------------------------------------
# Add BED information for the bootstrapped data (in 2c)
# Last Modified: 09/11/18
# -----------------------------------------------------------------------------
ensGene.rt.start <- setEnsGeneBED(ensGene.rt.start, bed.gc, chrs[1:22], isStartPosition=T)
ensGene.rt.end   <- setEnsGeneBED(ensGene.rt.end,   bed.gc, chrs[1:22], isStartPosition=F)

save(ensGene.rt.start, ensGene.rt.end, file=file.path(wd.rt.data, paste0("ensGene.rt_", base, "_bstrps", BSTRPS, ".RData")))   ## Store new ensGene.rt.start and ensGene.rt.end
