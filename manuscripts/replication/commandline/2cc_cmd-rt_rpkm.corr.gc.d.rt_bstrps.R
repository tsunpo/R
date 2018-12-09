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

load(file.path(wd.rt.data, paste0("ensGene.rt_", base, "_bstrps", BSTRPS, "_full.RData")))

ensGene.rt.start$RFD <- mapply(x = 1:nrow(ensGene.rt.start), function(x) as.numeric(getRFD(ensGene.rt.start[x, ])))
ensGene.rt.end$RFD   <- mapply(x = 1:nrow(ensGene.rt.end),   function(x) as.numeric(getRFD(ensGene.rt.end[x, ])))
save(ensGene.rt.start, ensGene.rt.end, file=file.path(wd.rt.data, paste0("ensGene.rt_", base, "_bstrps", BSTRPS, "_full.RData")))

# -----------------------------------------------------------------------------
# Add BED information for the bootstrapped data (in 2c)
# Last Modified: 09/11/18
# -----------------------------------------------------------------------------
ensGene.rt.start <- setEnsGeneBED(ensGene.rt.start, bed.gc, chrs[1:22], isStartPosition=T)
ensGene.rt.end   <- setEnsGeneBED(ensGene.rt.end,   bed.gc, chrs[1:22], isStartPosition=F)

save(ensGene.rt.start, ensGene.rt.end, file=file.path(wd.rt.data, paste0("ensGene.rt_", base, "_bstrps", BSTRPS, ".RData")))   ## Store new ensGene.rt.start and ensGene.rt.end
