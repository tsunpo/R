#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
BASE   <- args[1]   ## Cancer type
BSTRPS <- as.numeric(args[2])
COLUMN <- args[3]
#CHR    <- as.numeric(args[4])
base   <- tolower(BASE)
method <- "rpkm"

# =============================================================================
# Name: cmd-rt_3b_nrd.gc.cn.d.rt.RT_bstrps.R (commandline mode)
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 20/09/19; 01/11/18
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"         ## tyang2@cheops
#wd.src <- "/re/home/tyang2/dev/R"                ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"             ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "ReplicationTiming.R", "Bootstrapping.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.1kb.gc.RData"))

# -----------------------------------------------------------------------------
# Read in bootstrapped data (in 3a)
# Last Modified: 01/11/18
# -----------------------------------------------------------------------------
wd <- file.path("/projects/cangen/tyang2", BASE)   ## tyang2@cheops
wd.anlys   <- file.path(wd, "analysis")
wd.rt      <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data <- file.path(wd.rt, "data/bstrps")

load(file.path(wd.rt.data, 1, paste0(base, "_rpkm.gc.cn.d.rt.RT.RData")))
nrds.RT.BSTRPS <- getBSTRPS(nrds.RT, COLUMN, 1)
for (b in 2:BSTRPS) {
   load(file.path(wd.rt.data, b, paste0(base, "_rpkm.gc.cn.d.rt.RT.RData")))
   nrds.RT.BSTRP <- getBSTRPS(nrds.RT, COLUMN, b)
  
   #overlaps <- intersect(rownames(nrds.RT.BSTRPS), rownames(nrds.RT.BSTRP))
   #nrds.RT.BSTRPS <- cbind(nrds.RT.BSTRPS[overlaps,], nrds.RT.BSTRP[overlaps, -1])
   nrds.RT.BSTRPS <- cbind(nrds.RT.BSTRPS, nrds.RT.BSTRP[, -1])
   colnames(nrds.RT.BSTRPS)[1+b] <- colnames(nrds.RT.BSTRP)[2]
}
nrds.RT.BSTRPS <- nrds.RT.BSTRPS[,-1]
save(nrds.RT.BSTRPS, file=file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d.rt.RT.SPLINE.", COLUMN, ".RData")))
