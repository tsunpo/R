#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
BASE   <- args[1]   ## Cancer type
CHR    <- as.numeric(args[2])
BSTRPS <- 1000
#CV     <- args[3]
COLUMN <- "SLOPE"
#CHR    <- as.numeric(args[4])
base   <- tolower(BASE)
method <- "rpkm"
#cv     <- tolower(CV)

# =============================================================================
# Name: cmd-rt_3c_nrd.gc.cn.d.rt.RT_bstrps.R (commandline mode)
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 20/09/19; 01/11/18
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"         ## tyang2@cheops
#wd.src <- "/re/home/tyang2/dev/R"                ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"             ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "ReplicationForkDirectionality.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.bed.gc.icgc.RData"))

# -----------------------------------------------------------------------------
# 3b: Read in bootstrapped data 
# Last Modified: 05/03/22; 01/11/18
# -----------------------------------------------------------------------------
wd <- file.path("/projects/cangen/tyang2", "ICGC", BASE)   ## tyang2@cheops
wd.anlys   <- file.path(wd, "analysis")
wd.rt      <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data <- file.path(wd.rt, "data/resampling")
#if (CV != "NA")
#   wd.rt.data <- file.path(wd.rt, "data", cv)

load(file.path(wd.rt.data, 1, paste0(base, "_", method, ".gc.cn.m.rt.RT_chr", CHR, ".RData")))
nrds.RT.BSTRPS.chr <- getBSTRPS(nrds.RT.chr, COLUMN, 1)
rm(nrds.RT.chr)
for (b in 2:BSTRPS) {
   load(file.path(wd.rt.data, b, paste0(base, "_", method, ".gc.cn.m.rt.RT_chr", CHR, ".RData")))
   nrds.RT.BSTRP.chr <- getBSTRPS(nrds.RT.chr, COLUMN, b)
   rm(nrds.RT.chr)
   
   overlaps <- intersect(rownames(nrds.RT.BSTRPS.chr), rownames(nrds.RT.BSTRP.chr))
   nrds.RT.BSTRPS.chr <- cbind(nrds.RT.BSTRPS.chr[overlaps,], nrds.RT.BSTRP.chr[overlaps, -1])
   #nrds.RT.BSTRPS <- cbind(nrds.RT.BSTRPS, nrds.RT.BSTRP[, -1])
   colnames(nrds.RT.BSTRPS.chr)[1+b] <- colnames(nrds.RT.BSTRP.chr)[2]
   rm(nrds.RT.BSTRP.chr)   ## ADD 27/02/22
}
nrds.RT.BSTRPS.chr <- nrds.RT.BSTRPS.chr[,-1]
#save(nrds.RT.BSTRPS, file=file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.m.rt.RT.", COLUMN, ".BSTRPS.RData")))

# -----------------------------------------------------------------------------
# 3c: Determine replication timing from bootstrapping data 
# Last Modified: 05/03/22; 01/11/18
# ----------------------------------------------------------------------------
nrds.RT.BSTRPS.chr <- pipeBootstrap(nrds.RT.BSTRPS.chr, BSTRPS)

nrds.RT.BSTRPS.chr <- nrds.RT.BSTRPS.chr[,c("POS", "NEG")]
save(nrds.RT.BSTRPS.chr, file=file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.m.rt.RT.", COLUMN, "_chr", CHR, ".RData")))
