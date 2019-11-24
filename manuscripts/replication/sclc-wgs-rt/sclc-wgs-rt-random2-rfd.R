# =============================================================================
# Manuscript   : 
# Chapter II   : 
# Name         : manuscripts/replicaiton/sclc-wgs-rt-rfd.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 11/09/19; 12/11/18
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/re/home/tyang2/dev/R"                ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "ReplicationForkDirectionality.R", "ReplicationTiming.R", "DifferentialExpression.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.bed.gc.1kb.RData"))
load(file.path(wd.src.ref, "hg19.ensGene.bed.1kb.RData"))

# -----------------------------------------------------------------------------
# 
# Last Modified: 31/10/18
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
#wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE <- "SCLC"
PAIR1 <- "M2"
PAIR0 <- "M1"
base <- tolower(BASE)
method <- "rpkm"
bstrps        <- 1000
boundary.upper <- 520   ## 500-520 breaks
boundary.lower <- 480   ## 480-500 breaks
boundary.break <- 2     ## 1 breaks each centering 500
n1 <- 25
n0 <- 26

wd.ngs   <- file.path(wd, BASE, "ngs/WGS")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.rt    <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data/random2")
#wd.rt.plots <- file.path(wd.rt, "plots/bstrps")

# -----------------------------------------------------------------------------
# Bootstrap distribution
# Last Modified: 02/11/18
# -----------------------------------------------------------------------------
nrds.RT.BSTRPS <- getBootstrap(base, "SLOPE")
save(nrds.RT.BSTRPS, file=file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.d.rt.RT.SLOPE.RData")))
# > nrow(nrds.RT.BSTRPS)
# [1] 2650061

# -----------------------------------------------------------------------------
# Create RT + RFD data
# Last Modified: 22/10/1
# -----------------------------------------------------------------------------
#load(file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.d.rt.RT.SLOPE.RData")))
load(file.path(wd.anlys, "replication", paste0(base, "-wgs-rt-m2"), "data", paste0(base, "_", method, ".gc.cn.d.rt.log2s_", "m2-m1", ".RData")))
nrds.RT <- getRT(nrds, bed.gc)
# > nrow(nrds.RT)
# [1] 2650083

nrds.RT.RFD <- getRTRFD(nrds.RT, nrds.RT.BSTRPS)
save(nrds.RT.RFD, file=file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.d.rt.log2s.rfd_", "m2-m1", ".RData")))
writeTable(nrds.RT.RFD, gzfile(file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.d.rt.log2s.rfd_", "m2-m1", ".txt.gz"))), colnames=T, rownames=T, sep="\t")
nrds.RT.RFD.2 <- nrds.RT.RFD
# > nrow(nrds.RT.RFD.2)
# [1] 2650061
