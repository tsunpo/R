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
n1 <- 25
n0 <- 25

wd.ngs   <- file.path(wd, BASE, "ngs/WGS")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.rt    <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data/random1")
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
# nrow(nrds)
# [1] 

#install.packages("zoo", method="wget")
library("zoo")
kb <- 20
nrds.RT.NRFD <- getRTNRFD(nrds, nrds.RT.BSTRPS, bed.gc, kb)

save(nrds.RT.NRFD, file=file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.d.rt.log2s.nrfd.", kb, "kb_", "m2-m1", ".RData")))
writeTable(nrds.RT.NRFD, gzfile(file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.d.rt.log2s.nrfd.", kb, "kb_", "m2-m1", ".txt.gz"))), colnames=T, rownames=T, sep="\t")
nrds.RT.NRFD.sclc.1 <- nrds.RT.NRFD
# nrow(nrds.RT.NRFD.sclc.1)
# [1] 2650061













# -----------------------------------------------------------------------------
# |RFD| ≥ 0.9
# Last Modified: 24/09/19
# -----------------------------------------------------------------------------
boundary.upper <- 950   ## RFD > +0.9
boundary.lower <-  50   ## RFD < -0.9

nrds.RT.RFD.1.t <- getBootstrapTTR(nrds.RT.RFD.1, boundary.lower, boundary.upper)
nrds.RT.RFD.2.t <- getBootstrapTTR(nrds.RT.RFD.2, boundary.lower, boundary.upper)
nrow(nrds.RT.RFD.1.t)
# [1] 2031894
# > 2031894/2650061
# [1] 0.7667348
nrow(nrds.RT.RFD.2.t)
# [1] 1833114
# > 1833114/2650061
# [1] 0.6917252

# -----------------------------------------------------------------------------
# ALL TTR (between RANDOM1 and RANDOM2)
# Last Modified: 21/11/19; 27/10/19; 22/09/19
# -----------------------------------------------------------------------------
overlaps.t <- intersect(rownames(nrds.RT.RFD.1.t), rownames(nrds.RT.RFD.2.t))
length(overlaps.t)
# [1] 1706325

nrds.1.RT.o <- nrds.RT.RFD.1.t[overlaps.t,]
nrds.2.RT.o <- nrds.RT.RFD.2.t[overlaps.t,]

##
nrds.1.RT.o$SIGN <- nrds.1.RT.o$SLOPE * nrds.2.RT.o$SLOPE
length(which(nrds.1.RT.o$SIGN > 0))
# [1] 1833114

# -----------------------------------------------------------------------------
# |RFD| < 0.9
# Last Modified: 19/10/19
# -----------------------------------------------------------------------------
boundary.upper <- 950   ## RFD < +0.9
boundary.lower <-  50   ## RFD > -0.9

nrds.RT.RFD.1.c <- getBootstrapCTR(nrds.RT.RFD.1, boundary.lower, boundary.upper)
nrds.RT.RFD.2.c <- getBootstrapCTR(nrds.RT.RFD.2, boundary.lower, boundary.upper)
nrow(nrds.RT.RFD.1.c)
# [1] 618167
nrow(nrds.RT.RFD.2.c)
# [1] 816947

nrds.RT.RFD.1.c.e <- subset(nrds.RT.RFD.1.c, SPLINE >= 0)
nrds.RT.RFD.1.c.l <- subset(nrds.RT.RFD.1.c, SPLINE < 0)
nrow(nrds.RT.RFD.1.c.e)
# [1] 350743
nrow(nrds.RT.RFD.1.c.l)
# [1] 267424

nrds.RT.RFD.2.c.e <- subset(nrds.RT.RFD.2.c, SPLINE >= 0)
nrds.RT.RFD.2.c.l <- subset(nrds.RT.RFD.2.c, SPLINE < 0)
nrow(nrds.RT.RFD.2.c.e)
# [1] 431965
nrow(nrds.RT.RFD.2.c.l)
# [1] 384982

# -----------------------------------------------------------------------------
# CTR (E)
# Last Modified: 19/10/19
# -----------------------------------------------------------------------------
overlaps.c.e <- intersect(rownames(nrds.RT.RFD.1.c.e), rownames(nrds.RT.RFD.2.c.e))
length(overlaps.c.e)
# [1] 276271

nrds.1.RT.o <- nrds.RT.RFD.1.c.e[overlaps.c.e,]
nrds.2.RT.o <- nrds.RT.RFD.2.c.e[overlaps.c.e,]

##
nrds.1.RT.o$SIGN <- nrds.1.RT.o$SLOPE * nrds.2.RT.o$SLOPE
length(which(nrds.1.RT.o$SIGN > 0))
# [1] 276271

# -----------------------------------------------------------------------------
# CTR (L)
# Last Modified: 19/10/19
# -----------------------------------------------------------------------------
overlaps.c.l <- intersect(rownames(nrds.RT.RFD.1.c.l), rownames(nrds.RT.RFD.2.c.l))
length(overlaps.c.l)
# [1] 215107

nrds.1.RT.o <- nrds.RT.RFD.1.c.l[overlaps.c.l,]
nrds.2.RT.o <- nrds.RT.RFD.2.c.l[overlaps.c.l,]

##
nrds.1.RT.o$SIGN <- nrds.1.RT.o$SLOPE * nrds.2.RT.o$SLOPE
length(which(nrds.1.RT.o$SIGN > 0))
# [1] 215107
