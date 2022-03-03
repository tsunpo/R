# =============================================================================
# Manuscript   : 
# Chapter      : 
# Name         : manuscripts/replication/icgc/clle-es-wgs-rt.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 26/02/22; 14/10/19; 26/02/19
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "ReplicationTiming.R", "Transcription.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.bed.gc.icgc.RData"))

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 30/01/18
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
#wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE  <- "ESAD-UK"
PAIR1 <- "T"
PAIR0 <- "N"
base  <- tolower(BASE)
method <- "rpkm"

wd.ngs   <- file.path(wd, "ICGC", BASE, "ngs/WGS")
wd.anlys <- file.path(wd, "ICGC", BASE, "analysis")

wd.rt       <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

wd.ngs.data <- file.path(wd.ngs, "data")
samples1 <- readTable(file.path(wd.ngs, "esad-uk_wgs_n88.list"), header=F, rownames=F, sep="")$V3
samples0 <- readTable(file.path(wd.ngs, "esad-uk_wgs_n88.list"), header=F, rownames=F, sep="")$V3
n1 <- length(samples1)
n0 <- length(samples0)

# -----------------------------------------------------------------------------
# CLL T vs LCL S/G1
# http://pklab.med.harvard.edu/scw2014/subpop_tutorial.html
# Last Modified: 06/03/19
# -----------------------------------------------------------------------------
cors.samples <- getSAMPLEvsRT(wd.rt.data, samples1)
save(cors.samples, file=file.path(wd.rt.data, paste0("samples-vs-rt_", base, "-vs-lcl_spline_spearman.RData")))
# > min(cors.samples[,-c(1:4)])
# [1] -0.8424011
# > max(cors.samples[,-c(1:4)])
# [1] 0.7085663

# -----------------------------------------------------------------------------
# Overall correlation with LCL S/G1
# Last Modified: 19/11/19; 16/06/19; 04/06/19; 06/03/19
# -----------------------------------------------------------------------------
samples.esad.uk <- setSamplesQ4(wd.rt.data, samples1)
writeTable(samples.esad.uk, file.path(wd.ngs, "esad-uk_wgs_m2_n88.txt"), colnames=T, rownames=F, sep="\t")
#         0%        25%        50%        75%       100%
# -0.7472742 -0.7202768 -0.6393568 -0.3887980  0.2835392

writeTable(subset(samples.esad.uk, Q4 %in% c(4,1)), file.path(wd.ngs, "esad-uk_wgs_q4_n44.txt"), colnames=T, rownames=F, sep="\t")
