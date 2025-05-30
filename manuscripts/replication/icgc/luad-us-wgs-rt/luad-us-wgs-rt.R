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
BASE  <- "LUAD-US"
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
samples1 <- readTable(file.path(wd.ngs, "luad-us_wgs_n37.list"), header=F, rownames=F, sep="")$V3
samples0 <- readTable(file.path(wd.ngs, "luad-us_wgs_n37.list"), header=F, rownames=F, sep="")$V3
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
# [1] -0.8334016
# > max(cors.samples[,-c(1:4)])
# [1] 0.830849

# -----------------------------------------------------------------------------
# Overall correlation with LCL S/G1
# Last Modified: 19/11/19; 16/06/19; 04/06/19; 06/03/19
# -----------------------------------------------------------------------------
samples.luad.us <- setSamplesQ4(wd.rt.data, samples1)
writeTable(samples.luad.us, file.path(wd.ngs, "luad-us_wgs_m2_n37.txt"), colnames=T, rownames=F, sep="\t")
#         0%        25%        50%        75%       100% 
# -0.7497030 -0.7271070 -0.5994878  0.2939130  0.7477155 

writeTable(subset(samples.luad.us, Q4 %in% c(4,1)), file.path(wd.ngs, "luad-us_wgs_q4_n19.txt"), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# To compare between 1-kb and ICGC partionings
# Last Modified: 04/03/23
# -----------------------------------------------------------------------------
load("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/LUAD-US/analysis/replication/luad-us-wgs-rt/data/samples-vs-rt_luad-us-vs-lcl_spline_spearman.RData")
cors.samples.clle <- cors.samples
cors.samples.clle <- as.data.frame(t(cors.samples.clle)[-c(1:4),])
m.clle <- mapply(x = 1:nrow(cors.samples.clle), function(x) median(as.numeric(cors.samples.clle[x,])))
cors.samples.clle$MEDIAN <- m.clle

load("/Users/tpyang/Work/uni-koeln/tyang2/LUAD/analysis/replication/luad-wgs-rt/data/samples-vs-rt_luad-vs-lcl_spline_spearman.RData")
cors.samples.cll <- cors.samples
cors.samples.cll <- as.data.frame(t(cors.samples.cll)[-c(1:4),])
m.cll <- mapply(x = 1:nrow(cors.samples.cll), function(x) median(as.numeric(cors.samples.cll[x,])))
cors.samples.cll$MEDIAN <- m.cll

overlaps <- intersect(rownames(cors.samples.clle), rownames(cors.samples.cll))
cors.samples.clle.o <- cors.samples.clle[overlaps,]
cors.samples.cll.o  <- cors.samples.cll[overlaps,]

file.name <- "/Users/tpyang/Work/uni-koeln/tyang2/ICGC/LUAD-US/analysis/replication/luad-us-wgs-rt/plots/correlation_luad-us-vs-luad_spearman"
plotCorrelation(file.name, "LUAD-US", "ICGC partioning", "1-kb partioning", cors.samples.clle.o$MEDIAN, cors.samples.cll.o$MEDIAN, "bottomright")
