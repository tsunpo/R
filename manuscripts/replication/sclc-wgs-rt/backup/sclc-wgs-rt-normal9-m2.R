# =============================================================================
# Manuscript   : 
# Chapter      : Chromosome replication timing of the human genome
# Name         : manuscripts/replication/sclc-wgs-rt-m2.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 09/08/19; 24/04/19; 05/03/19; 25/02/19; 30/01/18
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "ReplicationTiming.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.bed.gc.1kb.RData"))
load(file.path(wd.src.ref, "hg19.lcl.koren.woodfine.RData"))

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 30/01/18
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
#wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE  <- "SCLC"
PAIR1 <- "N"
PAIR0 <- "N"
base <- tolower(BASE)
method <- "rpkm"

wd.ngs   <- file.path(wd, BASE, "ngs/WGS")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.rt       <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt-normal9"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

samples1 <- readTable(file.path(wd.ngs, "sclc_wgs_n9_WB.txt"), header=T, rownames=T, sep="")
samples1 <- subset(samples1, M2 == 1)[,1]
samples0 <- readTable(file.path(wd.ngs, "sclc_wgs_n9_WB.txt"), header=T, rownames=T, sep="")
samples0 <- subset(samples0, M2 == 0)[,1]
n1 <- length(samples1)
n0 <- length(samples0)

# -----------------------------------------------------------------------------
# Plot RD and RT (see ReplicationTiming.R)
# Last Modified: 09/08/19; 14/02/19; 10/01/19; 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
nrds <- getLog2ScaledRT(wd.rt.data, base, method, PAIR1, PAIR0, n1, n0, chrs, bed.gc)
save(nrds, file=file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d.rt.log2s_", "m2-m1", ".RData")))
#load(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d.rt.log2s_", "m2-m1", ".RData")))
# nrow(nrds)
# [1] 2658679
nrds.sclc.m2 <- nrds

ymax <- 0.6
ymin <- 0.14
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]
   lcl.rt.chr <- subset(lcl.rt, CHR == chr)   ## Koren 2012
 
   ## Plot RT
   main.text <- paste0(BASE, " (WB) M2/M1 read depth ratio between normal (n=", n1, ") and normal (n=", n0, ") whole bloods")  
   file.name <- file.path(wd.rt.plots, paste0("RT_", base, "_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ""))   
   plotRT(file.name, main.text, chr, NA, NA, nrds.chr, bed.gc.chr, c("red", "blue"), c("M2 normal", "M1 normal"), c("lightcoral", "lightskyblue3"), c("M2", "M1"), "png", width=10, peaks=c(), ylim=c(ymin, ymax), lcl.rt.chr)
}

# -----------------------------------------------------------------------------
# RD vs RT (RDS and SPR)
# Last Modified: 11/07/19; 31/05/19
# -----------------------------------------------------------------------------
sprs <- getSPR(nrds, bed.gc)
save(sprs, file=file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-m2-m1_spline_spearman.RData")))
writeTable(sprs, file=file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-m2-m1_spline_spearman.txt")), colnames=T, rownames=F, sep="\t")
#load(file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-m2-m1_spline_spearman.RData")))

for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   
   nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]
   nrds.chr.T  <- setSpline(nrds.chr, bed.gc.chr, "T")
   nrds.chr.N  <- setSpline(nrds.chr, bed.gc.chr, "N")
   nrds.chr.RT <- setSpline(nrds.chr, bed.gc.chr, "RT")

   main.text <- c(paste0("SCLC (WB) read depth correlation (", "Chr", c, ")"), paste0("rho = ", round0(sprs$cor[c], digits=2), " (M2 vs. M1)"))
   xlab.text <- "SCLC (WB) M2/M1"
   ylab.text <- "SCLC (WB) read depth [RPKM]"
   file.name <- file.path(wd.rt.plots, "chrs", paste0("RD-vs-RT_SCLC-M2-M1_chr", c, "_spline_spearman"))
   plotRD2vsRT(nrds.chr.T$SPLINE, nrds.chr.N$SPLINE, nrds.chr.RT$SPLINE, file.name, main.text, ylab.text, xlab.text, c("red", "blue"), c("M2", "M1"), method="spearman")
}

## S-phase progression rate (SPR)
ylab.text <- "SPR"
file.name <- file.path(wd.rt.plots, "SPR_SCLC-M2-M1_spline_spearman")
main.text <- c(paste0(BASE, " (WB) M2/M1 S-phase progression rate"), "SPR = (E-L)/(E+L)")
plotSPR(sprs, file.name, main.text, c(13, 17), digits=3, unit=5, ylab.text)

## SPR vs Read depth correlation
file.name <- file.path(wd.rt.plots, "SPR-RDC_SCLC-M2-M1_spline_spearman")
main.text <- c(paste0(BASE, " (WB) M2/M1 SPR vs. Read depths correlation"), "")
xlab.text <- "M2 vs. M1 [rho]"
plotSPRRDC(sprs$spr, sprs$cor, file.name, main.text, c(4, 13, 17, 19, 22), xlab.text, unit=5, ylab.text)

## SPR vs Woodfine 2004
file.name <- file.path(wd.rt.plots, "SPR-Woodfine_SCLC-M2-M1_spline_spearman")
main.text <- c(paste0(BASE, " (WB) M2/M1 SPR vs. Woodfine 2004"), "Mean replication timing ratio")
xlab.text <- "Woodfine et al. 2004"
plotSPRRDC(sprs$spr, lcl.mean$Mean, file.name, main.text, c(13, 17, 19, 22), xlab.text, unit=5, ylab.text)

# -----------------------------------------------------------------------------
# RT vs LCL S/G1
# Last Modified: 09/08/19; 06/06/19
# -----------------------------------------------------------------------------
nrds.tmp <- nrds
load(file.path(wd, "LCL/analysis/replication/lcl-wgs-rt/data/lcl_rpkm.gc.cn.d.rt.log2s_s-g1.RData"))
nrds.lcl <- nrds
nrds <- nrds.tmp

cors <- getRTvsRT(nrds, nrds.lcl, bed.gc)
save(cors, file=file.path(wd.rt.data, paste0("rt-vs-rt_", base, "-m2-m1-vs-lcl-s-g1_spline_spearman.RData")))
#load(file.path(wd.rt.data, paste0("rt-vs-rt_", base, "-m2-m1-vs-lcl-s-g1_spline_spearman.RData")))

ylab.text <- "Spearman's rho"
xlab.text <- "Chromosome"
file.name <- file.path(wd.rt.plots, "RT-vs-RT_SCLC-M2-M1-vs-LCL-S-G1_spline_spearman_test")
main.text <- paste0("SCLC (WB) M2/M1 vs. LCL S/G1")
ymin <- 0.25
ymax <- 1.05
plotRTvsRTALL(cors, file.name, main.text, ylab.text, xlab.text, ymin, ymax, col="black", c=2, pos=1)

##
file.name <- file.path(wd.rt.plots, "RTD-vs-RT_SCLC-M2-M1-vs-LCL-S-G1_spline_spearman")
ymin <- -1.1
ymax <- 1.1
plotRD3vsRTALL(cors, file.name, main.text, ymin, ymax, cols=c("red", "blue", "black"), c("M2", "M1", "M2/M1"), c=NA, isRT=T)
