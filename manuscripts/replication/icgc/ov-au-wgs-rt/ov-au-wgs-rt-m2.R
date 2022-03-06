# =============================================================================
# Manuscript   : 
# Chapter      : Chromosome replication timing of the human genome
# Name         : manuscripts/replication/cll-wgs-rt-m2.R
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
load(file.path(wd.src.ref, "hg19.bed.gc.icgc.RData"))
#load(file.path(wd.src.ref, "hg19.lcl.koren.woodfine.RData"))

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 30/01/18
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
#wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE  <- "OV-AU"
PAIR1 <- "T"
PAIR0 <- "T"
base <- tolower(BASE)
method <- "rpkm"

wd.ngs   <- file.path(wd, "ICGC", BASE, "ngs/WGS")
wd.anlys <- file.path(wd, "ICGC", BASE, "analysis")

wd.rt       <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt-m2"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

samples1 <- readTable(file.path(wd.ngs, "ov-au_wgs_m2_n70.txt"), header=T, rownames=T, sep="")   ## M2/M1
#samples1 <- readTable(file.path(wd.ngs, ""), header=T, rownames=T, sep="")    ## Q4/Q1
samples1 <- subset(samples1, M2 == 2)[,1]
samples0 <- readTable(file.path(wd.ngs, "ov-au_wgs_m2_n70.txt"), header=T, rownames=T, sep="")
#samples0 <- readTable(file.path(wd.ngs, ""), header=T, rownames=T, sep="")
samples0 <- subset(samples0, M2 == 1)[,1]
n1 <- length(samples1)
n0 <- length(samples0)

# -----------------------------------------------------------------------------
# Plot RD and RT (see ReplicationTiming.R)
# Last Modified: 09/08/19; 14/02/19; 10/01/19; 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
nrds <- getLog2ScaledRT(wd.rt.data, base, method, PAIR1, PAIR0, n1, n0, chrs, bed.gc)
save(nrds, file=file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.m.rt.log2s_", "m2-m1", ".RData")))
#load(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.m.rt.log2s_", "m2-m1", ".RData")))
nrow(nrds)
# [1] 
nrds.ov.au.m2 <- nrds










load(file.path(wd, "LCL/analysis/replication/lcl-wgs-rt/data/lcl_rpkm.gc.cn.d.rt.log2s_s-g1.RData"))
nrds.lcl <- nrds
load(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d.rt.log2s_", "m2-m1", ".RData")))

ymax <- 0.6
ymin <- 0.15
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   
   overlaps <- intersect(rownames(bed.gc.chr), rownames(nrds))   ## 29/11/19: Changed from intersect(rownames(nrds), rownames(bed.gc.chr))
   nrds.chr   <- nrds[overlaps,]
   bed.gc.chr <- bed.gc.chr[overlaps,]
   #lcl.rt.chr <- subset(lcl.rt, CHR == chr)   ## Koren 2012
   #nrds.lcl.chr <- nrds.lcl[intersect(nrds.lcl$BED, rownames(bed.gc.chr)),]
   
   ## Plot RT
   main.text <- paste0(BASE, " tumour replication timing")  
   file.name <- file.path(wd.rt.plots, paste0("RT_", BASE, "_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ""))   
   plotRT(file.name, main.text, chr, NA, NA, nrds.chr, bed.gc.chr, c(red, blue, green), c("M2 cells", "M1 cells"), c(red, blue), c("M2", "M1"), "png", width=10, peaks=c(), ylim=c(ymin, ymax), NULL, NULL)

   #file.name <- file.path(wd.rt.plots, "with-LCL", paste0("RT_", BASE, "_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ""))  
   #plotRT(file.name, main.text, chr, NA, NA, nrds.chr, bed.gc.chr, c(red, blue, green), c("M2 tumour", "M1 tumour"), c(red, blue), c("M2", "M1"), "png", width=10, peaks=c(), ylim=c(ymin, ymax), NULL, nrds.lcl.chr)
}

# > nrds.chr.o[270:280,]
#                 BED             T          N          RT       SPLINE1  SPLINE2  SPLINE3
# ICGC_PartNo_294 ICGC_PartNo_294 0.16045728 0.12269221 3.003422 2.030693 2.038254 1.778576
# ICGC_PartNo_295 ICGC_PartNo_295 0.15644159 0.11296179 3.669864 2.032079 2.040764 1.783627
# ICGC_PartNo_296 ICGC_PartNo_296 0.17987896 0.13692516 3.055784 2.032327       NA 1.788664
# ICGC_PartNo_297 ICGC_PartNo_297 0.14863725 0.11395546 2.972619 2.032616 2.041244 1.793689
# ICGC_PartNo_298 ICGC_PartNo_298 0.14912141 0.11700476 2.703235 2.033060 2.041499 1.798702
# ICGC_PartNo_299 ICGC_PartNo_299 0.14343473 0.10493446 3.517588 2.034666 2.042104 1.803701
# ICGC_PartNo_300 ICGC_PartNo_300 0.13820485 0.10113914 3.514047 2.036757 2.044610 1.808688
# ICGC_PartNo_301 ICGC_PartNo_301 0.13102423 0.08791241 4.523893 2.037865 2.046152 1.813662
# ICGC_PartNo_302 ICGC_PartNo_302 0.12249880 0.08535697 4.084337 2.039194 2.046755 1.818624
# ICGC_PartNo_303 ICGC_PartNo_303 0.12830215 0.08139079 5.176377 2.040367 2.048726 1.823572
# ICGC_PartNo_304 ICGC_PartNo_304 0.09346876 0.06252074 4.559813 2.040686 2.049028 1.828508

# SPLINE1: smooth.spline(x=bed.gc.chr.o$MEAN, y=nrds.chr.o[, column], tol=50)
# SPLINE2: spline <- smooth.spline(x=bed.gc.chr.o$START, y=nrds.chr.o[, column])
#          bed.gc.chr.o.spline <- subset(bed.gc.chr.o, START %in% spline$x)
# SPLINE3: smooth.spline(nrds.chr.o[, column])

# -----------------------------------------------------------------------------
# RD vs RT (RDS and SPR)
# Last Modified: 09/08/19; 31/05/19
# -----------------------------------------------------------------------------
sprs <- getSPR(nrds, bed.gc)
save(sprs, file=file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-m2-m1_spline_spearman.RData")))
writeTable(sprs, file=file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-m2-m1_spline_spearman.txt")), colnames=T, rownames=F, sep="\t")
#load(file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-m2-m1_spline_spearman.RData")))
sprs.cll <- sprs

for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   
   nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]
   nrds.chr.T  <- setSpline(nrds.chr, bed.gc.chr, "T")
   nrds.chr.N  <- setSpline(nrds.chr, bed.gc.chr, "N")
   nrds.chr.RT <- setSpline(nrds.chr, bed.gc.chr, "RT")
  
   main.text <- c(paste0("CLL read depth correlation (", "Chr", c, ")"), paste0("rho = ", round0(sprs$cor[c], digits=2), " (Q3 vs. Q1)"))
   xlab.text <- "CLL M2/M1"
   ylab.text <- "CLL read depth [RPKM]"
   file.name <- file.path(wd.rt.plots, "chrs", paste0("RD-vs-RT_CLL-M2-M1_chr", c, "_spline_spearman"))
   plotRD2vsRT(nrds.chr.T$SPLINE, nrds.chr.N$SPLINE, nrds.chr.RT$SPLINE, file.name, main.text, ylab.text, xlab.text, c("red", "blue"), c("M2", "M1"), method="spearman")
}

## Replication timing skew (RTS)
file.name <- file.path(wd.rt.plots, "RTS_CLL-M2-M1_spline_spearman_chr2")
main.text <- c("Replication timing skew", "RTS = (E-L)/(E+L)")
ylab.text <- "CLL M2/M1"
plotRTS(sprs.cll, file.name, main.text, c(4, 13, 17, 19), digits=3, unit=5, ylab.text, cex=1.2, chr2=-0.15, offset="            ")

### Figure 4D
## CLL vs. LCL
file.name <- file.path(wd.rt.plots, "RTS2_CLL-M2-M1_vs_LCL_spline_spearman")
main.text <- c("Replication timing skew", "")
xlab.text <- "LCL S/G1"
ylab.text <- "CLL M2/M1"
plotRTS2(sprs.cll$spr, sprs.lcl$spr, file.name, main.text, c(4, 13, 17, 19, 22), xlab.text, unit=5, ylab.text, cex=1.2)









## S-phase progression rate (SPR)
ylab.text <- "SPR"
file.name <- file.path(wd.rt.plots, "SPR_CLL-M2-M1_spline_spearman")
main.text <- c(paste0(BASE, " M2/M1 S-phase progression rate"), "SPR = (E-L)/(E+L)")
plotSPR(sprs, file.name, main.text, c(13, 17), digits=3, unit=5, ylab.text)

## SPR vs Read depth correlation
file.name <- file.path(wd.rt.plots, "SPR-RDC_CLL-M2-M1_spline_spearman")
main.text <- c(paste0(BASE, " M2/M1 SPR vs. Read depths correlation"), "")
xlab.text <- "M2 vs. M1 [rho]"
plotSPRRDC(sprs$spr, sprs$cor, file.name, main.text, c(4, 13, 17, 19, 22), xlab.text, unit=5, ylab.text)

## SPR vs Woodfine 2004
file.name <- file.path(wd.rt.plots, "SPR-Woodfine_CLL-M2-M1_spline_spearman")
main.text <- c(paste0(BASE, " M2/M1 SPR vs. Woodfine 2004"), "Mean replication timing ratio")
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

##
file.name <- file.path(wd.rt.plots, "RTD-vs-RT_CLL-M2-M1-vs-LCL-S-G1_spline_spearman")
main.text <- paste0("CLL M2/M1 vs. LCL S/G1")
ymin <- -1.1
ymax <- 1.1
plotRD3vsRTALL(cors, file.name, main.text, ymin, ymax, cols=c(red, blue, "black"), c("M2", "M1", "M2/M1"), c=NA, isRT=T)

ylab.text <- "Spearman's rho"
xlab.text <- "Chromosome"
file.name <- file.path(wd.rt.plots, "RT-vs-RT_CLL-M2-M1-vs-LCL-S-G1_spline_spearman")
main.text <- paste0("CLL M2/M1 vs. LCL S/G1")
ymin <- 0.25
ymax <- 1.05
plotRTvsRTALL(cors, file.name, main.text, ylab.text, xlab.text, ymin, ymax, col="black", c=2, pos=1)


# -----------------------------------------------------------------------------
# M2/M1 vs. Q4/Q4, LCL and SCLC-NL RTs
# Last Modified: 05/02/20
# -----------------------------------------------------------------------------
cors <- getRTvsRT3(nrds.cll.m2, nrds.cll.q4, nrds.sclc.nl.m2, nrds.lcl, bed.gc)
save(cors, file=file.path(wd.rt.data, paste0("rt-vs-rt3_", base, "-m2-m1-vs-ALL_spline_spearman.RData")))
#load(file.path(wd.rt.data, paste0("rt-vs-rt3_", base, "-m2-m1-vs-ALL_spline_spearman.RData")))

file.name <- file.path(wd.rt.plots, "RT-vs-RT2_CLL-M2-M1-vs-ALL")   ## gold (#f6c700)
main.text <- paste0(BASE, " RT")
ymin <- 0.85
ymax <- 1
plotRTvsRT2(cors, file.name, main.text, ymin, ymax, cols=c("black", yellow), c("M2/M1 vs. Q4/Q1", "CLL vs. SCLC-NL"))

file.name <- file.path(wd.rt.plots, "RT-vs-RT2_CLL-M2-M1-vs-ALL_0")   ## gold (#f6c700)
plotRTvsRT2(cors, file.name, main.text, ymin, ymax, cols=c("black", "white"), c("M2/M1 vs. Q4/Q1", "CLL vs. SCLC-NL"))

#file.name <- file.path(wd.rt.plots, "RT-vs-RT3_CLL-M2-M1-vs-ALL_spline_spearman")
#main.text <- paste0(BASE, " M2/M1")
#ymin <- 0.55
#ymax <- 1
#plotRTvsRT3(cors, file.name, main.text, ymin, ymax, cols=c("gold", "black", "#01DF01"), c("M2/M1 vs. Q4/Q1", "CLL vs. SCLC-NL", "CLL vs. LCL S/G1"))

# -----------------------------------------------------------------------------
# M2/M1 vs. NBL-CL RTs
# Last Modified: 14/11/20
# -----------------------------------------------------------------------------
cors <- getRTvsRT3(nrds.cll.m2, nrds.cll.m2, nrds.nbl.cl.m2, nrds.lcl, bed.gc)
save(cors, file=file.path(wd.rt.data, paste0("rt-vs-rt3_", base, "-m2-m1-vs-NBL-CL_spline_spearman.RData")))
#load(file.path(wd.rt.data, paste0("rt-vs-rt3_", base, "-m2-m1-vs-NBL-CL_spline_spearman.RData")))
# > median(cors$cor1)
# [1] 0.6621615

file.name <- file.path(wd.rt.plots, "RT-vs-RT2_CLL-M2-M1-vs-NBL-CL_spline_spearman")   ## gold (#f6c700)
main.text <- paste0(BASE, " M2/M1")
ymin <- 0.2
ymax <- 1
plotRTvsRT2(cors, file.name, main.text, ymin, ymax, cols=c("white", yellow), c("", "CLL vs. NBL-CL   "))   #00BA38 green for LCL RT??

# -----------------------------------------------------------------------------
# M2/M1 vs. LUAD RTs
# Last Modified: 26/11/20
# -----------------------------------------------------------------------------
cors <- getRTvsRT3(nrds.cll.m2, nrds.cll.m2, nrds.luad.m2, nrds.lcl, bed.gc)
save(cors, file=file.path(wd.rt.data, paste0("rt-vs-rt3_", base, "-m2-m1-vs-LUAD_spline_spearman.RData")))
#load(file.path(wd.rt.data, paste0("rt-vs-rt3_", base, "-m2-m1-vs-LUAD_spline_spearman.RData")))

file.name <- file.path(wd.rt.plots, "RT-vs-RT2_CLL-M2-M1-vs-LUAD_spline_spearman")   ## gold (#f6c700)
main.text <- paste0(BASE, " M2/M1")
ymin <- 0.55
ymax <- 0.9
plotRTvsRT2(cors, file.name, main.text, ymin, ymax, cols=c("white", yellow), c("", "CLL vs. LUAD "))   #00BA38 green for LCL RT??
# median(cors$cor1)
# [1] 0.8061612







# -----------------------------------------------------------------------------
# NBL M2/M1 vs NBL Q4/Q1
# Last Modified: 09/08/19; 11/07/19; 06/06/19
# -----------------------------------------------------------------------------
## NBL Q4/Q1
nrds.tmp <- nrds
load(file.path(wd, "CLL/analysis/replication/cll-wgs-rt-m2/data", paste0("nrds_cll-t-t_", method, ".RData")))
nrds.q4 <- nrds
nrds <- nrds.tmp

cors <- toTable(0, 3, 22, c("chr", "length", "cor"))
cors$chr <- 1:22
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]
   nrds.chr.RT <- setSpline(nrds.chr, bed.gc.chr, "RT")
 
   nrds.q4.chr <- nrds.q4[intersect(nrds.q4$BED, rownames(bed.gc.chr)),]  ## Reference LCL S/G1 ratio
   nrds.q4.chr.RT <- setSpline(nrds.q4.chr, bed.gc.chr, "RT")
 
   ## Keep only overlapping 1kb windows
   overlaps <- intersect(nrds.chr.RT$BED, nrds.q4.chr.RT$BED)
   cors$length[c] <- length(overlaps)
   cors$cor[c] <- getCor(nrds.chr.RT[overlaps,]$SPLINE, nrds.q4.chr.RT[overlaps,]$SPLINE, method="spearman")
}
#save(cors, file=file.path(wd.rt.data, paste0("rt-vs-rt_", base, "-m2-m1-vs-nbl-q4-q1_spline_spearman.RData")))
save(cors, file=file.path(wd.rt.data, paste0("rt-vs-rt_", base, "-m2-m1-vs-cll-m2-m1_spline_spearman.RData")))
#save(cors, file=file.path(wd.rt.data, paste0("rt-vs-rt_", base, "-m2-m1-vs-sclc-t-n_spline_spearman.RData")))

ylab.text <- "Spearman's rho"
xlab.text <- "Chromosome"
#file.name <- file.path(wd.rt.plots, "RT-vs-RT_NBL-M2-M1-vs-NBL-Q4-Q1_spline_spearman")
file.name <- file.path(wd.rt.plots, "RT-vs-RT_NBL-M2-M1-vs-CLL-M2-M1_spline_spearman")
#file.name <- file.path(wd.rt.plots, "RT-vs-RT_NBL-M2-M1-vs-SCLC-T-N_spline_spearman")
main.text <- paste0("NBL M2/M1 vs. CLL M2/M1")
ymin <- 0.25
ymax <- 1.05
plotRTvsRTALL(cors, file.name, main.text, ylab.text, xlab.text, ymin, ymax, col="black", c=2, pos=3)
