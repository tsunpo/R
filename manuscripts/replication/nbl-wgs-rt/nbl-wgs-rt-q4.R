# =============================================================================
# Manuscript   : 
# Chapter      : Chromosome replication timing of the human genome
# Name         : manuscripts/replication/nbl-wgs-rt-q4.R
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
BASE  <- "NBL"
PAIR1 <- "T"
PAIR0 <- "T"
base <- tolower(BASE)
method <- "rpkm"

wd.ngs   <- file.path(wd, BASE, "ngs/WGS")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.rt       <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt-q4"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

samples1 <- readTable(file.path(wd.ngs, "nbl_wgs_n57-1.txt"), header=T, rownames=T, sep="")   ## M2/M1
#samples1 <- readTable(file.path(wd.ngs, "nbl_wgs_q4_n28.txt"), header=T, rownames=T, sep="")    ## Q4/Q1
samples1 <- subset(samples1, Q4 == 2)[,1]
samples0 <- readTable(file.path(wd.ngs, "nbl_wgs_n57-1.txt"), header=T, rownames=T, sep="")
#samples0 <- readTable(file.path(wd.ngs, "nbl_wgs_q4_n28.txt"), header=T, rownames=T, sep="")
samples0 <- subset(samples0, Q4 == 1)[,1]
n1 <- length(samples1)
n0 <- length(samples0)

# -----------------------------------------------------------------------------
# Plot RD and RT (see ReplicationTiming.R)
# Last Modified: 09/08/19; 14/02/19; 10/01/19; 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
nrds <- getLog2ScaledRT(wd.rt.data, base, method, PAIR1, PAIR0, n1, n0, chrs, bed.gc)
save(nrds, file=file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d.rt.log2s_", "q4-q1", ".RData")))
#load(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d.rt.log2s_", "q4-q1", ".RData")))
# > nrow(nrds)
# [1] 2659570
nrds.nbl.q4 <- nrds

##
nrds.nbl.q4 <- setSpline(nrds.nbl.q4, bed.gc.chr, "RT", returnAll=T)
snr.q4$S[3] <- sd(nrds.nbl.q4$SPLINE)
snr.q4$N[3] <- sd(nrds.nbl.q4$RT - nrds.nbl.q4$SPLINE)

##
ymax <- 0.6
ymin <- 0.15
for (c in 2:2) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]
   #lcl.rt.chr <- subset(lcl.rt, CHR == chr)   ## Koren 2012
   #nrds.lcl.chr <- nrds.lcl[intersect(nrds.lcl$BED, rownames(bed.gc.chr)),]
   
   ## Plot RT
   main.text <- paste0("Replication timing (RT)")  
   file.name <- file.path(wd.rt.plots, paste0("RT_", BASE, "_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ""))   
   #plotRT(file.name, main.text, chr, NA, NA, nrds.chr, bed.gc.chr, c(red, blue), c("Q4 tumour", "Q1 tumour"), c(red, blue), c("Q4", "Q1"), "png", width=10, peaks=c(), ylim=c(ymin, ymax), NULL, NULL)
   
   ## Plot RT
   main.text <- paste0(BASE, "") 
   #plotRT(file.name, main.text, chr, 39000000, 89000000, nrds.chr, bed.gc.chr, c(red, blue, green), c("Q4", "Q1"), c(red, blue), c("Q4", "Q1"), "png", width=5, peaks=c(), ylim=c(ymin, ymax), NULL, NULL)
   
   ## chr2
   main.text <- paste0("M2/M1 replication timing")  
   file.name <- file.path(wd.rt.plots, paste0("RT_", BASE, "_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, "_example_size=6_M2-M1_cex=0.6"))   
   plotRT(file.name, main.text, chr,  13000000,  17000000, nrds.chr, bed.gc.chr, c(red, blue, green), c("M2", "M1"), c(red, blue), c("M2", "M1"), "png", width=6, peaks=c(), ylim=c(ymin, ymax), NULL, NULL)
   #plotRT(file.name, main.text, chr,  39000000,  89000000, nrds.chr, bed.gc.chr, c(red, blue, green), c("S", "G1"), c(red, blue), c("S", "G1"), "png", width=6, peaks=c(), ylim=c(ymin, ymax), NULL, NULL)
   #plotRTAbstract(file.name, main.text, chr,  69500000,  81000000, nrds.chr, bed.gc.chr, c(red, blue, green), c("Q4", "Q1"), c(red, blue), c("Q4", "Q1"), "png", width=5, peaks=c(), ylim=c(ymin, ymax), NULL, NULL, 0.05)
   #plotRTAbstract(file.name, main.text, chr,  70000000,  80000000, nrds.chr, bed.gc.chr, c(red, blue, green), c("Q4", "Q1"), c(red, blue), c("Q4", "Q1"), "png", width=5, peaks=c(), ylim=c(ymin, ymax), NULL, NULL, 0.05)
   #plotRTAbstract(file.name, main.text, chr, 160000000, 170000000, nrds.chr, bed.gc.chr, c(red, blue, green), c("Q4", "Q1"), c(red, blue), c("Q4", "Q1"), "png", width=5, peaks=c(), ylim=c(ymin, ymax), NULL, NULL, 0.05)
}

# -----------------------------------------------------------------------------
# RD vs RT (RDS and SPR)
# Last Modified: 09/08/19; 31/05/19
# -----------------------------------------------------------------------------
sprs <- getSPR(nrds, bed.gc)
save(sprs, file=file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-q4-q1_spline_spearman.RData")))
writeTable(sprs, file=file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-q4-q1_spline_spearman.txt")), colnames=T, rownames=F, sep="\t")
#load(file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-q4-q1_spline_spearman.RData")))

for (c in 2:2) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
 
   nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]
   nrds.chr.T  <- setSpline(nrds.chr, bed.gc.chr, "T")
   nrds.chr.N  <- setSpline(nrds.chr, bed.gc.chr, "N")
   nrds.chr.RT <- setSpline(nrds.chr, bed.gc.chr, "RT")
 
   ## Figure 1
   xlab.text <- "RT [log2]"
   ylab.text <- "Read depth [RPKM]"
   main.text <- c(paste0("Spearman's correlation (", "Chr", c, ")"), "")   #, paste0("rho = ", round0(sprs$cor[c], digits=2), " (S vs. G1)"))
   file.name <- file.path(wd.rt.plots, "chrs", paste0("RD-vs-RT_NBL-Q4-Q1_chr", c, "_test"))
   plotRD2vsRT(nrds.chr.T$SPLINE, nrds.chr.N$SPLINE, nrds.chr.RT$SPLINE, file.name, main.text, ylab.text, xlab.text, c(red, blue), c("Q4", "Q1"), method="spearman")
 
 ## SFigure 1
 #xlab.text <- "RT [log2]"
 #ylab.text <- "Read depth [RPKM]"
 #main.text <- c(paste0("Chr", c))   #, paste0("rho = ", round0(sprs$cor[c], digits=2), " (S vs. G1)"))
 #file.name <- file.path(wd.rt.plots, "chrs", paste0("RD-vs-RT_LCL-S-G1_chr", c, "_spline_spearman"))
 #plotRD2vsRT(nrds.chr.T$SPLINE, nrds.chr.N$SPLINE, nrds.chr.RT$SPLINE, file.name, main.text, ylab.text, xlab.text, c(red, blue), c("S", "G1"), method="spearman")
 
 #main.text <- c(paste0("Chr", c), "")
 #file.name <- file.path(wd.rt.plots, "chrs", paste0("RD-vs-RT_LCL-G1_chr", c, "_spline_spearman"))
 #plotRDvsRT(nrds.chr.N$SPLINE, nrds.chr.RT$SPLINE, file.name, main.text, ylab.text, xlab.text, c(blue, adjustcolor(lighterblue, alpha.f=0.01)), c("S", "G1"), method="spearman")
 
 #main.text <- c(paste0("Chr", c), "")
 #file.name <- file.path(wd.rt.plots, "chrs", paste0("RD-vs-RT_LCL-S_chr", c, "_spline_spearman"))
 #plotRDvsRT(nrds.chr.T$SPLINE, nrds.chr.RT$SPLINE, file.name, main.text, ylab.text, xlab.text, c(red, adjustcolor(lighterred, alpha.f=0.01)), c("S", "G1"), method="spearman")
}

###
##
file.name <- file.path(wd.rt.plots, "RD2_NBL-Q4-Q1-vs-NBL-Q4-Q1")
main.text <- c("Correlation (Chr1-22)", "Spearman's rho")
plotRD2(sprs, file.name, main.text, 0.74, 1.00005)

###
## 24/01/22
file.name <- file.path(wd.rt.plots, "RD2_NBL-Q4-Q1-vs-LCL-S-G1")
main.text <- c("Correlation with RT (Chr1-22)", "Spearman's rho")
plotRD2(cors, file.name, main.text, 0.3468236, 0.8311754)
# > max(cors$cor1, abs(cors$cor2))
# [1] 0.8311754
# > min(cors$cor1, abs(cors$cor2))
# [1] 0.3468236


## S-phase progression rate (SPR)
file.name <- file.path(wd.rt.plots, "SPR_NBL-Q4-Q1_spline_spearman")
main.text <- c("S-phase progression rate (SPR)")
ylab.text <- "SPR = (E-L)/(E+L)"
plotSPR(sprs, file.name, main.text, c(13, 17), digits=3, unit=5, ylab.text)

## SPR vs Woodfine 2004
file.name <- file.path(wd.rt.plots, "SPR-Woodfine_NBL-Q4-Q1_spline_spearman")
main.text <- c("Mean replication timing ratio", "")
xlab.text <- "Woodfine et al. 2004"
ylab.text <- "SPR"
plotSPRRDC(sprs$spr, lcl.mean$Mean, file.name, main.text, c(13, 17, 19, 22), xlab.text, unit=5, ylab.text)

## SPR vs Read depth correlation
#file.name <- file.path(wd.rt.plots, "SPR-RDC_NBL-Q4-Q1_spline_spearman")
#main.text <- c(paste0(BASE, " Q4/Q1 SPR vs. Read depths correlation"), "")
#xlab.text <- "Q4 vs. Q1 [rho]"
#plotSPRRDC(sprs$spr, sprs$cor, file.name, main.text, c(4, 13, 17, 19, 22), xlab.text, unit=5, ylab.text)

# -----------------------------------------------------------------------------
# RT vs LCL S/G1
# Last Modified: 09/08/19; 06/06/19
# -----------------------------------------------------------------------------
nrds.tmp <- nrds
load(file.path(wd, "LCL/analysis/replication/lcl-wgs-rt/data/lcl_rpkm.gc.cn.d.rt.log2s_s-g1.RData"))
nrds.lcl <- nrds
nrds <- nrds.tmp

cors <- getRTvsRT(nrds, nrds.lcl, bed.gc)
save(cors, file=file.path(wd.rt.data, paste0("rt-vs-rt_", base, "-q4-q1-vs-lcl-s-g1_spline_spearman.RData")))
#load(file.path(wd.rt.data, paste0("rt-vs-rt_", base, "-q4-q1-vs-lcl-s-g1_spline_spearman.RData")))

##
file.name <- file.path(wd.rt.plots, "RTD-vs-RT_NBL-Q4-Q1-vs-LCL-S-G1_spline_spearman")
main.text <- paste0("NBL Q4/Q1 vs. LCL RT")
ylab.text <- "Correlation [rho]"
xlab.text <- "Chromosome"
ymin <- -1.1
ymax <- 1.1
plotRD3vsRTALL(cors, file.name, main.text, ymin, ymax, cols=c(red, blue, "black"), c("Q4", "Q1", "LCL RT"), c=NA, isRT=T)

#ylab.text <- "Spearman's rho"
#xlab.text <- "Chromosome"
#file.name <- file.path(wd.rt.plots, "RT-vs-RT_NBL-Q4-Q1-vs-LCL-S-G1_spline_spearman")
#main.text <- paste0("NBL Q4/Q1 vs. LCL S/G1")
#ymin <- 0.25
#ymax <- 1.05
#plotRTvsRTALL(cors, file.name, main.text, ylab.text, xlab.text, ymin, ymax, col="black", c=2, pos=1)

###
##
file.name <- file.path(wd.rt.plots, "RD2_NBL-Q4-Q1-vs-LCL-S-G1_spline_spearman_nasa.blue_google.red")
main.text <- "NBL Q4/Q1"
plotRD2(cors, file.name, main.text, 0.3, 0.85)







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
