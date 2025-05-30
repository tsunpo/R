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
#wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/ty2/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE  <- "SCLC"
PAIR1 <- "N"
PAIR0 <- "N"
base <- tolower(BASE)
method <- "rpkm"

wd.ngs   <- file.path(wd, BASE, "ngs/WGS")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.rt       <- file.path(wd.anlys, "replication", paste0(base, "-nl-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

samples1 <- readTable(file.path(wd.ngs, "sclc_wgs_n92_NL.txt"), header=T, rownames=T, sep="")
samples1 <- subset(samples1, M2 == 2)[,1]
samples0 <- readTable(file.path(wd.ngs, "sclc_wgs_n92_NL.txt"), header=T, rownames=T, sep="")
samples0 <- subset(samples0, M2 == 1)[,1]
n1 <- length(samples1)
n0 <- length(samples0)

# -----------------------------------------------------------------------------
# Plot RD and RT (see ReplicationTiming.R)
# Last Modified: 09/08/19; 14/02/19; 10/01/19; 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
nrds <- getLog2ScaledRT(wd.rt.data, base, method, PAIR1, PAIR0, n1, n0, chrs, bed.gc)
save(nrds, file=file.path(wd.rt.data, paste0(base, "_", method, ".cn.m.rt.log2s_", "m2-m1", ".RData")))
#load(file.path(wd.rt.data, paste0(base, "_", method, ".cn.m.rt.log2s_", "m2-m1", ".RData")))
nrow(nrds)
# [1] 2658679
# [1] 2673826
nrds.sclc.nl.m2 <- nrds

load(file.path(wd, "LCL/analysis/replication/lcl-wgs-rt/data/lcl_rpkm.cn.m.rt.log2s_s-g1.RData"))
nrds.lcl <- nrds
load(file.path(wd.rt.data, paste0(base, "_", method, ".cn.m.rt.log2s_", "m2-m1", ".RData")))

ymax <- 0.6
ymin <- 0.15
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]
   #lcl.rt.chr <- subset(lcl.rt, CHR == chr)   ## Koren 2012
   nrds.lcl.chr <- nrds.lcl[intersect(nrds.lcl$BED, rownames(bed.gc.chr)),]
   
   ## Plot RT
   main.text <- paste0(BASE, "\u2212", "NL M2/M1 read depth ratio RT")  
   file.name <- file.path(wd.rt.plots, paste0("RT_", BASE, "-NL_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ""))   
   #plotRT(file.name, main.text, chr, NA, NA, nrds.chr, bed.gc.chr, c(red, blue, green), c("M2 normal", "M1 normal"), c(red, blue), c("M2", "M1"), "png", width=10, peaks=c(), ylim=c(ymin, ymax), NULL, NULL)
   
   file.name <- file.path(wd.rt.plots, "with-LCL", paste0("RT_", BASE, "_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ""))  
   plotRT(file.name, main.text, chr, NA, NA, nrds.chr, bed.gc.chr, c(red, blue, green), c("M2 normal", "M1 normal"), c(red, blue), c("M2", "M1"), "png", width=13, peaks=c(), ylim=c(ymin, ymax), NULL, nrds.lcl.chr)
}

# -----------------------------------------------------------------------------
# RD vs RT (RDS and SPR)
# Last Modified: 11/07/19; 31/05/19
# -----------------------------------------------------------------------------
sprs <- getSPR(nrds, bed.gc)
save(sprs, file=file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-m2-m1_spline_spearman.RData")))
writeTable(sprs, file=file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-m2-m1_spline_spearman.txt")), colnames=T, rownames=F, sep="\t")
#load(file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-m2-m1_spline_spearman.RData")))
sprs.sclc.nl <- sprs

for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   
   nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]
   nrds.chr.T  <- setSpline(nrds.chr, bed.gc.chr, "T")
   nrds.chr.N  <- setSpline(nrds.chr, bed.gc.chr, "N")
   nrds.chr.RT <- setSpline(nrds.chr, bed.gc.chr, "RT")

   main.text <- c(paste0("SCLC (TN) read depth correlation (", "Chr", c, ")"), paste0("rho = ", round0(sprs$cor[c], digits=2), " (M2 vs. M1)"))
   xlab.text <- "SCLC (TN) M2/M1"
   ylab.text <- "SCLC (TN) read depth [RPKM]"
   file.name <- file.path(wd.rt.plots, "chrs", paste0("RD-vs-RT_SCLC-M2-M1_chr", c, "_spline_spearman"))
   plotRD2vsRT(nrds.chr.T$SPLINE, nrds.chr.N$SPLINE, nrds.chr.RT$SPLINE, file.name, main.text, ylab.text, xlab.text, c("red", "blue"), c("M2", "M1"), method="spearman")
}

## Replication timing skew (RTS)
file.name <- file.path(wd.rt.plots, "RTS_SCLC-NL-M2-M1_spline_spearman_chr2")
main.text <- c("Replication timing skew", "RTS = (E-L)/(E+L)")
ylab.text <- "SCLC-NL M2/M1"
plotRTS(sprs.sclc.nl, file.name, main.text, c(4, 13, 17, 19), digits=3, unit=5, ylab.text, cex=2, size=6)

### Figure 4D
## SCLC vs. LCL
file.name <- file.path(wd.rt.plots, "RTS2_SCLC-NL-M2-M1_vs_LCL_spline_spearman")
main.text <- c("Replication timing skew", "")
xlab.text <- "LCL S/G1"
ylab.text <- "SCLC-NL M2/M1"
plotRTS2(sprs.sclc.nl$spr, sprs.lcl$spr, file.name, main.text, c(4, 13, 17, 19), xlab.text, unit=5, ylab.text, cex=2, size=6)


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
file.name <- file.path(wd.rt.plots, "RTD-vs-RT_SCLC-NL-M2-M1-vs-LCL-S-G1_spline_spearman")
main.text <- paste0("SCLC-NL M2/M1 vs. LCL S/G1")
ymin <- -1.1
ymax <- 1.1
plotRD3vsRTALL(cors, file.name, main.text, ymin, ymax, cols=c(red, blue, "black"), c("M2", "M1", "M2/M1"), c=NA, isRT=T)

# -----------------------------------------------------------------------------
# Plot signal-to-noice
# Last Modified: 10/12/19
# -----------------------------------------------------------------------------
plotSNR2 <- function(n1, snr1, n2, snr2, file.name, main.text, xlab.text, ylab.text, col, pos) {
   n <- c(n1, n2)
   snr <- c(snr1, snr2)
 
   unit <- (max(snr) - min(snr))/10
   xlim <- c(min(snr) - unit, max(snr) + unit)
   unit <- (max(n) - min(n))/10
   ylim <- c(min(n) - unit, max(n) + unit)
 
   pdf(paste0(file.name, ".pdf"), height=6, width=6)
   plot(n1 ~ snr1, ylim=ylim, xlim=xlim, ylab="", xlab=xlab.text, main=main.text[1], col=col[1], pch=19, cex=2, cex.axis=1.8, cex.lab=1.9, cex.main=2.1)
   points(n2 ~ snr2, col=col[2], pch=19, cex=2)
 
   samples <- c("SCLC-NL", "SCLC", "NBL", "CLL")
   text(snr1, n1, samples, col=col[1], pos=c(1,2,2,2), cex=1.9)
   text(snr2, n2, samples, col=col[2], pos=c(1,4,4,4), cex=1.9)
 
   lm.fit <- lm(n2 ~ snr2)
   abline(lm.fit, col=col[2], lwd=4.5)
 
   lm.fit <- lm(n1 ~ snr1)
   abline(lm.fit, col=col[1], lwd=4.5)
 
   cor <- cor.test(n1, snr1, method="spearman", exact=F)
   #legend(pos[1], c(paste0("rho = ", round0(cor[[4]], digits=1)), paste0("p-value = ", scientific(cor[[3]], digits=1))), text.col=col[1], bty="n", cex=1.75)
   legend(pos[1], "M2/M1", text.col=col[1], pch=c(19, NA), cex=1.9)
 
   cor <- cor.test(n2, snr2, method="spearman", exact=F)
   #legend(pos[2], c(paste0("rho = ", round0(cor[[4]], digits=1)), paste0("p-value = ", scientific(cor[[3]], digits=1))), text.col=col[2], bty="n", cex=1.75)
   legend(pos[2], "Q4/Q1", text.col=col[2], pch=c(19, NA), col=col[2], box.col=col[2], cex=1.9)
 
   mtext(ylab.text, side=2, line=2.75, cex=1.9)
   #mtext(main.text[2], cex=1.2, line=0.3)
   dev.off()
}

## Sample size
snr    <- readTable(file.path(wd.rt.data, paste0("SNR_ALL.txt")), header=T, rownames=T, sep="\t")
snr.q4 <- readTable(file.path(wd.rt.data, paste0("SNR.Q4_ALL.txt")), header=T, rownames=T, sep="\t")

file.name <- file.path(wd.rt.plots, "STN2_ALL_SIZE")
main.text <- c("Signal-to-noise", "")
xlab.text <- "Signal-to-noise ratio"
ylab.text <- "Sample size"
plotSNR2(c(92, 101, 56, 96), snr$SNR, c(46, 51, 28, 48), snr.q4$SNR, file.name, main.text, xlab.text, ylab.text, c("black", grey), c("topleft", "bottomright"))

# -----------------------------------------------------------------------------
# M2/M1 vs. NBL-CL RTs
# Last Modified: 14/11/20
# -----------------------------------------------------------------------------
cors <- getRTvsRT3(nrds.sclc.nl.m2, nrds.sclc.nl.m2, nrds.nbl.cl.m2, nrds.lcl, bed.gc)
save(cors, file=file.path(wd.rt.data, paste0("rt-vs-rt3_", base, "-m2-m1-vs-NBL-CL_spline_spearman.RData")))
#load(file.path(wd.rt.data, paste0("rt-vs-rt3_", base, "-m2-m1-vs-NBL-CL_spline_spearman.RData")))

file.name <- file.path(wd.rt.plots, "RT-vs-RT2_SCLC-NL-M2-M1-vs-NBL-CL_spline_spearman")   ## gold (#f6c700)
main.text <- paste0(BASE, "-NL M2/M1")
ymin <- 0.2
ymax <- 1
plotRTvsRT2(cors, file.name, main.text, ymin, ymax, cols=c("white", yellow), c("", "SCLC-NL vs. NBL-CL "))   #00BA38 green for LCL RT??
# > median(cors$cor1)
# [1] 0.7754852

# -----------------------------------------------------------------------------
# M2/M1 vs. LUAD RTs
# Last Modified: 26/11/20
# -----------------------------------------------------------------------------
cors <- getRTvsRT3(nrds.sclc.nl.m2, nrds.sclc.nl.m2, nrds.luad.m2, nrds.lcl, bed.gc)
save(cors, file=file.path(wd.rt.data, paste0("rt-vs-rt3_", base, "-m2-m1-vs-LUAD_spline_spearman.RData")))
#load(file.path(wd.rt.data, paste0("rt-vs-rt3_", base, "-m2-m1-vs-LUAD_spline_spearman.RData")))

file.name <- file.path(wd.rt.plots, "RT-vs-RT2_SCLC-NL-M2-M1-vs-LUAD_spline_spearman")   ## gold (#f6c700)
main.text <- paste0(BASE, "-NL M2/M1")
ymin <- 0.55
ymax <- 0.9
plotRTvsRT2(cors, file.name, main.text, ymin, ymax, cols=c("white", yellow), c("", "SCLC-NL vs. LUAD "))   #00BA38 green for LCL RT??
# > median(cors$cor1)
# [1] 0.7838435
