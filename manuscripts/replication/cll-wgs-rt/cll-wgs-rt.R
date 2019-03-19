# =============================================================================
# Manuscript   : 
# Chapter      : Reconstruction of replication timing profile in tumour cells
# Name         : manuscripts/replication/cll-wgs-rt.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 26/02/19
# =============================================================================
#wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Common.R", "DifferentialExpression.R", "ReplicationTiming.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.1kb.gc.RData"))

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 30/01/18
# -----------------------------------------------------------------------------
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE1 <- "CLL"
PAIR1 <- "T"
BASE0 <- "CLL"
PAIR0 <- "T"
base1 <- tolower(BASE1)
base0 <- tolower(BASE0)

wd.anlys <- file.path(wd, BASE1, "analysis")
wd.rt       <- file.path(wd.anlys, "replication", paste0(base1, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

wd.ngs <- file.path(wd, BASE1, "ngs/WGS")
#samples1 <- readTable(file.path(wd.ngs, "cll_wgs_n96.list"), header=F, rownames=F, sep="")
#samples0 <- readTable(file.path(wd.ngs, "cll_wgs_n96.list"), header=F, rownames=F, sep="")
samples  <- readTable(file.path(wd.ngs, "cll_wgs_n93.txt"), header=T, rownames=T, sep="")
samples1 <- readTable(file.path(wd.ngs, "cll_rna_n93-rt29.list"), header=F, rownames=F, sep="")
samples0 <- readTable(file.path(wd.ngs, "cll_rna_n93-wt33.list"), header=F, rownames=F, sep="")
n1 <- length(samples1)
n0 <- length(samples0)

# -----------------------------------------------------------------------------
# Plot RD and RT (see ReplicationTiming.R)
# Last Modified: 14/02/19; 10/01/19; 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
for (c in 1:22) {
   chr <- chrs[c]
   
   rpkms.chr.rt <- readTable(file.path(wd.rt.data, paste0(base1, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, flipRT=F, scaledRT=T) 
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   bed.gc.chr <- bed.gc.chr[rpkms.chr.rt$BED,]

   ## Colours (was "lightcoral", "skyblue3")
   adjustcolor.red  <- adjustcolor("lightcoral", alpha.f=0.3)
   adjustcolor.blue <- adjustcolor("skyblue3", alpha.f=0.3)
   adjustcolor.gray <- adjustcolor("gray", alpha.f=0.3)
   
   ## RD 
   ylab.text <- "Read depth"
   file.name <- file.path(wd.rt.plots, paste0("RD_", base1, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "_n", n1))
   main.text <- paste0("Read depth in ", BASE1, " tumour cells (n=", n1, ")")
   #plotRD(file.name, main.text, ylab.text, chr, NA, NA, rpkms.chr.rt$T, bed.gc.chr, c(adjustcolor.red, "red"), "png", 7.75, 9.25)   #(7.75, 9.25) for chr2
   
   file.name <- file.path(wd.rt.plots, paste0("RD_", base0, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR0, "_n", n0))
   main.text <- paste0("Read depth in ", BASE0, " normal cells (n=", n0, ")")
   #plotRD(file.name, main.text, ylab.text, chr, NA, NA, rpkms.chr.rt$N, bed.gc.chr, c(adjustcolor.blue, "blue"), "png", 7.75, 9.25)
   
   file.name <- file.path(wd.rt.plots, paste0("RD_", base0, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "+", PAIR0, "_n", n1, "-", n0))
   main.text <- paste0("Read depth in ", BASE1, " tumour (n=", n1, ") and normal (n=", n0, ") cells")
   #plotRD2(file.name, main.text, ylab.text, chr, NA, NA, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("Tumour", "Normal"), "png", 7.75, 9)

   ## RD & RT (N/T)
   #main.text <- paste0(BASE1, " N/T read depth ratio between normal (n=", n1, ") and tumour (n=", n0, ") cells")   
   #file.name <- file.path(wd.rt.plots, paste0("RTD_", base0, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0))   
   #plotRD3(file.name, main.text, chr, NA, NA, rpkms.chr.rt, bed.gc.chr, c("blue", "red"), c("Normal", "Tumour"), c(adjustcolor.gray, adjustcolor.gray), c("N", "T"), "png", width=10, peaks=c(74353001, 85951001), 7.5, 9.25, 3, 3)
   #plotRD3(file.name, paste0(BASE1, " N/T read depth ratio"), chr, 71500000, 90500000, rpkms.chr.rt, bed.gc.chr, c("blue", "red"), c("Normal", "Tumour"), c(adjustcolor.blue, adjustcolor.red), c("N", "T"), "png", width=5, peaks=c(74353001, 85951001), 7.5, 9.25, 3, 3)
   
   ## RD & RT (T29/T33)
   main.text <- paste0(BASE1, " T29/T33 read depth ratio between tumour (n=", n1, ") and tumour (n=", n0, ") cells")
   
   file.name <- file.path(wd.rt.plots, paste0("RTD_", base0, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0))   
   #plotRD3(file.name, main.text, chr, NA, NA, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("Tumour (n=29)", "Tumour (n=33)"), c(adjustcolor.gray, adjustcolor.gray), c("T29", "T33"), "png", width=10, peaks=c(), 7.5, 9.25, 3, 3)
   plotRD3(file.name, main.text, chr, NA, NA, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("Tumour (n=29)", "Tumour (n=33)"), c(adjustcolor.gray, adjustcolor.gray), c("T29", "T33"), "png", width=10, peaks=c(74353001, 85951001), 7.5, 9.25, 3, 3)
   plotRD3(file.name, paste0(BASE1, " T29/T33 read depth ratio"), chr, 71500000, 90500000, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("Tumour (n=29)", "Tumour (n=33)"), c(adjustcolor.red, adjustcolor.blue), c("T29", "T33"), "png", width=5, peaks=c(74353001, 85951001), 7.5, 9.25, 3, 3)
   
   ## HLA-DRB5 (chr6)
   file.name <- file.path(wd.rt.plots, paste0("RTD_", base0, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, "_HLA-DRB5"))
   plotRD3(file.name, main.text, chr, NA, NA, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("Tumour (n=29)", "Tumour (n=33)"), c(adjustcolor.gray, adjustcolor.gray), c("T29", "T33"), "png", width=10, peaks=c(32485120, 32498064), 7.5, 9.25, 3, 3)
   plotRD3(file.name, paste0(BASE1, " T29/T33 read depth ratio"), chr, 33000000, 34000000, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("Tumour (n=29)", "Tumour (n=33)"), c(adjustcolor.red, adjustcolor.blue), c("T29", "T33"), "png", width=5, peaks=c(32485120, 32498064), 7.5, 9.25, 3, 3)
   
   ## CCDC142 & TGOLN2 (chr2)
   file.name <- file.path(wd.rt.plots, paste0("RTD_", base0, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, "_CCDC142&TGOLN2"))
   plotRD3(file.name, main.text, chr, NA, NA, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("Tumour (n=29)", "Tumour (n=33)"), c(adjustcolor.gray, adjustcolor.gray), c("T29", "T33"), "png", width=10, peaks=c(74699113, 74710535, 85545147,	85555548), 7.5, 9.25, 3, 3)
   plotRD3(file.name, paste0(BASE1, " T29/T33 read depth ratio"), chr, 71500000, 90500000, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("Tumour (n=29)", "Tumour (n=33)"), c(adjustcolor.red, adjustcolor.blue), c("T29", "T33"), "png", width=5, peaks=c(74699113,	74710535, 85545147,	85555548), 7.5, 9.25, 3, 3)
   ## B2M (chr15)
   #plotRD3(file.name, main.text, chr, NA, NA, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("Tumour", "Tumour"), c(adjustcolor.gray, adjustcolor.gray), c("T", "T"), "png", width=10, peaks=c(45003675,	45011075), 7.5, 9.25, 3, 3)
   #plotRD3(file.name, paste0(BASE1, " T/T read depth ratio"), chr, 45000000,	46000000, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("Tumour", "Tumour"), c(adjustcolor.red, adjustcolor.blue), c("T", "T"), "png", width=5, peaks=c(45003675,	45011075), 7.5, 9.25, 3, 3)
   ## KLRC3 (chr12)
   #plotRD3(file.name, main.text, chr, NA, NA, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("Tumour", "Tumour"), c(adjustcolor.gray, adjustcolor.gray), c("T", "T"), "png", width=10, peaks=c(10564911,	10573194), 7.5, 9.25, 3, 3)
   #plotRD3(file.name, paste0(BASE1, " T/T read depth ratio"), chr, 10000000,	11000000, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("Tumour", "Tumour"), c(adjustcolor.red, adjustcolor.blue), c("T", "T"), "png", width=5, peaks=c(10564911,	10573194), 7.5, 9.25, 3, 3)
   ## RUNDC1 (chr17)
   #plotRD3(file.name, main.text, chr, NA, NA, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("Tumour", "Tumour"), c(adjustcolor.gray, adjustcolor.gray), c("T", "T"), "png", width=10, peaks=c(41132582,	41145707), 7.5, 9.25, 3, 3)
   #plotRD3(file.name, paste0(BASE1, " T/T read depth ratio"), chr, 40500000,	41500000, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("Tumour", "Tumour"), c(adjustcolor.red, adjustcolor.blue), c("T", "T"), "png", width=5, peaks=c(41132582,	41145707), 7.5, 9.25, 3, 3)
   
   ## RT
   ylab.text <- "Replication timing"
   file.name <- file.path(wd.rt.plots, paste0("RT_", base0, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0))
   main.text <- paste0(BASE1, " T29/T33 read depth ratio between tumour (n=", n1, ") and tumour (n=", n0, ") cells")
   plotRT(file.name, main.text, ylab.text, chr, NA, NA, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c(adjustcolor.gray, adjustcolor.gray), c("T", "T"), "png", 3, 3)
}

##
bed.gc.chr.t <- subset(subset(bed.gc.chr, START > 74000000), END < 75000000)
rpkms.chr.rt.t <- rpkms.chr.rt[rownames(bed.gc.chr.t),]
spline <- smooth.spline(x=bed.gc.chr.t$START, y=rpkms.chr.rt.t$RT)
# > max(spline$y)
# [1] 2.38142
# > rpkms.chr.rt.t[which(spline$y >= 2.3),]
# BED        T        N       RT
# P302861 P302861 7.837823 7.709673 1.622024
# P302862 P302862 7.937029 7.796177 1.782966
# P302863 P302863 7.835020 7.655530 2.272577
# P302864 P302864 8.006141 7.894438 1.413616
# P302865 P302865 7.843089 7.588346 3.226154
# P302866 P302866 7.395729 7.032869 4.596170
# > spline$y[which(spline$y >= 2.3)]
# [1] 2.311091 2.359962 2.381420 2.380529 2.362357 2.331969
# > bed.gc.chr.t["P302863",]
# CHR    START      END       GC
# P302863 chr2 74689001 74690000 0.599401
# > bed.gc.chr.t[c("P302863","P302864"),]
# CHR    START      END       GC
# P302863 chr2 74689001 74690000 0.599401
# P302864 chr2 74690001 74691000 0.503497

# > max(rpkms.chr.rt.t$RT)
# [1] 5.723453
# > rpkms.chr.rt.t[which(as.numeric(rpkms.chr.rt.t$RT) >= 5.72345),]
# BED        T        N       RT
# P302386 P302386 6.909331 6.457509 5.723453
# > bed.gc.chr.t["P302386",]
# CHR    START      END       GC
# P302386 chr2 74212001 74213000 0.778222

# -----------------------------------------------------------------------------
# CLL RD vs RT
# Last Modified: 18/02/19
# -----------------------------------------------------------------------------
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
 
   rpkms.chr.rt <- readTable(file.path(wd.rt.data, paste0(base1, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, flipRT=F, scaledRT=T) 
 
   rpkms.chr.rt.T  <- setSlope(rpkms.chr.rt, bed.gc.chr, "T")
   rpkms.chr.rt.N  <- setSlope(rpkms.chr.rt, bed.gc.chr, "N")
   rpkms.chr.rt.RT <- setSlope(rpkms.chr.rt, bed.gc.chr, "RT")
 
   main.text <- paste0("CLL RT (T29/T33) vs. Read depths (", "Chr", c, ")")
   xlab.text <- "CLL RT (T29/T33)"
   ylab.text <- "CLL read depth"
   file.name <- file.path(wd.rt.plots, paste0("plot_RT-vs-RD_CLL-T-T_pearson_chr", c, "_n"))
   xmin <- min(rpkms.chr.rt.RT$SLOPE)
   xmax <- max(rpkms.chr.rt.RT$SLOPE)
   ymin <- min(c(rpkms.chr.rt.T$SLOPE, rpkms.chr.rt.N$SLOPE))
   ymax <- max(c(rpkms.chr.rt.T$SLOPE, rpkms.chr.rt.N$SLOPE))
   plotRD2vsRT(rpkms.chr.rt.T$SLOPE, rpkms.chr.rt.N$SLOPE, rpkms.chr.rt.RT$SLOPE, file.name, main.text, ylab.text, xlab.text, c("red", "blue"), c("T29", "T33"), xmin, xmax, ymin, ymax)
}

# -----------------------------------------------------------------------------
# CLL RT vs LCL RT
# Last Modified: 05/03/19; 18/02/19
# -----------------------------------------------------------------------------
cors <- toTable(0, 2, 22, c("chr", "cor"))
cors$chr <- 1:22
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
 
   ## SCLC and LCL RTs
   rpkms.chr.rt <- readTable(file.path(wd.rt.data, paste0(base1, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, flipRT=F, scaledRT=T)
   rpkms.chr.rt.RT <- setSlope(rpkms.chr.rt, bed.gc.chr, "RT")
 
   rpkms.chr.rt.lcl <-readTable(paste0("/Users/tpyang/Work/uni-koeln/tyang2/LCL/analysis/replication/lcl-wgs-rt/data/lcl_rpkm.corr.gc.d.rt_", chr, "_LCL-LCL_n7-7.txt.gz"), header=T, rownames=T, sep="\t")
   rpkms.chr.rt.lcl <- setScaledRT(rpkms.chr.rt.lcl, pseudocount=0.01, recaliRT=T, flipRT=F, scaledRT=T) 
   rpkms.chr.rt.lcl.RT <- setSlope(rpkms.chr.rt.lcl, bed.gc.chr, "RT")
 
   ## Keep 1kb slopes based on overlapping windows
   overlaps <- intersect(rpkms.chr.rt.RT$BED, rpkms.chr.rt.lcl.RT$BED)
   #rpkms.chr.rt.RT     <- rpkms.chr.rt.RT[overlaps,]       ## Too slow
   #rpkms.chr.rt.lcl.RT <- rpkms.chr.rt.lcl.RT[overlaps,]   ## Too slow
 
   ##
   main.text <- paste0("CLL RT (T29/T33) vs. LCL RT (", "Chr", c, ")")
   xlab.text <- "CLL RT (T29/T33)"
   ylab.text <- "LCL RT (S/G1)"
   file.name <- file.path(wd.rt.plots, paste0("plot_RT-vs-RT_CLL-N-T-vs-LCL_pearson_chr", c, "_n"))
   plotRTvsRT(rpkms.chr.rt.lcl.RT[overlaps,]$SLOPE, rpkms.chr.rt.RT[overlaps,]$SLOPE, file.name, main.text, ylab.text, xlab.text, c(adjustcolor.gray, "black"))
 
   cors$cor[c] <- getCor(rpkms.chr.rt.lcl.RT[overlaps,]$SLOPE, rpkms.chr.rt.RT[overlaps,]$SLOPE)
}
save(cors, file=file.path(wd.rt.data, paste0("rt-vs-rt_", "cll-T-T", "-vs-lcl_cors-pearson.RData")))

ylab.text <- "Pearson's r"
xlab.text <- "Chromosome"
file.name <- file.path(wd.rt.plots, "plot_RT-vs-RT_CLL-T-T-vs-LCL_pearson")
main.text <- paste0("CLL RT (T29/T33) vs. LCL RT")
ymin <- 0
ymax <- 0.8
plotRTvsRTALL(cors, file.name, main.text, ylab.text, xlab.text, ymin, ymax, line0=T)






# -----------------------------------------------------------------------------
# SCLC-NBL RT vs LCL RT
# Last Modified: 18/02/19
# -----------------------------------------------------------------------------
BASE0 <- "LCL"
base0 <- tolower(BASE0)
n0 <- 7

cors <- toTable(0, 2, 22, c("chr", "cor"))
cors$chr <- 1:22
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   
   rpkms.chr.rt <- readTable(file.path(wd.rt.data, "hybrid", paste0(base1, "_rpkm.corr.gc.d.rt_", chr, "_", BASE1, "-", BASE0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, flipRT=F, scaledRT=T) 
   rpkms.chr.rt.RT <- setSlope(rpkms.chr.rt, bed.gc.chr, "RT")
   
   rpkms.chr.rt.lcl <-readTable(paste0("/Users/tpyang/Work/uni-koeln/tyang2/LCL/analysis/replication/lcl-wgs-rt/data/lcl_rpkm.corr.gc.d.rt_", chr, "_LCL-LCL_n7-7.txt.gz"), header=T, rownames=T, sep="\t")
   rpkms.chr.rt.lcl <- setScaledRT(rpkms.chr.rt.lcl, pseudocount=0.01, recaliRT=T, flipRT=F, scaledRT=T) 
   rpkms.chr.rt.lcl.RT <- setSlope(rpkms.chr.rt.lcl, bed.gc.chr, "RT")
   
   overlaps <- intersect(rpkms.chr.rt.RT$BED, rpkms.chr.rt.lcl.RT$BED)
   #rpkms.chr.rt.RT     <- rpkms.chr.rt.RT[overlaps,]
   #rpkms.chr.rt.lcl.RT <- rpkms.chr.rt.lcl.RT[overlaps,]
   
   cors$cor[c] <- getCor(rpkms.chr.rt.lcl.RT[overlaps,]$SLOPE, rpkms.chr.rt.RT[overlaps,]$SLOPE)
}
save(cors, file=file.path(wd.rt.data, paste0("rt-vs-rt_", base1, "&", base0, "-vs-lcl_cors-pearson.RData")))

ylab.text <- "Pearson's r"
xlab.text <- "Chromosome"
file.name <- file.path(wd.rt.plots, "plot_RT-vs-RT_SCLC&LCL-vs-LCL_pearson_All")
main.text <- paste0("SCLC–LCL RT (T/G1) vs. LCL RT")
ymin <- -0.2
ymax <- 0.85
plotRTvsRTALL(cors, file.name, main.text, ylab.text, xlab.text, ymin, ymax)

# -----------------------------------------------------------------------------
# 
# http://pklab.med.harvard.edu/scw2014/subpop_tutorial.html
# Last Modified: 06/03/19
# -----------------------------------------------------------------------------
plotSAMPLEvsRTALL <- function(cors.samples, samples1, file.name, main.text=NA, ylab.text=NA, xlab.text=NA, ymin=NA, ymax=NA, line0=NA) {
   cors.samples.plot <- toTable(0, 2, 22*length(samples1), c("chr", "cor"))
   n <- length(samples1)
   cnt <- 0
   for (c in 1:22) {
      start <- n * cnt + 1
      end   <- n * (cnt + 1)
      cors.samples.plot[start:end, 1] <- c
      cors.samples.plot[start:end, 2] <- as.numeric(cors.samples[c, samples1])
  
      cnt <- cnt + 1
   }
 
   png(paste0(file.name, "_boxplot.png"), height=5, width=5.5, units="in", res=300)
   boxplot(cor ~ chr, data=cors.samples.plot, ylab="Pearson's r (n=96)", outline=T, xaxt="n")
   axis(side=1, at=seq(2, 22, by=2))
   abline(h=0, lty=5)
   dev.off()
 
   png(paste0(file.name, "_scatter.png"), height=5, width=5.5, units="in", res=300)
   plot(cors.samples$mean, log(cors.samples$cv2), ylab="log(cv2)", xlab="mean")
   text(cors.samples$mean[2], log(cors.samples$cv2[2]), "chr2", cex=1.2, pos=3)
   text(cors.samples$mean[4], log(cors.samples$cv2[4]), "chr4", cex=1.2, pos=3)
   dev.off()
 
   png(paste0(file.name, "_var.png"), height=5, width=5.5, units="in", res=300)
   #plot(cors.samples$var ~ cors$chr, ylim=c(ymin, ymax), ylab=ylab.text, xlab=xlab.text, main=main.text, col="black", xaxt="n", yaxt="n", pch=19)
   plot(cors.samples$var ~ cors.samples$chr, ylab="var", xlab="Chromosome", col="black", xaxt="n", pch=19)
   lines(cors.samples$var, y=NULL, type="l", lwd=3)
   axis(side=1, at=seq(2, 22, by=2))
   dev.off()
}

##
cors.samples <- toTable(0, length(samples1)+4, 22, c("chr", "mean", "var", "cv2", samples1))
cors.samples$chr <- 1:22
for (s in 1:length(samples1)) {
   sample <- samples1[s]
   load(file.path(wd.rt.data, "samples", paste0("rt-vs-rt_", sample, "-vs-lcl_cors-pearson.RData")))

   cors.samples[, sample] <- cors$cor
}

for (c in 1:22) {
   cors.samples$mean[c] <- mean(as.numeric(cors.samples[c, samples1]))
   cors.samples$var[c]  <- var(as.numeric(cors.samples[c, samples1]))
   cors.samples$cv2[c]  <- cors.samples$var[c]/cors.samples$mean[c]^2
}
save(cors.samples, file=file.path(wd.rt.data, paste0("rt-vs-rt_samples-vs-lcl_cors-pearson.RData")))

ylab.text <- "Pearson's r (n=96)"
xlab.text <- "Chromosome"
file.name <- file.path(wd.rt.plots, "plot_RT-vs-RT_SAMPLES-vs-LCL_pearson")
main.text <- paste0("SCLC–LCL RT (T/G1) vs. LCL RT")   ## TO-DO
ymin <- -0.2
ymax <- 0.85
plotSAMPLEvsRTALL(cors.samples, samples1, file.name, main.text, ylab.text, xlab.text, ymin, ymax, line0=T)

# -----------------------------------------------------------------------------
# 
# Last Modified: 06/03/19
# -----------------------------------------------------------------------------
# > length(which(as.numeric(cors.samples[2, -c(1:4)]) > 0))
# [1] 40
# > length(which(as.numeric(cors.samples[2, -c(1:4)]) < 0))
# [1] 56

idx.pos.chrs <- which(as.numeric(cors.samples[1, -c(1:4)]) > 0) + 4
idx.neg.chrs <- which(as.numeric(cors.samples[1, -c(1:4)]) < 0) + 4
for (c in 2:22) {
   idx.pos.chr <- which(as.numeric(cors.samples[c, -c(1:4)]) > 0) + 4
   idx.neg.chr <- which(as.numeric(cors.samples[c, -c(1:4)]) < 0) + 4
   
   idx.pos.chrs <- intersect(idx.pos.chrs, idx.pos.chr)
   idx.neg.chrs <- intersect(idx.neg.chrs, idx.neg.chr)
}
# > length(idx.pos.chrs)
# [1] 29
# > length(idx.neg.chrs)
# [1] 33
idx.var.chrs <- c(5:100)
idx.var.chrs <- setdiff(idx.var.chrs, c(idx.pos.chrs, idx.neg.chrs))
# > length(idx.var.chrs)
# [1] 34
# > 29+33+34
# [1] 96

##
phenos <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/CLL/ngs/WGS/peiflyne/cll_cn_summary.txt", header=F, rownames=T, sep="")
colnames(phenos) <- c("V1", "purity", "ploidy", "V4", "sex", "V6", "V7")

samples.pos <- colnames(cors.samples[, idx.pos.chrs])
samples.neg <- colnames(cors.samples[, idx.neg.chrs])
samples.var <- colnames(cors.samples[, idx.var.chrs])
# > testU(phenos[samples.pos,]$purity, phenos[samples.neg,]$purity)
# [1] 0.8709254
# > testU(phenos[samples.pos,]$ploidy, phenos[samples.neg,]$ploidy)
# [1] 0.6065899

##
wd.rna <- file.path(wd, "CLL", "ngs/RNA")
ega <- readTable(file.path(wd.rna, "EGAD00001002130", "EGAD00001002130.list"), header=F, rownames=T, sep="")

# -----------------------------------------------------------------------------
# Mutation burdens between T29 and T33
# Last Modified: 15/03/19
# -----------------------------------------------------------------------------
for (s in 1:22) {
 
}









# -----------------------------------------------------------------------------
# 
# Last Modified: 18/02/19
# -----------------------------------------------------------------------------
setSlopes <- function(rpkms.chr.rt, bed.gc.chr, column) {
   spline <- smooth.spline(x=bed.gc.chr$START, y=rpkms.chr.rt[, column])
   slopes <- diff(spline$y)/diff(bed.gc.chr$START/1E6)   ## ADD 31/10/18
   rpkms.chr.rt <- rpkms.chr.rt[-nrow(rpkms.chr.rt),]
   rpkms.chr.rt$SLOPE <- slopes   ## length(slopes) is 1 less than nrow(bed.gc.chr), as no slope for the last 1kb window 

   sizes <- diff(bed.gc.chr$START)
   gaps <- which(sizes != 1000)
   
   return(rpkms.chr.rt[-gaps, c("BED", column, "SLOPE")])
}

plotRDvsRT <- function(reads, timings, file.name, main.text, ylab.text, xlab.text, colours) {
   lm.fit <- lm(reads ~ timings)
   #r2 <- summary(lm.fit)$r.squared
 
   #pdf(paste0(file.name, ".pdf"), height=6, width=6)
   png(paste0(file.name, ".png"), height=5, width=5, units="in", res=300)
   plot(reads ~ timings, ylab=ylab.text, xlab=xlab.text, main=main.text, col=colours[1])
   #plot(reads ~ timings, ylab=ylab.text, xlab=xlab.text, main=main.text, col="white")
   abline(lm.fit, col=colours[2], lwd=3)
   #mtext(paste0("R^2=", paste0(round0(r2*100, digits=2), "%")), cex=1.2, line=0.3)
   
   #rho <- cor.test(reads, timings, method="spearman", exact=F)[[4]]
   #mtext(paste0("SRC's rho = ", round0(rho, digits=2)), cex=1.2, line=0.3)
   
   cor <- cor.test(reads, timings, method="pearson")$estimate
   mtext(paste0("Pearson's r = ", round0(cor, digits=2)), cex=1.2, line=0.3) 
   dev.off()
   
   return(cor)
}

plotRD2vsRT <- function(reads1, reads2, timings, file.name, main.text, ylab.text, xlab.text, colours, legends, xmin, xmax, ymin, ymax) {
   png(paste0(file.name, ".png"), height=5, width=5, units="in", res=300)
   plot(NULL, xlim=c(xmin, xmax), ylim=c(ymin, ymax), ylab=ylab.text, xlab=xlab.text, main=main.text)

   cors <- c()
   lm.fit <- lm(reads1 ~ timings)
   abline(lm.fit, col=colours[1], lwd=3)
   cor <- cor.test(reads1, timings, method="pearson")$estimate
   cors <- c(cors, round0(abs(cor), digits=2))
   
   lm.fit <- lm(reads2 ~ timings)
   abline(lm.fit, col=colours[2], lwd=3)
   cor <- cor.test(reads2, timings, method="pearson")$estimate
   cors <- c(cors, round0(abs(cor), digits=2))
   
   legend("bottomright", c(paste0("r = –", cors[1], " (RT vs. ", legends[1], ")"), paste0("r = –", cors[2], " (RT vs. ", legends[2], ")")), text.col=colours, bty="n", cex=1.2)
   #legend("topright", paste0("r = ", cors[1], " (RT vs. ", legends[1], ")"), text.col=colours[1], bty="n", cex=1.2)
   #legend("bottomright", paste0("r = –", cors[2], " (RT vs. ", legends[2], ")"), text.col=colours[2], bty="n", cex=1.2)
   
   mtext("Pearson correlation", cex=1.2, line=0.3)
   dev.off()
}

plotRTvsRTALL <- function(cors, file.name, main.text, ylab.text, xlab.text, ymin, ymax) {
   png(paste0(file.name, ".png"), height=5, width=5, units="in", res=300)
   plot(cors$cor ~ cors$chr, ylim=c(ymin, ymax), ylab=ylab.text, xlab=xlab.text, main=main.text, col="black", xaxt="n", yaxt="n", pch=19)
   lines(cors$cor, y=NULL, type="l", lwd=3)
   abline(h=0, lty=5)
 
   text(cors$chr[2], cors$cor[2], round0(cors$cor[2], digits=2), cex=1.2, pos=1)
   axis(side=1, at=seq(2, 22, by=2))
   axis(side=2, at=seq(-0.6, 0.6, by=0.2), labels=c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6))
   dev.off()
}

###
## NBL RD vs RT
for (c in 1:22) {
   chr <- chrs[c]
   
   #rpkms.chr.rt <-readTable(paste0("/Users/tpyang/Work/uni-koeln/tyang2/LCL/analysis/replication/lcl-wgs-rt/data/lcl_rpkm.corr.gc.d.rt_", chr, "_LCL-LCL_n7-7.txt.gz"), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- readTable(file.path(wd.rt.data, paste0(base1, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, scaledRT=T) 
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   bed.gc.chr <- bed.gc.chr[rpkms.chr.rt$BED,]

   rpkms.chr.rt.N  <- setSlopes(rpkms.chr.rt, bed.gc.chr, "T")
   rpkms.chr.rt.T  <- setSlopes(rpkms.chr.rt, bed.gc.chr, "N")
   rpkms.chr.rt.RT <- setSlopes(rpkms.chr.rt, bed.gc.chr, "RT")
   
   main.text <- paste0("CLL RT vs. Read depths (", "Chr2", ")")
   xlab.text <- "NBL RT (N/T)"
   ylab.text <- "NBL read depth"
   file.name <- file.path(wd.rt.plots, paste0("plot_RT-vs-RD_CLL_pearson_chr2"))
   xmin <- min(rpkms.chr.rt.RT$SLOPE)
   xmax <- max(rpkms.chr.rt.RT$SLOPE)
   ymin <- min(c(rpkms.chr.rt.T$SLOPE, rpkms.chr.rt.N$SLOPE))
   ymax <- max(c(rpkms.chr.rt.T$SLOPE, rpkms.chr.rt.N$SLOPE))
   plotRD2vsRT(rpkms.chr.rt.N$SLOPE, rpkms.chr.rt.T$SLOPE, rpkms.chr.rt.RT$SLOPE, file.name, main.text, ylab.text, xlab.text, c("blue", "red"), c("N", "T"), xmin, xmax, ymin, ymax)
}

BASE0 <- "NBL"
base0 <- tolower(BASE0)
n0 <- 56

cors <- toTable(0, 2, 22, c("chr", "cor"))
cors$chr <- 1:22
for (c in 1:22) {
   chr <- chrs[c]
 
   rpkms.chr.rt <- readTable(file.path(wd.rt.data, paste0(base1, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, scaledRT=T) 
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   bed.gc.chr <- bed.gc.chr[rpkms.chr.rt$BED,]
   rpkms.chr.rt.RT     <- setSlopes(rpkms.chr.rt, bed.gc.chr, "RT")
 
   rpkms.chr.rt.lcl <-readTable(paste0("/Users/tpyang/Work/uni-koeln/tyang2/LCL/analysis/replication/lcl-wgs-rt/data/lcl_rpkm.corr.gc.d.rt_", chr, "_LCL-LCL_n7-7.txt.gz"), header=T, rownames=T, sep="\t")
   rpkms.chr.rt.lcl <- setScaledRT(rpkms.chr.rt.lcl, pseudocount=0.01, recaliRT=T, scaledRT=T) 
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   bed.gc.chr <- bed.gc.chr[rpkms.chr.rt.lcl$BED,]
   rpkms.chr.rt.lcl.RT <- setSlopes(rpkms.chr.rt.lcl, bed.gc.chr, "RT")
 
   overlaps <- intersect(rpkms.chr.rt.RT$BED, rpkms.chr.rt.lcl.RT$BED)
   rpkms.chr.rt.RT     <- rpkms.chr.rt.RT[overlaps,]
   rpkms.chr.rt.lcl.RT <- rpkms.chr.rt.lcl.RT[overlaps,]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   bed.gc.chr <- bed.gc.chr[overlaps,]
 
   cors$cor[c] <- getCor(rpkms.chr.rt.lcl.RT$SLOPE, rpkms.chr.rt.RT$SLOPE)
}
ylab.text <- "Pearson's r"
xlab.text <- "Chromosome"
file.name <- file.path(wd.rt.plots, "plot_RT-vs-RT_NBL&SCLC-vs-LCL_pearson_All")
main.text <- paste0("NBL-T/SCLC-N RT vs. LCL RT (All)")
ymin <- 0
ymax <- 0.8
plotRTvsRTALL(cors, file.name, main.text, ylab.text, xlab.text, ymin, ymax)
save(cors, file=file.path(wd.rt.data, paste0("rt-vs-rt_", base1, "&", base0, "-vs-lcl_cors-pearson.RData")))

###
## CLL RT vs LCL RT
cors <- toTable(0, 2, 22, c("chr", "cor"))
cors$chr <- 1:22
for (c in 1:22) {
   chr <- chrs[c]

   rpkms.chr.rt <- readTable(file.path(wd.rt.data, paste0(base1, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, scaledRT=T) 
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   bed.gc.chr <- bed.gc.chr[rpkms.chr.rt$BED,]
   rpkms.chr.rt.RT     <- setSlopes(rpkms.chr.rt, bed.gc.chr, "RT")
   
   rpkms.chr.rt.lcl <-readTable(paste0("/Users/tpyang/Work/uni-koeln/tyang2/LCL/analysis/replication/lcl-wgs-rt/data/lcl_rpkm.corr.gc.d.rt_", chr, "_LCL-LCL_n7-7.txt.gz"), header=T, rownames=T, sep="\t")
   rpkms.chr.rt.lcl <- setScaledRT(rpkms.chr.rt.lcl, pseudocount=0.01, recaliRT=T, scaledRT=T) 
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   bed.gc.chr <- bed.gc.chr[rpkms.chr.rt.lcl$BED,]
   rpkms.chr.rt.lcl.RT <- setSlopes(rpkms.chr.rt.lcl, bed.gc.chr, "RT")

   overlaps <- intersect(rpkms.chr.rt.RT$BED, rpkms.chr.rt.lcl.RT$BED)
   rpkms.chr.rt.RT     <- rpkms.chr.rt.RT[overlaps,]
   rpkms.chr.rt.lcl.RT <- rpkms.chr.rt.lcl.RT[overlaps,]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   bed.gc.chr <- bed.gc.chr[overlaps,]
   
   ylab.text <- "LCL RT (S/G1)"
   xlab.text <- "CLL RT (N/T)"
   file.name <- file.path(wd.rt.plots, paste0("plot_RT-vs-RT_CLL-vs-LCL_pearson_chr", c))
   main.text <- paste0("CLL RT vs. LCL RT (", "Chr", c, ")")
   cors$cor[c] <- plotRDvsRT(rpkms.chr.rt.lcl.RT$SLOPE, rpkms.chr.rt.RT$SLOPE, file.name, main.text, ylab.text, xlab.text, c(adjustcolor.gray, "black"))
   #cors$cor[c] <- getCor(rpkms.chr.rt.lcl.RT$SLOPE, rpkms.chr.rt.RT$SLOPE)
}
ylab.text <- "Pearson's r"
xlab.text <- "Chromosome"
file.name <- file.path(wd.rt.plots, "plot_RT-vs-RT_CLL-T-N-vs-LCL_pearson_All")
main.text <- paste0("CLL RT (T/N) vs. LCL RT (All)")
ymin <- -0.6
ymax <- 0.2
plotRTvsRTALL(cors, file.name, main.text, ylab.text, xlab.text, ymin, ymax)


save(cors, file=file.path(wd.rt.data, paste0("rt-vs-rt_", base1, "-N", "&", base0, "-T-vs-lcl_cors-pearson.RData")))

