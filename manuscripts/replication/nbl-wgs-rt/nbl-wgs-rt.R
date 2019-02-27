# =============================================================================
# Manuscript   : 
# Chapter      : Reconstruction of replication timing profile in tumour cells
# Name         : manuscripts/replication/nbl-wgs-rt.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 25/02/19
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
BASE1 <- "NBL"
PAIR1 <- "T"
BASE0 <- "NBL"
PAIR0 <- "N"
base1 <- tolower(BASE1)
base0 <- tolower(BASE0)

wd.anlys <- file.path(wd, BASE1, "analysis")
wd.rt       <- file.path(wd.anlys, "replication", paste0(base1, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

wd.ngs <- file.path(wd, BASE1, "ngs/WGS")
samples1 <- readTable(file.path(wd.ngs, "nbl_wgs_n57-1.list"), header=F, rownames=F, sep="")
samples0 <- readTable(file.path(wd.ngs, "nbl_wgs_n57-1.list"), header=F, rownames=F, sep="")
n1 <- length(samples1)
n0 <- length(samples0)

# -----------------------------------------------------------------------------
# Plot RD and RT (see ReplicationTiming.R)
# Last Modified: 14/02/19; 10/01/19; 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
for (c in 1:22) {
   chr <- chrs[c]
   
   rpkms.chr.rt <- readTable(file.path(wd.rt.data, paste0(base1, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, scaledRT=T) 
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
   plotRD(file.name, main.text, ylab.text, chr, NA, NA, rpkms.chr.rt$T, bed.gc.chr, c(adjustcolor.red, "red"), "png", 7.75, 9.25)   #(7.75, 9.25) for chr2
   
   file.name <- file.path(wd.rt.plots, paste0("RD_", base0, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR0, "_n", n0))
   main.text <- paste0("Read depth in ", BASE0, " normal cells (n=", n0, ")")
   plotRD(file.name, main.text, ylab.text, chr, NA, NA, rpkms.chr.rt$N, bed.gc.chr, c(adjustcolor.blue, "blue"), "png", 7.75, 9.25)
   
   file.name <- file.path(wd.rt.plots, paste0("RD_", base0, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "+", PAIR0, "_n", n1, "-", n0))
   main.text <- paste0("Read depth in ", BASE1, " tumour (n=", n1, ") and normal (n=", n0, ") cells")
   plotRD2(file.name, main.text, ylab.text, chr, NA, NA, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("Tumour", "Normal"), "png", 7.75, 9)

   ## RD & RT 
   main.text <- paste0(BASE1, " T/N read depth ratio between tumour (n=", n1, ") and normal (n=", n0, ") cells")   
   file.name <- file.path(wd.rt.plots, paste0("RTD_", base0, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0))   
   plotRD3(file.name, main.text, chr, NA, NA, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("Tumour", "Normal"), c(adjustcolor.gray, adjustcolor.gray), c("T", "N"), "png", width=10, peaks=c(74353001, 85951001), 7.75, 9, 3, 3)
   plotRD3(file.name, paste0(BASE1, " T/N read depth ratio"), chr, 71500000, 90500000, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("Tumour", "Normal"), c(adjustcolor.red, adjustcolor.blue), c("T", "N"), "png", width=5, peaks=c(74353001, 85951001), 7.75, 9, 3, 3)
   
   ## RT
   ylab.text <- "Replication timing"
   file.name <- file.path(wd.rt.plots, paste0("RT_", base0, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0))
   main.text <- paste0(BASE1, " T/N read depth ratio between tumour (n=", n1, ") and normal (n=", n0, ") cells")
   plotRT(file.name, main.text, ylab.text, chr, NA, NA, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c(adjustcolor.gray, adjustcolor.gray), c("T", "N"), "png", 3, 3)
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

plotChrs <- function(cors, file.name, main.text, ylab.text, xlab.text, ymin, ymax) {
   png(paste0(file.name, ".png"), height=5, width=5, units="in", res=300)
   plot(cors$cor ~ cors$chr, ylim=c(ymin, ymax), ylab=ylab.text, xlab=xlab.text, main=main.text, col="black", xaxt="n", pch=19)
   lines(cors$cor, y=NULL, type="l", lwd=3)

   text(cors$chr[2], cors$cor[2], round0(cors$cor[2], digits=2), cex=1.2, pos=3)   
   axis(side=1, at=seq(2, 22, by=2))
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

   rpkms.chr.rt.T  <- setSlopes(rpkms.chr.rt, bed.gc.chr, "T")
   rpkms.chr.rt.N  <- setSlopes(rpkms.chr.rt, bed.gc.chr, "N")
   rpkms.chr.rt.RT <- setSlopes(rpkms.chr.rt, bed.gc.chr, "RT")
   
   main.text <- paste0("NBL RT vs. Read depths (", "Chr2", ")")
   xlab.text <- "NBL RT (T/N)"
   ylab.text <- "NBL read depth"
   file.name <- file.path(wd.rt.plots, paste0("plot_RT-vs-RD_NBL_pearson_chr2"))
   xmin <- min(rpkms.chr.rt.RT$SLOPE)
   xmax <- max(rpkms.chr.rt.RT$SLOPE)
   ymin <- min(c(rpkms.chr.rt.T$SLOPE, rpkms.chr.rt.N$SLOPE))
   ymax <- max(c(rpkms.chr.rt.T$SLOPE, rpkms.chr.rt.N$SLOPE))
   plotRD2vsRT(rpkms.chr.rt.T$SLOPE, rpkms.chr.rt.N$SLOPE, rpkms.chr.rt.RT$SLOPE, file.name, main.text, ylab.text, xlab.text, c("red", "blue"), c("T", "N"), xmin, xmax, ymin, ymax)
}

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
## NBL RT vs LCL RT
BASE0 <- "CLL"
base0 <- tolower(BASE0)
n0 <- 96

cors <- toTable(0, 2, 22, c("chr", "cor"))
cors$chr <- 1:22
for (c in 1:22) {
   chr <- chrs[c]

   rpkms.chr.rt <- readTable(file.path(wd.rt.data, "hybrid", paste0(base1, "_rpkm.corr.gc.d.rt_", chr, "_", BASE1, "-", BASE0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
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
   
   #ylab.text <- "LCL RT (S/G1)"
   #xlab.text <- "NBL RT (T/N)"
   #file.name <- file.path(wd.rt.plots, paste0("plot_RT-vs-RT_NBL-vs-LCL_pearson_chr", c))
   #main.text <- paste0("NBL RT vs. LCL RT (", "Chr", c, ")")
   #cors$cor[c] <- plotRDvsRT(rpkms.chr.rt.lcl.RT$SLOPE, rpkms.chr.rt.RT$SLOPE, file.name, main.text, ylab.text, xlab.text, c(adjustcolor.gray, "black"))
   cors$cor[c] <- getCor(rpkms.chr.rt.lcl.RT$SLOPE, rpkms.chr.rt.RT$SLOPE)
}
ylab.text <- "Pearson's r"
xlab.text <- "Chromosome"
file.name <- file.path(wd.rt.plots, "plot_RT-vs-RT_NBL&CLL-T-vs-LCL_pearson_All")
main.text <- paste0("NBL–CLL (T/T) vs. LCL RT")
ymin <- -0.2
ymax <- 0.85
plotRTvsRTALL(cors, file.name, main.text, ylab.text, xlab.text, ymin, ymax)
save(cors, file=file.path(wd.rt.data, paste0("rt-vs-rt_", base1, "&", base0, "-T-vs-lcl_cors-pearson.RData")))



ylab.text <- "LCL (S/G1)"
xlab.text <- "SCLC (T/N)"
file.name <- file.path(wd.rt.plots, paste0("plot_RT-LCL_vs_RT-SCLC_slope_1kb_chr2"))
main.text <- paste0("Correlation between RT profiles (", "Chr2", ")")
cor <- plotRDvsRT(rpkms.chr.rt.lcl.RT$SLOPE, rpkms.chr.rt.RT$SLOPE, file.name, main.text, ylab.text, xlab.text, c("grey", "black"))




   ## LCL RD vs RT
   rpkms.chr.rt.lcl <-readTable(paste0("/Users/tpyang/Work/uni-koeln/tyang2/LCL/analysis/replication/lcl-wgs-rt/data/lcl_rpkm.corr.gc.d.rt_", chr, "_LCL-LCL_n7-7.txt.gz"), header=T, rownames=T, sep="\t")
   rpkms.chr.rt.lcl$RT <- scale(rpkms.chr.rt.lcl$RT)   ## ADD 19/02/19

bed.gc.chr <- subset(bed.gc, CHR == chr)
bed.gc.chr <- bed.gc.chr[rpkms.chr.rt.lcl$BED,]
overlaps <- intersect(rpkms.chr.rt.lcl$BED, rownames(bed.gc.chr))
idxes <- seq(1, length(overlaps), 1)
bed.gc.chr <- bed.gc.chr[overlaps[idxes],]

rpkms.chr.rt.lcl.RT <- setSlopes(rpkms.chr.rt.lcl, bed.gc.chr)
rpkms.chr.rt.lcl.RD  <- setSlopes2(rpkms.chr.rt.lcl, bed.gc.chr)

ylab.text <- "Read depth (G1)"
xlab.text <- "RT (LCL)"
file.name <- file.path(wd.rt.plots, "plot_RD_vs_RT_LCL_G1-LCL_slope_rho_1kb_log2")
main.text <- paste0("Slopes of each 1kb window (Chr2)")
plotRDvsRT(rpkms.chr.rt.lcl.RD$SLOPE_N, rpkms.chr.rt.lcl.RT$SLOPE, file.name, main.text, ylab.text, xlab.text, c("lightskyblue1", "blue"))

ylab.text <- "Read depth (S)"
file.name <- file.path(wd.rt.plots, "plot_RD_vs_RT_LCL_S-LCL_slope_rho_1kb_log2")
plotRDvsRT(rpkms.chr.rt.lcl.RD$SLOPE_T, rpkms.chr.rt.lcl.RT$SLOPE, file.name, main.text, ylab.text, xlab.text, c("pink", "red"))










# -----------------------------------------------------------------------------
# Insert size (each sample)
# Last Modified: 20/09/18
# -----------------------------------------------------------------------------
wd.meta <- file.path(wd.ngs, "insert_size")

for (s in 1:length(samples)) {
   table <- readTable(file.path(wd.meta, paste0(samples[s], "_qc.txt")), header=F, rownames=F, sep="")   
}



# -----------------------------------------------------------------------------
# Insert size
# Last Modified: 10/09/18
# -----------------------------------------------------------------------------
wd.meta <- file.path(wd, BASE, "metadata/George 2015")
table <- readTable(file.path(wd.meta, "nature14664-s1_ST2.txt"), header=T, rownames=T, sep="\t")
table <- table[samples,]

table1 <- table[,c("Insert.Size", "Mean.Coverage")]
table1$Group <- 1
table0 <- table[,c("Insert.Size.1", "Mean.Coverage.1")]
colnames(table0) <- c("Insert.Size", "Mean.Coverage")
table0$Group <- 0
table <- rbind(table1, table0)
table$Group <- as.factor(table$Group)

file.name <- file.path(wd.rt.plots, paste0("boxplot_", base, "_wgs_insert-size.pdf"))
pdf(file.name, height=6, width=3)
boxplot(Insert.Size ~ Group, data=table, outline=F, names=c("Normal", "Tumour"), col=c("dodgerblue", "red"), ylab="Insert size", main=BASE)
dev.off()

median(table1$Insert.Size)
# [1] 312.312
median(table0$Insert.Size)
# [1] 308.56
median(table0[normals,]$Insert.Size)

sd(table1$Insert.Size)
# [1] 16.19834
sd(table0$Insert.Size)
# [1] 13.34431

## Tumour vs. Normal
testW(table1$Insert.Size, table0$Insert.Size)
# [1] 0.5489026


# -----------------------------------------------------------------------------
# Step 6.1: Define replicaiton timing direction for expressed genes (Following Step 4 in "asym-sclc-tx.R" and Step 5 from rt-sclc-wgs.R)
# Link(s):  http://www.mun.ca/biology/scarr/2250_DNA_replication_&_transcription.html
#           https://stackoverflow.com/questions/43615469/how-to-calculate-the-slope-of-a-smoothed-curve-in-r
# Last Modified: 29/01/18
# -----------------------------------------------------------------------------
plotRT0 <- function(wd.rt.plots, BASE, chr, n1, n0, xmin, xmax, rpkms.chr.rt, bed.gc.chr.rt, pair1, pair0, ext, spline) {
   file.name  <- file.path(wd.rt.plots, paste0(tolower(BASE), "_wgs_rt_bstrp1000_", chr, "_", pair1, "-", pair0, "_n", n1, "-", n0, "_spline"))
   main.text <- paste0("Bootstrapped read depth (CN-, GC-corrected) ratio (", pair1, "/", pair0, ") in ", BASE)
   xlab.text <- paste0("Chromosome ", gsub("chr", "", chr), " coordinate (Mb)")
   ylab.text <- "Replication time (log2 FC)"
 
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   ymin <- min(rpkms.chr.rt$MEDIAN)
   ymax <- max(rpkms.chr.rt$MEDIAN)
   if (!is.na(xmin) && !is.na(xmax)) file.name <- paste0(file.name, "_", xmin/1E6, "-", xmax/1E6, "Mb")
   if (is.na(xmin)) xmin <- 0
   if (is.na(xmax)) xmax <- subset(chromInfo, chrom == chr)$size

   if (ext == "pdf") {
      pdf(paste0(file.name, ".pdf"), height=4, width=10)
   } else if (ext == "png")
      png(paste0(file.name, ".png"), height=4, width=10, units="in", res=300)   ## ADD 16/05/17: res=300
   
   plot(NULL, ylim=c(ymin, ymax), xlim=c(xmin/1E6, xmax/1E6), xlab=xlab.text, ylab=ylab.text, main=main.text)
   points(bed.gc.chr.rt$START/1E6, rpkms.chr.rt$MEDIAN, col="red", cex=0.3)
   abline(h=0, lwd=0.5, col="grey")
   lines(bed.gc.chr.rt$START/1E6, smooth.spline(rpkms.chr.rt$MEDIAN)$y)
   lines(spline$x/1E6, spline$y, col="blue")
   
   #slopes <- diff(smooth.spline(rpkms.chr)$y)/diff((bed.gc.chr$START)/1E6)
   #slopes2 <- diff(smooth.spline(rpkms.chr)$y)/diff(smooth.spline(rpkms.chr)$x)
   
   #temp <- loess.smooth(bed.gc.chr$START, rpkms.chr)
   #slopes3 <- diff(temp$y)/diff(temp$x)
   
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.4, col="lightgrey")
 
   dev.off()
}

getEnsGeneBED <- function(pos, bed.gc.chr) {
   bed.gc.chr.start <- subset(bed.gc.chr, pos >= START)
   bed.gc.chr.start.end <- subset(bed.gc.chr.start, pos <= END)
 
   return(rownames(bed.gc.chr.start.end))
}

BASE  <- "SCLC"
PAIR1 <- "T"
PAIR0 <- "N"
PAIR  <- paste0(PAIR1, "-", PAIR0)
CHR   <- 6
CUTOFF <- 0.15

###
##
bed.gc <- bed[which(bed$GC > 0),]   ## Only keep partitions (in the BED file) with a GC content
ensGene.tx <- ensGene[rownames(tpm.gene.input),]

ensGene.tx.rt <- ensGene.tx[1,]
ensGene.tx.rt$SLOPE_START <- 0
ensGene.tx.rt$SLOPE_END <- 0
ensGene.tx.rt <- ensGene.tx.rt[-1,]
for (c in 1:22) {
   #chr <- chrs[CHR]
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)

   ## Replication timing
   #rpkms.chr <- readTable(file.path(wd.rt.data, paste0(base, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR, "_n", length(samples), "-", length(normals), ".txt.gz")), header=T, rownames=T, sep="\t") 
   rpkms.chr <- readTable(file.path(wd.rt.data, paste0(base, "_rpkm.corr.gc.d.rt_bstrp1000_", chr, "_T-N_n101-92.txt.gz")), header=T, rownames=T, sep="\t") 
   rpkms.chr$MEDIAN <- rpkms.chr$RT   ## REMOVED 07/10/18
   
   ##
   rpkms.chr.rt <- rpkms.chr[which(rpkms.chr$MEDIAN > -CUTOFF),]
   rpkms.chr.rt <- rpkms.chr.rt[which(rpkms.chr.rt$MEDIAN < CUTOFF),]
   overlaps <- intersect(rownames(rpkms.chr.rt), rownames(bed.gc.chr))
   bed.gc.chr.rt <- bed.gc.chr[overlaps,]
   
   #plotRT0(wd.rt.plots, BASE, chr, 101, 92, NA, NA, rpkms.chr.rt, bed.gc.chr.rt, PAIR1, PAIR0, "png")
   plotRT0(wd.rt.plots, BASE, chr, 101, 92, 87730060,98505459, rpkms.chr.rt, bed.gc.chr.rt, PAIR1, PAIR0, "png")   ## CCDC67
   #plotRT0(wd.rt.plots, BASE, chr, 101, 92, 96754840, 131794839, rpkms.chr.rt, bed.gc.chr.rt, PAIR1, PAIR0, "png")   ## HDAC2
   #plotRT0(wd.rt.plots, BASE, chr, length(samples), 140813453, 153118090, rpkms.chr.rt, bed.gc.chr, PAIR1, PAIR0, "png")   ## CNTNAP2
   #plotRT0(wd.rt.plots, BASE, chr, length(samples), 43877887, 54056122, rpkms.chr.rt, bed.gc.chr, PAIR1, PAIR0, "png")   ## RB1
   #plotRT0(wd.rt.plots, BASE, chr, length(samples), 147504475, 149581413, rpkms.chr.rt, bed.gc.chr, PAIR1, PAIR0, "png")   ## EZH2
   #plotRT0(wd.rt.plots, BASE, chr, length(samples), 7314246, 11612723, rpkms.chr.rt, bed.gc.chr, PAIR1, PAIR0, "png")   ## PTPRD
   #plotRT0(wd.rt.plots, BASE, chr, length(samples), 114521235, 118716095, rpkms.chr.rt, bed.gc.chr, PAIR1, PAIR0, "png")   ## LSAMP
   
   ## Determin replication direction for each expressed gene
   slopes <- diff(smooth.spline(rpkms.chr.rt$MEDIAN)$y)/diff((bed.gc.chr$START)/1E7)   ## WHY?
 
   ensGene.tx.chr <- subset(ensGene.tx, chromosome_name == chr)
   ensGene.tx.chr$SLOPE_START <- NA
   ensGene.tx.chr$SLOPE_START <- NA
   for (g in 1:nrow(ensGene.tx.chr)) {
      gene <- ensGene.tx.chr[g,]
      bed.s <- getEnsGeneBED(gene$start_position, bed.gc.chr)
      bed.e <- getEnsGeneBED(gene$end_position, bed.gc.chr)
      
      if (length(bed.s) != 0) ensGene.tx.chr$SLOPE_START[g] <- slopes[which(rownames(bed.gc.chr) == bed.s[1])]
      if (length(bed.e) != 0) ensGene.tx.chr$SLOPE_END[g] <- slopes[which(rownames(bed.gc.chr) == bed.e[1])]
   }
   ensGene.tx.rt <- rbind(ensGene.tx.rt, ensGene.tx.chr)
}
save(ensGene.tx.rt, file=file.path(wd.asym.data, paste0(base, "_asym_tx_rt_bstrp1000.RData")))
# > nrow(ensGene.tx.rt)   ## All genes
# [1] 18440
# > nrow(ensGene.tx.rt)   ## All pcg genes
# [1] 16410
# > nrow(ensGene.tx.rt)   ## All non-pcg genes
# [1] 10604

