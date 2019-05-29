# =============================================================================
# Library      : Replication Timing
# Name         : handbook-of/ReplicationTiming.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 05/03/19; 15/11/18
# =============================================================================

# -----------------------------------------------------------------------------
# Methods: Calculate ratio of raw sequencing reads between tumour and matched normal (in 1a_cmd-rt_bam.rpkm.R)
# Last Modified: 15/11/18
# -----------------------------------------------------------------------------
initRPKM <- function(bam, bed, pair) {
   overlaps <- intersect(rownames(bam), rownames(bed))
   bam.o <- bam[overlaps,]   ## ADD 23/06/17: Autosomal-only RPKMs
   bed.o <- bed[overlaps,]
   bed.o$LENGTH <- bed.o$END - bed.o$START + 1   ## ADD 18/11/18: E.g. P2 chr1 10001 11000 0.646000
   
   rpkm <- bam.o[,c("BED", paste0("BAM_", pair))]
   TOTAL <- sum(as.numeric(rpkm$BAM))   ## BUG FIX 17/05/17: Warning message: In sum(rpkm$BAM) : integer overflow - use sum(as.numeric(.))
   rpkm$RPKM <- rpkm$BAM / TOTAL * 1E6 * bed.o$LENGTH   ## ADD 15/11/18: Replace to "bed.o$LENGTH" from "1000"
 
   return(rpkm)
}

# -----------------------------------------------------------------------------
# Methods: Calculate copy number-, GC-corrected read counts (in 1b_cmd-rt_bam.rpkm.corr.gc.R)
# Last Modified: 15/05/17
# -----------------------------------------------------------------------------
initRPKMCORRGC <- function(rpkm) {
   colnames <- c("BED", "BAM", "RPKM", "Ratio", "RPKM_CORR", "RPKM_CORR_GC")
   rpkm.corr.gc <- toTable(NA, length(colnames), nrow(rpkm), colnames)
   rpkm.corr.gc[,1:3] <- rpkm
 
   return(rpkm.corr.gc)
}

getBEDFromSegment <- function(bed.gc, seg.gc) {
   #bed.gc.chr <- subset(bed.gc, CHR == seg.gc$CHR)                       ## Using subset() takes ~15 MINS on 1,040 segments over 2,861,558 entries (S00022)
   #bed.gc.chr.end <- subset(bed.gc.chr, END > seg.gc$START)              ## Important! Not >= seg$START if e.g. P2 chr1 10000 11000 0.646000
   #bed.gc.chr.end.start <- subset(bed.gc.chr.end, START <= seg.gc$END)
   bed.gc.chr <- bed.gc[which(bed.gc$CHR == seg.gc$CHR),]                 ## Using which() takes ~5 MINS on 1,040 segments over 2,861,558 entries (S00022)
   bed.gc.chr.end <- bed.gc.chr[which(bed.gc.chr$END >= seg.gc$START),]   ## Important! Not >= seg$START if e.g. P2 chr1 10000 11000 0.646000
   bed.gc.chr.end.start <- bed.gc.chr.end[which(bed.gc.chr.end$START <= seg.gc$END),]
 
   return(bed.gc.chr.end.start)
}

setRPKMCORR <- function(rpkm.corr.gc, segs.gc, bed.gc) {
   for (s in 1:nrow(segs.gc)) {
      seg.gc <- segs.gc[s,]
      bed.gc.seg <- getBEDFromSegment(bed.gc, seg.gc)
      bed.gc.seg.idx <- which(rpkm.corr.gc$BED %in% bed.gc.seg$BED)
  
      rpkm.corr.gc$Ratio[bed.gc.seg.idx] <- seg.gc$Ratio
   }

   return(rpkm.corr.gc[which(rpkm.corr.gc$Ratio != 0),])   ## ADD 15/05/17: No information from chrY in female (e.g. S00035)
}                                                          ## 2852532 P2868190    0    0.0000000     0       NaN          NaN
                                                           ## 2852533 P2868191    1    0.9901062     0       Inf          Inf
## Main function (in 1b_cmd-rt_bam.rpkm.corr.gc.R)
getRPKMCORRGC <- function(rpkm, segs, bed, PAIR, CORR) {
   bed$BED  <- rownames(bed)                        ## To seed up?
   overlaps <- intersect(rownames(rpkm), bed$BED)   ## ADD 24/06/17: Double-check if they have the same rows
   rpkm.gc <- rpkm[overlaps,]
   bed.gc  <- bed[overlaps,]
   
   gc.mean <- mean(bed.gc$GC)
   rpkm.corr.gc <- initRPKMCORRGC(rpkm.gc)
   if (CORR) {   ## ADD 02/07/17: No need to correct for copy numbers in the normals (e.g. LCLs)
      segs.gc <- subset(segs, CHR %in% unique(bed.gc$CHR))   ## ADD 24/06/17: If bed.gc.au, remove chrXY from segs accordingly
    
      rpkm.corr.gc <- setRPKMCORR(rpkm.corr.gc, segs.gc, bed.gc)
      rpkm.corr.gc$RPKM_CORR    <- rpkm.corr.gc$RPKM / (rpkm.corr.gc$Ratio / 2)  
      rpkm.corr.gc$RPKM_CORR_GC <- rpkm.corr.gc$RPKM_CORR / (bed.gc$GC / gc.mean)
   } else {
      rpkm.corr.gc$RPKM_CORR_GC <- rpkm.corr.gc$RPKM / (bed.gc$GC / gc.mean)
   }
    
   return(rpkm.corr.gc)
}

# -----------------------------------------------------------------------------
# Methods: Gather coverages (per chromsome) from all samples (in 1c_cmd-rt_bam.rpkm.corr.gc_chr.R)
# Last Modified: 15/05/17
# -----------------------------------------------------------------------------
initReadDepthPerChromosome <- function(samples, bed.gc.chr) {
   colnames <- samples
   rt.chr <- toTable(NA, length(colnames), length(bed.gc.chr), colnames)
   rownames(rt.chr) <- bed.gc.chr
 
   return(rt.chr)
}

# -----------------------------------------------------------------------------
# Methods: Bootstraps (in 2a_cmd-rt_rpkm.corr.gc.d_bstrp.R)
# Last Modified: 15/05/17
# -----------------------------------------------------------------------------
getDetectedRD <- function(rpkms) {   ## Find not dected (RPKM_CORR_GC = 0) windows in any of the samples 
   return(mapply(x = 1:nrow(rpkms), function(x) !any(as.numeric(rpkms[x, -1]) == 0)))   ## ADD 15/02/19; To skip the first column "BED"
}

## Read in rpkm.corr.gc.txt.gz and getDetectedRD()
pipeGetDetectedRD <- function(wd.ngs.data, BASE, chr, PAIR) {
   rpkms.chr <- readTable(file.path(wd.ngs.data, paste0(tolower(BASE), "_rpkm.corr.gc_", chr, "_", PAIR, ".txt.gz")), header=T, rownames=T, sep="")##[, samples]   ## REMOVED 15/02/19; if length(samples) == 1
   rpkms.chr.d <- rpkms.chr[getDetectedRD(rpkms.chr),]   ## ADD 13/06/17; getDetectedRD()
 
   return(rpkms.chr.d)
}

getLog2RDRatio <- function(rpkms.T.chr, rpkms.N.chr, pseudocount) {
   return(log2(as.numeric(rpkms.T.chr) + pseudocount) - log2(as.numeric(rpkms.N.chr) + pseudocount))
}

outputRT <- function(rpkms.chr) {
   samples <- colnames(rpkms.chr)
   rpkms.chr$BED <- rownames(rpkms.chr)
 
   return(rpkms.chr[,c("BED", samples)])
}

# -----------------------------------------------------------------------------
# Plot RD and RT in sclc-wgs-rt.R (also refer to plotBootstrapsRT)
# Link(s): http://www.mun.ca/biology/scarr/2250_DNA_replication_&_transcription.html
#          https://stackoverflow.com/questions/43615469/how-to-calculate-the-slope-of-a-smoothed-curve-in-r
# Last Modified: 14/02/19
# -----------------------------------------------------------------------------
## Colours (was "lightcoral", "skyblue3")
## http://r.789695.n4.nabble.com/plot-function-color-transparency-td4682424.html
adjustcolor.red  <- adjustcolor("lightcoral", alpha.f=0.3)
adjustcolor.blue <- adjustcolor("skyblue3", alpha.f=0.3)
adjustcolor.gray <- adjustcolor("gray", alpha.f=0.3)

setScaledRT <- function(rpkms.chr.rt, pseudocount, recaliRT, scaledRT) {
   rpkms.chr.rt$T <- log2(rpkms.chr.rt$T + pseudocount)
   rpkms.chr.rt$N <- log2(rpkms.chr.rt$N + pseudocount)

   if (recaliRT == T)
      rpkms.chr.rt$RT <- rpkms.chr.rt$T - rpkms.chr.rt$N   ## ADD 24/02/19

   if (scaledRT == T)
      rpkms.chr.rt$RT <- scale(rpkms.chr.rt$RT)               ## ADD 19/02/19
   
   return(rpkms.chr.rt)
}

setSpline <- function(rpkms.chr.rt, bed.gc.chr, column) {
   bed.gc.chr <- bed.gc.chr[rpkms.chr.rt$BED,]
 
   spline <- smooth.spline(x=bed.gc.chr$START, y=rpkms.chr.rt[, column])
   rpkms.chr.rt$SPLINE <- spline$y
   #slopes <- diff(spline$y)/diff(bed.gc.chr$START/1E6)   ## ADD 31/10/18
   #rpkms.chr.rt <- rpkms.chr.rt[-nrow(rpkms.chr.rt),]    ## length(slopes) is 1 less than nrow(bed.gc.chr), as no slope for the last 1kb window
   #rpkms.chr.rt$SLOPE <- slopes
 
   sizes <- diff(bed.gc.chr$START)
   gaps <- which(sizes != 1000)
   return(rpkms.chr.rt[-gaps, c("BED", column, "SPLINE")])
}

plotRT <- function(file.name, main.text, chr, xmin, xmax, rpkms.chr.rt, bed.gc.chr, colours, legends, colours2, legends2, ext, width, peaks, ymin=NA, ymax=NA, cutoff, scale) {
   xlab.text <- paste0("Chromosome ", gsub("chr", "", chr), " coordinate [Mb]")
   if (!is.na(xmin) && !is.na(xmax)) file.name <- paste0(file.name, "_", xmin/1E6, "-", xmax/1E6, "Mb")
   if (is.na(xmin)) xmin <- 0
   if (is.na(xmax)) xmax <- subset(chromInfo, chrom == chr)$size

   if (ext == "pdf") {
      pdf(paste0(file.name, ".pdf"), height=5, width=width)
   } else if (ext == "png")
      png(paste0(file.name, ".png"), height=5, width=width, units="in", res=300)   ## ADD 16/05/17: res=300
   
   ###
   ## Initiate RD plot
   layout(matrix(c(1,2), 2, 1), widths=1, heights=c(1,1))   ## One figure each in row 1 and row 2   ## See plotBootstrapsHist()
   par(mar=c(1,4,4,1))
   ylab.text <- "Read depth [log2]"
   if (is.na(ymin) || is.na(ymax)) {
      plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=c(min(rpkms.chr.rt$T), max(rpkms.chr.rt$T)), xlab=xlab.text, ylab=ylab.text, main=main.text, xaxt="n")
   } else
      plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=c(ymin, ymax), xlab=xlab.text, ylab=ylab.text, main=main.text, xaxt="n")
   #points(bed.gc.chr$START/1E6, rpkms.chr.rt, col=colours[1], cex=0.3)
   abline(h=0, lwd=0.5, col="lightgrey")
 
   ## Plot cytobands (before smoothing spline)
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.4, col="lightgrey") 
   
   ## Plot smoothing splines
   rpkms.chr.rt.T  <- setSpline(rpkms.chr.rt, bed.gc.chr, "T")
   rpkms.chr.rt.N  <- setSpline(rpkms.chr.rt, bed.gc.chr, "N")
   rpkms.chr.rt.RT <- setSpline(rpkms.chr.rt, bed.gc.chr, "RT")
   bed.gc.chr <- bed.gc.chr[rownames(rpkms.chr.rt.RT),]
   
   #spline <- smooth.spline(x=bed.gc.chr$START, y=rpkms.chr.rt$T)
   points(bed.gc.chr$START/1E6, rpkms.chr.rt.T$SPLINE, col=colours[1], pch=16, cex=0.2)
   #spline <- smooth.spline(x=bed.gc.chr$START, y=rpkms.chr.rt$N)   ## TO-DO: Change it back to T and N
   points(bed.gc.chr$START/1E6, rpkms.chr.rt.N$SPLINE, col=colours[2], pch=16, cex=0.2)
   
   ## Plot legend and peaks
   if (xmin == 0) {
      legend("bottomright", legends, col=colours, lty=1, lwd=2, bty="n", horiz=T)
   } else
      legend("bottomright", legends2, col=colours, lty=1, lwd=2, bty="n", horiz=T)
    
   if (length(peaks) != 0)
      for (p in 1:length(peaks))
         abline(v=peaks[p]/1E6, lty=5, lwd=1, col="black")

   ### 
   ## Initiate RT plot
   par(mar=c(5.5,4,0,1))
   ylab.text <- "Replication timing"
   
   plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=c(-cutoff, cutoff), xlab=xlab.text, ylab=ylab.text, main="", yaxt="n")
   idx <- which(rpkms.chr.rt.RT$RT == 0)
   points(bed.gc.chr[idx,]$START/1E6, rpkms.chr.rt.RT[idx,]$RT, col="lightgrey", cex=0.3)
   idx <- which(rpkms.chr.rt.RT$RT > 0)
   points(bed.gc.chr[idx,]$START/1E6, rpkms.chr.rt.RT[idx,]$RT, col=colours2[1], cex=0.3)
   idx <- which(rpkms.chr.rt.RT$RT < 0)
   points(bed.gc.chr[idx,]$START/1E6, rpkms.chr.rt.RT[idx,]$RT, col=colours2[2], cex=0.3)
   
   abline(h=0, lwd=0.5, col="lightgrey")
   axis(side=2, at=seq(-scale, scale, by=scale), labels=c(-scale, 0, scale))
  
   ## Plot cytobands (before smoothing spline)
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.4, col="lightgrey") 
  
   ## Plot smoothing spline
   #spline <- smooth.spline(x=bed.gc.chr$START, y=rpkms.chr.rt$RT)
   points(bed.gc.chr$START/1E6, rpkms.chr.rt.RT$SPLINE, col="black", pch=16, cex=0.2)
   
   ## Plot legend and peaks
   if (xmin == 0) {
      legend("topright", paste0(legends2[1], "/", legends2[2], " read depth ratio"), col="black", lty=1, lwd=2, bty="n", horiz=T)
      legend("topleft", "Early", bty="n", text.col="black")   
      legend("bottomleft", "Late", bty="n", text.col="black")
   } else {
      legend("topright", paste0(legends2[1], "/", legends2[2], " ratio"), col="black", lty=1, lwd=2, bty="n", horiz=T)
      legend("topleft",    paste0("Early: ", legends2[1], " > ", legends2[2]), bty="n", text.col=colours[1])   
      legend("bottomleft", paste0("Late:  ", legends2[1], " < ", legends2[2]), bty="n", text.col=colours[2])
   }
   
   if (length(peaks) != 0)
      for (p in 1:length(peaks))
         abline(v=peaks[p]/1E6, lty=5, lwd=1, col="black")

   dev.off()
}

# -----------------------------------------------------------------------------
# Compare RD and RT in sclc-wgs-rt.R
# Last Modified: 05/03/19
# -----------------------------------------------------------------------------
#setSlope <- function(rpkms.chr.rt, bed.gc.chr, column) {
#   bed.gc.chr <- bed.gc.chr[rpkms.chr.rt$BED,]
# 
#   spline <- smooth.spline(x=bed.gc.chr$START, y=rpkms.chr.rt[, column])
#   slopes <- diff(spline$y)/diff(bed.gc.chr$START/1E6)   ## ADD 31/10/18
#   rpkms.chr.rt <- rpkms.chr.rt[-nrow(rpkms.chr.rt),]    ## length(slopes) is 1 less than nrow(bed.gc.chr), as no slope for the last 1kb window
#   rpkms.chr.rt$SLOPE <- slopes
# 
#   sizes <- diff(bed.gc.chr$START)
#   gaps <- which(sizes != 1000)
#   return(rpkms.chr.rt[-gaps, c("BED", column, "SLOPE")])
#}

plotRD2vsRT <- function(reads1, reads0, timings, file.name, main.text, ylab.text, xlab.text, colours, legends, xmin, xmax, ymin, ymax, line0) {
   #png(paste0(file.name, ".png"), height=5, width=5, units="in", res=300)
   pdf(paste0(file.name, ".pdf"), height=4.5, width=4.5)
   plot(NULL, xlim=c(xmin, xmax), ylim=c(ymin, ymax), ylab=ylab.text, xlab=xlab.text, main=main.text)
   if (line0)
      abline(v=0, lty=5)
   
   cors <- c()
   lm.fit <- lm(reads1 ~ timings)
   abline(lm.fit, col=colours[1], lwd=3)
   cor <- cor.test(reads1, timings, method="pearson")$estimate
   cors <- c(cors, cor)
 
   lm.fit <- lm(reads0 ~ timings)
   abline(lm.fit, col=colours[2], lwd=3)
   cor <- cor.test(reads0, timings, method="pearson")$estimate
   cors <- c(cors, cor)
 
   if (cors[1] < 0 && cors[2] < 0) {
      legend("bottomright", c(paste0("r = -", round0(abs(cors[1]), digits=2), " (", legends[1], " vs. RT)"), paste0("r = -", round0(abs(cors[2]), digits=2), " (", legends[2], " vs. RT)")), text.col=colours, bty="n", cex=1.1)
   } else if (as.numeric(cors[1]) > 0 && as.numeric(cors[2]) < 0) {
      legend("topright", paste0("r = ", round0(cors[1], digits=2), " (", legends[1], " vs. RT)"), text.col=colours[1], bty="n", cex=1.1)        
      legend("bottomright", paste0("r = -", round0(abs(cors[2]), digits=2), " (", legends[2], " vs. RT)"), text.col=colours[2], bty="n", cex=1.1)
   }
   
   mtext("Pearson correlation", cex=1.2, line=0.3)
   dev.off()
}

getCor <- function(reads, timings) {
   cor <- cor.test(reads, timings, method="pearson")$estimate
 
   return(cor)
}

# -----------------------------------------------------------------------------
# Compare betweeen RT and LCL RT in sclc-wgs-rt.R
# Last Modified: 05/03/19
# -----------------------------------------------------------------------------
plotRTvsRT <- function(reads, timings, file.name, main.text, ylab.text, xlab.text, colours) {
   lm.fit <- lm(reads ~ timings)
   #r2 <- summary(lm.fit)$r.squared
 
   png(paste0(file.name, ".png"), height=4.9, width=4.9, units="in", res=300)
   #pdf(paste0(file.name, ".pdf"), height=5, width=5)
   plot(reads ~ timings, ylab=ylab.text, xlab=xlab.text, main=main.text, col=colours[1])
   #plot(reads ~ timings, ylab=ylab.text, xlab=xlab.text, main=main.text, col="white")
   abline(lm.fit, col=colours[2], lwd=3)
   abline(h=0, lty=5)
   abline(v=0, lty=5)
   #mtext(paste0("R^2=", paste0(round0(r2*100, digits=2), "%")), cex=1.2, line=0.3)
 
   #rho <- cor.test(reads, timings, method="spearman", exact=F)[[4]]
   #mtext(paste0("SRC's rho = ", round0(rho, digits=2)), cex=1.2, line=0.3)
 
   cor <- cor.test(reads, timings, method="pearson")$estimate
   mtext(paste0("Pearson's r = ", round0(cor, digits=2)), cex=1.1, line=0.3) 
   dev.off()
}

plotRTvsRTALL <- function(cors, file.name, main.text, ylab.text, xlab.text, ymin, ymax, col, c=NA) {
   #png(paste0(file.name, ".png"), height=5, width=5, units="in", res=300)
   pdf(paste0(file.name, ".pdf"), height=4.5, width=4.5)
   plot(cors$cor ~ cors$chr, ylim=c(ymin, ymax), ylab=ylab.text, xlab=xlab.text, main=main.text, col=col, xaxt="n", yaxt="n", pch=19)
   lines(cors$cor, y=NULL, type="l", lwd=3, col=col)
   abline(h=0, lty=5)
 
   if (!is.na(c))
      text(c+2, cors$cor[c], paste0(round0(cors$cor[c], digits=2), " (Chr", c, ")"), cex=1.1, col=col, pos=3, offset=1)
   axis(side=1, at=seq(2, 22, by=2))
   axis(side=2, at=seq(-1, 1, by=0.2), labels=c(-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1))
   dev.off()
}

plotRD2vsRTALL <- function(cors, file.name, main.text, ylab.text, xlab.text, ymin, ymax, cols, c=NA) {
   #png(paste0(file.name, ".png"), height=5, width=5, units="in", res=300)
   pdf(paste0(file.name, ".pdf"), height=4.5, width=4.5)
   plot(cors$cor1 ~ cors$chr, ylim=c(ymin, ymax), ylab=ylab.text, xlab=xlab.text, main=main.text, col=cols[1], xaxt="n", yaxt="n", pch=19)
   lines(cors$cor1, y=NULL, type="l", lwd=3, col=cols[1])
   points(cors$chr, cors$cor2, col=cols[2], pch=19)
   lines(cors$cor2, y=NULL, type="l", lwd=3, col=cols[2])
   abline(h=0, lty=5)
 
   if (!is.na(c)) {
      text(c+2, cors$cor1[c], paste0(round0(cors$cor1[c], digits=2), " (Chr", c, ")"), cex=1.1, col=cols[1], pos=3, offset=1)
      text(c+2, cors$cor2[c], paste0(round0(cors$cor2[c], digits=2), " (Chr", c, ")"), cex=1.1, col=cols[2], pos=3)
   }
   axis(side=1, at=seq(2, 22, by=2))
   axis(side=2, at=seq(-1, 1, by=0.2), labels=c(-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1))
   dev.off()
}

plotSAMPLEvsRTALL <- function(cors.samples, samples1, file.name, main.text=NA, ymin=NA, ymax=NA) {
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
 
   pdf(paste0(file.name, "_boxplot_spline.pdf"), height=5.1, width=5.1)
   boxplot(cor ~ chr, data=cors.samples.plot, ylim=c(ymin, ymax), ylab="Pearson's r", xlab="Chromosome", outline=T, xaxt="n", main=main.text[1])
   axis(side=1, at=seq(2, 22, by=2))
   abline(h=0, lty=5)
   dev.off()
 
   png(paste0(file.name, "_scatter_spline.png"), height=5, width=5.5, units="in", res=300)
   plot(cors.samples$mean, log(cors.samples$cv2), ylab="log(cv2)", xlab="mean", main=main.text[2])
   text(cors.samples$mean[2], log(cors.samples$cv2[2]), "chr2", cex=1.2, pos=3)
   text(cors.samples$mean[4], log(cors.samples$cv2[4]), "chr4", cex=1.2, pos=3)
   text(cors.samples$mean[20], log(cors.samples$cv2[20]), "chr20", cex=1.2, pos=3)
   text(cors.samples$mean[12], log(cors.samples$cv2[12]), "chr12", cex=1.2, pos=3)
   dev.off()
 
   png(paste0(file.name, "_var_spline.png"), height=5, width=5.5, units="in", res=300)
   #plot(cors.samples$var ~ cors$chr, ylim=c(ymin, ymax), ylab=ylab.text, xlab=xlab.text, main=main.text, col="black", xaxt="n", yaxt="n", pch=19)
   plot(cors.samples$var ~ cors.samples$chr, ylab="var", xlab="Chromosome", col="black", xaxt="n", pch=19, main=main.text[2])
   lines(cors.samples$var, y=NULL, type="l", lwd=3)
   axis(side=1, at=seq(2, 22, by=2))
   dev.off()
}

# =============================================================================
# Inner Class  : PeifLyne File Reader
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 01/05/17
# =============================================================================
read.peiflyne.cn.txt <- function(cn.file) {
   cn <- readTable(cn.file, header=F, rownames=T, sep="")
   colnames(cn) <- c("BED", "CHR", "START", "END", "V5", "V6", "V7", "BAM_T", "BAM_N", "V10", "INSERT_SIZE_T", "INSERT_SIZE_N")
 
   return(cn)
}

read.peiflyne.cn.seg <- function(cn.file) {
   cn <- readTable(cn.file, header=F, rownames=F, sep="")
   colnames(cn) <- c("SAMPLE", "CHR", "START", "END", "V5", "Ratio")
 
   return(cn)
}

# =============================================================================
# Inner Class  : Bedtools File Reader
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 10/05/17
# =============================================================================
read.bedtools.multicov.cov <- function(cov.file) {
   cov <- readTable(cov.file, header=F, rownames=F, sep="")
   colnames(cov) <- c("CHR", "START", "END", "BAM_T", "BAM_N")
   cov$BED <- paste0("P", 1:nrow(cov))
   rownames(cov) <- cov$BED
   
   return(cov)
}
