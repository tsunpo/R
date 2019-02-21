# =============================================================================
# Library      : Replication Timing
# Name         : handbook-of/ReplicationTiming.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 15/11/18
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
setScaledRT <- function(rpkms.chr.rt, pseudocount, scaled) {
   rpkms.chr.rt$T <- log2(rpkms.chr.rt$T + pseudocount)
   rpkms.chr.rt$N <- log2(rpkms.chr.rt$N + pseudocount)
   if (scaled == T)
      rpkms.chr.rt$RT <- scale(rpkms.chr.rt$RT)   ## ADD 19/02/19
   
   return(rpkms.chr.rt)
}

# https://www.r-bloggers.com/cubic-and-smoothing-splines-in-r/
plotRD <- function(file.name, main.text, ylab.text, chr, xmin, xmax, rpkms.chr.rt, bed.gc.chr, colors, ext, ymin=NA, ymax=NA) {
   xlab.text <- paste0("Chromosome ", gsub("chr", "", chr), " coordinate (Mb)")
   if (!is.na(xmin) && !is.na(xmax)) file.name <- paste0(file.name, "_", xmin/1E6, "-", xmax/1E6, "Mb")
   if (is.na(xmin)) xmin <- 0
   if (is.na(xmax)) xmax <- subset(chromInfo, chrom == chr)$size
 
   ## Initiate plot
   if (ext == "pdf") {
      pdf(paste0(file.name, ".pdf"), height=4, width=10)
   } else if (ext == "png")
      png(paste0(file.name, ".png"), height=4, width=10, units="in", res=300)   ## ADD 16/05/17: res=300
 
   if (is.na(ymin) || is.na(ymax))
      plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=c(min(rpkms.chr.rt), max(rpkms.chr.rt)), xlab=xlab.text, ylab=ylab.text, main=main.text)
   else
      plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=c(ymin, ymax), xlab=xlab.text, ylab=ylab.text, main=main.text)
   points(bed.gc.chr$START/1E6, rpkms.chr.rt, col=colors[1], cex=0.3)
   abline(h=0, lwd=0.5, col="lightgrey")
 
   ## Plot smoothing spline
   spline <- smooth.spline(x=bed.gc.chr$START, y=rpkms.chr.rt, cv=T)
   points(bed.gc.chr$START/1E6, spline$y, col=colors[2], pch=16, cex=0.2)   ## CHANGED 21/02/19; From smooth.spline(rpkms.chr.rt$RT)$y
 
   ## Plot cytobands and legend
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.4, col="lightgrey")
 
   dev.off()
}

plotRD2 <- function(file.name, main.text, ylab.text, chr, xmin, xmax, rpkms.chr.rt, bed.gc.chr, colors, legends, ext, ymin=NA, ymax=NA) {
   xlab.text <- paste0("Chromosome ", gsub("chr", "", chr), " coordinate (Mb)")
   if (!is.na(xmin) && !is.na(xmax)) file.name <- paste0(file.name, "_", xmin/1E6, "-", xmax/1E6, "Mb")
   if (is.na(xmin)) xmin <- 0
   if (is.na(xmax)) xmax <- subset(chromInfo, chrom == chr)$size

   ## Initiate plot
   if (ext == "pdf") {
      pdf(paste0(file.name, ".pdf"), height=4, width=10)
   } else if (ext == "png")
      png(paste0(file.name, ".png"), height=4, width=10, units="in", res=300)   ## ADD 16/05/17: res=300
 
   if (is.na(ymin) || is.na(ymax))
      plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=c(min(rpkms.chr.rt$T), max(rpkms.chr.rt$T)), xlab=xlab.text, ylab=ylab.text, main=main.text)
   else
      plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=c(ymin, ymax), xlab=xlab.text, ylab=ylab.text, main=main.text)
   #points(bed.gc.chr$START/1E6, rpkms.chr.rt, col=colors[1], cex=0.3)
   abline(h=0, lwd=0.5, col="lightgrey")
 
   ## Plot cytobands (before smoothing spline)
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.4, col="lightgrey") 
   
   ## Plot smoothing spline
   spline <- smooth.spline(x=bed.gc.chr$START, y=rpkms.chr.rt$T, cv=T)
   points(bed.gc.chr$START/1E6, spline$y, col=colors[1], pch=16, cex=0.2)
 
   spline <- smooth.spline(x=bed.gc.chr$START, y=rpkms.chr.rt$N, cv=T)
   points(bed.gc.chr$START/1E6, spline$y, col=colors[2], pch=16, cex=0.2)
 
   ## Plot legend
   legend("bottomright", legends, col=colors, lty=1, lwd=2, bty="n", horiz=T)
   
   dev.off()
}

plotRD3 <- function(file.name, main.text, chr, xmin, xmax, rpkms.chr.rt, bed.gc.chr, colors, legends, colors2, legends2, ext, width, peaks, ymin=NA, ymax=NA, cutoff, scale) {
   xlab.text <- paste0("Chromosome ", gsub("chr", "", chr), " coordinate (Mb)")
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
   ylab.text <- "Read depth"
   if (is.na(ymin) || is.na(ymax))
      plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=c(min(rpkms.chr.rt$T), max(rpkms.chr.rt$T)), xlab=xlab.text, ylab=ylab.text, main=main.text, xaxt="n")
   else
      plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=c(ymin, ymax), xlab=xlab.text, ylab=ylab.text, main=main.text, xaxt="n")
   #points(bed.gc.chr$START/1E6, rpkms.chr.rt, col=colors[1], cex=0.3)
   abline(h=0, lwd=0.5, col="lightgrey")
 
   ## Plot cytobands (before smoothing spline)
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.4, col="lightgrey") 
   
   ## Plot smoothing splines
   spline <- smooth.spline(x=bed.gc.chr$START, y=rpkms.chr.rt$T, cv=T)
   points(bed.gc.chr$START/1E6, spline$y, col=colors[1], pch=16, cex=0.2)
   
   spline <- smooth.spline(x=bed.gc.chr$START, y=rpkms.chr.rt$N, cv=T)
   points(bed.gc.chr$START/1E6, spline$y, col=colors[2], pch=16, cex=0.2)

   ## Plot legend and peaks
   legend("bottomright", legends, col=colors, lty=1, lwd=2, bty="n", horiz=T)
   if (length(peaks) != 0)
      for (p in 1:length(peaks))
         abline(v=peaks[p]/1E6, lty=5, lwd=1, col="black")
 
   ### 
   ## Initiate RT plot
   par(mar=c(5.5,4,0,1))
   ylab.text <- "Replication timing"
   
   plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=c(-cutoff, cutoff), xlab=xlab.text, ylab=ylab.text, main="", yaxt="n")
   points(bed.gc.chr$START/1E6, rpkms.chr.rt$RT, col=colors2[1], cex=0.3)
   idx <- which(rpkms.chr.rt$RT < 0)
   points(bed.gc.chr[idx,]$START/1E6, rpkms.chr.rt[idx,]$RT, col=colors2[2], cex=0.3)
   
   abline(h=0, lwd=0.5, col="lightgrey")
   axis(side=2, at=seq(-scale, scale, by=scale), labels=c(-scale, 0, scale))
  
   ## Plot cytobands (before smoothing spline)
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.4, col="lightgrey") 
  
   ## Plot smoothing spline
   spline <- smooth.spline(x=bed.gc.chr$START, y=rpkms.chr.rt$RT, cv=T)
   points(bed.gc.chr$START/1E6, spline$y, col="black", pch=16, cex=0.2)
   
   ## Plot legend and peaks
   if (length(legends2) != 0) {
      legend("topright",    paste0(legends2[1], " > ", legends2[2]), bty="n")   
      legend("bottomright", paste0(legends2[1], " < ", legends2[2]), bty="n")
   }
   if (length(peaks) != 0)
      for (p in 1:length(peaks))
         abline(v=peaks[p]/1E6, lty=5, lwd=1, col="black")
   
   dev.off()
}

plotRT <- function(file.name, main.text, ylab.text, chr, xmin, xmax, rpkms.chr.rt, bed.gc.chr, colors, ext, cutoff, scale) {
   xlab.text <- paste0("Chromosome ", gsub("chr", "", chr), " coordinate (Mb)")
   if (!is.na(xmin) && !is.na(xmax)) file.name <- paste0(file.name, "_", xmin/1E6, "-", xmax/1E6, "Mb")
   if (is.na(xmin)) xmin <- 0
   if (is.na(xmax)) xmax <- subset(chromInfo, chrom == chr)$size
 
   ## Initiate plot
   if (ext == "pdf") {
      pdf(paste0(file.name, ".pdf"), height=4, width=10)
   } else if (ext == "png")
      png(paste0(file.name, ".png"), height=4, width=10, units="in", res=300)   ## ADD 16/05/17: res=300
 
   plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=c(-cutoff, cutoff), xlab=xlab.text, ylab=ylab.text, main=main.text, yaxt="n")
   points(bed.gc.chr$START/1E6, rpkms.chr.rt$RT, col=colors[1], cex=0.3)
   idx <- which(rpkms.chr.rt$RT < 0)
   points(bed.gc.chr[idx,]$START/1E6, rpkms.chr.rt[idx,]$RT, col=colors[2], cex=0.3)
   
   abline(h=0, lwd=0.5, col="lightgrey")
   axis(side=2, at=seq(-scale, scale, by=scale), labels=c(-scale, 0, scale))
   
   ## Plot cytobands (before smoothing spline)
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.4, col="lightgrey") 
   
   ## Plot smoothing spline
   spline <- smooth.spline(x=bed.gc.chr$START, y=rpkms.chr.rt$RT, cv=T)
   points(bed.gc.chr$START/1E6, spline$y, col="black", pch=16, cex=0.2)

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










# -----------------------------------------------------------------------------
# Methods: File manipulation and plotting
# Last Modified: 13/06/17
# -----------------------------------------------------------------------------
## This method is for plotRT()
getBEDFromPosition <- function(bed.gc, chr, start, end) {
 bed.gc.chr <- bed.gc[which(bed.gc$CHR == chr),]
 bed.gc.chr.end <- bed.gc.chr[which(bed.gc.chr$END >= start),]
 bed.gc.chr.end.start <- bed.gc.chr.end[which(bed.gc.chr.end$START <= end),]
 
 return(bed.gc.chr.end.start)
}

plotRTPerChromPerSample <- function(wd.rt.plots, chr, sample, rpkm.chr, bed.gc.chr, dao, ext) {
 main.text <- " read depth (copy number-, GC-corrected RPKM"
 xlab.text <- paste0("Chromosome ", gsub("chr", "", chr), " coordinate (Mb)")
 ylab.text <- "Replication timing"
 
 if (ext == "pdf")
  pdf(paste0(wd.rt.plots, "pdf/", chr, "/", sample, "_rpkm.corr.gc.q_", chr, ".pdf"), height=4, width=10)
 else if (ext == "png")
  png(paste0(wd.rt.plots, "png/", chr, "/", sample, "_rpkm.corr.gc.q_", chr, ".png"), height=4, width=10, units="in", res=300)   ## ADD 16/05/17: res=300
 
 plot(NULL, ylim=c(dao$ymin, dao$ymax), xlim=c(dao$xmin/1E6, dao$xmax/1E6), xlab=xlab.text, ylab=ylab.text, main=paste0(sample, main.text, ")"))
 points(bed.gc.chr$START/1E6, rpkm.chr, col="red", cex=0.3)
 lines(bed.gc.chr$START/1E6, smooth.spline(rpkm.chr)$y)
 
 dev.off()
}

removeNAinAnyRows <- function(rpkms.chr) {   ## Remove windows if FPKM = NA in any of the samples 
 rpkms.chr$KEEP <- mapply(x = 1:nrow(rpkms.chr), function(x) !any(is.na(rpkms.chr[x,])))
 rpkms.chr <- subset(rpkms.chr, KEEP == T)[,-ncol(rpkms.chr)]
 
 return(rpkms.chr)
}

plotRT_old <- function(wd.rt.plots, chr, sample, rpkms.chr, bed.gc, cytoBand.chr, xmin, xmax, ext) {
 main.text <- " read depth (copy number-, GC-corrected RPKM"
 xlab.text <- paste0("Chromosome ", gsub("chr", "", chr), " coordinate (Mb)")
 ylab.text <- "Read depth"
 
 rpkm.chr <- log2(rpkms.chr$MEDIAN_RPKM + 0.01)
 ymin <- min(as.numeric(rpkm.chr))
 ymax <- max(as.numeric(rpkm.chr))
 bed.gc.chr <- bed.gc[rownames(rpkms.chr),]
 
 if (ext == "pdf")
  pdf(paste0(wd.rt.plots, "pdf/sclc_rpkm.corr.gc.q_", chr, "_N.pdf"), height=4, width=10)
 else if (ext == "png")
  png(paste0(wd.rt.plots, "png/sclc_rpkm.corr.gc.q_", chr, "-", xmin, ":", xmax, "_N.png"), height=4, width=10, units="in", res=300)   ## ADD 16/05/17: res=300
 
 plot(NULL, ylim=c(ymin, ymax), xlim=c(xmin/1E6, xmax/1E6), xlab=xlab.text, ylab=ylab.text, main=paste0("Normal SCLC median", main.text, ")"))
 points(bed.gc.chr$START/1E6, rpkm.chr, col="red", cex=0.3)
 lines(bed.gc.chr$START/1E6, smooth.spline(rpkm.chr)$y)
 
 for (c in 1:nrow(cytoBand.chr))
  abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.5, col="lightgrey")
 
 dev.off()
}

plotRT2 <- function(wd.rt.plots, chr, sample, rpkm.chr, bed.gc.chr, dao, ext) {
 main.text <- " read depth (copy number-, GC-corrected RPKM"
 xlab.text <- paste0("Chromosome ", gsub("chr", "", chr), " coordinate (Mb)")
 ylab.text <- "Replication timing"
 
 if (ext == "pdf")
  pdf(paste0(wd.rt.plots, "pdf/sclc_rpkm.corr.gc.q.sd_", chr, "_N.pdf"), height=4, width=10)
 else if (ext == "png")
  png(paste0(wd.rt.plots, "png/sclc_rpkm.corr.gc.q.sd_", chr, "_N.png"), height=4, width=10, units="in", res=300)   ## ADD 16/05/17: res=300
 
 plot(NULL, ylim=c(dao$ymin, dao$ymax), xlim=c(dao$xmin/1E6, dao$xmax/1E6), xlab=xlab.text, ylab=ylab.text, main=paste0("Normal SCLC median", main.text, ")"))
 points(bed.gc.chr$START/1E6, rpkm.chr, col="red", cex=0.3)
 lines(bed.gc.chr$START/1E6, smooth.spline(rpkm.chr)$y)
 
 dev.off()
}

# =============================================================================
# Inner Class  : Collections of test/obsolete/deprecated methods
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified:
# =============================================================================
initRPKMRatio <- function(cn) {
   colnames <- c("RPKM_T", "RPKM_N", "RPKM_RATIO", "RPKM_LOG2")
   rr <- toTable(NA, length(colnames), nrow(cn), colnames)
   rt.rr <- cbind(cn[,c("Partition", "BAM_T", "BAM_N")], rr)
 
   TOTAL_T <- sum(cn$BAM_T)
   TOTAL_N <- sum(cn$BAM_N)
   rt.rr$RPKM_T <- cn$BAM_T / TOTAL_T * 1E6 * 1000
   rt.rr$RPKM_N <- cn$BAM_N / TOTAL_N * 1E6 * 1000
   rt.rr$RPKM_RATIO <- rt.rr$RPKM_T / rt.rr$RPKM_N
   rt.rr$RPKM_LOG2 <- log2(rt.rr$RPKM_T + 0.01) - log2(rt.rr$RPKM_N + 0.01)
 
   return(rt.rr)
}

removeZeroReads <- function(cn) {
   cn.cln <- subset(cn, BAM_T != 0)
   cn.cln <- subset(cn.cln, BAM_N != 0)
 
   return(cn.cln)
}
