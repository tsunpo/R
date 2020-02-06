# =============================================================================
# Library      : Replication Timing
# Name         : handbook-of/ReplicationTiming.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 01/08/19; 11/07/19; 05/03/19; 15/11/18
# =============================================================================

# -----------------------------------------------------------------------------
# Methods: Normalising raw sequencing reads (in cmd-rt_1a_rpkb.R)
# Last Modified: 01/08/19
# -----------------------------------------------------------------------------
initNRD <- function(rd, bed.gc, pair, method) {
   overlaps <- intersect(rownames(rd), rownames(bed.gc))
   rd.o  <- rd[overlaps,]
   bed.o <- bed.gc[overlaps,]
   bed.o$LENGTH <- bed.o$END - bed.o$START + 1   ## ADD 18/11/18: E.g. P2 chr1 10001 11000 0.646000
 
   nrd <- rd.o[,c("BED", paste0("RD_", pair))]    ## ADD 05/08/19 Normalised read depth (NRD)
   colnames(nrd) <- c("BED", "RD")
   
   if (method == "RPKM") {               ## ADD 05/08/19: Reads per kilobase billion
      TOTAL <- sum(as.numeric(nrd$RD))
      nrd$NRD <- nrd$RD / (TOTAL / 1E6) / (bed.o$LENGTH / 1E3)   ## ADD 15/11/18: Replace to "bed.o$LENGTH" from "1000"
   } else if (method == "MEAN") {
      MEAN <- mean(as.numeric(nrd$RD))   ## ADD 03/08/19
      nrd$NRD <- nrd$RD / MEAN
   }
   
   return(nrd)
}

# -----------------------------------------------------------------------------
# Methods: Calculate copy number-, GC-corrected read counts (in 1b_cmd-rt_bam.nrd.corr.gc.R)
# Last Modified: 15/05/17
# -----------------------------------------------------------------------------
initNRDGCCN <- function(nrd) {
   colnames <- c("BED", "RD", "NRD", "NRD_GC", "CN", "NRD_GC_CN")
   nrd.gc.cn <- toTable(NA, length(colnames), nrow(nrd), colnames)
   nrd.gc.cn[,1:3]  <- nrd
 
   return(nrd.gc.cn)
}

getBEDFromSegment <- function(bed.gc.o, seg.gc) {
   #bed.gc.chr <- subset(bed.gc, CHR == seg.gc$CHR)                       ## Using subset() takes ~15 MINS on 1,040 segments over 2,861,558 entries (S00022)
   #bed.gc.chr.end <- subset(bed.gc.chr, END > seg.gc$START)              ## Important! Not >= seg$START if e.g. P2 chr1 10000 11000 0.646000
   #bed.gc.chr.end.start <- subset(bed.gc.chr.end, START <= seg.gc$END)
   bed.gc.chr <- bed.gc.o[which(bed.gc.o$CHR == seg.gc$CHR),]                 ## Using which() takes ~5 MINS on 1,040 segments over 2,861,558 entries (S00022)
   bed.gc.chr.end <- bed.gc.chr[which(bed.gc.chr$END >= seg.gc$START),]   ## Important! Not >= seg$START if e.g. P2 chr1 10000 11000 0.646000
   bed.gc.chr.end.start <- bed.gc.chr.end[which(bed.gc.chr.end$START <= seg.gc$END),]
 
   return(bed.gc.chr.end.start)
}

setNRDCN <- function(nrd.gc.cn, segs.gc, bed.gc.o) {
   bed.gc.o$BED <- rownames(bed.gc.o)
 
   for (s in 1:nrow(segs.gc)) {
      seg.gc <- segs.gc[s,]
      bed.gc.seg <- getBEDFromSegment(bed.gc.o, seg.gc)
      bed.gc.seg.idx <- which(nrd.gc.cn$BED %in% bed.gc.seg$BED)
  
      nrd.gc.cn$CN[bed.gc.seg.idx] <- seg.gc$CN
   }

   return(nrd.gc.cn[which(nrd.gc.cn$CN != 0),])   ## ADD 15/05/17: No information from chrY in female (e.g. S00035)
}                                                        ## 2852532 P2868190    0    0.0000000     0       NaN          NaN
                                                         ## 2852533 P2868191    1    0.9901062     0       Inf          Inf
## Main function (in 1b_cmd-rt_rd.nrd.cn.R)
getNRDGCCN <- function(nrd, segs, bed.gc, PAIR, CN) {
   overlaps <- intersect(rownames(nrd), rownames(bed.gc))   ## ADD 24/06/17: Double-check if they have the same rows
   nrd.o <- nrd[overlaps,]
   bed.gc.o <- bed.gc[overlaps,]
   
   nrd.gc.cn <- initNRDGCCN(nrd.o)
   MEAN  <- mean(bed.gc.o$GC)
   nrd.gc.cn$NRD_GC    <- nrd.gc.cn$NRD / (bed.gc.o$GC / MEAN)
   nrd.gc.cn$NRD_GC_CN <- nrd.gc.cn$NRD_GC   ## ADD 02/07/17: No need to correct for copy numbers in the normals (e.g. LCLs)
   
   if (CN) {
      segs.gc <- subset(segs, CHR %in% unique(bed.gc.o$CHR))   ## ADD 24/06/17: If bed.gc.au, remove chrXY from segs accordingly
    
      nrd.gc.cn <- setNRDCN(nrd.gc.cn, segs.gc, bed.gc.o)
      nrd.gc.cn$NRD_GC_CN <- nrd.gc.cn$NRD_GC / (nrd.gc.cn$CN / 2) 
   }
   
   return(nrd.gc.cn)
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
getDetectedRD <- function(rds) {   ## Find not dected (RD_CN = 0) windows in any of the samples 
   return(mapply(x = 1:nrow(rds), function(x) !any(as.numeric(rds[x, -1]) == 0)))   ## ADD 15/02/19; To skip the first column "BED"
}

## Read in rpkm.corr.gc.txt.gz and getDetectedRD()
pipeGetDetectedRD <- function(wd.ngs.data, BASE, chr, PAIR, method) {
   nrds.chr <- readTable(file.path(wd.ngs.data, paste0(tolower(BASE), "_", method, ".gc.cn_", chr, "_", PAIR, ".txt.gz")), header=T, rownames=T, sep="")##[, samples]   ## REMOVED 15/02/19; if length(samples) == 1
   nrds.chr.d <- nrds.chr[getDetectedRD(nrds.chr),]   ## ADD 13/06/17; getDetectedRD()
 
   return(nrds.chr.d)
}

outputRT <- function(nrds.chr) {
   samples <- colnames(nrds.chr)
   nrds.chr$BED <- rownames(nrds.chr)
 
   return(nrds.chr[,c("BED", samples)])
}

# -----------------------------------------------------------------------------
# Plot RD and RT in sclc-wgs-rt.R (also refer to plotBootstrapsRT)
# Link(s): http://www.mun.ca/biology/scarr/2250_DNA_replication_&_transcription.html
#          https://stackoverflow.com/questions/43615469/how-to-calculate-the-slope-of-a-smoothed-curve-in-r
# Last Modified: 14/02/19
# -----------------------------------------------------------------------------
getLog2ScaledRT <- function(wd.rt.data, base, method, PAIR1, PAIR0, n1, n0, chrs, bed.gc, isFlip=F) {
   nrds <- toTable(NA, 4, 0, c("BED", "T", "N", "RT"))
   for (c in 1:22) {
      chr <- chrs[c]
      bed.gc.chr <- subset(bed.gc, CHR == chr)
      nrds.chr <- readTable(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
  
      nrds <- rbind(nrds, nrds.chr)
   }
   if (isFlip) {
      nrds$RT <- nrds$N / nrds$T
   }
   nrds$RT <- log2(nrds$RT)   ## MUY MUY IMPORTANTE!! 2019/10/10
   nrds$RT <- scale(nrds$RT)
   
   return(nrds)
}

setSpline <- function(nrds.chr, bed.gc.chr, column, kb=1, returnAll=F) {
   overlaps <- intersect(rownames(bed.gc.chr), nrds.chr$BED)
   bed.gc.chr.o <- bed.gc.chr[overlaps,]
   nrds.chr.o <- nrds.chr[overlaps,]
 
   spline <- smooth.spline(x=bed.gc.chr.o$START, y=nrds.chr.o[, column])
   nrds.chr.o$SPLINE <- spline$y
 
   ## https://stackoverflow.com/questions/43615469/how-to-calculate-the-slope-of-a-smoothed-curve-in-r
   slopes <- diff(spline$y)/diff(bed.gc.chr.o$START/1E6)   ## ADD 31/10/18
   nrds.chr.o$SLOPE <- NA
   nrds.chr.o$SLOPE[1:length(slopes)] <- slopes   ## length(slopes) is 1 less than nrow(bed.gc.chr.o), as no slope for the last 1kb window
   #nrds.chr.o <- nrds.chr.o[1:length(slopes),]   ## REMOVED 24/11/19
   #nrds.chr.o$SLOPE <- slopes
 
   sizes <- diff(bed.gc.chr.o$START/1E6)
   gaps <- which(sizes > kb)   ## gaps <- which(sizes != 1000)
   if (length(gaps) != 0)       ## For chr14
      nrds.chr.o <- nrds.chr.o[-gaps, ]
 
   if (returnAll)
      return(nrds.chr.o)
   return(nrds.chr.o[, c("BED", column, "SPLINE", "SLOPE")])
}

#setSpline <- function(nrds.chr, bed.gc.chr, column) {
#   overlaps <- intersect(rownames(bed.gc.chr), nrds.chr$BED)
#   bed.gc.chr.o <- bed.gc.chr[overlaps,]
#   nrds.chr.o <- nrds.chr[overlaps,]
#  
#   spline <- smooth.spline(x=bed.gc.chr.o$START, y=nrds.chr.o[, column])
#   nrds.chr.o$SPLINE <- spline$y
   
   ## https://stackoverflow.com/questions/43615469/how-to-calculate-the-slope-of-a-smoothed-curve-in-r
#   slopes <- diff(spline$y)/diff(bed.gc.chr.o$START/1E6)   ## ADD 31/10/18
#   nrds.chr.o$SLOPE <- NA
#   nrds.chr.o$SLOPE[1:length(slopes)] <- slopes   ## length(slopes) is 1 less than nrow(bed.gc.chr.o), as no slope for the last 1kb window
   #nrds.chr.o <- nrds.chr.o[1:length(slopes),]   ## REMOVED 24/11/19
   #nrds.chr.o$SLOPE <- slopes
   
#   sizes <- diff(bed.gc.chr.o$START)
#   gaps <- which(sizes != 1000)
#   return(nrds.chr.o[-gaps, c("BED", column, "SPLINE", "SLOPE")])
#}

getRT <- function(nrds, bed.gc) {
   nrds.RT <- NULL
   
   for (c in 1:22) {
      chr <- chrs[c]
      bed.gc.chr <- subset(bed.gc, CHR == chr)
      nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]
      
      nrds.chr.RT <- setSpline(nrds.chr, bed.gc.chr, "RT")
      nrds.chr.RT <- cbind(nrds.chr[rownames(nrds.chr.RT),], nrds.chr.RT[, c("SPLINE", "SLOPE")])
      if (is.null(nrds.RT))
         nrds.RT <- nrds.chr.RT
      else
         nrds.RT <- rbind(nrds.RT, nrds.chr.RT)
   }
   
   return(nrds.RT)
}

plotRT <- function(file.name, main.text, chr, xmin, xmax, nrds.chr, bed.gc.chr, colours, legends, colours2, legends2, ext, width, peaks, ylim=NULL, lcl.rt.chr=NULL, nrds.lcl.chr=NULL) {
   ## Colours (was "lightcoral", "skyblue3")
   ## http://r.789695.n4.nabble.com/plot-function-color-transparency-td4682424.html
   adjustcolor.red  <- adjustcolor(colours2[1], alpha.f=0.08)
   adjustcolor.blue <- adjustcolor(colours2[2], alpha.f=0.08)

   nrds.chr.T  <- setSpline(nrds.chr, bed.gc.chr, "T")
   nrds.chr.N  <- setSpline(nrds.chr, bed.gc.chr, "N")
   nrds.chr.RT <- setSpline(nrds.chr, bed.gc.chr, "RT")
   nrds.lcl.chr.RT <- setSpline(nrds.lcl.chr, bed.gc.chr, "RT")
   if (!is.null(nrds.lcl.chr))
    
   #nrds.chr.RT$SPLINE <- scale(nrds.chr.RT$SPLINE)
   #bed.gc.chr <- bed.gc.chr[rownames(nrds.chr.RT),]   ## NOT HERE?
 
   xlab.text <- paste0("Chromosome ", gsub("chr", "", chr), " coordinate [Mb]")
   if (!is.na(xmin) && !is.na(xmax)) file.name <- paste0(file.name, "_", xmin/1E6, "-", xmax/1E6, "Mb")
   if (is.na(xmin)) {
      start <- bed.gc.chr[rownames(nrds.chr)[1],]$START
      if (start < 5000000) xmin <- 0
      else xmin <- start
   }
   if (is.na(xmax)) xmax <- subset(chromInfo, chrom == chr)$size

   ## Initiation plot
   if (!is.null(lcl.rt.chr))
      file.name <- paste0(file.name, "_with-koren")
   if (!is.null(nrds.lcl.chr))
      file.name <- paste0(file.name, "_with-LCL")
   if (ext == "pdf") {
      pdf(paste0(file.name, ".pdf"), height=5, width=width)
   } else if (ext == "png")
      png(paste0(file.name, ".png"), height=5, width=width, units="in", res=300)   ## ADD 16/05/17: res=300
   
   ###
   ## Initiate RD plot
   layout(matrix(c(1,2), 2, 1), widths=1, heights=c(1,1))   ## One figure each in row 1 and row 2   ## See plotBootstrapsHist()
   par(mar=c(1,4,4,1))
   ylab.text <- "Read depth [RPKM]"
   if (is.null(ylim)) {
      rds <- c(nrds.chr.T$SPLINE, nrds.chr.N$SPLINE)
      plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=c(min(rds), max(rds)), xlab="", ylab=ylab.text, main=main.text, xaxt="n", cex.axis=1.05, cex.lab=1.08)
   } else
      plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=ylim, xlab="", ylab=ylab.text, main=main.text, xaxt="n", cex.axis=1.05, cex.lab=1.08)
   #points(bed.gc.chr$START/1E6, nrds.chr, col=colours[1], cex=0.3)
   abline(h=0, lwd=0.5, col="lightgrey")
 
   ## Plot cytobands (before smoothing spline)
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.4, col="lightgrey") 

   ## Plot smoothing splines
   bed.gc.chr.lcl <- NULL
   if (!is.null(nrds.lcl.chr))
      bed.gc.chr.lcl <- bed.gc.chr[rownames(nrds.lcl.chr.RT),]
   bed.gc.chr <- bed.gc.chr[rownames(nrds.chr.RT),]   ## BUT HERE?
   points(bed.gc.chr$START/1E6, nrds.chr.T$SPLINE, col=colours[1], pch=16, cex=0.2)
   points(bed.gc.chr$START/1E6, nrds.chr.N$SPLINE, col=colours[2], pch=16, cex=0.2)
   
   ## Plot legend and peaks
   legend("topright", legends, col=colours, lty=1, lwd=2, bty="n", horiz=T, cex=1.15)
   if (length(peaks) != 0)
      for (p in 1:length(peaks))
         abline(v=peaks[p]/1E6, lty=5, lwd=1, col="black")
   
   ### 
   ## Initiate RT plot
   par(mar=c(5.5,4,0,1))
   ylab.text <- "RT [log2]"
   
   plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=c(-2, 2), xlab=xlab.text, ylab=ylab.text, main="", yaxt="n", cex.axis=1.05, cex.lab=1.08, cex.main=1.25)
   idx <- which(nrds.chr.RT$RT == 0)
   points(bed.gc.chr[idx,]$START/1E6, nrds.chr.RT[idx,]$RT, col="lightgrey", cex=0.3)
   idx <- which(nrds.chr.RT$RT < 0)
   points(bed.gc.chr[idx,]$START/1E6, nrds.chr.RT[idx,]$RT, col=adjustcolor.blue, cex=0.3)
   idx <- which(nrds.chr.RT$RT > 0)
   points(bed.gc.chr[idx,]$START/1E6, nrds.chr.RT[idx,]$RT, col=adjustcolor.red, cex=0.3)
   
   ## Plot cytobands (before smoothing spline)
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.4, col="lightgrey") 
  
   ## Plot Koren 2012 (before smoothing spline)
   if (!is.null(lcl.rt.chr))
      points(lcl.rt.chr$POS/1E6, lcl.rt.chr$RT, col=colours[3], pch=16, cex=0.2)
   if (!is.null(nrds.lcl.chr)) {
      points(bed.gc.chr.lcl$START/1E6, nrds.lcl.chr.RT$SPLINE, col=colours[3], pch=16, cex=0.2)
   }
   ## Plot smoothing spline
   points(bed.gc.chr$START/1E6, nrds.chr.RT$SPLINE, col="black", pch=16, cex=0.2)
   abline(h=0, lty=5, lwd=1, col="black")
   axis(side=2, at=seq(-2, 2, by=4), labels=c("\u22122", 2), cex.axis=1.08)
   axis(side=2, at=seq(-1, 1, by=1), labels=c("\u22121", 0, 1), cex.axis=1.08)
   
   ## Plot legend and peaks
   legend("topleft", paste0("Early (", legends2[1], " > ", legends2[2], ")"), bty="n", text.col="black", pt.cex=0.9, pch=1, col=colours2[1], cex=1.15)   
   legend("bottomleft", paste0("Late (", legends2[1], " < ", legends2[2], ")"), bty="n", text.col="black", pt.cex=0.9, pch=1, col=colours2[2], cex=1.15)
   if (!is.na(xmin) && !is.na(xmax))
      legend("topright", paste0(legends2[1], "/", legends2[2], " read depth ratio"), col="black", lty=1, lwd=2, bty="n", horiz=T, cex=1.15)
   else
      legend("topright", paste0(legends2[1], "/", legends2[2], " ratio"), col="black", lty=1, lwd=2, bty="n", horiz=T, cex=1.15)
   if (!is.null(lcl.rt.chr))
      legend("bottomright", "Koren et al. 2012", col=colours[3], lty=1, lwd=2, bty="n", horiz=T, cex=1.15)
   if (!is.null(nrds.lcl.chr))
      legend("bottomright", "S/G1 read depth ratio", col=colours[3], lty=1, lwd=2, bty="n", horiz=T, cex=1.15)
   if (length(peaks) != 0)
      for (p in 1:length(peaks))
         abline(v=peaks[p]/1E6, lty=5, lwd=1, col="black")

   dev.off()
}

# -----------------------------------------------------------------------------
# Compare RD and RT in sclc-wgs-rt.R
# Last Modified: 05/03/19
# -----------------------------------------------------------------------------
getCor <- function(reads, timings, method) {
   if (method == "pearson") {
      return(cor.test(reads, timings, method="pearson")$estimate)
   } else if (method == "spearman") {
      return(cor.test(reads, timings, method="spearman", exact=F)[[4]])
   } else if (method == "linear") {
      lm.fit <- lm(reads ~ timings)
      return(summary(lm.fit)$r.squared)
   }
}

plotRDvsRT <- function(reads1, timings, file.name, main.text, ylab.text, xlab.text, colours, legends, method) {
   xlim <- c(min(timings), max(timings))
   ylim <- c(min(reads1), max(reads1))
   cor <- "rho"
   main.text2 <- "Spearman correlation"
   if (method == "pearson") {
      cor <- "r" 
      main.text2 <- "Pearson correlation"
   }
 
   png(paste0(file.name, ".png"), height=5, width=5, units="in", res=300)
   plot(reads1 ~ timings, xlim=xlim, ylim=ylim, ylab=ylab.text, xlab=xlab.text, main=main.text[1], xaxt="n", col=colours[2], cex.axis=1.1, cex.lab=1.2, cex.main=1.25)
   #pdf(paste0(file.name, ".pdf"), height=5, width=5)
   #plot(NULL, xlim=xlim, ylim=ylim, ylab=ylab.text, xlab=xlab.text, main=main.text[1], xaxt="n", cex.axis=1.1, cex.lab=1.2, cex.main=1.25)
   abline(v=0, lty=5)
 
   lm.fit1 <- lm(reads1 ~ timings)
   abline(lm.fit1, col=colours[1], lwd=3)
   cor1 <- getCor(reads1, timings, method)

   RT <- paste0(legends[1], "/", legends[2])
   if (cor1 > 0) {
      legend("topright", paste0(cor, " = ", round0(cor1, digits=2), " (", legends[1], " vs. ", RT, ")"), text.col=colours[1], bty="n", cex=1.2)        
   } else {
      legend("bottomright", paste0(cor, " = \u2212", round0(abs(cor1), digits=2), " (", legends[2], " vs. ", RT, ")"), text.col=colours[1], bty="n", cex=1.2)
   }
   axis(side=1, at=seq(-3, 3, by=0.5), labels=c(-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3), cex.axis=1.1)
   mtext(main.text[2], line=0.3, cex=1.2)
   dev.off()
}

plotRD2vsRT <- function(reads1, reads2, timings, file.name, main.text, ylab.text, xlab.text, colours, legends, method) {
   xlim <- c(min(timings), max(timings))
   ylim <- c(min(c(reads1, reads2)), max(c(reads1, reads2)))
   cor <- "rho"
   main.text2 <- "Spearman correlation"
   if (method == "pearson") {
      cor <- "r" 
      main.text2 <- "Pearson correlation"
   }

   #png(paste0(file.name, ".png"), height=5, width=5, units="in", res=300)
   pdf(paste0(file.name, ".pdf"), height=5, width=5)
   plot(NULL, xlim=xlim, ylim=ylim, ylab=ylab.text, xlab=xlab.text, main=main.text[1], xaxt="n", cex.axis=1.1, cex.lab=1.2, cex.main=1.25)
   abline(v=0, lty=5)
   
   lm.fit1 <- lm(reads1 ~ timings)
   abline(lm.fit1, col=colours[1], lwd=3)
   cor1 <- getCor(reads1, timings, method)

   lm.fit2 <- lm(reads2 ~ timings)
   abline(lm.fit2, col=colours[2], lwd=3)
   cor2 <- getCor(reads2, timings, method)

   RT <- paste0(legends[1], "/", legends[2])
   if (cor1 < 0 && cor2 < 0) {
      legend("bottomright", c(paste0(cor, " = ", round0(cor1, digits=2), " (", legends[1], " vs. ", RT, ")"), paste0(cor, " = ", round0(cor2, digits=2), " (", legends[2], " vs. ", RT, ")")), text.col=colours, bty="n", cex=1.2)
   } else if (as.numeric(cor1) > 0 && as.numeric(cor2) < 0) {
      legend("topright", paste0(cor, " = ", round0(cor1, digits=2), " (", legends[1], " vs. ", RT, ")"), text.col=colours[1], bty="n", cex=1.2)        
      legend("bottomright", paste0(cor, " = ", round0(cor2, digits=2), " (", legends[2], " vs. ", RT, ")"), text.col=colours[2], bty="n", cex=1.2)
   }
   axis(side=1, at=seq(-3, 3, by=0.5), labels=c(-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3), cex.axis=1.1)
   mtext(main.text[2], line=0.3, cex=1.2)
   dev.off()
}

plotRD2vsRTALL <- function(cors, file.name, main.text, ymin, ymax, cols, legends, c=NA) {
   ylab.text <- "Spearman's rho"
   xlab.text <- "Chromosome"
 
   #png(paste0(file.name, ".png"), height=5, width=5, units="in", res=300)
   pdf(paste0(file.name, ".pdf"), height=6, width=6)
   plot(cors$cor1 ~ cors$chr, ylim=c(ymin, ymax), ylab="", xlab=xlab.text, main=main.text, col=cols[1], xaxt="n", yaxt="n", pch=19, cex.lab=1.7, cex.main=1.8)
   lines(cors$cor1, y=NULL, type="l", lwd=3, col=cols[1])
   points(cors$chr, cors$cor2, col=cols[2], pch=19)
   lines(cors$cor2, y=NULL, type="l", lwd=3, col=cols[2])
   abline(h=0, lty=5)
   
   RT <- paste0(legends[1], "/", legends[2])
   legend("topright", paste0(legends[1], " vs. ", RT), text.col=cols[1], bty="n", cex=1.75)        
   legend("bottomright", paste0(legends[2], " vs. ", RT), text.col=cols[2], bty="n", cex=1.75)
   
   if (!is.na(c)) {
      text(c, cors$cor1[c], round0(cors$cor1[c], digits=2), col=cols[1], pos=3, cex=1.75)   ##, offset=1.3)
      text(c, cors$cor2[c], round0(cors$cor2[c], digits=2), col=cols[2], pos=3, cex=1.75)
   }
   axis(side=1, at=seq(2, 22, by=4), cex.axis=1.5)
   axis(side=1, at=seq(4, 20, by=4), cex.axis=1.5)
   axis(side=2, at=seq(-0.8, 0.8, by=0.4), labels=c(-0.8, -0.4, 0, 0.4, 0.8), cex.axis=1.5)
   mtext(ylab.text, side=2, line=2.85, cex=1.7)
   dev.off()
}

plotRD3vsRTALL <- function(cors, file.name, main.text, ymin, ymax, cols, legends, c=NA, isRT=F) {
   ylab.text <- "Spearman's rho"
   xlab.text <- "Chromosome"
 
   #png(paste0(file.name, ".png"), height=5, width=5, units="in", res=300)
   pdf(paste0(file.name, ".pdf"), height=5, width=5)
   plot(cors$cor1 ~ cors$chr, ylim=c(ymin, ymax), ylab=ylab.text, xlab=xlab.text, main=main.text, col=cols[1], xaxt="n", yaxt="n", pch=19, cex.lab=1.2, cex.main=1.25)
   lines(cors$cor1, y=NULL, type="l", lwd=3, col=cols[1])
   points(cors$chr, cors$cor2, col=cols[2], pch=19)
   lines(cors$cor2, y=NULL, type="l", lwd=3, col=cols[2])
   if (isRT) {
      points(cors$chr, cors$cor, col="black", pch=19)
      lines(cors$cor, y=NULL, type="l", lwd=3, col="black")
   }
   abline(h=0, lty=5)
 
   RT <- "S/G1"
   #legend("topright", c(paste0(legends[3], " vs. ", RT), paste0("      ", legends[1], " vs. ", RT)), text.col=c(cols[3], cols[1]), bty="n", cex=1.2) 
   #legend("bottomright", paste0(legends[2], " vs. ", RT), text.col=cols[2], bty="n", cex=1.2)
   legend("topright", paste0(legends[3], " vs. ", RT), text.col=cols[3], bty="n", cex=1.2) 
   legend("bottomright", c(paste0(legends[1], " vs. ", RT), paste0(legends[2], " vs. ", RT)), text.col=c(cols[1], cols[2]), bty="n", cex=1.2)        
   ##legend("bottomright", paste0(legends[2], " vs. ", RT), text.col=cols[2], bty="n", cex=1.1)
 
   if (!is.na(c)) {
      text(c, cors$cor1[c], round0(cors$cor1[c], digits=2), col=cols[1], pos=3, cex=1.1)   ##, offset=1.3)
      text(c, cors$cor2[c], round0(cors$cor2[c], digits=2), col=cols[2], pos=3, cex=1.1)
   }
   axis(side=1, at=seq(2, 22, by=4), cex.axis=1.1)
   axis(side=1, at=seq(4, 20, by=4), cex.axis=1.1)
   axis(side=2, at=seq(-0.8, 0.8, by=0.4), labels=c(-0.8, -0.4, 0, 0.4, 0.8), cex.axis=1.1)
   dev.off()
}

# -----------------------------------------------------------------------------
# Compare betweeen RT and LCL RT in *-wgs-rt.R
# Last Modified: 05/03/19
# -----------------------------------------------------------------------------
getRTvsRT <- function(nrds, nrds.lcl, bed.gc) {
   cors <- toTable(0, 5, 22, c("chr", "length", "cor", "cor1", "cor2"))
   cors$chr <- 1:22
   
   for (c in 1:22) {
      chr <- chrs[c]
      bed.gc.chr <- subset(bed.gc, CHR == chr)
      nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]
      nrds.chr.RT <- setSpline(nrds.chr, bed.gc.chr, "RT")
      nrds.chr.T  <- setSpline(nrds.chr, bed.gc.chr, "T")
      nrds.chr.N  <- setSpline(nrds.chr, bed.gc.chr, "N")
  
      nrds.lcl.chr <- nrds.lcl[intersect(nrds.lcl$BED, rownames(bed.gc.chr)),]  ## Reference LCL S/G1 ratio
      nrds.lcl.chr.RT <- setSpline(nrds.lcl.chr, bed.gc.chr, "RT")
  
      ## Keep only overlapping 1kb windows
      overlaps <- intersect(nrds.chr.RT$BED, nrds.lcl.chr.RT$BED)
      cors$length[c] <- length(overlaps)
      cors$cor[c]  <- getCor(nrds.chr.RT[overlaps,]$SPLINE, nrds.lcl.chr.RT[overlaps,]$SPLINE, method="spearman")
      cors$cor1[c] <- getCor(nrds.chr.T[overlaps,]$SPLINE,  nrds.lcl.chr.RT[overlaps,]$SPLINE, method="spearman")
      cors$cor2[c] <- getCor(nrds.chr.N[overlaps,]$SPLINE,  nrds.lcl.chr.RT[overlaps,]$SPLINE, method="spearman")
   }
   
   return(cors)
}

plotRTvsRTALL <- function(cors, file.name, main.text, ylab.text, xlab.text, ymin, ymax, col, c=NA, pos) {
   #png(paste0(file.name, ".png"), height=5, width=5, units="in", res=300)
   pdf(paste0(file.name, ".pdf"), height=5, width=5)
   plot(cors$cor ~ cors$chr, ylim=c(ymin, ymax), ylab=ylab.text, xlab=xlab.text, main=main.text, col=col, xaxt="n", yaxt="n", pch=19, cex.lab=1.2)
   lines(cors$cor, y=NULL, type="l", lwd=3, col=col)
   abline(h=0, lty=5)
 
   if (!is.na(c))
      text(cors[c,]$chr, cors[c,]$cor, round0(cors[c,]$cor, digits=2), col=col, pos=pos, cex=1.2)
   axis(side=1, at=seq(2, 22, by=2))
   axis(side=2, at=seq(-1, 1, by=0.2), labels=c(-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1))
   dev.off()
}

# -----------------------------------------------------------------------------
# Compare betweeen M2/M1, Q4/Q4, LCL and SCLC-NL RTs in *-wgs-rt.R
# Last Modified: 05/02/20
# -----------------------------------------------------------------------------
getRTvsRT3 <- function(nrds.m2, nrds.q4, nrds.sclc.nl.m2, nrds.lcl, bed.gc) {
   cors <- toTable(0, 5, 22, c("chr", "length", "cor", "cor1", "cor2"))
   cors$chr <- 1:22
 
   for (c in 1:22) {
      chr <- chrs[c]
      bed.gc.chr <- subset(bed.gc, CHR == chr)
      nrds.m2.chr <- nrds.m2[intersect(nrds.m2$BED, rownames(bed.gc.chr)),]
      nrds.m2.chr.RT <- setSpline(nrds.m2.chr, bed.gc.chr, "RT")
  
      nrds.q4.chr <- nrds.q4[intersect(nrds.q4$BED, rownames(bed.gc.chr)),]
      nrds.q4.chr.RT <- setSpline(nrds.q4.chr, bed.gc.chr, "RT")
      
      nrds.sclc.nl.m2.chr <- nrds.sclc.nl.m2[intersect(nrds.sclc.nl.m2$BED, rownames(bed.gc.chr)),]
      nrds.sclc.nl.m2.chr.RT <- setSpline(nrds.sclc.nl.m2.chr, bed.gc.chr, "RT")
      
      nrds.lcl.chr <- nrds.lcl[intersect(nrds.lcl$BED, rownames(bed.gc.chr)),]  ## Reference LCL S/G1 ratio
      nrds.lcl.chr.RT <- setSpline(nrds.lcl.chr, bed.gc.chr, "RT")
  
      ## Keep only overlapping 1kb windows
      overlaps <- intersect(intersect(intersect(nrds.m2.chr.RT$BED, nrds.q4.chr.RT$BED), nrds.sclc.nl.m2.chr.RT$BED), nrds.lcl.chr.RT$BED)
      cors$length[c] <- length(overlaps)
      cors$cor[c]  <- getCor(nrds.m2.chr.RT[overlaps,]$SPLINE, nrds.q4.chr.RT[overlaps,]$SPLINE, method="spearman")
      cors$cor1[c] <- getCor(nrds.m2.chr.RT[overlaps,]$SPLINE, nrds.sclc.nl.m2.chr.RT[overlaps,]$SPLINE, method="spearman")
      cors$cor2[c] <- getCor(nrds.m2.chr.RT[overlaps,]$SPLINE, nrds.lcl.chr.RT[overlaps,]$SPLINE, method="spearman")
   }
 
   return(cors)
}

plotRTvsRT3 <- function(cors, file.name, main.text, ymin, ymax, cols, legends) {
   ylab.text <- "Spearman's rho"
   xlab.text <- "Chromosome"
 
   #png(paste0(file.name, ".png"), height=5, width=5, units="in", res=300)
   pdf(paste0(file.name, ".pdf"), height=6, width=6)
   plot(cors$cor ~ cors$chr, ylim=c(ymin, ymax), ylab="", xlab=xlab.text, main=main.text, col=cols[1], xaxt="n", yaxt="n", pch=19, cex.axis=1.7, cex.lab=1.7, cex.main=1.8)
   lines(cors$cor, y=NULL, type="l", lwd=3, col=cols[1])
   points(cors$chr, cors$cor1, col=cols[2], pch=19)
   lines(cors$cor1, y=NULL, type="l", lwd=3, col=cols[2])
   points(cors$chr, cors$cor2, col=cols[3], pch=19)
   lines(cors$cor2, y=NULL, type="l", lwd=3, col=cols[3])

   legend("bottomright", c(legends[1], legends[2], legends[3]), text.col=c(cols[1], cols[2], cols[3]), bty="n", cex=1.75)        

   axis(side=1, at=seq(2, 22, by=4), cex.axis=1.6)
   axis(side=1, at=seq(4, 20, by=4), cex.axis=1.6)
   mtext(ylab.text, side=2, line=2.85, cex=1.7)
   axis(side=2, at=seq(0.5, 1, by=0.1), labels=c(0.5, 0.6, 0.7, 0.8, 0.9, 1), cex.axis=1.7)
   dev.off()
}

plotRTvsRT2 <- function(cors, file.name, main.text, ymin, ymax, cols, legends) {
   ylab.text <- "Spearman's rho"
   xlab.text <- "Chromosome"
 
   #png(paste0(file.name, ".png"), height=5, width=5, units="in", res=300)
   pdf(paste0(file.name, ".pdf"), height=6, width=6)
   plot(cors$cor ~ cors$chr, ylim=c(ymin, ymax), ylab="", xlab=xlab.text, main=main.text, col=cols[1], xaxt="n", yaxt="n", pch=19, cex.axis=1.7, cex.lab=1.7, cex.main=1.8)
   lines(cors$cor, y=NULL, type="l", lwd=3, col=cols[1])
   points(cors$chr, cors$cor1, col=cols[2], pch=19)
   lines(cors$cor1, y=NULL, type="l", lwd=3, col=cols[2])

   legend("bottomright", c(legends[1], legends[2]), text.col=c(cols[1], cols[2]), bty="n", cex=1.75)        
 
   axis(side=1, at=seq(2, 22, by=4), cex.axis=1.6)
   axis(side=1, at=seq(4, 20, by=4), cex.axis=1.6)
   mtext(ylab.text, side=2, line=2.85, cex=1.7)
   axis(side=2, at=seq(0.8, 1, by=0.05), labels=c(0.8, 0.85, 0.9, 0.95, 1), cex.axis=1.7)
   dev.off()
}

###
##
plotRD2 <- function(cors, file.name, main.text) {
   ylab.text <- "NBL Q4 vs. LCL S/G1 [rho]"
   xlab.text <- "NBL Q1 vs. LCL S/G1 [-rho]"
 
   #png(paste0(file.name, ".png"), height=5, width=5, units="in", res=300)
   pdf(paste0(file.name, ".pdf"), height=6, width=6)
   plot(cors$cor1 ~ abs(cors$cor2), ylab="", xlab=xlab.text, main=main.text, col="purple", pch=19, cex.axis=1.7, cex.lab=1.7, cex.main=1.8)
   lines(abs(cors$cor2), abs(cors$cor2), pch=NULL, col="purple", type="l", lty=2, lwd=2)
   #lm.fit <- lm(abs(cors$cor2) ~ abs(cors$cor2))
   #abline(lm.fit, col="purple", lty=2, lwd=2)
   
   for (c in 1:22)
      if (cors$cor1[c] > abs(cors$cor2[c]))
         text(abs(cors$cor2[c]), cors$cor1[c], paste0("chr", c), col="black", pos=3, cex=1.5)

   mtext(ylab.text, side=2, line=2.85, cex=1.7)
   dev.off()
}

# -----------------------------------------------------------------------------
# Compare betweeen RT and LCL RT in sclc-wgs-rt.R
# Last Modified: 03/06/19
# -----------------------------------------------------------------------------
plotSAMPLEvsRTALL <- function(cors.samples, samples, file.name, main.text=NA, ymin=NA, ymax=NA) {
   cors.samples.plot <- toTable(0, 2, 22*length(samples1), c("chr", "cor"))
   n <- length(samples)
   cnt <- 0
   for (c in 1:22) {
      start <- n * cnt + 1
      end   <- n * (cnt + 1)
      cors.samples.plot[start:end, 1] <- c
      cors.samples.plot[start:end, 2] <- as.numeric(cors.samples[c, samples])
  
      cnt <- cnt + 1
   }
 
   pdf(paste0(file.name, ".pdf"), height=6, width=6)
   boxplot(cor ~ chr, data=cors.samples.plot, ylim=c(ymin, ymax), ylab="", xlab="Chromosome", outline=T, xaxt="n", yaxt="n", main=main.text[1], cex.lab=1.7, cex.main=1.8)#, medcol="red")
   
   axis(side=1, at=seq(2, 22, by=4), cex.axis=1.5)
   axis(side=1, at=seq(4, 20, by=4), cex.axis=1.5)
   axis(side=2, at=seq(-0.8, 0.8, by=0.4), labels=c(-0.8, -0.4, 0, 0.4, 0.8), cex.axis=1.5)
   abline(h=0, lty=5)
   mtext("Spearman's rho", side=2, line=2.8, cex=1.7)
   dev.off()
}

# -----------------------------------------------------------------------------
# Mean replication timing ratio (MRTR)
# Last Modified: 06/10/19
# -----------------------------------------------------------------------------
getSPR <- function(nrds, bed.gc) {
   sprs <- toTable(0, 7, 22, c("chr", "cor", "cor1", "cor2", "e", "l", "spr"))
   sprs$chr <- 1:22
   
   for (c in 1:22) {
      chr <- chrs[c]
      bed.gc.chr <- subset(bed.gc, CHR == chr)
      nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]
      nrds.chr.T  <- setSpline(nrds.chr, bed.gc.chr, "T")
      nrds.chr.N  <- setSpline(nrds.chr, bed.gc.chr, "N")
      nrds.chr.RT <- setSpline(nrds.chr, bed.gc.chr, "RT")
  
      sprs$cor[c]  <- getCor(nrds.chr.T$SPLINE, nrds.chr.N$SPLINE,  method="spearman")
      sprs$cor1[c] <- getCor(nrds.chr.T$SPLINE, nrds.chr.RT$SPLINE, method="spearman")
      sprs$cor2[c] <- getCor(nrds.chr.N$SPLINE, nrds.chr.RT$SPLINE, method="spearman")
  
      e <- nrow(subset(nrds.chr.RT, SPLINE > 0))
      l <- nrow(subset(nrds.chr.RT, SPLINE < 0))
      sprs$e[c] <- e
      sprs$l[c] <- l
      sprs$spr[c] <- (e - l)/(e + l)
   }
   
   return(sprs)
}

plotSPR <- function(sprs, file.name, main.text, cs=NULL, digits, unit, ylab.text) {
   xlab.text <- "Chromosome"
   cols <- c("red", "blue", "black")
   #ylim <- getYlim(sprs$spr, unit)
   ylim <- c(-1, 1)
    
   pdf(paste0(file.name, ".pdf"), height=5, width=5)
   #plot(sprs$skew ~ sprs$chr, ylim=c(ymin, ymax), ylab=ylab.text, xlab=xlab.text, main=main.text, col=cols[3], xaxt="n", pch=19)   ## yaxt="n",
   plot(NULL, xlim=c(1, 22), ylim=ylim, xlab=xlab.text, ylab=ylab.text, main=main.text[1], col=cols[3], xaxt="n", yaxt="n", pch=19, cex.axis=1.1, cex.lab=1.2, cex.main=1.25)
   
   abline(h=sprs$spr[2], lty=5)
   lines(sprs$spr, y=NULL, type="l", lwd=3, col=cols[3])
   
   idx <- which(sprs$spr > sprs$spr[2])
   points(sprs$spr[idx] ~ sprs$chr[idx], col=cols[1], pch=19)
   idx <- which(sprs$spr < sprs$spr[2])
   points(sprs$spr[idx] ~ sprs$chr[idx], col=cols[2], pch=19)
   points(sprs$spr[2] ~ sprs$chr[2], col=cols[3], pch=19)
 
   text(sprs$chr[2], sprs$spr[2], paste0("Chr", 2), col=cols[3], pos=3, cex=1.2)
   #text(sprs$chr[2]+1.8, sprs$spr[2], paste0("Chr2 (", round0(sprs$spr[2], digits=digits), ")"), cex=1.1, col=cols[3], pos=3)
   if (!is.null(cs))
      for (c in 1:length(cs)) {
         c <- cs[c]
         if (sprs$spr[c] > sprs$spr[2])
            text(sprs$chr[c], sprs$spr[c], paste0("Chr", c), col=cols[1], pos=3, cex=1.2)
            #text(sprs$chr[c]+1.8, sprs$spr[c], paste0("Chr", c, " (", round0(sprs$spr[c], digits=digits), ")"), cex=1.1, col=cols[1], pos=3)
         else
            text(sprs$chr[c], sprs$spr[c], paste0("Chr", c), col=cols[2], pos=1, cex=1.2)
            #text(sprs$chr[c]+1.8, sprs$spr[c], paste0("Chr", c, " (", round0(sprs$spr[c], digits=digits), ")"), cex=1.1, col=cols[2], pos=1)
      }
   legend("topleft", "Earlier than chr2", text.col=cols[1], pch=16, col=cols[1], cex=1.05)   
   legend("bottomleft", "Later than chr2", text.col=cols[2], pch=16, col=cols[2], cex=1.05)

   axis(side=1, at=seq(2, 22, by=4), cex.axis=1.1)
   axis(side=1, at=seq(4, 20, by=4), cex.axis=1.1)
   axis(side=2, at=seq(-1, 1, by=0.5), labels=c(-1, -0.5, 0, 0.5, 1), cex.axis=1.1)
   mtext(main.text[2], line=0.3, cex=1.2)
   dev.off()
}

plotSPRRDC <- function(sprs, means, file.name, main.text, cs, xlab.text, unit, ylab.text) {
   cols <- c("red", "blue", "black", "purple")
   ylim <- c(-1, 1)
   
   unit <- (max(means) - min(means))/15
   xlim <- c(min(means) - unit, max(means) + unit)
   
   pdf(paste0(file.name, ".pdf"), height=5, width=5)
   if (is.null(xlim))
      plot(sprs ~ means, ylim=ylim, ylab=ylab.text, xlab=xlab.text, main=main.text[1], col="black", yaxt="n", pch=19, cex.axis=1.1, cex.lab=1.2, cex.main=1.25)
   else
      plot(sprs ~ means, ylim=ylim, xlim=xlim, ylab=ylab.text, xlab=xlab.text, main=main.text[1], col="black", yaxt="n", pch=19, cex.axis=1.1, cex.lab=1.2, cex.main=1.25)

   abline(h=sprs[2], lty=5)
   lm.fit <- lm(sprs ~ means)
   abline(lm.fit, col=cols[4], lwd=3)

   idx <- which(sprs > sprs[2])
   points(sprs[idx] ~ means[idx], col=cols[1], pch=19)
   idx <- which(sprs < sprs[2])
   points(sprs[idx] ~ means[idx], col=cols[2], pch=19)
   points(sprs[2] ~ means[2], col=cols[3], pch=19)
 
   text(means[2], sprs[2], paste0("Chr", 2), col=cols[3], pos=3, cex=1.2)
   if (!is.null(cs))
      for (c in 1:length(cs)) {
         c <- cs[c]
         if (sprs[c] > sprs[2])
            text(means[c], sprs[c], paste0("Chr", c), col=cols[1], pos=3, cex=1.2)
         else
            text(means[c], sprs[c], paste0("Chr", c), col=cols[2], pos=1, cex=1.2)
      }
   
   cor <- cor.test(sprs, means, method="spearman", exact=F)
   legends <- c("topright", "bottomleft")
   if (cor[[4]] > 0) legends[1] <- "topleft"
   legend(legends[1], "Earlier than chr2", text.col=cols[1], pch=16, col=cols[1], cex=1.05)   ## bty="n"
   legend(legends[2], "Later than chr2", text.col=cols[2], pch=16, col=cols[2], cex=1.05)
   
   legend("bottomright", c(paste0("rho = ", round0(cor[[4]], digits=2)), paste0("p-value = ", scientific(cor[[3]]))), text.col=cols[4], bty="n", cex=1.2)
   #legend("bottomright", c(paste0("R^2 = ", round0(summary(lm.fit)$r.squared, digits=2)), paste0("p-value = ", scientific(summary(lm.fit)$coefficients[2, 4]))), text.col=cols[4], bty="n", cex=1.1)
   #axis(side=1, at=seq(2, 22, by=2))
   axis(side=2, at=seq(-1, 1, by=0.5), labels=c(-1, -0.5, 0, 0.5, 1), cex.axis=1.1)
   mtext(main.text[2], cex=1.2, line=0.3)
   dev.off()
}

# -----------------------------------------------------------------------------
# Find S-like and G1-like tumour samples
# Last Modified: 04/06/19
# -----------------------------------------------------------------------------
setSamplesQ4 <- function(wd.rt.data, samples1) {
   cors.all <- toTable(0, 4, length(samples1), c("SAMPLE_ID", "COR", "Q4", "M2"))
   rownames(cors.all) <- samples1
   cors.all$SAMPLE_ID <- samples1
   for (s in 1:length(samples1)) {
      sample <- samples1[s]
      load(file.path(wd.rt.data, "samples", paste0("rd-vs-rt_", sample, "-vs-lcl_spline_spearman.RData")))
  
      cors.all$COR[s] <- cor
   }
 
   q <- quantile(as.numeric(cors.all$COR))
   print(q)

   samples.q4 <- list()
   samples.q4[[4]] <- rownames(subset(cors.all, COR > as.numeric(q[4])))
   samples.q4[[3]] <- rownames(subset(subset(cors.all, COR > as.numeric(q[3])), COR <= as.numeric(q[4])))
   samples.q4[[2]] <- rownames(subset(subset(cors.all, COR > as.numeric(q[2])), COR <= as.numeric(q[3])))
   samples.q4[[1]] <- rownames(subset(cors.all, COR <= as.numeric(q[2])))
 
   cors.all$Q4[which(cors.all$SAMPLE_ID %in% samples.q4[[4]])] <- 4
   cors.all$Q4[which(cors.all$SAMPLE_ID %in% samples.q4[[3]])] <- 3
   cors.all$Q4[which(cors.all$SAMPLE_ID %in% samples.q4[[2]])] <- 2
   cors.all$Q4[which(cors.all$SAMPLE_ID %in% samples.q4[[1]])] <- 1
 
   cors.all$M2[which(cors.all$Q4 %in% c(3, 4))] <- 1
   cors.all$M2[which(cors.all$Q4 %in% c(1, 2))] <- 0

   return(cors.all)
}

# =============================================================================
# Inner Class  : PeifLyne File Reader
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 01/05/17
# =============================================================================
read.peiflyne.cn.txt <- function(cn.file) {
   cn <- readTable(cn.file, header=F, rownames=T, sep="")
   colnames(cn) <- c("BED", "CHR", "START", "END", "V5", "V6", "V7", "RD_T", "RD_N", "V10", "INSERT_SIZE_T", "INSERT_SIZE_N")
 
   return(cn)
}

read.peiflyne.cn.seg <- function(cn.file) {
   cn <- readTable(cn.file, header=F, rownames=F, sep="")
   colnames(cn) <- c("SAMPLE", "CHR", "START", "END", "V5", "CN")
 
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

# =============================================================================
# Inner Class: Collections of test/obsolete/deprecated methods
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 01/02/18
# =============================================================================
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

setSamplesSG1 <- function(wd.rt.data, samples1, cors.samples) {
   cors.all <- toTable(0, 4, length(samples1), c("SAMPLE_ID", "COR", "SG1", "M2"))
   rownames(cors.all) <- samples1
   cors.all$SAMPLE_ID <- samples1
   for (s in 1:length(samples1)) {
      sample <- samples1[s]
      load(file.path(wd.rt.data, "samples", paste0("rd-vs-rt_", sample, "-vs-lcl_spline_spearman.RData")))
  
      cors.all$COR[s] <- cor
   }
   #print(quantile(as.numeric(cors.all$COR)))
 
   ##
   s_likes <- c()
   g1_likes <- c()
   for (c in 1:22) {
      cors <- cors.samples[c, -c(1:4)]
      median <- median(as.numeric(cors))
      #median <- 0
      
      if (c == 1) {
         s_likes  <- colnames(cors)[which(as.numeric(cors) > median)]
         g1_likes <- colnames(cors)[which(as.numeric(cors) <= median)]
      } else {
         s_likes  <- intersect(s_likes, colnames(cors)[which(as.numeric(cors) > median)])
         g1_likes <- intersect(g1_likes, colnames(cors)[which(as.numeric(cors) <= median)])
      }
   }
   print(length(s_likes))
   print(length(g1_likes))
   print(s_likes)
   print(g1_likes)
   
   cors.all$SG1 <- NA
   cors.all[s_likes, ]$SG1 <- "SL"
   cors.all[g1_likes,]$SG1 <- "G1L"
   
   cors.all$M2 <- NA
   cors.all$M2[which(cors.all$SG1 == "SL") ] <- 1
   cors.all$M2[which(cors.all$SG1 == "G1L")] <- 0
   return(cors.all)
}
