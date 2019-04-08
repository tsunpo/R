# =============================================================================
# Library      : Replication Fork Directionality
# Name         : handbook-of/ReplicationForkDirectionality.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 05/04/19; 21/02/19
# =============================================================================

# -----------------------------------------------------------------------------
# OK-seq (in lcl-ok-rfd.R)
# -----------------------------------------------------------------------------
setScaledOK <- function(rpkms.chr.rt, pseudocount) {
   rpkms.chr.rt$C  <- log2(rpkms.chr.rt$T + pseudocount)
   rpkms.chr.rt$W  <- log2(rpkms.chr.rt$N + pseudocount)
   rpkms.chr.rt$RFD <- (rpkms.chr.rt$C - rpkms.chr.rt$W) / (rpkms.chr.rt$C + rpkms.chr.rt$W)
   #rpkms.chr.rt$RFD <- scale(rpkms.chr.rt$RFD)   ## REMOVE 05/04/19; RFD is between -1 to 1
 
   rpkms.chr.rt$CW <- rpkms.chr.rt$C + rpkms.chr.rt$W
   rpkms.chr.rt$CW <- scale(rpkms.chr.rt$CW)

   return(rpkms.chr.rt[,c("BED", "C", "W", "RFD", "CW")])
}

plotCW <- function(file.name, main.text, chr, xmin, xmax, rpkms.chr.rt, bed.gc.chr, colours, legends, colours2, legends2, ext, width, peaks, ymin=NA, ymax=NA, cutoff, scale) {
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
   #par(mfrow=c(3,1)) 
   par(mar=c(1,4,4,1))
   ylab.text <- "Read depth"
   if (is.na(ymin) || is.na(ymax)) {
      plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=c(-max(rpkms.chr.rt$W), max(rpkms.chr.rt$C)), xlab=xlab.text, ylab=ylab.text, main=main.text, xaxt="n")
   } else
      plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=c(ymin, ymax), xlab=xlab.text, ylab=ylab.text, main=main.text, xaxt="n")
   #points(bed.gc.chr$START/1E6, rpkms.chr.rt, col=colours[1], cex=0.3)
   abline(h=0, lwd=0.5, col="lightgrey")
 
   ## Plot cytobands (before smoothing spline)
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.4, col="lightgrey") 
 
   ## Plot smoothing splines
   lines(bed.gc.chr$START/1E6, rpkms.chr.rt$C, col=colours[1], pch=16, cex=0.2, type="S")
   lines(bed.gc.chr$START/1E6, -rpkms.chr.rt$W, col=colours[2], pch=16, cex=0.2, type="S")
 
   ## Plot legend and peaks
   if (xmin == 0) {
      legend("bottomright", legends, col=colours, lty=1, lwd=2, bty="n", horiz=T)
   } else {
      legend("bottomright", legends2, col=colours, lty=1, lwd=2, bty="n", horiz=T)
   }
   if (length(peaks) != 0)
      for (p in 1:length(peaks))
         abline(v=peaks[p]/1E6, lty=5, lwd=1, col="black")

   ### 
   ## Initiate C+W plot
   par(mar=c(5.5,4,0,1))
   ylab.text <- "C+W"
   
   plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=c(-cutoff, cutoff), xlab=xlab.text, ylab=ylab.text, main="", yaxt="n")
   points(bed.gc.chr$START/1E6, rpkms.chr.rt$CW, col="lightgrey", cex=0.3)
   #idx <- which(rpkms.chr.rt$RT == 0)
   #points(bed.gc.chr[idx,]$START/1E6, rpkms.chr.rt[idx,]$CW, col="lightgrey", cex=0.3)
   #idx <- which(rpkms.chr.rt$RT > 0)
   #points(bed.gc.chr[idx,]$START/1E6, rpkms.chr.rt[idx,]$CW, col=colours2[1], cex=0.3)
   #idx <- which(rpkms.chr.rt$RT < 0)
   #points(bed.gc.chr[idx,]$START/1E6, rpkms.chr.rt[idx,]$CW, col=colours2[2], cex=0.3)
   
   abline(h=0, lwd=0.5, col="lightgrey")
   axis(side=2, at=seq(-scale, scale, by=scale), labels=c(-scale, 0, scale))
   
   ## Plot cytobands (before smoothing spline)
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.4, col="lightgrey") 
   
   ## Plot smoothing spline
   spline <- smooth.spline(x=bed.gc.chr$START, y=rpkms.chr.rt$CW)
   points(bed.gc.chr$START/1E6, spline$y, col="black", pch=16, cex=0.2)
   
   ## Plot legend and peaks
   if (xmin == 0) {
      legend("topright", paste0(legends2[1], "+", legends2[2], " read depth"), col="black", lty=1, lwd=2, bty="n", horiz=T)
      #legend("topleft", "Early", bty="n", text.col="black")   
      #legend("bottomleft", "Late", bty="n", text.col="black")
   } else {
      legend("topright", paste0(legends2[1], "+", legends2[2]), col="black", lty=1, lwd=2, bty="n", horiz=T)
      #legend("topleft",    paste0("Early: ", legends2[1], " > ", legends2[2]), bty="n", text.col=colours[1])   
      #legend("bottomleft", paste0("Late:  ", legends2[1], " < ", legends2[2]), bty="n", text.col=colours[2])
   }
   if (length(peaks) != 0)
      for (p in 1:length(peaks))
         abline(v=peaks[p]/1E6, lty=5, lwd=1, col="black")
   
   dev.off()
}


## Change to plotRFD
plotCW <- function(file.name, main.text, chr, xmin, xmax, rpkms.chr.rt, bed.gc.chr, colours, legends, colours2, legends2, ext, width, peaks, ymin=NA, ymax=NA, cutoff, scale) {
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
 layout(atrix(c(1,2), 2, 1), widths=1, heights=c(1,1))   ## One figure each in row 1 and row 2   ## See plotBootstrapsHist()
 #par(mfrow=c(3,1)) 
 par(mar=c(1,4,4,1))
 ylab.text <- "Read depth"
 if (is.na(ymin) || is.na(ymax)) {
  plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=c(-max(rpkms.chr.rt$W), max(rpkms.chr.rt$C)), xlab=xlab.text, ylab=ylab.text, main=main.text, xaxt="n")
 } else
  plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=c(ymin, ymax), xlab=xlab.text, ylab=ylab.text, main=main.text, xaxt="n")
 #points(bed.gc.chr$START/1E6, rpkms.chr.rt, col=colours[1], cex=0.3)
 abline(h=0, lwd=0.5, col="lightgrey")
 
 ## Plot cytobands (before smoothing spline)
 cytoBand.chr <- subset(cytoBand, chrom == chr)
 for (c in 1:nrow(cytoBand.chr))
  abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.4, col="lightgrey") 
 
 ## Plot smoothing splines
 lines(bed.gc.chr$START/1E6, rpkms.chr.rt$C, col=colours[1], pch=16, cex=0.2, type="S")
 lines(bed.gc.chr$START/1E6, -rpkms.chr.rt$W, col=colours[2], pch=16, cex=0.2, type="S")
 
 ## Plot legend and peaks
 if (xmin == 0) {
  legend("bottomright", legends, col=colours, lty=1, lwd=2, bty="n", horiz=T)
 } else {
  legend("bottomright", legends2, col=colours, lty=1, lwd=2, bty="n", horiz=T)
 }
 if (length(peaks) != 0)
  for (p in 1:length(peaks))
   abline(v=peaks[p]/1E6, lty=5, lwd=1, col="black")
 
 ### 
 ## Initiate RFD plot
 par(mar=c(5.5,4,0,1))
 ylab.text <- "RFD"
 
 plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=c(-1, 1), xlab=xlab.text, ylab=ylab.text, main="", yaxt="n")
 idx <- which(rpkms.chr.rt$RT == 0)
 points(bed.gc.chr[idx,]$START/1E6, rpkms.chr.rt[idx,]$RFD, col="lightgrey", cex=0.3)
 idx <- which(rpkms.chr.rt$RT > 0)
 points(bed.gc.chr[idx,]$START/1E6, rpkms.chr.rt[idx,]$RFD, col=colours2[1], cex=0.3)
 idx <- which(rpkms.chr.rt$RT < 0)
 points(bed.gc.chr[idx,]$START/1E6, rpkms.chr.rt[idx,]$RFD, col=colours2[2], cex=0.3)
 
 abline(h=0, lwd=0.5, col="lightgrey")
 axis(side=2, at=seq(-scale, scale, by=scale), labels=c(-scale, 0, scale))
 
 ## Plot cytobands (before smoothing spline)
 cytoBand.chr <- subset(cytoBand, chrom == chr)
 for (c in 1:nrow(cytoBand.chr))
  abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.4, col="lightgrey") 
 
 ## Plot legend and peaks
 if (xmin == 0) {
  legend("topright", "(C-W)/(C+W)", col="black", lty=1, lwd=2, bty="n", horiz=T)
  #legend("topleft", "Early", bty="n", text.col="black")   
  #legend("bottomleft", "Late", bty="n", text.col="black")
 } else {
  legend("topright", "(C-W)/(C+W)", col="black", lty=1, lwd=2, bty="n", horiz=T)
  #legend("topleft",    paste0("Early: ", legends2[1], " > ", legends2[2]), bty="n", text.col=colours[1])   
  #legend("bottomleft", paste0("Late:  ", legends2[1], " < ", legends2[2]), bty="n", text.col=colours[2])
 }
 if (length(peaks) != 0)
  for (p in 1:length(peaks))
   abline(v=peaks[p]/1E6, lty=5, lwd=1, col="black")
 
 ### 
 ## Initiate C+W plot
 par(mar=c(5.5,4,0,1))
 ylab.text <- "C+W"
 
 plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=c(-cutoff, cutoff), xlab=xlab.text, ylab=ylab.text, main="", yaxt="n")
 idx <- which(rpkms.chr.rt$RT == 0)
 points(bed.gc.chr[idx,]$START/1E6, rpkms.chr.rt[idx,]$CW, col="lightgrey", cex=0.3)
 idx <- which(rpkms.chr.rt$RT > 0)
 points(bed.gc.chr[idx,]$START/1E6, rpkms.chr.rt[idx,]$CW, col=colours2[1], cex=0.3)
 idx <- which(rpkms.chr.rt$RT < 0)
 points(bed.gc.chr[idx,]$START/1E6, rpkms.chr.rt[idx,]$CW, col=colours2[2], cex=0.3)
 
 abline(h=0, lwd=0.5, col="lightgrey")
 axis(side=2, at=seq(-scale, scale, by=scale), labels=c(-scale, 0, scale))
 
 ## Plot cytobands (before smoothing spline)
 cytoBand.chr <- subset(cytoBand, chrom == chr)
 for (c in 1:nrow(cytoBand.chr))
  abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.4, col="lightgrey") 
 
 ## Plot smoothing spline
 spline <- smooth.spline(x=bed.gc.chr$START, y=rpkms.chr.rt$CW)
 points(bed.gc.chr$START/1E6, spline$y, col="black", pch=16, cex=0.2)
 
 ## Plot legend and peaks
 if (xmin == 0) {
  legend("topright", paste0(legends2[1], "+", legends2[2], " read depth ratio"), col="black", lty=1, lwd=2, bty="n", horiz=T)
  #legend("topleft", "Early", bty="n", text.col="black")   
  #legend("bottomleft", "Late", bty="n", text.col="black")
 } else {
  legend("topright", paste0(legends2[1], "+", legends2[2], " ratio"), col="black", lty=1, lwd=2, bty="n", horiz=T)
  #legend("topleft",    paste0("Early: ", legends2[1], " > ", legends2[2]), bty="n", text.col=colours[1])   
  #legend("bottomleft", paste0("Late:  ", legends2[1], " < ", legends2[2]), bty="n", text.col=colours[2])
 }
 if (length(peaks) != 0)
  for (p in 1:length(peaks))
   abline(v=peaks[p]/1E6, lty=5, lwd=1, col="black")
 
 dev.off()
}




plotOK <- function(file.name, main.text, ylab.text, chr, xmin, xmax, rpkms.chr.rt, bed.gc.chr, colors, ext, ymin=NA, ymax=NA) {
   xlab.text <- paste0("Chromosome ", gsub("chr", "", chr), " coordinate (Mb)")
   if (!is.na(xmin) && !is.na(xmax)) file.name <- paste0(file.name, "_", xmin/1E6, "-", xmax/1E6, "Mb")
   if (is.na(xmin)) xmin <- 0
   if (is.na(xmax)) xmax <- subset(chromInfo, chrom == chr)$size
 
   ## Initiate plot
   if (ext == "pdf") {
      pdf(paste0(file.name, ".pdf"), height=4, width=10)
   } else if (ext == "png")
      png(paste0(file.name, ".png"), height=4, width=10, units="in", res=300)   ## ADD 16/05/17: res=300
 
   plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=c(ymin, ymax), xlab=xlab.text, ylab=ylab.text, main=main.text)   #, yaxt="n")
   points(bed.gc.chr$START/1E6, rpkms.chr.rt$RT, col=colors[1], cex=0.3)
   idx <- which(rpkms.chr.rt$RT < 0)
   points(bed.gc.chr[idx,]$START/1E6, rpkms.chr.rt[idx,]$RT, col=colors[2], cex=0.3)
 
   #abline(h=0, lwd=0.5, col="lightgrey")
   #axis(side=2, at=seq(-scale, scale, by=scale), labels=c(-scale, 0, scale))
 
   ## Plot cytobands (before smoothing spline)
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.4, col="lightgrey") 
 
   ## Plot smoothing spline
   spline <- smooth.spline(x=bed.gc.chr$START, y=rpkms.chr.rt$RT)
   points(bed.gc.chr$START/1E6, spline$y, col="black", pch=16, cex=0.2)
 
   dev.off()
}

# -----------------------------------------------------------------------------
# Read in bootstrapped data (in 2b)
# -----------------------------------------------------------------------------
getEnsGeneBED <- function(pos, bed.gc.chr) {
   bed.gc.chr.start <- subset(bed.gc.chr, pos >= START)
   bed.gc.chr.start.end <- subset(bed.gc.chr.start, pos <= END)    ## ADD "=" 19/11/18
 
   return(rownames(bed.gc.chr.start.end))
}

# -----------------------------------------------------------------------------
# Read in bootstrapped data (in 2c and 2d)
# Last Modified: 31/10/18
# -----------------------------------------------------------------------------
loadEnsGeneRT <- function(wd1.rt.data, b, type) {
   load(file.path(wd1.rt.data, b, paste0(type, "_ensGene.rt.RData")))
   ensGene.rt <- ensGene.rt[!is.na(ensGene.rt$SLOPE_START),]
   ensGene.rt <- ensGene.rt[!is.na(ensGene.rt$SLOPE_END),]
 
   return(ensGene.rt)
}

getEnsGeneRT <- function(ensGene.rt, colname, b) {
   ensGene.rt$ensembl_gene_id <- rownames(ensGene.rt)
   ensGene.rt.position <- ensGene.rt[, c("ensembl_gene_id", colname)]
   colnames(ensGene.rt.position)[1+b] <- paste0(colname, "_", b)
 
   return(ensGene.rt.position)
}

loadBEDGCRT <- function(wd1.rt.data, b, type) {
   load(file.path(wd1.rt.data, b, paste0(type, "_ensGene.rt.RData")))
   bed.gc.rt <- bed.gc.rt[!is.na(bed.gc.rt$SLOPE),]
 
   return(bed.gc.rt)
}

getBEDGCRT <- function(bed.gc.rt, colname, b) {
   bed.gc.rt$BED <- rownames(bed.gc.rt)
   bed.gc.rt <- bed.gc.rt[, c("BED", colname)]
   colnames(bed.gc.rt)[1+b] <- paste0(colname, "_", b)
 
   return(bed.gc.rt)
}

# -----------------------------------------------------------------------------
# Determine RFD from bootstrapped data (in 2c and 2d)
# Last Modified: 31/10/18
# -----------------------------------------------------------------------------
getRightLeading <- function(ensGene.rt.bstrp) {
   return(as.numeric(length(which(ensGene.rt.bstrp < 0))))
}

getLeftLeading <- function(ensGene.rt.bstrp) {
   return(as.numeric(length(which(ensGene.rt.bstrp > 0))))
}

#getLeading <- function(ensGene.rt.bstrp) {
#   right <- ensGene.rt.bstrp$RIGHT_LEADING
#   left  <- ensGene.rt.bstrp$LEFT_LEADING
# 
#   if (right > left) return(1)
#   else if (right < left) return(-1)
#   else return(0)
#}

getRFD <- function(bed.gc.rt) {
   return((bed.gc.rt$RIGHT_LEADING - bed.gc.rt$LEFT_LEADING)/(bed.gc.rt$RIGHT_LEADING + bed.gc.rt$LEFT_LEADING))
}

pipeLeading <- function(ensGene.rt.bstrp, bstrps) {
   ensGene.rt.bstrp$RIGHT_LEADING <- mapply(x = 1:nrow(ensGene.rt.bstrp), function(x) as.numeric(getRightLeading(ensGene.rt.bstrp[x, 1:bstrps])))   ## BUG FIX: 01/11/18
   ensGene.rt.bstrp$LEFT_LEADING  <- mapply(x = 1:nrow(ensGene.rt.bstrp), function(x) as.numeric(getLeftLeading(ensGene.rt.bstrp[x, 1:bstrps])))    ## BUG FIX: 01/11/18
   ensGene.rt.bstrp$RFD           <- mapply(x = 1:nrow(ensGene.rt.bstrp), function(x) as.numeric(getRFD(ensGene.rt.bstrp[x, ])))
 
   return(ensGene.rt.bstrp)
}

# -----------------------------------------------------------------------------
# Add BED information for the bootstrapped data (in 2c)
# Last Modified: 09/11/18
# -----------------------------------------------------------------------------
setEnsGeneBED <- function(ensGene.rt, bed.gc, chrs, isStartPosition) {
   ensGene.rt <- cbind(ensGene[rownames(ensGene.rt), ], ensGene.rt[, c("RIGHT_LEADING", "LEFT_LEADING", "RFD")])
   ensGene.rt.chrs <- ensGene.rt[1,]
   ensGene.rt.chrs$BED <- 0
   ensGene.rt.chrs <- ensGene.rt.chrs[-1,]
   for (c in 1:length(chrs)) {
      chr <- chrs[c]
      bed.gc.chr <- subset(bed.gc, CHR == chr)
      ensGene.rt.chr <- subset(ensGene.rt, chromosome_name == chr)
  
      if (isStartPosition)
         ensGene.rt.chr$BED <- mapply(x = 1:nrow(ensGene.rt.chr), function(x) getEnsGeneBED(ensGene.rt.chr$start_position[x], bed.gc.chr))
      else
         ensGene.rt.chr$BED <- mapply(x = 1:nrow(ensGene.rt.chr), function(x) getEnsGeneBED(ensGene.rt.chr$end_position[x],   bed.gc.chr))   ## BUG FIXED: 19/11/18
      ensGene.rt.chrs <- rbind(ensGene.rt.chrs, ensGene.rt.chr)
   }
 
   return(ensGene.rt.chrs[, c("BED", "RIGHT_LEADING", "LEFT_LEADING", "RFD")])
}

# -----------------------------------------------------------------------------
# Visualisation of bootstrapping data (Histogram, RFD, and RT)
# Last Modified: 13/11/18
# -----------------------------------------------------------------------------
## https://www.r-graph-gallery.com/190-mirrored-histogram/
## http://www.r-graph-gallery.com/72-set-margin-size-with-par-mar-function/
## https://www.statmethods.net/advgraphs/layout.html
## https://stackoverflow.com/questions/6461209/how-to-round-up-to-the-nearest-10-or-100-or-x
## https://stackoverflow.com/questions/15717545/set-the-intervals-of-x-axis-using-r
plotBootstrapsHist <- function(bed.gc.rt.chr, file.name, main.text, BASE, breaks, origin.break) {
   main.text <- paste0(main.text, " (", BASE, ")")
   xlab.text <- c("Number of right-leading counts", "(out of 1,000 bootstraps)")
   cols <- rep("steelblue1", breaks)
   cols[(breaks/2 + origin.break):breaks] <- "sandybrown"
   cols[(breaks/2 - origin.break + 1):(breaks/2 + origin.break)] <- "red"
   h <- hist(bed.gc.rt.chr$RIGHT_LEADING, breaks=breaks) 
   ymax <- max(c(h$counts[2:4], h$counts[(breaks-3):(breaks-1)]))   ## Calculatte max frequency in row 2 before next line
   h$counts <- h$counts/1000                                        ## Change frequency scale to x1000 in row 1
   
   pdf(file.name, height=6, width=6)
   #par(mfrow=c(2,1))
   layout(matrix(c(1,2), 2, 1), widths=1, heights=c(1,2))   ## One figure each in row 1 and row 2; row 1 is 1/3 the height of row 2
   ylim <- sort(c(h$counts[1], h$counts[breaks]), decreasing=F)     ## Min and max frequencies in row 1 (in x1000 scale)
   if (ylim[1] > 500)                                               ## Round up to the nearest power of 10 (in x1000 scale)
      ylim[1] <- floor(ylim[1]/100)*100
   else if (ylim[1] > 10)
      ylim[1] <- floor(ylim[1]/10)*10
   else
      ylim[1] <- floor(ylim[1])
   par(mar=c(1,4,4,1))
   plot(h, main=main.text, ylab="Frequency (x1000)", xlab="", ylim=ylim, col=cols, xaxt="n")
   
   par(mar=c(5.5,4,0,1))
   hist(bed.gc.rt.chr$RIGHT_LEADING, main="" , ylab="Frequency", xlab=xlab.text, ylim=c(0, ymax), breaks=breaks, col=cols, las=1, axes=F)
   if (ymax < 1000)
      axis(side=2, at=seq(0, ymax, by=250))
   else if (ymax < 3000)
      axis(side=2, at=seq(0, ymax, by=500))
   else if (ymax > 50000)
      axis(side=2, at=seq(0, ymax, by=10000))
   else
      axis(side=2, at=seq(0, ymax, by=1000))
   axis(side=1, at=seq(0, 1000, by=250))
   mtext(paste0("Distribution of 1kb windows (n=", separator(nrow(bed.gc.rt.chr)), ")"), cex=1.2, line=6.3)
   dev.off()
}

plotBootstrapsRFD <- function(file.name, BASE, chr, xmin, xmax, bed.gc.chr, bed.gc.rt.chr, right.idx, left.idx, origin.idx, ext) {
   main.text <- paste0("Bootstrapped replication fork directionality (RFD) in ", BASE)
   xlab.text <- paste0("Chromosome ", gsub("chr", "", chr), " coordinate (Mb)")
   ylab.text <- "RFD = (R-L)/(R+L)"
   if (!is.na(xmin) && !is.na(xmax)) file.name <- paste0(file.name, "_", xmin/1E6, "-", xmax/1E6, "Mb")
   if (is.na(xmin)) xmin <- 0
   if (is.na(xmax)) xmax <- subset(chromInfo, chrom == chr)$size
   ymin <- -2
   ymax <- 1.5
 
   ## Initiation plot
   if (ext == "pdf") {
      pdf(paste0(file.name, ".pdf"), height=3.5, width=10)
   } else if (ext == "png")
      png(paste0(file.name, ".png"), height=3.5, width=10, units="in", res=300)   ## ADD 16/05/17: res=300
   plot(NULL, ylim=c(ymin, ymax), xlim=c(xmin/1E6, xmax/1E6), xlab=xlab.text, ylab=ylab.text, main=main.text, yaxt="n")
   axis(side=2, at=seq(-1, 1, by=1))
   
   ## Plot right-/left-leading positions
   points(bed.gc.chr$START[left.idx]/1E6,   bed.gc.rt.chr$RFD[left.idx],   col="steelblue1", cex=0.2)
   points(bed.gc.chr$START[right.idx]/1E6,  bed.gc.rt.chr$RFD[right.idx],  col="sandybrown", cex=0.2)
   points(bed.gc.chr$START[origin.idx]/1E6, bed.gc.rt.chr$RFD[origin.idx], col="red", cex=0.5)
 
   #spline <- smooth.spline(x=bed.gc.chr$START, y=bed.gc.rt.chr$RATIO)   ## The smooth line is a bit confusing
   #lines(spline$x/1E6, spline$y)
 
   ## Plot cytobands and legend
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.4, col="lightgrey")
 
   legend("bottomright", c("Left-leading", "Left / Right", "Right-leading"), col=c("steelblue1", "red", "sandybrown"), pch=16, cex=0.8, horiz=T)
   dev.off()
}

plotBootstrapsRT <- function(file.name, BASE, chr, xmin, xmax, rt.chr, bed.gc.chr, right.idx, left.idx, origin.idx, ymax, ext, gene="") {
   main.text <- paste0("Read depth (CN-, GC-corrected) ratio (T/N) in ", BASE)
   xlab.text <- paste0("Chromosome ", gsub("chr", "", chr), " coordinate (Mb)")
   if (gene != "")
      xlab.text <- paste0(xlab.text, " near ", gene)
   ylab.text <- "Replication timing (log2)"
   if (!is.na(xmin) && !is.na(xmax)) file.name <- paste0(file.name, "_", xmin/1E6, "-", xmax/1E6, "Mb")
   if (is.na(xmin)) xmin <- 0
   if (is.na(xmax)) xmax <- subset(chromInfo, chrom == chr)$size

   ## Initiate plot
   if (ext == "pdf") {
      pdf(paste0(file.name, ".pdf"), height=4, width=10)
   } else if (ext == "png")
      png(paste0(file.name, ".png"), height=4, width=10, units="in", res=300)   ## ADD 16/05/17: res=300
   plot(NULL, ylim=c(-ymax, ymax), xlim=c(xmin/1E6, xmax/1E6), xlab=xlab.text, ylab=ylab.text, main=main.text)
   points(bed.gc.chr$START/1E6, rt.chr$RT, col="gray", cex=0.3)
   abline(h=0, lwd=0.5, col="lightgrey")
 
   ## Plot right-/left-leading positions
   spline <- smooth.spline(x=bed.gc.chr$START, y=rt.chr$RT)
   points(spline$x[left.idx]/1E6,   spline$y[left.idx], col="steelblue1", pch=16, cex=0.4)
   points(spline$x[right.idx]/1E6,  spline$y[right.idx], col="sandybrown", pch=16, cex=0.4)
   points(spline$x[origin.idx]/1E6, spline$y[origin.idx], col="red", pch=16, cex=0.4)
 
   ## Plot cytobands and legend
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.4, col="lightgrey")
 
   legend("bottomright", c("Left-leading", "Left / Right", "Right-leading"), col=c("steelblue1", "red", "sandybrown"), pch=16, cex=0.8, horiz=T)
   dev.off()
}

# -----------------------------------------------------------------------------
# Visualisation of bootstrapping data (Leading-count ratio for expressed genes)
# Last Modified: 13/11/18
# -----------------------------------------------------------------------------
#getLeadingRatio <- function(bed.gc.rt) {
#   if (bed.gc.rt$RT == 1) return(log2(bed.gc.rt$RIGHT_LEADING/500))
#   else if (bed.gc.rt$RT == -1) return(log2(bed.gc.rt$LEFT_LEADING/500) * -1)
#   else return(0)
#}

getEnsGeneTxRFD <- function(ensGene, ensGene.rt.start, ensGene.rt.end, tpm.gene.log2) {
   txs <- intersect(rownames(ensGene.rt.start), rownames(tpm.gene.log2))
   ensGene.rt.tx <- cbind(ensGene[txs, -1], ensGene.rt.start[txs,], ensGene.rt.end[txs,])
   ensGene.rt.tx$LENGTH <- abs(ensGene.rt.tx$start_position - ensGene.rt.tx$end_position)
   ensGene.rt.tx$MEDIAN <- tpm.gene.log2[txs,]$MEDIAN
   
   txs.fw.idx <- which(ensGene.rt.tx$strand == 1)
   txs.re.idx <- which(ensGene.rt.tx$strand == -1)   
   ensGene.rt.tx$RFD_TSS <- NA
   ensGene.rt.tx$RFD_TTS <- NA
   ensGene.rt.tx$RFD_TSS[txs.fw.idx] <- mapply(x = 1:length(txs.fw.idx), function(x) as.numeric(getRFD(ensGene.rt.tx[txs.fw.idx[x],  7:10])))   ## FW: TSS = start_position
   ensGene.rt.tx$RFD_TTS[txs.fw.idx] <- mapply(x = 1:length(txs.fw.idx), function(x) as.numeric(getRFD(ensGene.rt.tx[txs.fw.idx[x], 11:14])))   ##     TES = end_position
   ensGene.rt.tx$RFD_TSS[txs.re.idx] <- mapply(x = 1:length(txs.re.idx), function(x) as.numeric(getRFD(ensGene.rt.tx[txs.re.idx[x], 11:14])))   ## RE: TSS = end_position
   ensGene.rt.tx$RFD_TTS[txs.re.idx] <- mapply(x = 1:length(txs.re.idx), function(x) as.numeric(getRFD(ensGene.rt.tx[txs.re.idx[x],  7:10])))   ##     TES = start_position
   
   colnames(ensGene.rt.tx)[ 7:10] <- paste0(colnames(ensGene.rt.tx)[ 7:10], "_1")
   colnames(ensGene.rt.tx)[11:14] <- paste0(colnames(ensGene.rt.tx)[11:14], "_2")
   ensGene.rt.tx$CONSIST <- ensGene.rt.tx$RT_1 * ensGene.rt.tx$RT_2   ## This may include "50/50" if cutoff for leading ratio change
   
   ensGene.rt.tx$RT <- NA
   ensGene.rt.tx$RT[txs.fw.idx] <- ensGene.rt.tx[txs.fw.idx,]$RT_1   ## FW: TSS = start_position
   ensGene.rt.tx$RT[txs.re.idx] <- ensGene.rt.tx[txs.re.idx,]$RT_2   ## RE: TSS = end_position
   
   ensGene.rt.tx$CD <- ensGene.rt.tx$strand * ensGene.rt.tx$RT
   return(ensGene.rt.tx)
}

plotEnsGeneTxRFD <- function(file.name, BASE, ensGene.rt.tx, origin.upper, ext) {
   main.text  <- paste0("Transcription level and RFD in ", BASE)
   mtext.text <- paste0("Expressed genes (n=", separator(nrow(ensGene.rt.tx)), ")")
   xlab.text <- c("RFD = (R-L)/(R+L)", "TSS")
   ylab.text <- "log2(TPM+0.01)"
   rfd <- (origin.upper - (1000 - origin.upper))/1000
   ensGene.rt.tx.50.50 <- subset(subset(ensGene.rt.tx, RFD_TSS <= rfd), RFD_TSS >= -rfd)
   
   if (ext == "pdf") {
      pdf(paste0(file.name, ".pdf"), height=6, width=6)
   } else if (ext == "png")
      png(paste0(file.name, ".png"), height=6, width=6, units="in", res=300)   ## ADD 16/05/17: res=300
   cols <- rep("steelblue1", nrow(ensGene.rt.tx))
   plot(  MEDIAN ~ RFD_TSS, data=ensGene.rt.tx, xlab=xlab.text, ylab=ylab.text, main=main.text, col=cols)
   points(MEDIAN ~ RFD_TSS, data=subset(ensGene.rt.tx, RFD_TSS > leading.ratio), col="sandybrown")
   points(MEDIAN ~ RFD_TSS, data=subset(ensGene.rt.tx, CONSIST < 0), col="purple1")
   points(MEDIAN ~ RFD_TSS, data=ensGene.rt.tx.50.50, col="red")
   
   legend("top", c("Left-leading", paste0("L/R (n=", separator(nrow(ensGene.rt.tx.50.50)), ")"), "Right-leading", "Inconsistent"), col=c("steelblue1", "red", "sandybrown", "purple1"), pch=1, cex=0.75, horiz=T)
   mtext(mtext.text, cex=1, line=0.4)
   dev.off()
}

# -----------------------------------------------------------------------------
# Visualisation of bootstrapping data (Expression level and gene length)
# Last Modified: 13/11/18
# -----------------------------------------------------------------------------
testWbyEnsGeneRTTx <- function(ensGene.rt.tx, sclc.tx.right, sclc.tx.left) {
   tx.right <- ensGene.rt.tx[sclc.tx.right,]$MEDIAN
   tx.left  <- ensGene.rt.tx[sclc.tx.left, ]$MEDIAN
   
   return(testW(tx.right, tx.left))
}

plotQ4SS <- function(q4ss, file.name, file.main, mtext, ylim, ylab, isLog10) {
   counts <- c()
   quantiles <- c()
   names1 <- c()
   names2 <- c()
   for (q in 1:4) {
      counts    <- c(counts, q4ss[[q]])
      quantiles <- c(quantiles, mapply(x = 1:length(q4ss[[q]]), function(x) paste0("q", q)))
      names1 <- c(names1, paste0("Q", q, "\n"))
      names2 <- c(names2, paste0("", "\nn=", length(q4ss[[q]])))
   }
 
   pdf(file.name, height=6, width=4)
   #par(mgp=c(3,2,0))
   if (isLog10)
      boxplot(log10(as.numeric(counts))~quantiles, ylab=ylab, xlab=c("", "Expression"), names=c("", "", "", ""), main=file.main, ylim=log10(ylim))
   else
      boxplot(as.numeric(counts)~quantiles, ylab=ylab, xlab=c("", "Expression"), names=c("", "", "", ""), main=file.main, ylim=ylim)
   mtext(mtext, cex=0.8, line=0.5)
   #mtext("2                                   ", cex=0.6, line=0.25)
   axis(side=1, at=1:4, names1, cex.axis=1, line=1, lwd=0)
   axis(side=1, at=1:4, names2, cex.axis=0.8, line=1, lwd=0)
   dev.off()
}
