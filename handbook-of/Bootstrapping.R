# =============================================================================
# Library      : Bootstrapping
# Name         : handbook-of/Bootstrapping.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 13/09/19; 05/04/19; 21/02/19
# =============================================================================

# -----------------------------------------------------------------------------
# Read in bootstrapping data (in 3b)
# Last Modified: 20/09/19; 31/10/18
# -----------------------------------------------------------------------------
getBSTRPS <- function(nrds.RT, colname, b) {
   nrds.RT <- nrds.RT[, c("BED", colname)]
   colnames(nrds.RT)[2] <- paste0(colname, "_", b)
 
   return(nrds.RT)
}

# -----------------------------------------------------------------------------
# Determine RFD from bootstrapping data (in 3b)
# Last Modified: 20/09/19; 31/10/18
# -----------------------------------------------------------------------------
getPOS <- function(nrds.RT) {
   return(as.numeric(length(which(nrds.RT > 0))))
}

getNEG <- function(nrds.RT) {
 return(as.numeric(length(which(nrds.RT < 0))))
}

getRFD <- function(nrds.RT) {
   return((nrds.RT$POS - nrds.RT$NEG)/(nrds.RT$POS + nrds.RT$NEG))
}

pipeBootstrapping <- function(nrds.RT, bstrps) {
   nrds.RT$POS <- mapply(x = 1:nrow(nrds.RT), function(x) as.numeric(getPOS(nrds.RT[x, 1:bstrps])))   ## BUG FIX: 01/11/18
   nrds.RT$NEG <- mapply(x = 1:nrow(nrds.RT), function(x) as.numeric(getNEG(nrds.RT[x, 1:bstrps])))   ## BUG FIX: 01/11/18
   nrds.RT$RFD <- mapply(x = 1:nrow(nrds.RT), function(x) as.numeric(getRFD(nrds.RT[x,])))
 
   return(nrds.RT)
}

getBootstrapping <- function(nrds.RT, origin.lower, origin.upper) {
   nrds.RT.l   <- subset(nrds.RT,   POS >= origin.lower)
   nrds.RT.l.u <- subset(nrds.RT.l, POS <= origin.upper)
   
   diff <- setdiff(rownames(nrds.RT), rownames(nrds.RT.l.u))
   return(nrds.RT[diff,])
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
plotBootstrapsHist <- function(bed.gc.rt.chr, file.name, main.text, xlab.text, breaks, origin.break) {
   cols <- rep("steelblue1", breaks)
   cols[(breaks/2 + origin.break):breaks] <- "sandybrown"
   cols[(breaks/2 - origin.break + 1):(breaks/2 + origin.break)] <- "red"
   h <- hist(bed.gc.rt.chr$NEG, breaks=breaks) 
   ymax <- max(c(h$counts[2:4], h$counts[(breaks-3):(breaks-1)]))   ## Calculatte max frequency in row 2 before next line
   h$counts <- h$counts/1000                                        ## Change frequency scale to x1000 in row 1
   
   pdf(file.name, height=6, width=6)
   #par(mfrow=c(2,1))
   layout(matrix(c(1,2), 2, 1), widths=1, heights=c(1,2))           ## One figure each in row 1 and row 2; row 1 is 1/3 the height of row 2
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
   hist(bed.gc.rt.chr$NEG, main="" , ylab="Frequency", xlab=xlab.text, ylim=c(0, ymax), breaks=breaks, col=cols, las=1, axes=F)
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
