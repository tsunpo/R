# =============================================================================
# Library      : Bootstrap
# Name         : handbook-of/Bootstrap.R
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
getPOS <- function(nrds.RT) {   ## Left leadings are with a positive slope
   return(as.numeric(length(which(nrds.RT > 0))))
}

getNEG <- function(nrds.RT) {   ## Right leadings are with a negative slope
   return(as.numeric(length(which(nrds.RT < 0))))
}

getRFD <- function(nrds.RT) {   ## (R-L)/(R+L)
   return((nrds.RT$NEG - nrds.RT$POS)/(nrds.RT$NEG + nrds.RT$POS))
}

pipeBootstrap <- function(nrds.RT, bstrps) {
   nrds.RT$POS <- mapply(x = 1:nrow(nrds.RT), function(x) as.numeric(getPOS(nrds.RT[x, 1:bstrps])))   ## BUG FIX: 01/11/18
   nrds.RT$NEG <- mapply(x = 1:nrow(nrds.RT), function(x) as.numeric(getNEG(nrds.RT[x, 1:bstrps])))   ## BUG FIX: 01/11/18
   nrds.RT$RFD <- mapply(x = 1:nrow(nrds.RT), function(x) as.numeric(getRFD(nrds.RT[x,])))
 
   return(nrds.RT)
}

getBootstrapRFD <- function(nrds.RT, boundary.lower, boundary.upper) {
   nrds.RT.l   <- subset(nrds.RT,   NEG > boundary.lower)
   nrds.RT.l.u <- subset(nrds.RT.l, NEG < boundary.upper)
   
   diff <- setdiff(rownames(nrds.RT), rownames(nrds.RT.l.u))
   return(nrds.RT[diff,])
}

# -----------------------------------------------------------------------------
# Visualisation of bootstrap re-sampling data (Histogram, RFD, and RT)
# Last Modified: 13/11/18
# -----------------------------------------------------------------------------
## https://www.r-graph-gallery.com/190-mirrored-histogram/
## http://www.r-graph-gallery.com/72-set-margin-size-with-par-mar-function/
## https://www.statmethods.net/advgraphs/layout.html
## https://stackoverflow.com/questions/6461209/how-to-round-up-to-the-nearest-10-or-100-or-x
## https://stackoverflow.com/questions/15717545/set-the-intervals-of-x-axis-using-r
plotBootstrapHist <- function(nrds.RT.BSTRPS, file.name, main.text, xlab.text, breaks, boundary.break) {
   cols <- rep(lightblue, breaks)
   cols[(breaks/2 + boundary.break + 1):breaks] <- orange
   #cols[(breaks/2 - boundary.break + 1):(breaks/2 + boundary.break)] <- "gray"
   #cols[51:breaks] <- "sandybrown"
   
   h <- hist(nrds.RT.BSTRPS$NEG, breaks=breaks)
   ymax <- max(c(h$counts[2:4], h$counts[(breaks-3):(breaks-1)]))   ## Calculatte max frequency in row 2 before next line
   h$counts <- h$counts/1000                                        ## Change frequency scale to x1000 in row 1
   
   pdf(file.name, height=5, width=5)
   #par(mfrow=c(2,1))
   layout(matrix(c(1,2), 2, 1), widths=1, heights=c(1,2))           ## One figure each in row 1 and row 2; row 1 is 1/3 the height of row 2
   ylim <- sort(c(h$counts[1], h$counts[breaks]), decreasing=F)     ## Min and max frequencies in row 1 (in x1000 scale)
   if (ylim[1] > 500) {                                             ## Round up to the nearest power of 10 (in x1000 scale)
      ylim[1] <- floor(ylim[1]/100)*100
   } else if (ylim[1] > 10) {
      ylim[1] <- floor(ylim[1]/10)*10
   } else
      ylim[1] <- floor(ylim[1])
   par(mar=c(1,4,3.6,1))
   plot(h, main=main.text[1], ylab="Fr. [x1000]", xlab="", ylim=ylim, border="darkgray", col=cols, xaxt="n", cex.axis=1.2, cex.lab=1.3, cex.main=1.35)
   abline(v=500, lty=5, lwd=1.5, col="black")
   #abline(v=950, lty=5, lwd=1, col="red")
   #abline(v=50,  lty=5, lwd=1, col="red")
   #text(950, max(ylim)-10, "950", cex=1.2, col="red") 
   #text(50, max(ylim)-10, "50", cex=1.2, col="blue") 
   #mtext("RFD = (R\u2212L)/(R+L)", cex=1.3, line=4.7)   ## separator(nrow(nrds.RT.BSTRPS))

   ##
   par(mar=c(5,4,0,1))
   hist(nrds.RT.BSTRPS$NEG, main="", ylab="Frequency", xlab=xlab.text, ylim=c(0, ymax), breaks=breaks, border="darkgray", col=cols, las=1, axes=F, cex.axis=1.2, cex.lab=1.3, cex.main=1.4)
   if (ymax < 1000) {
      axis(side=2, at=seq(0, ymax, by=250), cex.axis=1.2)
   } else if (ymax < 3000) {
      axis(side=2, at=seq(0, ymax, by=500), cex.axis=1.2)
   } else if (ymax > 20000) {
      axis(side=2, at=seq(0, ymax, by=10000), cex.axis=1.2)
   } else
      axis(side=2, at=seq(0, ymax, by=1000), cex.axis=1.2)
   axis(side=1, at=seq(0, 1000, by=250), cex.axis=1.2)
   abline(v=500, lty=5, lwd=1.5, col="black")
   #abline(v=950, lty=5, lwd=1, col="red")
   #abline(v=50,  lty=5, lwd=1, col="red")
   text(500, ymax*4/5, "RFD = 0", cex=1.3, col="black", font=2) 
   
   dev.off()
}

boundary.upper <- 500   ## 500-520 breaks
boundary.lower <- 500   ## 480-500 breaks
boundary.break <- 0   ## 1 breaks each centering 500
file.name <- file.path(wd.rt.plots, paste0("hist_", base, "_rpkm_SLOPE_RFD=0.pdf"))
main.text <- c(paste0(BASE, " bootstrap distribution"), paste0(""))
xlab.text <- "Number of rightward forks per kb window"
plotBootstrapHist(nrds.RT.BSTRPS, file.name, main.text, xlab.text, 100, boundary.break)

plotBootstrapRFD0 <- function(file.name, BASE, chr, xmin, xmax, nrds.RT.NRFD, bed.gc.chr, boundary.upper, boundary.lower, ext, width, kb, withUnclassified=F) {
   overlaps <- intersect(rownames(bed.gc.chr), nrds.RT.NRFD$BED)
   nrds.RT.NRFD.chr <- nrds.RT.NRFD[overlaps,]
   bed.gc.chr <- bed.gc.chr[overlaps,]

   adjustcolor.gray <- adjustcolor("darkgray", alpha.f=0.08)
 
   rights <- rownames(subset(nrds.RT.NRFD.chr, RFD >= 0))
   lefts  <- rownames(subset(nrds.RT.NRFD.chr, RFD <= 0))
   
   if (width == 10) main.text <- paste0(BASE, " bootstrap replication fork directionality (RFD)")
   else main.text <- paste0(BASE, " bootstrap-based RFD")
   if (withUnclassified)
      main.text <- paste0(main.text, " (", kb, " kb)")
   
   if (!is.na(xmin) && !is.na(xmax)) file.name <- paste0(file.name, "_", xmin/1E6, "-", xmax/1E6, "Mb_", kb, "kb")
   if (withUnclassified) file.name <- paste0(file.name, "_with-un")
   if (is.na(xmin)) {
      start <- bed.gc.chr[rownames(nrds.RT.NRFD.chr)[1],]$START
      if (start < 5000000) xmin <- 0
      else xmin <- start
   }
   if (is.na(xmax)) xmax <- subset(chromInfo, chrom == chr)$size
   
   ## Initiation plot
   if (ext == "pdf") {
      pdf(paste0(file.name, ".pdf"), height=5, width=width)
   } else if (ext == "png")
      png(paste0(file.name, ".png"), height=5, width=width, units="in", res=300)   ## ADD 16/05/17: res=300
   
   ###
   ## Initiate RT plot
   layout(matrix(c(1,2), 2, 1), widths=1, heights=c(1,1))   ## One figure each in row 1 and row 2   ## See plotBootstrapsHist()
   par(mar=c(1,4,4,1))
   ylab.text <- "RT [log2]"
   
   plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=c(-2, 2), xlab="", ylab=ylab.text, main=main.text, xaxt="n", yaxt="n", cex.axis=1.2, cex.lab=1.3, cex.main=1.35)
   points(bed.gc.chr$START/1E6, nrds.RT.NRFD.chr$RT, col=adjustcolor.gray, pch=16, cex=0.3)
   
   axis(side=2, at=seq(-2, 2, by=4), labels=c("\u22122", 2), cex.axis=1.15)
   axis(side=2, at=seq(-1, 1, by=1), labels=c("\u22121", 0, 1), cex.axis=1.15)
   abline(h=0, lty=5, lwd=1.5, col="black")
   
   ## Plot cytobands (before smoothing spline)
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.4, col="lightgrey") 
   
   ## Plot smoothing spline with right-/left-leading positions
   points(bed.gc.chr[lefts,]$START/1E6,  nrds.RT.NRFD.chr[lefts,]$SPLINE, col=lightblue, pch=19, cex=0.5)
   points(bed.gc.chr[rights,]$START/1E6, nrds.RT.NRFD.chr[rights,]$SPLINE, col=orange, pch=19, cex=0.5)

   ## Plot legend
   #legend("topright", c("Leftward", "Rightward   "), col=c(blue, orange), bty="n", pt.cex=1, lty=c(1, 1), lwd=c(3, 3), pch=c(NA, NA), horiz=T, cex=1.2)
   legend("topright", "Rightward ", col=orange, bty="n", pt.cex=1, lty=1, lwd=3, pch=NA, horiz=T, cex=1.3)
   legend("topleft",  "Leftward",  col=lightblue, bty="n", pt.cex=1, lty=1, lwd=3, pch=NA, horiz=T, cex=1.3, inset=c(0.03, 0))
   mtext("RFD = (R\u2212L)/(R+L)", line=0.25, cex=1.3)   ## separator(nrow(nrds.RT.BSTRPS)),
   
   ###
   ## Initiate RFD plot
   par(mar=c(5.5,4,0,1))
   xlab.text <- paste0("Chromosome ", gsub("chr", "", chr), " position [Mb]")
   ylab.text <- "RFD"
   
   plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=c(-1.8, 1.8), xlab=xlab.text, ylab=ylab.text, main="", yaxt="n", cex.axis=1.2, cex.lab=1.3)
   axis(side=2, at=seq(-1, 1, by=1), labels=c("\u22121", 0, 1), cex.axis=1.15)
   abline(h=0, lty=5, lwd=1.5, col="black")
   #abline(h=0.9, lty=5, lwd=1, col="black")
   #abline(h=-0.9, lty=5, lwd=1, col="black")
   
   ## Plot cytobands (before points)
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.4, col="lightgrey")
   
   points(bed.gc.chr[lefts,]$START/1E6,  nrds.RT.NRFD.chr[lefts,]$RFD,  col=lightblue, pch=19, cex=0.5)
   points(bed.gc.chr[rights,]$START/1E6, nrds.RT.NRFD.chr[rights,]$RFD, col=orange, pch=19, cex=0.5)

   ## Plot legend
   legend("topright", "TTR (R)", col=orange, bty="n", pt.cex=1, lty=1, lwd=3, pch=NA, horiz=T, cex=1.3)
   legend("bottomright", "TTR (L)", col=lightblue, bty="n", pt.cex=1, lty=1, lwd=3, pch=NA, horiz=T, cex=1.3)
   dev.off()
}






boundary.upper <- 500   ## 500-520 breaks
boundary.lower <- 500   ## 480-500 breaks
boundary.break <- 0   ## 1 breaks each centering 500

## Chr2
c <- 2
chr <- chrs[c]
bed.gc.chr <- subset(bed.gc, CHR == chr)
nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]

## RFD
load(file=file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.d.rt.RT.SLOPE_", chr, ".RData")))
nrds.RT.BSTRPS.chr$RFD <- getRFD(nrds.RT.BSTRPS.chr)
file.name <- file.path(wd.rt.plots, paste0("RFD_", base, "_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, "_RFD=0"))
plotBootstrapRFD(file.name, BASE, chr,  70000000,  80000000, nrds.chr, bed.gc.chr, nrds.RT.BSTRPS.chr, boundary.upper, boundary.lower, "png", width=5)
plotBootstrapRFD(file.name, BASE, chr, 110000000, 130000000, nrds.chr, bed.gc.chr, nrds.RT.BSTRPS.chr, boundary.upper, boundary.lower, "png", width=5)

## Chr12
c <- 12
chr <- chrs[c]
bed.gc.chr <- subset(bed.gc, CHR == chr)
nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]

## RFD   
load(file=file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.d.rt.RT.SLOPE_", chr, ".RData")))
nrds.RT.BSTRPS.chr$RFD <- getRFD(nrds.RT.BSTRPS.chr)
file.name <- file.path(wd.rt.plots, paste0("RFD_", base, "_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, "_RFD=0"))
plotBootstrapRFD(file.name, BASE, chr,  97500000, 105000000, nrds.chr, bed.gc.chr, nrds.RT.BSTRPS.chr, boundary.upper, boundary.lower, "png", width=5)





plotBootstrapsRT <- function(file.name, BASE, chr, xmin, xmax, rt.chr, bed.gc.chr, right.idx, left.idx, boundary.idx, ymax, ext, gene="") {
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
   points(spline$x[boundary.idx]/1E6, spline$y[boundary.idx], col="red", pch=16, cex=0.4)
 
   ## Plot cytobands and legend
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.4, col="lightgrey")
 
   legend("bottomright", c("Left-leading", "Left / Right", "Right-leading"), col=c("steelblue1", "red", "sandybrown"), pch=16, cex=0.8, horiz=T)
   dev.off()
}
