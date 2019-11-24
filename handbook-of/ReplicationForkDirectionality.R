# =============================================================================
# Library      : Replication Fork Directionality
# Name         : handbook-of/ReplicationForkDirectionality.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 21/10/19; 13/09/19; 05/04/19; 21/02/19
# =============================================================================

# -----------------------------------------------------------------------------
# Read in bootstrap data (in 3b)
# Last Modified: 20/09/19; 31/10/18
# -----------------------------------------------------------------------------
getBSTRPS <- function(nrds.RT, colname, b) {
   nrds.RT <- nrds.RT[, c("BED", colname)]
   colnames(nrds.RT)[2] <- paste0(colname, "_", b)
 
   return(nrds.RT)
}

# -----------------------------------------------------------------------------
# Determine RFD from bootstrap data (in 3b)
# Last Modified: 20/09/19; 31/10/18
# -----------------------------------------------------------------------------
getPOS <- function(nrds.RT.BSTRPS) {   ## Left leadings are with a positive slope
   return(as.numeric(length(which(nrds.RT.BSTRPS > 0))))
}

getNEG <- function(nrds.RT.BSTRPS) {   ## Right leadings are with a negative slope
   return(as.numeric(length(which(nrds.RT.BSTRPS < 0))))
}

getRFD <- function(nrds.RT.BSTRPS) {   ## (R-L)/(R+L)
   return((nrds.RT.BSTRPS$NEG - nrds.RT.BSTRPS$POS)/(nrds.RT.BSTRPS$NEG + nrds.RT.BSTRPS$POS))
}

pipeBootstrap <- function(nrds.RT.BSTRPS, bstrps) {
   nrds.RT.BSTRPS$POS <- mapply(x = 1:nrow(nrds.RT.BSTRPS), function(x) as.numeric(getPOS(nrds.RT.BSTRPS[x, 1:bstrps])))   ## BUG FIX: 01/11/18
   nrds.RT.BSTRPS$NEG <- mapply(x = 1:nrow(nrds.RT.BSTRPS), function(x) as.numeric(getNEG(nrds.RT.BSTRPS[x, 1:bstrps])))   ## BUG FIX: 01/11/18

   return(nrds.RT.BSTRPS)
}

getBootstrap <- function(base, column) {
   load(file=file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.d.rt.RT.", column, "_", "chr1", ".RData")))
   nrds.RT.BSTRPS <- nrds.RT.BSTRPS.chr
   for (c in 2:22) {
      chr <- chrs[c]

      load(file=file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.d.rt.RT.", column, "_", chr, ".RData")))
      nrds.RT.BSTRPS <- rbind(nrds.RT.BSTRPS, nrds.RT.BSTRPS.chr)
   }
   
   return(nrds.RT.BSTRPS)
}

getRTRFD <- function(nrds.RT, nrds.RT.BSTRPS) {
   nrds.RT.BSTRPS$RFD <- getRFD(nrds.RT.BSTRPS)
   colnames(nrds.RT.BSTRPS) <- c("L", "R", "RFD")
   overlaps <- intersect(rownames(nrds.RT), rownames(nrds.RT.BSTRPS))
   
   return(cbind(nrds.RT[overlaps,], nrds.RT.BSTRPS[overlaps,]))   ## CHANGE FROM nrds.RT.BSTRPS[rownames(nrds.RT),]
}

# -----------------------------------------------------------------------------
# Determine TTR and CTR
# Last Modified: 20/09/19; 31/10/18
# -----------------------------------------------------------------------------
getBootstrapCTR <- function(nrds.RFD, boundary.lower, boundary.upper) {
   nrds.RFD.l   <- subset(nrds.RFD,   R > boundary.lower)
   nrds.RFD.l.u <- subset(nrds.RFD.l, R < boundary.upper)
 
   return(nrds.RFD.l.u)
}

getBootstrapTTR <- function(nrds.RFD, boundary.lower, boundary.upper) {
   nrds.RFD.l.u <- getBootstrapCTR(nrds.RFD, boundary.lower, boundary.upper) 
   
   diff <- setdiff(rownames(nrds.RFD), rownames(nrds.RFD.l.u))
   return(nrds.RFD[diff,])
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
   adjustcolor.blue <- adjustcolor("steelblue1", alpha.f=0.15)
   adjustcolor.orange <- adjustcolor("sandybrown", alpha.f=0.15)
 
   cols <- rep("steelblue1", breaks)
   cols[(breaks/2 + boundary.break):breaks] <- "sandybrown"
   cols[(breaks/2 - boundary.break + 1):50] <- "white"
   cols[51:(breaks/2 + boundary.break)] <- "white"
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
   plot(h, main=main.text[1], ylab="Freq. (x1000)", xlab="", ylim=ylim, col=cols, xaxt="n", cex.axis=1.1, cex.lab=1.2, cex.main=1.25)
   #abline(v=500, lty=5, lwd=1, col="black")
   abline(v=950, lty=5, lwd=1, col="black")
   abline(v=50,  lty=5, lwd=1, col="black")

   ##
   par(mar=c(5,4,0,1))
   hist(nrds.RT.BSTRPS$NEG, main="", ylab="Frequency", xlab=xlab.text, ylim=c(0, ymax), breaks=breaks, col=cols, las=1, axes=F, cex.axis=1.1, cex.lab=1.2, cex.main=1.25)
   if (ymax < 1000) {
      axis(side=2, at=seq(0, ymax, by=250))
   } else if (ymax < 3000) {
      axis(side=2, at=seq(0, ymax, by=500))
   } else if (ymax > 20000) {
      axis(side=2, at=seq(0, ymax, by=10000))
   } else
      axis(side=2, at=seq(0, ymax, by=1000))
   axis(side=1, at=seq(0, 1000, by=250))
   #abline(v=500, lty=5, lwd=1, col="black")
   abline(v=950, lty=5, lwd=1, col="black")
   abline(v=50,  lty=5, lwd=1, col="black")
   text(500, ymax*4/5, "|RFD| > 0.9 (80.7%)", cex=1.2, col="black") 
   
   mtext(main.text[2], line=4.8, cex=1.2)   ## separator(nrow(nrds.RT.BSTRPS)),
   dev.off()
}

boundary.upper <- 950   ## 500-520 breaks
boundary.lower <-  50   ## 480-500 breaks
boundary.break <-  45   ## 1 breaks each centering 500
file.name <- file.path(wd.rt.plots, paste0("hist_", base, "_rpkm_SLOPE_RFD>0.9_white.pdf"))
main.text <- c(paste0(BASE, " bootstrap distribution"), "Chr1-22 (1-kbs)")   #paste0("Chr1-22 (1-kbs)"))
xlab.text <- "Number of right-leading resamplings"
plotBootstrapHist(nrds.RT.BSTRPS, file.name, main.text, xlab.text, 100, boundary.break)



plotBootstrapRFD <- function(file.name, BASE, chr, xmin, xmax, nrds.chr, bed.gc.chr, nrds.RT.BSTRPS.chr, boundary.upper, boundary.lower, ext, width) {
   adjustcolor.gray <- adjustcolor("darkgray", alpha.f=0.08)

   nrds.chr.RT <- setSpline(nrds.chr, bed.gc.chr, "RT")
   overlaps <- intersect(rownames(nrds.RT.BSTRPS.chr), rownames(nrds.chr.RT))
   nrds.RT.BSTRPS.chr <- nrds.RT.BSTRPS.chr[overlaps,]
   nrds.chr.RT <- nrds.chr.RT[overlaps,]
   bed.gc.chr  <- bed.gc.chr[overlaps,]
   
   rights <- rownames(subset(nrds.RT.BSTRPS.chr, RFD >= 0.9))
   lefts  <- rownames(subset(nrds.RT.BSTRPS.chr, RFD <= -0.9))
   boundary.rights <- rownames(subset(subset(nrds.RT.BSTRPS.chr, RFD < 0.9), RFD >= 0))
   boundary.lefts <- rownames(subset(subset(nrds.RT.BSTRPS.chr, RFD < 0), RFD > -0.9))
   
   boundaries <- c(boundary.rights, boundary.lefts)
   earlies <- rownames(subset(nrds.chr.RT[boundaries,], SPLINE >= 0))
   lates <- rownames(subset(nrds.chr.RT[boundaries,], SPLINE < 0))
   
   if (width == 10) main.text <- paste0(BASE, " bootstrap replication fork directionality (RFD)")
   else main.text <- paste0(BASE, " (WB) bootstrap RFD")

   if (!is.na(xmin) && !is.na(xmax)) file.name <- paste0(file.name, "_", xmin/1E6, "-", xmax/1E6, "Mb")
   if (is.na(xmin)) {
      start <- bed.gc.chr[rownames(nrds.chr.RT)[1],]$START
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

   plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=c(-2, 2), xlab="", ylab=ylab.text, main=main.text, xaxt="n", yaxt="n", cex.axis=1.1, cex.lab=1.2, cex.main=1.25)
   points(bed.gc.chr$START/1E6, nrds.chr.RT$RT, col=adjustcolor.gray, pch=16, cex=0.35)
   
   axis(side=2, at=seq(-2, 2, by=4), labels=c("\u22122", 2), cex.axis=1.1)
   axis(side=2, at=seq(-1, 1, by=1), labels=c("\u22121", 0, 1), cex.axis=1.1)
   abline(h=0, lty=5, lwd=1, col="black")
   
   ## Plot cytobands (before smoothing spline)
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.4, col="lightgrey") 
   
   ## Plot smoothing spline with right-/left-leading positions
   points(bed.gc.chr[lefts,]$START/1E6,  nrds.chr.RT[lefts,]$SPLINE, col="steelblue1", cex=0.4)
   points(bed.gc.chr[rights,]$START/1E6, nrds.chr.RT[rights,]$SPLINE, col="sandybrown", cex=0.4)
   points(bed.gc.chr[lates,]$START/1E6,   nrds.chr.RT[lates,]$SPLINE,   col="blue", pch=19, cex=0.4)
   points(bed.gc.chr[earlies,]$START/1E6, nrds.chr.RT[earlies,]$SPLINE, col="red", pch=19, cex=0.4)
   
   ## Plot legend
   legend("topright", "CTR (E)", col="red", bty="n", pt.cex=1, lty=1, lwd=3, pch=NA, horiz=T, cex=1.2)
   legend("bottomright", "CTR (L)", col="blue", bty="n", pt.cex=1, lty=1, lwd=3, pch=NA, horiz=T, cex=1.2)
   mtext("", line=0.25, cex=1.2)   ## separator(nrow(nrds.RT.BSTRPS)),
      
   ###
   ## Initiate RFD plot
   par(mar=c(5.5,4,0,1))
   xlab.text <- paste0("Chromosome ", gsub("chr", "", chr), " coordinate [Mb]")
   ylab.text <- "RFD"
   
   plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=c(-1.8, 1.8), xlab=xlab.text, ylab=ylab.text, main="", yaxt="n", cex.axis=1.1, cex.lab=1.2)
   axis(side=2, at=seq(-1, 1, by=1), labels=c("\u22121", 0, 1), cex.axis=1.1)
   #abline(h=0, lty=5, lwd=1, col="black")
   abline(h=0.9, lty=5, lwd=1, col="black")
   abline(h=-0.9, lty=5, lwd=1, col="black")
   
   ## Plot cytobands (before points)
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.4, col="lightgrey") 

   points(bed.gc.chr[lates,]$START/1E6,   nrds.RT.BSTRPS.chr[lates,]$RFD,   col="blue", cex=0.4)
   points(bed.gc.chr[earlies,]$START/1E6, nrds.RT.BSTRPS.chr[earlies,]$RFD, col="red", cex=0.4)
   points(bed.gc.chr[lefts,]$START/1E6,  nrds.RT.BSTRPS.chr[lefts,]$RFD,  col="steelblue1", cex=0.4)
   points(bed.gc.chr[rights,]$START/1E6, nrds.RT.BSTRPS.chr[rights,]$RFD, col="sandybrown", cex=0.4)
   
   ## Plot legend
   legend("topright", "TTR (R)", col="sandybrown", bty="n", pt.cex=1, lty=1, lwd=3, pch=NA, horiz=T, cex=1.2)
   legend("bottomright", "TTR (L)", col="steelblue1", bty="n", pt.cex=1, lty=1, lwd=3, pch=NA, horiz=T, cex=1.2)
   dev.off()
}

boundary.upper <- 950   ## 500-520 breaks
boundary.lower <-  50   ## 480-500 breaks
boundary.break <-  45   ## 45 breaks each centering 500

## Chr2
c <- 2
chr <- chrs[c]
bed.gc.chr <- subset(bed.gc, CHR == chr)
nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]

## RFD
load(file=file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.d.rt.RT.SLOPE_", chr, ".RData")))
nrds.RT.BSTRPS.chr$RFD <- getRFD(nrds.RT.BSTRPS.chr)
file.name <- file.path(wd.rt.plots, paste0("RFD_", base, "_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, "_TTR"))
plotBootstrapRFD(file.name, BASE, chr, 110000000, 130000000, nrds.chr, bed.gc.chr, nrds.RT.BSTRPS.chr, boundary.upper, boundary.lower, "png", width=5)

## Chr12
c <- 12
chr <- chrs[c]
bed.gc.chr <- subset(bed.gc, CHR == chr)
nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]

## RFD   
load(file=file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.d.rt.RT.SLOPE_", chr, ".RData")))
nrds.RT.BSTRPS.chr$RFD <- getRFD(nrds.RT.BSTRPS.chr)
file.name <- file.path(wd.rt.plots, paste0("RFD_", base, "_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, "_TTR"))
plotBootstrapRFD(file.name, BASE, chr,  97500000, 105000000, nrds.chr, bed.gc.chr, nrds.RT.BSTRPS.chr, boundary.upper, boundary.lower, "png", width=5)
