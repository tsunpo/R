# =============================================================================
# Library      : Genomic Properties
# Name         : handbook-of/GenomicProperty.R
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 29/06/23
# =============================================================================
getRandomBreakpoint1 <- function() {
	  return(sample(c("start1", "end1"), 1, replace=F))
}

getRandomBreakpoint2 <- function() {
	  return(sample(c("start2", "end2"), 1, replace=F))
}

getGenomicProperty <- function(chr, bp, rt) {
	  rt.chr <- subset(rt, CHR == chr)
	
	  rt.chr.start <- subset(rt.chr, START <= bp)
	  rt.chr.start.end <- subset(rt.chr.start, END >= bp)
	
	  if (nrow(rt.chr.start.end) == 0)
	  	  return(NA)
	  else if (nrow(rt.chr.start.end) == 1)
	  	  return(rt.chr.start.end$BED)
	  else
		    return(paste(rt.chr.start.end$BED, collapse=","))
}

getGenomicPropertyCount <- function(chr, bp, rt) {
	  rt.chr <- subset(rt, CHR == chr)
	
	  rt.chr.start <- subset(rt.chr, START <= bp)
	  rt.chr.start.end <- subset(rt.chr.start, END >= bp)
	
	  return(nrow(rt.chr.start.end))
}

getDELRTNRFD <- function(del, nrds.RT.NRFD, bed.gc) {
	  chr <- paste0("chr", del$chrom1)
	  bed.gc.chr <- subset(bed.gc, CHR == chr)
	
	  bed.gc.chr.start <- subset(bed.gc.chr, START <= del$end2)
	  bed.gc.chr.start.end <- subset(bed.gc.chr.start, END >= del$start1)
	
  	return(intersect(rownames(bed.gc.chr.start.end), rownames(nrds.RT.NRFD)))
}

#getRandomBreakpointBED <- function(del, nrds.RT.NRFD, bed.gc) {
#	  chr <- paste0("chr", del$chrom1)
#	  bed.gc.chr <- subset(bed.gc, CHR == chr)
#	  breakpoints <- which(colnames(del) %in% c("start1", "end1", "start2", "end2"))
#	  breakpoint <- sample(breakpoints, 1, replace=F)
#	
#	  bed.gc.chr.start <- subset(bed.gc.chr, START <= del[, breakpoint])
#	  bed.gc.chr.start.end <- subset(bed.gc.chr.start, END >= del[, breakpoint])
#	
#	  overlaps <- intersect(rownames(bed.gc.chr.start.end), nrds.RT.NRFD$BED)
#	  if (length(overlaps) == 0)
#	     return(NA)
#	  else
#   	  return(overlaps)
#}

getPixel <- function(dels, nrds.RT.NRFD, bed.gc) {
	  dels$BED  
	
	for (d in 1:nrow(dels)) {
		    overlaps <- getRandomBreakpoint(dels[d,], nrds.RT.NRFD, bed.gc)
	    	pixels <- c(pixels, overlaps)
	  }
	  pixels <- pixels[order(pixels)]
	
	  return(nrds.RT.NRFD.nona[pixels, ])
}

toFrequencyTable <- function(partitions) {
	  test <- as.data.frame(table(partitions))
	  test <- test[order(test$Freq, decreasing=T),]
	  colnames(test) <- c("PIXEL", "Y")
	
	  return(test)
}

# =============================================================================
# Methods: Density plot
# Last Modified: 30/05/23
# =============================================================================
plotDensity <- function(reals, file.name, col, main.text, xlab.text="", showMedian=F, min=NA, max=NA, rt=NA, rev=F) {
	  ylab.text <- "Density"
	  d <- density(reals)
	  xlim <- c(min(reals), max(reals))
	  if (!is.na(max))
	  	  xlim <- c(-max, max)
	  if (!is.na(min))
	    	xlim <- c( min, max)
	  if (rev)
	  	  xlim <- rev(range(reals))
	
	  pdf(file.name, height=6, width=6)
	  par(mar=c(5.1, 4.7, 4.1, 1.4))
	  plot(d, xlab=xlab.text, ylab=ylab.text, main=main.text, xlim=xlim, col=col, lwd=3, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
	
	  if (!is.na(rt))
	  	  abline(v=rt, lty=5, lwd=2)
	  if (showMedian)
	  	  abline(v=median(reals), col=col, lty=5, lwd=3)

	  dev.off()
}

plotDensity2 <- function(reals, randoms, file.name, cols, legends, main.text, xlab.text="", showMedian=F, min=NA, max=NA, rt=NA, rev=F) {
	  ylab.text <- "Density"
	  d <- density(reals)
	  d2 <- density(randoms)
	  ylim <- c(min(c(d$y, d2$y)), max(c(d$y, d2$y)))
	  xlim <- c(min(c(reals, randoms)), max(c(reals, randoms)))
	  if (!is.na(max))
		    xlim <- c(-max, max)
	  if (!is.na(min))
	  	  xlim <- c( min, max)
	  if (rev)
	  	  xlim <- rev(range(c(reals, randoms)))
	  
	  pdf(file.name, height=6, width=6)
	  par(mar=c(5.1, 4.7, 4.1, 1.4))
	  plot(d, xlab=xlab.text, ylab=ylab.text, main=main.text, xlim=xlim, ylim=ylim, col=cols[1], lwd=3, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
	  lines(d2, col=cols[2], lwd=3)
	  
	  if (!is.na(rt))
	     abline(v=rt, lty=5, lwd=2)
	  if (showMedian)
		    abline(v=median(reals), col=cols[1], lty=5, lwd=3)
	
	  legend("topright", legend=legends, col=cols, lty=1, lwd=5, pt.cex=1.5, cex=1.8)
	  dev.off()
}

plotDensity3 <- function(reals, randoms, randoms3, file.name, cols, legends, main.text, xlab.text="", showMedian=F, min=NA, max=NA, rt=NA, rev=F) {
	  ylab.text <- "Density"
	  d <- density(reals)
	  d2 <- density(randoms)
	  d3 <- density(randoms3)
	  ylim <- c(min(c(d$y, d2$y, d3$y)), max(c(d$y, d2$y, d3$y)))
	  xlim <- c(min(c(reals, randoms, randoms3)), max(c(reals, randoms, randoms3)))
	  if (!is.na(max))
		    xlim <- c(-max, max)
	  if (!is.na(min))
		    xlim <- c( min, max)
	  if (rev)
		    xlim <- rev(range(c(reals, randoms, randoms3)))
	
	  pdf(file.name, height=6, width=6)
	  par(mar=c(5.1, 4.7, 4.1, 1.4))
	  plot(d, xlab=xlab.text, ylab=ylab.text, main=main.text, xlim=xlim, ylim=ylim, col=cols[1], lwd=3, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
	  lines(d2, col=cols[2], lwd=3)
	  lines(d3, col=cols[3], lwd=3)
	
	  if (!is.na(rt))
		    abline(v=rt, lty=5, lwd=2)
	  if (showMedian)
		    abline(v=median(reals), col=cols[1], lty=5, lwd=3)
	
	  legend("topright", legend=legends, col=cols, lty=1, lwd=5, pt.cex=1.5, cex=1.8)
	  dev.off()
}

plotDensityWilcox <- function(reals, randoms, file.name, col, main.text, xlab.text, showMedian=F, max=2) {
	  ylab.text <- "Density"
	  pval <- testU(reals, randoms)
	
	  d <- density(reals)
	  xlim <- c(min(reals), max(reals))
	  if (!is.na(max))
		    xlim <- c(min(-max), max(max))
	  pdf(file.name, height=6, width=6)
	  par(mar=c(5.1, 4.7, 4.1, 1.4))
	  plot(d, xlab=xlab.text, ylab=ylab.text, main=main.text, xlim=xlim, col=col, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
	
	  if (showMedian) {
		    abline(v=median(randoms), col="black", lty=5, lwd=3)
		    abline(v=median(reals), col=col, lty=5, lwd=3)
	  }
	
	  text(0, (max(d$y) + min(d$y))/2, getPvalueSignificanceLevel(pval), col="black", cex=5)
	  dev.off()
}

getMonteCarloSimulations <- function(nrds.RT.2, size) {
	  random.idx <- sort(sample(1:nrow(nrds.RT.2), size, replace=T))
	
	  return(median(nrds.RT.2[random.idx, "RT"]))
}

plotDensityMonteCarlo <- function(reals, nrds.RT.2, file.name, col, main.text, xlab.text, showMedian=F, max=2) {
	  randoms <- replicate(1000, getMonteCarloSimulations(nrds.RT.2, length(reals)))
	
	  ylab.text <- "Density"
	  ranks <- c(randoms, median(reals))
	  pval <- sum(ranks >= median(reals)) / length(ranks)
	  if (median(reals) < median(randoms))
	  	  pval <- sum(ranks <= median(reals)) / length(ranks)
	  
	  d <- density(reals)
	  xlim <- c(min(reals), max(reals))
	  if (!is.na(max))
	  	  xlim <- c(min(-max), max(max))
	  pdf(file.name, height=6, width=6)
	  par(mar=c(5.1, 4.7, 4.1, 1.4))
	  plot(d, xlab=xlab.text, ylab=ylab.text, main=main.text, xlim=xlim, col=col, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
	
	  if (showMedian) {
	     abline(v=median(randoms), col="black", lty=5, lwd=3)
	     abline(v=median(reals), col=col, lty=5, lwd=3)
   }

	  text(0, (max(d$y) + min(d$y))/2, getPvalueSignificanceLevel(pval), col="black", cex=5)
	  dev.off()
}





plotMonteCarloSimulation <- function(reals, randoms, file.name, cols, main.text, xlab.text) {
	  ylab.text <- "Density"
	  ranks <- c(randoms, median(reals))
	  d <- density(ranks)
	
	  col <- cols[1]
	  if (median(reals) < 0)
	  	  col <- cols[2]
	  
	  pdf(file.name, height=6, width=6)
	  par(mar=c(5.1, 4.7, 4.1, 1.4))
	  plot(d, xlab=xlab.text, ylab=ylab.text, main=main.text, xlim=c(min(ranks), max(ranks)), col="black", cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
	
	  abline(v=0, col="black", lty=5, lwd=3)
	  abline(v=median(reals), col=col, lty=5, lwd=3)
	  #text(median(reals), max(d$y), paste0("Median breakpoints"), cex=1.8, col="black")
			
			dev.off()
}

plotPropertyDensity0 <- function(reals, file.name, cols, main.text, xlab.text) {
	  ylab.text <- "Density"
	  d <- density(reals)

	  pdf(file.name, height=6, width=6)
	  par(mar=c(5.1, 4.7, 4.1, 1.4))
	  plot(d, xlab=xlab.text, ylab=ylab.text, main=main.text, xlim=c(min(reals), max(reals)), col="black", cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
	
	  #abline(v=0, col="black", lty=5, lwd=3)

	  dev.off()
}

plotPropertyDensity <- function(reals, randoms, file.name, cols, main.text, xlab.text) {
	  ylab.text <- "Density"
   d <- density(reals)
   
   ranks <- c(randoms, median(reals))
   pval <- sum(ranks >= median(reals)) / length(ranks)
   
   col <- cols[1]
   q <- as.numeric(quantile(tests))
   if (q[3] < 0)
   	  col <- cols[2]

   pdf(file.name, height=6, width=6)
   par(mar=c(5.1, 4.7, 4.1, 1.4))
   plot(d, xlab=xlab.text, ylab=ylab.text, main=main.text, xlim=c(min(reals), max(reals)), col="black", cex.axis=1.7, cex.lab=1.8, cex.main=1.9)

   #abline(v=median(reals), col=col, lty=5, lwd=3)
   
   text(0, (max(d$y) + min(d$y))/2, getPvalueSignificanceLevel(pval), col="black", cex=5)
   dev.off()
}

plotPropertyDensity2 <- function(tests, randoms, file.name, cols, main.text, xlab.text) {
	  # <- "Size [kb]"
	  ylab.text <- "Density"
	  d <- density(c(tests, randoms))
	
	  col <- cols[1]
	  q <- as.numeric(quantile(tests))
	  if (q[3] < 0)
		    col <- cols[2]
	
	  pdf(file.name, height=6, width=6)
	  par(mar=c(5.1, 4.7, 4.1, 1.4))
	  plot(d, xlab=xlab.text, ylab=ylab.text, main=main.text, col="black", cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
	
	  if (length(randoms) > 0) {
		    p <- ks.test(c(tests, randoms), randoms)$p.value
		    text(0, (max(d$y) + min(d$y))/2, getPvalueSignificanceLevel(p), col="black", cex=2.5)
	  }
	
	  #rug(jitter(sizes))
	  dev.off()
}


# =============================================================================
# Methods: Quantile skew plot
# Last Modified: 04/06/23
# =============================================================================
plotQuantileSkew <- function(pixels, file.name, cols, main.text) {
	  pdf(file.name, height=6, width=6)
	  par(mar=c(5.1, 4.7, 4.1, 1.4))
	  
  	plot(pixels$X, pixels$Y, xlim=c(-1, 1), ylab="", xlab="", main=main.text, col=cols[1], cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
  	lines(x=pixels$X, y=pixels$X, type="l", lwd=3, col=cols[1])

	  dev.off()
}
