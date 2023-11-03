# =============================================================================
# Library      : Genomic Properties
# Name         : handbook-of/GenomicProperty.R
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 02/11/23
# =============================================================================

# -----------------------------------------------------------------------------
# Methods:
# Last Modified: 02/11/23
# -----------------------------------------------------------------------------
getArm <- function(chr, START) {
	  centromeres <- subset(cytoBand, gieStain == "acen")
	  colnames(centromeres)[1:3] <- c("CHR", "START", "END")
	  centromeres.chr <- subset(centromeres, CHR == chr)
	  centro <- centromeres.chr[1,]$END
	
	  if (START > centro) {
		    return("q")
	  } else {
		    return("p")
	  }
}

getTelo <- function(CHR, START, Arm) {
	  chromInfo.chr <- subset(chromInfo, chrom == CHR)
	
	  if (Arm == "q") {
		    return(chromInfo.chr$size - START)
	  } else {
		    return(START)
	  }
}

getCentro <- function(chr, START, Arm) {
	  centromeres.chr <- subset(centromeres, CHR == chr)
	  centro <- centromeres.chr[1,]$END
	
	  if (Arm == "q") {
		    return(START - centro)
	  } else {
		    return(centro - START )
	  }
}

# -----------------------------------------------------------------------------
# Methods: 
# Last Modified: 02/11/23
# -----------------------------------------------------------------------------
getDNS <- function(sv.del.nona, dist, phase) {
	  dns <- as.data.frame(table(subset(sv.del.nona, Telo < dist)$CHR))
	  rownames(dns) <- dns$Var1
	
	  dns.chr <- toTable(0, 2, 22, c("CHR", "Freq"))
	  for (c in 1:22) {
		    chr <- chrs[c]
		    dns.chr$CHR[c] <- c
		
		    if (length(which(rownames(dns) == chr)) != 0)
			      dns.chr$Freq[c] <- dns[chr,]$Freq
	  }
	
	  return(dns.chr)
}

plotDNS <- function(dns.f, file.name, main.text, ylab.text, ylim, legend, legends, cex=2.5, size=5) {
	  xlab.text <- "Chromosome"
  	cols <- c("black", blue.lighter, red.lighter)
	
  	pdf(paste0(file.name, ".pdf"), height=size, width=size)
  	par(mar=c(5.1, 4.6, 4.1, 1.5))
	  plot(NULL, xlim=c(1, 22), ylim=ylim, xlab=xlab.text, ylab=ylab.text, main=main.text, col=cols[3], xaxt="n", pch=19, cex.axis=1.8, cex.lab=1.9, cex.main=2)
	
	  points(dns.f$Freq ~ dns.f$CHR, col=cols[1], pch=19, cex=cex)
	  lines(dns.f$Freq, y=NULL, lty=5, lwd=2.5, col=cols[1])
	
	  axis(side=1, at=seq(2, 22, by=4), cex.axis=1.8)
	  axis(side=1, at=seq(4, 20, by=4), cex.axis=1.8)
	  legend(legend, legend=legends, col=cols[1], lty=2, lwd=5, pt.cex=2, cex=1.9)
	  dev.off()
}

# -----------------------------------------------------------------------------
# Methods: 
# Last Modified: 02/11/23
# -----------------------------------------------------------------------------
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

getIZ <- function(chr, start, end, orm.iz) {
	  rt.chr <- subset(orm.iz, CHR == chr)
	
	  rt.chr.start <- subset(rt.chr, START <= end)
	  rt.chr.start.end <- subset(rt.chr.start, END >= start)
	
	  if (nrow(rt.chr.start.end) != 0)
	  	  return(1)
	  else
	  	  return(0)
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
plotDensity <- function(reals, file.name, col, main.text, xlab.text="", showMedian=F, min=NA, max=NA, rt=NA, rev=F, ymax=NA) {
	  ylab.text <- "Density"
	  d <- density(reals)
	  
	  ylim <- c(min(d$y), max(d$y))
	  if (!is.na(ymax))
	  	  ylim <- c(0, ymax)
	  xlim <- c(min(reals), max(reals))
	  if (!is.na(max))
	  	  xlim <- c(-max, max)
	  if (!is.na(min))
	    	xlim <- c( min, max)
	  if (rev)
	  	  xlim <- rev(range(reals))
	
	  pdf(file.name, height=6, width=6)
	  par(mar=c(5.1, 4.7, 4.1, 1.4))
	  plot(d, xlab=xlab.text, ylab=ylab.text, main=main.text, xlim=xlim, ylim=ylim, col=col, lwd=3, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
	
	  if (!is.na(rt))
	  	  abline(v=rt, lty=5, lwd=2)
	  if (showMedian)
	  	  abline(v=median(reals), col=col, lty=5, lwd=3)

	  dev.off()
}

plotDensity2 <- function(reals, randoms, file.name, cols, legends, legend, main.text, xlab.text="", showMedian=F, min=NA, max=NA, rt=NA, rev=F, ymax=NA) {
	  ylab.text <- "Density"
	  d <- density(reals)
	  d2 <- density(randoms)
	  d$y <- d$y * length(reals) / (length(reals) + length(randoms))
	  d2$y <- d2$y * length(randoms) / (length(reals) + length(randoms))
	  
	  ylim <- c(min(c(d$y, d2$y)), max(c(d$y, d2$y)))
	  if (!is.na(ymax))
	  	  ylim <- c(0, ymax)
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
	
	  #abline(v=4, lty=5, lwd=3)
	  
	  legend(legend, legend=legends, col=cols, lty=1, lwd=5, pt.cex=1.5, cex=1.8)
	  dev.off()
}

plotDensity3 <- function(reals, randoms, randoms3, file.name, cols, legends, legend, main.text, xlab.text="", showMedian=F, min=NA, max=NA, rt=NA, rev=F) {
	  ylab.text <- "Density"
	  d <- density(reals)
	  d2 <- density(randoms)
	  d3 <- density(randoms3)
	  d$y <- d$y * length(reals) / (length(reals) + length(randoms) + length(randoms3))
	  d2$y <- d2$y * length(randoms) / (length(reals) + length(randoms) + length(randoms3))
	  d3$y <- d3$y * length(randoms3) / (length(reals) + length(randoms) + length(randoms3))
	  
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
	
	  legend(legend, legend=legends, col=cols, lty=1, lwd=5, pt.cex=1.5, cex=1.8)
	  dev.off()
}

plotBox2 <- function(wd.de.plots, file.name, tpm.1, tpm.2, main, names, cols, ylab, height=6, width=4) {
	  trait <- rep(0, length(tpm.1))
	  trait <- c(trait, rep(1, length(tpm.2)))
	  trait <- as.factor(trait)
	  expr <- as.numeric(c(tpm.1, tpm.2))
	  ylim <- c(min(expr), max(expr))
	  #ylim <- c(-3, 3)
	  
	  pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=height, width=width)
	  par(mar=c(5.1, 4.7, 4.1, 1.4))
	  boxplot(expr ~ trait, outline=F, xaxt="n", xlab="", ylab=ylab, ylim=ylim, main=main, col=cols, cex.axis=1.8, cex.lab=1.9, cex.main=2)
	  #abline(h=0, col="black", lty=5, lwd=3)
	  
	  p <- testU(tpm.1, tpm.2)
	  #text(1.5, ylim[2]-1.4, getPvalueSignificanceLevel(p), cex=2.5)
	  #lines(c(1, 2), y=c(ylim[2]-2, ylim[2]-2), type="l", lwd=2)

	  #text(1.5, 39, expression(italic('P')~"                   "), cex=1.9)
	  #text(1.5, 39, paste0("   = ", scientific(p)), cex=1.9)
	  text(1.5, 51, expression(italic('P')~"                   "), cex=1.9)
	  text(1.5, 51, paste0("   = ", scientific(p)), cex=1.9)
	  #text(1.5, ylim[2]-1, expression(italic('P')~"                   "), cex=1.9)
	  #text(1.5, ylim[2]-1, paste0("   = ", scientific(p)), cex=1.9)
	  #text(1.5, ylim[1]+0.125, expression(italic('P')~"                   "), cex=1.9)
	  #text(1.5, ylim[1]+0.125, paste0("   = ", scientific(p)), cex=1.9)

	  axis(side=2, at=0, labels=0, font=1, cex.axis=1.8)
	  axis(side=1, at=1, labels=names[1], font=1, cex.axis=1.8)
	  axis(side=1, at=2, labels=names[2], font=1, cex.axis=1.8)
	  axis(side=1, at=1, labels=paste0("n=", separator(length(tpm.1))), line=1.8, cex.axis=1.8, col.ticks="white")
	  axis(side=1, at=2, labels=paste0("n=", separator(length(tpm.2))), line=1.8, cex.axis=1.8, col.ticks="white")
	  dev.off()
}

plotCorrelationRTS <- function(file.name, main.text, xlab.text, ylab.text, dns, rts, pos="topright", cols=c("dimgray", "black"), size=5, pch=1, cex=1.9, chrs=c(4, 16)) {
	  x=dns$Freq
	  y=rts[dns$Var1,]$spr
	  sprs=rts[dns$Var1,]$spr
	  idxes <- which(dns$Var1 %in% chrs)
	
	  pdf(paste0(file.name, ".pdf"), height=size, width=size)
	  par(mar=c(5.1, 4.7, 4.1, 1.4))
	  plot(y ~ x, ylab=ylab.text, xlab=xlab.text, yaxt="n", main=main.text, pch=pch, cex=cex, col=cols[1], cex.axis=1.9, cex.lab=2, cex.main=2.1)
	
	  idx <- which(sprs >= 0)
	  points(y[idx] ~ x[idx], col=red, pch=19, cex=cex)
	  idx <- which(sprs < 0)
	  points(y[idx] ~ x[idx], col=blue, pch=19, cex=cex)
	  
	  lm.fit <- lm(y ~ x)
	  abline(lm.fit, lwd=3, col=cols[2])
	
	  cor <- cor.test(y, x, method="spearman", exact=F)
	  legend(pos, c(paste0("rho = ", round0(cor[[4]], digits=2)), paste0("P = 1.00E-00")), text.col=c(cols[2], "white"), text.font=1, bty="n", cex=2)
	  legend(pos, c("", expression(italic('P')~"                   ")), text.col=cols[2], text.font=1, bty="n", cex=2)
	  legend(pos, c("", paste0("   = ", scientific(cor[[3]]))), text.col=cols[2], text.font=1, bty="n", cex=2)
	
	  legend("topleft", legend=c("E > L", "L > E"), col=c(red, blue), pch=19, pt.cex=2.5, cex=1.8)
	  
	  par(xpd=T)
	  for (c in 1:length(idxes)) {
	  	  idx <- idxes[c]
	  	  text(x[idx], y[idx], paste0("Chr", chrs[c]), col="black", pos=3, cex=1.9)
	  }
	  axis(side=1, at=100, labels=100, cex.axis=1.9)
	  axis(side=2, at=seq(-0.5, 0.5, by=0.5), labels=c(-0.5, 0, 0.5), cex.axis=1.9)
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
