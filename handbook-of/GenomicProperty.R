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
# 
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

plotDNS <- function(dns.f, file.name, main.text, ylab.text, ylim, legend, legends, cex=2.5, size=5, cols="black") {
	  xlab.text <- "Chromosome"
	
  	pdf(paste0(file.name, ".pdf"), height=size, width=size)
  	par(mar=c(5.1, 4.6, 4.1, 1.5))
	  plot(NULL, xlim=c(1, 22), ylim=ylim, xlab=xlab.text, ylab=ylab.text, main=main.text, col=cols[3], xaxt="n", pch=19, cex.axis=1.8, cex.lab=1.9, cex.main=2)
	
	  points(dns.f$Freq ~ dns.f$CHR, col=cols, pch=19, cex=cex)
	  lines(dns.f$Freq, y=NULL, lty=5, lwd=2.5, col=cols[1])
	
	  axis(side=1, at=seq(2, 22, by=4), cex.axis=1.8)
	  axis(side=1, at=seq(4, 20, by=4), cex.axis=1.8)
	  legend(legend, legend=legends, col=cols[1], lty=2, lwd=5, pt.cex=2, cex=1.9)
	  dev.off()
}

# -----------------------------------------------------------------------------
# 
# Last Modified: 04/10/23
# -----------------------------------------------------------------------------
getDNS2 <- function(gel.del.nona, dist, phase) {
	  dns <- as.data.frame(table(subset(subset(gel.del.nona, Telo < dist), V21 == phase)$CHR))
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

plotDNS2Freq <- function(dns.f, dns.m, file.name, main.text, ylab.text, ylim, legend, legends, cex=2.5, size=5, cols=c("black", "gray", "black")) {
	  xlab.text <- "Chromosome"
	
	  pdf(paste0(file.name, ".pdf"), height=5, width=6)
	  par(mar=c(5.1, 4.6, 4.1, 1.5))
	  plot(NULL, xlim=c(1, 22), ylim=ylim, xlab=xlab.text, ylab=ylab.text, yaxt="n", main=main.text, col=cols[3], xaxt="n", pch=19, cex.axis=1.9, cex.lab=2, cex.main=2.2)
	
	  points(dns.f$Freq/total.f*100 ~ dns.f$CHR, col=cols[2], pch=19, cex=cex)
	  lines(dns.f$Freq/total.f*100, y=NULL, lty=5, lwd=2.5, col=cols[2])
	
	  points(dns.m$Freq/total.m*100 ~ dns.m$CHR, col=cols[1], pch=19, cex=cex)
	  lines(dns.m$Freq/total.m*100, y=NULL, lty=5, lwd=2.5, col=cols[1])
	
	  #axis(side=2, at=12, cex.axis=1.8)	  
	  axis(side=2, at=seq(0, 12, by=2), labels=c(0,"",4,"",8,"",12), cex.axis=1.9)
	  #axis(side=2, at=seq(0, 7, by=1), labels=c(0,"",2,"",4,"",6, ""), cex.axis=1.9)
	  axis(side=1, at=seq(2, 22, by=4), cex.axis=1.9)
	  axis(side=1, at=seq(4, 20, by=4), cex.axis=1.9)
	  #legend("topleft", legend=legends[1:2], col=cols[1:2], pch=19, lty=2, lwd=5, pt.cex=2, cex=2, bty="n")
	  dev.off()
}

plotDNS2 <- function(dns.f, dns.m, file.name, main.text, ylab.text, ylim, legend, legends, cex=2.5, size=5, cols=c(blue.lighter, red.lighter, "black")) {
	  xlab.text <- ""
	
	  pdf(paste0(file.name, ".pdf"), height=size, width=size)
	  par(mar=c(5.1, 4.6, 4.1, 1.5))
	  plot(NULL, xlim=c(1, 22), ylim=ylim, xlab=xlab.text, ylab=ylab.text, main=main.text, col=cols[3], xaxt="n", pch=19, cex.axis=1.8, cex.lab=1.9, cex.main=2)
	
	  points(dns.f$Freq ~ dns.f$CHR, col=cols[1], pch=19, cex=cex)
	  lines(dns.f$Freq, y=NULL, lty=5, lwd=2.5, col=cols[1])
	
	  points(dns.m$Freq ~ dns.m$CHR, col=cols[2], pch=19, cex=cex)
	  lines(dns.m$Freq, y=NULL, lty=5, lwd=2.5, col=cols[2])
	
	  axis(side=2, at=12, cex.axis=1.8)	  
	  axis(side=1, at=seq(2, 22, by=4), cex.axis=1.8)
	  axis(side=1, at=seq(4, 20, by=4), cex.axis=1.8)
	  #legend(legend, legend=legends[2:1], col=cols[1:2], lty=2, lwd=5, pt.cex=2, cex=1.9)
	  dev.off()
}

# -----------------------------------------------------------------------------
# 
# Last Modified: 04/10/23
# -----------------------------------------------------------------------------
getCounts <- function(sv.del.nona.chr16, sv.del.nona.others, histology) {
	  colnames <- names(table(c(sv.del.nona.chr16[, histology], sv.del.nona.others[, histology])))
	
	  size.1 <- toTable(0, length(colnames), 2, colnames)
	  rownames(size.1) <- c("Others", "Chr16")
	
	  gel.del.m <- sv.del.nona.others
	  size.1[1, ] <- table(gel.del.m[, histology])[colnames]
	  gel.del.m <- sv.del.nona.chr16
	  size.1[2, ] <- table(gel.del.m[, histology])[colnames]
	  size.1 <- as.data.frame(t(size.1))
	  
	  idx <- as.numeric(as.vector(which(is.na(size.1[,1]))))
	  size.1$Others[idx] <- 0
	  idx <- as.numeric(as.vector(which(is.na(size.1[,2]))))
	  size.1$Chr16[idx] <- 0
	  
	  size.1 <- size.1[order(size.1[,2]),]
	  return(as.matrix(size.1))
}

getProportions <- function(size.1) {
	  size.2 <- size.1
	  size.2[,1] <- size.2[,1] / sum(size.2[,1]) * 100
	  size.2[,2] <- size.2[,2] / sum(size.2[,2]) * 100
	  
	  return(size.2)
}

plotProportions <- function(file.name, main.text, xlab.text, ylab.text, labels, counts, counts.prop, outs=c(36), outs.col=c(red), height=6.1, cutoff=20) {
	  grays <- c("lightgray", "darkgray", "dimgray")
	  #grays <- c("dimgray", "darkgray", "lightgray")
	  cols <- grays[rep_len(1:3, nrow(counts))]
	
	  blues <- c(blue, blue.light)
	  idx <- as.numeric(which(counts.prop[,2] == 0))
	  cols[idx] <- blues[rep_len(1:2, length(idx))]
	  
	  if (length(outs) != 0) {
	  	  for (o in 1:length(outs))
	  	  	  cols[outs[o]] <- outs.col[o]
	  }
	  
	  par(xpd=T)
	  pdf(paste0(file.name, ".pdf"), height=height, width=6.8)
	  par(mar=c(5.1, 4.6, 4.2, 18), xpd=TRUE)
	  barplot(counts.prop, col=cols, ylim=c(0, 100), ylab=ylab.text, xaxt="n", main=main.text, cex.names=1.8, cex.axis=1.8, cex.lab=1.9, cex.main=2)
	  #text(labels=c("Others ", "Chr16 "), x=c(0.8, 2), y=par("usr")[3] - 4, srt=45, adj=0.965, xpd=NA, cex=1.8)
	  axis(side=1, at=1-0.3,     labels=labels[1], font=1, cex.axis=1.9)
	  axis(side=1, at=2-0.3/2/2, labels=labels[2], font=1, cex.axis=1.9)
	  axis(side=1, at=1-0.3,     labels=paste0("n=", sum(counts[,1])), line=1.8, cex.axis=1.9, col.ticks="white")
	  axis(side=1, at=2-0.3/2/2, labels=paste0("n=", sum(counts[,2])), line=1.8, cex.axis=1.9, col.ticks="white")
	
	  idx <- as.numeric(which(counts[,2] >= cutoff))
	  for (c in 1:ncol(counts))
	  	  #text(c - 0.3/c/c, counts.prop[1, c]/2, counts[1, c], cex=1.8)
		    for (r in 1:nrow(counts))
		    	  if (r %in% idx)
		    	  	  if (counts[r, c] > 0)
			   	        text(c - 0.3/c/c, sum(counts.prop[r-1:r, c]) + (counts.prop[r, c]/2), counts[r, c], cex=1.8)
	
	  #text(4.3, 86, expression(italic('P')~"                   "), cex=2)
	  #text(4.3, 86, paste0("   = ", scientific(fisher.test(counts)[[1]])), cex=2)
	  legend("right", c(rev(rownames(counts.prop)[idx]), paste0("Not in ", labels[2])), text.col="black", pch=15, col=c(rev(cols[idx]), blue), pt.cex=3, cex=1.9, horiz=F, bty="n", inset=c(-1.6, 0))
	  dev.off()
}

plotMH2 <- function(file.name, main.text, mh, gel.del.nona.mother.others, gel.del.nona.mother.chr16, legends, cols) {
	  mh.others <- as.data.frame(table(gel.del.nona.mother.others$Microhomolgy_length[which(!is.na(gel.del.nona.mother.others$Microhomolgy_length))]))
	  mh.chr16  <- as.data.frame(table(gel.del.nona.mother.chr16$Microhomolgy_length[which(!is.na(gel.del.nona.mother.chr16$Microhomolgy_length))]))
	
	  pdf(file.path(wd.rt.plots, paste0(file.name, ".pdf")), height=4, width=8)
	  par(mar=c(5.1, 4.8, 4.1, 1.3))
	  if (!is.null(mh))
	  	  plot(mh$Var1, mh$Freq, bty="n", pch="", xaxt="n", ylab="Frequency", xlab="Bases of microhomology [bp]", main=main.text, cex.axis=1.8, cex.lab=1.9, cex.main=2)
	  else
	  	  plot(mh.others$Var1, mh.others$Freq, bty="n", pch="", xaxt="n", ylab="Frequency", xlab="Bases of microhomology [bp]", main=main.text, cex.axis=1.8, cex.lab=1.9, cex.main=2)
	 
	  points(mh.others$Var1, mh.others$Freq, col=cols[1], pch=16, cex=2)
	  points(mh.chr16$Var1, mh.chr16$Freq, col=cols[2], pch=16, cex=2)
	
	  legend("topright", legend=c(paste0(legends[1], " (n=", nrow(gel.del.nona.mother.others), ")"), paste0(legends[2], " (n=", nrow(gel.del.nona.mother.chr16), ")")), col=cols, pch=19, pt.cex=2.5, cex=1.8)
	  axis(side=1, at=seq(1, 71, by=5), labels=c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70), cex.axis=1.8)
	  dev.off()
}

# -----------------------------------------------------------------------------
# 
# Last Modified: 09/11/23
# -----------------------------------------------------------------------------
getPublicRFD <- function(file.name, rfd) {
	  orm.iz <- readTable(file.name, header=F, rownames=F, sep="\t")
	  colnames(orm.iz) <- c("CHR", "START", "END")
	  orm.iz$BED <- paste0("B", rownames(orm.iz))
	
	  orm.iz$SIZE <- orm.iz$END - orm.iz$START
	  orm.iz$TTR <- 0
	  orm.iz$IZ  <- 0
	  orm.iz$TZ  <- 0
	  for (r in 1:nrow(orm.iz)) {
		    chr <- orm.iz$CHR[r]
		    rfd.chr <- subset(rfd, CHR == chr)
		
		    rfd.chr.start <- subset(rfd.chr, START <= orm.iz$END[r])
		    rfd.chr.start.end <- subset(rfd.chr.start, END >= orm.iz$START[r])
		
		    if (nrow(rfd.chr.start.end) != 0) {
			      freq <- as.data.frame(table(rfd.chr.start.end$BRFD))
			
			      for (f in 1:nrow(freq)) {
				        orm.iz[r, as.vector(freq$Var1[f])] <- as.numeric(freq$Freq[f])
			      }
		    }
	  }
	
	  return(orm.iz)
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

getGenomicProperty2 <- function(chr, start, end, rt) {
	  rt.chr <- subset(rt, CHR == chr)
	
	  rt.chr.start <- subset(rt.chr, START <= end)
	  rt.chr.start.end <- subset(rt.chr.start, END >= start)
	
	  if (nrow(rt.chr.start.end) == 0)
		    return(NA)
	  else
	  	  return(rt.chr.start.end)
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

# -----------------------------------------------------------------------------
# Methods: 
# Last Modified: 02/11/23
# -----------------------------------------------------------------------------
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
plotDensity0_1 <- function(reals, file.name, col, main.text, xlab.text="", showMedian=F, min=NA, max=NA, rt=NA, rev=F, ymax=NA, scale=F) {
	  ylab.text <- "Density"
	  reals.2 <- reals[!is.na(reals)]
	  if (scale) {
	     reals.2 <- scale(reals)
	     reals.2 <- (reals.2-min(reals.2))/(max(reals.2)-min(reals.2))
	  }
	  d <- density(reals.2)
	
	  pdf(file.name, height=6, width=6)
	  par(mar=c(5.1, 4.7, 4.1, 1.4))
	  plot(d, xlab=xlab.text, ylab=ylab.text, main=main.text, col=col, lwd=3, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)

	  dev.off()
}

plotDensity <- function(reals, file.name, col, main.text, xlab.text="", showMedian=F, min=NA, max=NA, rt=NA, rev=F, ymax=NA) {
	  ylab.text <- "Density"
	  reals <- reals[!is.na(reals)]
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
	
	  xlim <- c(0, 1)
	  
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

plotBox2 <- function(wd.de.plots, file.name, tpm.1, tpm.2, main, names, cols, ylab, height=4, width=4) {
	  trait <- rep(0, length(tpm.1))
	  trait <- c(trait, rep(1, length(tpm.2)))
	  trait <- as.factor(trait)
	  expr <- as.numeric(c(tpm.1, tpm.2))
	  ylim <- c(min(expr), max(expr))
	  #ylim <- c(-3, 3)
	  
	  pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=height, width=width)
	  par(mar=c(5.1, 4.7, 4.1, 1.4))
	  boxplot(expr ~ trait, outline=F, xaxt="n", yaxt="n", xlab="", ylab=ylab, ylim=ylim, main=main, col=cols, cex.axis=1.8, cex.lab=1.9, cex.main=2)
	  #abline(h=0, col="black", lty=5, lwd=3)
	  
	  p <- testU(tpm.1, tpm.2)
	  #text(1.5, ylim[2]-1.4, getPvalueSignificanceLevel(p), cex=2.5)
	  #lines(c(1, 2), y=c(ylim[2]-2, ylim[2]-2), type="l", lwd=2)

	  text(1.5, -0.1, expression(italic('P')~"                   "), cex=1.9)
	  text(1.5, -0.1, paste0("   = ", scientific(p)), cex=1.9)
	  if (main[1] == "Microhomolgy") {
	     text(1.5, 28, expression(italic('P')~"                   "), cex=1.9)
	     text(1.5, 28, paste0("   = ", scientific(p)), cex=1.9)
	  } else if (main[1] == "Mother's age at birth") {
	     text(1.5, 39, expression(italic('P')~"                   "), cex=1.9)
	     text(1.5, 39, paste0("   = ", scientific(p)), cex=1.9)
	  } else if (main[1] == "Father's age at birth") {
	     text(1.5, 51, expression(italic('P')~"                   "), cex=1.9)
	     text(1.5, 51, paste0("   = ", scientific(p)), cex=1.9)
	  } else if (main[1] == "Dist. to telomere") {
	     text(1.5, ylim[2]-1, expression(italic('P')~"                   "), cex=1.9)
	     text(1.5, ylim[2]-1, paste0("   = ", scientific(p)), cex=1.9)
	  } else if (main[1] == "Replication timing") {
	  	  #text(1.5, ylim[1]+0.125, expression(italic('P')~"                   "), cex=1.9)
	  	  #text(1.5, ylim[1]+0.125, paste0("   = ", scientific(p)), cex=1.9)
   }
	  axis(side=2, at=0, labels=0, font=1, cex.axis=1.8)
	  axis(side=2, at=seq(-0.2, 0.2, by=0.05), labels=c(-0.2, -0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2), cex.axis=1.8)
	  axis(side=2, at=-0.15, labels=-0.15, font=1, cex.axis=1.8)
	  axis(side=1, at=1, labels=names[1], font=1, cex.axis=1.8)
	  axis(side=1, at=2, labels=names[2], font=1, cex.axis=1.8)
	  axis(side=1, at=1, labels=paste0("n=", separator(length(tpm.1)/2)), line=1.8, cex.axis=1.8, col.ticks="white")
	  axis(side=1, at=2, labels=paste0("n=", separator(length(tpm.2)/2)), line=1.8, cex.axis=1.8, col.ticks="white")
	  dev.off()
}

plotBox2RT <- function(wd.de.plots, file.name, tpm.1, tpm.2, main, names, cols, ylab, height=5, width=4) {
	trait <- rep(0, length(tpm.1))
	trait <- c(trait, rep(1, length(tpm.2)))
	trait <- as.factor(trait)
	expr <- as.numeric(c(tpm.1, tpm.2))
	ylim <- c(min(expr), max(expr))

	pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=height, width=width)
	par(mar=c(5.1, 4.7, 4.1, 1.4))
	boxplot(expr ~ trait, outline=F, xaxt="n", yaxt="n", xlab="", ylab=ylab, ylim=ylim, main=main, col=cols, cex.axis=1.9, cex.lab=2, cex.main=2.2)

	p <- testU(tpm.1, tpm.2)
	text(1.5, -0.1, expression(italic('P')~"                   "), cex=2)
	text(1.5, -0.1, paste0("   = ", scientific(p)), cex=2)

	#axis(side=2, at=0, labels=0, font=1, cex.axis=1.9)
	axis(side=2, at=seq(-0.2, 0.2, by=0.05), labels=c(-0.2, -0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2), cex.axis=1.9)
	axis(side=2, at=-0.15, labels=-0.15, font=1, cex.axis=1.9)
	axis(side=1, at=1, labels=names[1], font=1, cex.axis=1.9)
	axis(side=1, at=2, labels=names[2], font=1, cex.axis=1.9)
	#axis(side=1, at=1, labels=paste0("n=", separator(length(tpm.1)/2)), line=1.8, cex.axis=1.8, col.ticks="white")
	#axis(side=1, at=2, labels=paste0("n=", separator(length(tpm.2)/2)), line=1.8, cex.axis=1.8, col.ticks="white")
	dev.off()
}


plotBox2Telo <- function(wd.de.plots, file.name, tpm.1, tpm.2, main, names, cols, ylab, height=5, width=4) {
	trait <- rep(0, length(tpm.1))
	trait <- c(trait, rep(1, length(tpm.2)))
	trait <- as.factor(trait)
	expr <- as.numeric(c(tpm.1, tpm.2))
	ylim <- c(min(expr), max(expr))

	pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=height, width=width)
	par(mar=c(5.1, 4.7, 4.1, 1.4))
	boxplot(expr ~ trait, outline=F, xaxt="n", xlab="", ylab=ylab, ylim=ylim, main=main, col=cols, cex.axis=1.9, cex.lab=2, cex.main=2.2)

	p <- testU(tpm.1, tpm.2)
		text(1.5, ylim[2]-1, expression(italic('P')~"                   "), cex=2)
		text(1.5, ylim[2]-1, paste0("   = ", scientific(p)), cex=2)

	axis(side=1, at=1, labels=names[1], font=1, cex.axis=1.9)
	axis(side=1, at=2, labels=names[2], font=1, cex.axis=1.9)
	axis(side=1, at=1, labels=paste0("n=", separator(length(tpm.1)/2)), line=1.8, cex.axis=1.9, col.ticks="white")
	axis(side=1, at=2, labels=paste0("n=", separator(length(tpm.2)/2)), line=1.8, cex.axis=1.9, col.ticks="white")
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

# =============================================================================
# Methods: Density plot
# Last Modified: 15/01/24
# =============================================================================
getMonteCarloSimulations <- function(nrds.RT.2, column, size) {
	random.idx <- sort(sample(1:nrow(nrds.RT.2), size, replace=T))
	
	return(median(nrds.RT.2[random.idx, column]))
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

getMonteCarloSimulations <- function(nrds.RT.2, column, size) {
	  random.idx <- sort(sample(1:nrow(nrds.RT.2), size, replace=T))
	
	  return(median(nrds.RT.2[random.idx, column]))
}

plotDensityMonteCarlo2 <- function(reals, nrds.RT.2, column, file.name, col, main.text, xlab.text, ylab.text, showMedian=F, max=NA, rev=F) {
	  reals <- reals[!is.na(reals)]
	  nrds.RT.2.nona <- nrds.RT.2[!is.na(nrds.RT.2$VALUE),]
	  randoms <- replicate(1000, getMonteCarloSimulations(nrds.RT.2.nona, column, length(reals)))
	
	  ranks <- c(randoms, median(reals))
	  pval <- sum(ranks >= median(reals)) / length(ranks)
	  if (median(reals) < median(randoms))
		    pval <- sum(ranks <= median(reals)) / length(ranks)
	
	  d <- density(reals)
	  xlim <- c(0, 1)
	  if (rev)
		     xlim <- rev(xlim)
	  #xlim <- c(0.2, 0.7)   ## GC contents
	  pdf(file.name, height=4, width=4)
	  par(mar=c(5.1, 4.7, 4.1, 1.4))
	  plot(d, xlab=xlab.text, ylab=ylab.text, main=main.text, xlim=xlim, col=col, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
	
	  text((xlim[1] + xlim[2])/2, (max(d$y) + min(d$y))/2, getPvalueSignificanceLevel(pval), col="black", cex=5)
	  dev.off()
}

scale_values <- function(x){(x-min(x))/(max(x)-min(x))}

plotDensityMonteCarlo <- function(reals, vs, column, file.name, col, main.text, xlab.text, ylab.text, showMedian=F, max=NA, rev=F) {
	  reals <- reals[!is.na(reals)]
	  vs.nona <- vs[!is.na(vs$VALUE),]
	  randoms <- replicate(1000, getMonteCarloSimulations(vs.nona, column, length(reals)))
	
	  ranks <- c(randoms, median(reals))
	  pval <- sum(ranks >= median(reals)) / length(ranks)
	  if (median(reals) < median(randoms))
	  	  pval <- sum(ranks <= median(reals)) / length(ranks)
	  
	  d <- density(reals)
	  xlim <- c(min(reals), max(reals))
	  if (!is.na(max))
	  	  xlim <- c(min(-max), max(max))
	  if (rev)
	  	  xlim <- rev(xlim)
	  #xlim <- c(0.2, 0.7)   ## GC contents
	  
	  #reals.2 <- scale(rank(reals), center=F)
	  reals.2 <- scale(reals, center=F)
	  reals.2 <- scale_values(reals.2)
	  median_normalized <- median(reals.2)
	  median_shift_normalized <- median_normalized - 0.5  # Assuming a uniform distribution has a median of 0.5
	  col.idx <- ceiling(median_shift_normalized/0.05)
	  if (median_shift_normalized < 0)
	  	  col.idx <- floor(median_shift_normalized/0.05)
	  #if (median(reals) < median(randoms))
	  #	  col.idx <- col.idx * -1
	  if (rev)
	  	  col.idx <- col.idx * -1
	  
	  pdf(file.name, height=4, width=4)
	  par(mar=c(5.1, 4.7, 4.1, 1.4))
	  plot(d, xlab=xlab.text, ylab=ylab.text, main=main.text, xlim=xlim, col=col, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
	
	  if (col.idx > 0)
	  	  polygon(d,	col=reds[col.idx])
	  else 
	     polygon(d,	col=blues[col.idx])
	  
	  text((xlim[1] + xlim[2])/2, (max(d$y) + min(d$y))/2, getPvalueSignificanceLevel(pval), col="black", cex=5)
	  #text((xlim[1] + xlim[2])/2, (max(d$y) + min(d$y))/2, median_shift_normalized, col="black", cex=1)
	  dev.off()
}

plotDensityKS <- function(reals, vs, column, file.name, col, main.text, xlab.text, ylab.text, showMedian=F, max=NA, rev=F) {
	  real_observations <- reals[!is.na(reals)]
	  vs <- vs[!is.na(vs$VALUE),]
	  
	  random.idx <- sort(sample(1:nrow(vs), 500000, replace=F))
	  random_observations <- vs[random.idx, column]

	  ## Pool real observations with random ones, rank-transform them, and normalize the values on a scale from 0 to 1
	  all_observations <- c(real_observations, random_observations)
	  all_observations <- unique(all_observations)
	  #ranked_observations <- rank(all_observations)
	  normalized_observations <- scale(all_observations, center=F)
	  normalized_observations <- scale_values(normalized_observations)
	  
	  ks_result <- ks.test(normalized_observations, "punif")
	  pval <- ks_result$p.value
	
	  d <- density(real_observations)
	  xlim <- c(min(real_observations), max(real_observations))
	  if (!is.na(max))
	  	  xlim <- c(min(-max), max(max))
	  if (rev)
	  	  xlim <- rev(xlim)
	  #xlim <- c(0.2, 0.7)   ## GC contents
	  
	  # Calculate the median shift using normalized and ranked observations on a scale from 0 to 1
	  reals.2 <- scale(real_observations, center=F)
	  reals.2 <- scale_values(reals.2)
	  median_normalized <- median(reals.2)
	  median_shift_normalized <- median_normalized - 0.5  # Assuming a uniform distribution has a median of 0.5
	  col.idx <- ceiling(median_shift_normalized/0.05)
	  if (median_shift_normalized < 0)
	  	  col.idx <- floor(median_shift_normalized/0.05)
	  #if (median(reals) < median(randoms))
	  #	  col.idx <- col.idx * -1
	  if (rev)
	  	  col.idx <- col.idx * -1
	  
	  pdf(file.name, height=4, width=4)
	  par(mar=c(5.1, 4.7, 4.1, 1.4))
	  plot(d, xlab=xlab.text, ylab=ylab.text, main=main.text, xlim=xlim, col=col, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)

	  if (col.idx > 0)
	  	polygon(d,	col=reds[col.idx])
	  else 
	  	polygon(d,	col=blues[col.idx])
	  
	  text((xlim[1] + xlim[2])/2, (max(d$y) + min(d$y))/2, getPvalueSignificanceLevel(pval), col="black", cex=5)
	  #text((xlim[1] + xlim[2])/2, (max(d$y) + min(d$y))/2, median_shift_normalized, col="black", cex=0.5)
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
