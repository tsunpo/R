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

getLog2RDRatio <- function(nrds.T.chr, nrds.N.chr, pseudocount) {
   return(log2(as.numeric(nrds.T.chr) + pseudocount) - log2(as.numeric(nrds.N.chr) + pseudocount))
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
setScaledRT <- function(nrds.chr, pseudocount, recaliRT, scaledRT) {
   nrds.chr$T <- log2(nrds.chr$T + pseudocount)
   nrds.chr$N <- log2(nrds.chr$N + pseudocount)
   
   if (recaliRT == T)
      nrds.chr$RT <- nrds.chr$T - nrds.chr$N   ## ADD 24/02/19

   if (scaledRT == T)
      nrds.chr$RT <- scale(nrds.chr$RT)            ## ADD 19/02/19
   
   return(nrds.chr)
}

setSpline <- function(nrds.chr, bed.gc.chr, column) {
   overlaps <- intersect(rownames(bed.gc.chr), nrds.chr$BED)
   bed.gc.chr.o <- bed.gc.chr[overlaps,]
   nrds.chr.o <- nrds.chr[overlaps,]
 
   spline <- smooth.spline(x=bed.gc.chr.o$START, y=nrds.chr.o[, column])
   nrds.chr.o$SPLINE <- spline$y
   #slopes <- diff(spline$y)/diff(bed.gc.chr$START/1E6)   ## ADD 31/10/18
   #nrds.chr <- nrds.chr[-nrow(nrds.chr),]    ## length(slopes) is 1 less than nrow(bed.gc.chr), as no slope for the last 1kb window
   #nrds.chr$SLOPE <- slopes
 
   sizes <- diff(bed.gc.chr.o$START)
   gaps <- which(sizes != 1000)
   return(nrds.chr.o[-gaps, c("BED", column, "SPLINE")])
}

plotRT <- function(file.name, main.text, chr, xmin, xmax, nrds.chr, bed.gc.chr, colours, legends, colours2, legends2, ext, width, peaks, ylim=NULL, lcl.rt.chr=NULL) {
   ## Colours (was "lightcoral", "skyblue3")
   ## http://r.789695.n4.nabble.com/plot-function-color-transparency-td4682424.html
   adjustcolor.red  <- adjustcolor(colours2[1], alpha.f=0.08)
   adjustcolor.blue <- adjustcolor(colours2[2], alpha.f=0.08)
   adjustcolor.gray <- adjustcolor("gray", alpha.f=0.08)
 
   nrds.chr.T  <- setSpline(nrds.chr, bed.gc.chr, "T")
   nrds.chr.N  <- setSpline(nrds.chr, bed.gc.chr, "N")
   nrds.chr.RT <- setSpline(nrds.chr, bed.gc.chr, "RT")
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
      plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=c(min(rds), max(rds)), xlab=xlab.text, ylab=ylab.text, main=main.text, xaxt="n")
   } else
      plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=ylim, xlab=xlab.text, ylab=ylab.text, main=main.text, xaxt="n")
   #points(bed.gc.chr$START/1E6, nrds.chr, col=colours[1], cex=0.3)
   abline(h=0, lwd=0.5, col="lightgrey")
 
   ## Plot cytobands (before smoothing spline)
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.4, col="lightgrey") 

   ## Plot smoothing splines
   bed.gc.chr <- bed.gc.chr[rownames(nrds.chr.RT),]   ## BUT HERE?
   points(bed.gc.chr$START/1E6, nrds.chr.T$SPLINE, col=colours[1], pch=16, cex=0.2)
   points(bed.gc.chr$START/1E6, nrds.chr.N$SPLINE, col=colours[2], pch=16, cex=0.2)
   
   ## Plot legend and peaks
   legend("bottomright", legends, col=colours, lty=1, lwd=2, bty="n", horiz=T)
   if (length(peaks) != 0)
      for (p in 1:length(peaks))
         abline(v=peaks[p]/1E6, lty=5, lwd=1, col="black")
   
   ### 
   ## Initiate RT plot
   par(mar=c(5.5,4,0,1))
   ylab.text <- "Replication timing"
   
   plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=c(-2, 2), xlab=xlab.text, ylab=ylab.text, main="", yaxt="n")
   idx <- which(nrds.chr.RT$RT == 0)
   points(bed.gc.chr[idx,]$START/1E6, nrds.chr.RT[idx,]$RT, col="lightgrey", cex=0.3)
   idx <- which(nrds.chr.RT$RT < 0)
   points(bed.gc.chr[idx,]$START/1E6, nrds.chr.RT[idx,]$RT, col=adjustcolor.blue, cex=0.3)
   idx <- which(nrds.chr.RT$RT > 0)
   points(bed.gc.chr[idx,]$START/1E6, nrds.chr.RT[idx,]$RT, col=adjustcolor.red, cex=0.3)
   
   abline(h=0, lwd=0.5, col="lightgrey")
   axis(side=2, at=seq(-2, 2, by=1), labels=c(-2, -1, 0, 1, 2))
  
   ## Plot cytobands (before smoothing spline)
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.4, col="lightgrey") 
  
   ## Plot Koren 2012
   if (!is.null(lcl.rt.chr))
      points(lcl.rt.chr$POS/1E6, lcl.rt.chr$RT, col="forestgreen", pch=16, cex=0.2)
   
   ## Plot smoothing spline
   points(bed.gc.chr$START/1E6, nrds.chr.RT$SPLINE, col="black", pch=16, cex=0.2)
   abline(h=0, lty=5, lwd=1, col="black")
   
   ## Plot legend and peaks
   legend("topleft", paste0("Early (", legends2[1], " > ", legends2[2], ")"), bty="n", text.col="black", pt.cex=0.7, pch=1, col=colours2[1])   
   legend("bottomleft", paste0("Late (", legends2[1], " < ", legends2[2], ")"), bty="n", text.col="black", pt.cex=0.7, pch=1, col=colours2[2])
   if (!is.na(xmin) && !is.na(xmax))
      legend("topright", paste0(legends2[1], "/", legends2[2], " read depth ratio"), col="black", lty=1, lwd=2, bty="n", horiz=T)
   else
      legend("topright", paste0(legends2[1], "/", legends2[2], " ratio"), col="black", lty=1, lwd=2, bty="n", horiz=T)
   if (!is.null(lcl.rt.chr))
      legend("bottomright", "Koren et al. 2012", col="forestgreen", lty=1, lwd=2, bty="n", horiz=T)
    
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

plotRD2vsRT <- function(reads1, reads2, timings, file.name, main.text, ylab.text, xlab.text, colours, legends, method, intercept=F) {
   #xmin <- min(nrds.chr.RT$SPLINE)
   #xmax <- max(nrds.chr.RT$SPLINE)
   #ymin <- min(c(nrds.chr.T$SPLINE, nrds.chr.N$SPLINE))
   #ymax <- max(c(nrds.chr.T$SPLINE, nrds.chr.N$SPLINE))
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
   plot(NULL, xlim=xlim, ylim=ylim, ylab=ylab.text, xlab=xlab.text, main=main.text[1])
   abline(v=0, lty=5)
   
   lm.fit1 <- lm(reads1 ~ timings)
   abline(lm.fit1, col=colours[1], lwd=3)
   cor1 <- getCor(reads1, timings, method)
   intercept1 <- lm.fit1[[1]][1]

   lm.fit2 <- lm(reads2 ~ timings)
   abline(lm.fit2, col=colours[2], lwd=3)
   cor2 <- getCor(reads2, timings, method)
   intercept2 <- lm.fit2[[1]][1]
   
   if (intercept) {
      points(0, intercept1, col=colours[1], pch=19)
      points(0, intercept2, col=colours[2], pch=19)
      if (intercept1 > intercept2) {
         text(0, intercept1, round0(intercept1, digits=3), cex=1.1, col=colours[1], pos=3)
         text(0, intercept2, round0(intercept2, digits=3), cex=1.1, col=colours[2], pos=1)
      } else {
         text(0, intercept1, round0(intercept1, digits=3), cex=1.1, col=colours[1], pos=1)
         text(0, intercept2, round0(intercept2, digits=3), cex=1.1, col=colours[2], pos=3)
      }
   }
   
   RT <- paste0(legends[1], "/", legends[2])
   if (cor1 < 0 && cor2 < 0) {
      legend("bottomright", c(paste0(cor, " = ", round0(cor1, digits=2), " (", legends[1], " vs. ", RT, ")"), paste0(cor, " = ", round0(cor2, digits=2), " (", legends[2], " vs. ", RT, ")")), text.col=colours, bty="n", cex=1.1)
   } else if (as.numeric(cor1) > 0 && as.numeric(cor2) < 0) {
      legend("topright", paste0(cor, " = ", round0(cor1, digits=2), " (", legends[1], " vs. ", RT, ")"), text.col=colours[1], bty="n", cex=1.1)        
      legend("bottomright", paste0(cor, " = ", round0(cor2, digits=2), " (", legends[2], " vs. ", RT, ")"), text.col=colours[2], bty="n", cex=1.1)
   }
   mtext(main.text[2], cex=1.2, line=0.3)
   dev.off()
}

plotRD2vsRTALL <- function(cors, file.name, main.text, ymin, ymax, cols, legends, c=NA) {
   ylab.text <- "Spearman's rho"
   xlab.text <- "Chromosome"
 
   #png(paste0(file.name, ".png"), height=5, width=5, units="in", res=300)
   pdf(paste0(file.name, ".pdf"), height=6, width=6)
   plot(cors$cor1 ~ cors$chr, ylim=c(ymin, ymax), ylab=ylab.text, xlab=xlab.text, main=main.text, col=cols[1], xaxt="n", yaxt="n", pch=19, cex.lab=1.5, cex.main=1.6)
   lines(cors$cor1, y=NULL, type="l", lwd=3, col=cols[1])
   points(cors$chr, cors$cor2, col=cols[2], pch=19)
   lines(cors$cor2, y=NULL, type="l", lwd=3, col=cols[2])
   abline(h=0, lty=5)
   
   RT <- paste0(legends[1], "/", legends[2])
   legend("topright", paste0(legends[1], " vs. ", RT), text.col=cols[1], bty="n", cex=1.6)        
   legend("bottomright", paste0(legends[2], " vs. ", RT), text.col=cols[2], bty="n", cex=1.6)
   
   if (!is.na(c)) {
      text(c, cors$cor1[c], round0(cors$cor1[c], digits=2), cex=1.6, col=cols[1], pos=3)   ##, offset=1.3)
      text(c, cors$cor2[c], round0(cors$cor2[c], digits=2), cex=1.6, col=cols[2], pos=3)
   }
   axis(side=1, at=seq(2, 22, by=2), cex.axis=1.25)
   axis(side=2, at=seq(-0.8, 0.8, by=0.4), labels=c(-0.8, -0.4, 0, 0.4, 0.8), cex.axis=1.4)
   dev.off()
}

plotRD3vsRTALL <- function(cors, file.name, main.text, ymin, ymax, cols, legends, c=NA, isRT=F) {
   ylab.text <- "Spearman's rho"
   xlab.text <- "Chromosome"
 
   #png(paste0(file.name, ".png"), height=5, width=5, units="in", res=300)
   pdf(paste0(file.name, ".pdf"), height=5, width=5)
   plot(cors$cor1 ~ cors$chr, ylim=c(ymin, ymax), ylab=ylab.text, xlab=xlab.text, main=main.text, col=cols[1], xaxt="n", yaxt="n", pch=19, cex.lab=1, cex.main=1.2)
   lines(cors$cor1, y=NULL, type="l", lwd=3, col=cols[1])
   points(cors$chr, cors$cor2, col=cols[2], pch=19)
   lines(cors$cor2, y=NULL, type="l", lwd=3, col=cols[2])
   if (isRT) {
      points(cors$chr, cors$cor, col="black", pch=19)
      lines(cors$cor, y=NULL, type="l", lwd=3, col="black")
   }
   abline(h=0, lty=5)
 
   RT <- "S/G1"
   #legend("topright", c(paste0(legends[3], " vs. ", RT), paste0("      ", legends[1], " vs. ", RT)), text.col=c(cols[3], cols[1]), bty="n", cex=1.1) 
   #legend("bottomright", paste0(legends[2], " vs. ", RT), text.col=cols[2], bty="n", cex=1.1) 
   legend("topright", paste0(legends[3], " vs. ", RT), text.col=cols[3], bty="n", cex=1.1) 
   legend("bottomright", c(paste0(legends[1], " vs. ", RT), paste0(legends[2], " vs. ", RT)), text.col=c(cols[1], cols[2]), bty="n", cex=1.1)        
   ##legend("bottomright", paste0(legends[2], " vs. ", RT), text.col=cols[2], bty="n", cex=1.1)
 
   if (!is.na(c)) {
      text(c, cors$cor1[c], round0(cors$cor1[c], digits=2), cex=1.1, col=cols[1], pos=3)   ##, offset=1.3)
      text(c, cors$cor2[c], round0(cors$cor2[c], digits=2), cex=1.1, col=cols[2], pos=3)
   }
   axis(side=1, at=seq(2, 22, by=2), cex.axis=1)
   axis(side=2, at=seq(-0.8, 0.8, by=0.4), labels=c(-0.8, -0.4, 0, 0.4, 0.8), cex.axis=1)
   dev.off()
}

plotRTvsRT <- function(timings1, timings2, file.name, main.text, xlab.text, ylab.text, method) {
   png(paste0(file.name, ".png"), height=5, width=5, units="in", res=300)
   #pdf(paste0(file.name, ".pdf"), height=5, width=5)
   plot(x=timings1, y=timings2, xlab=xlab.text, ylab=ylab.text, main=main.text[1], col=adjustcolor("black", alpha.f=0.20))
   abline(v=0, lty=5)
   abline(h=0, lty=5)
   
   lm.fit <- lm(timings1 ~ timings2)
   abline(lm.fit, col="black", lwd=3)
   #cor <- cor.test(timings1, timings2, method="spearman", exact=F)

   #legend("bottomright", c(paste0("rho = ", round0(cor[[4]], digits=2)), paste0("p-value = ", scientific(cor[[3]]))), text.col="black", bty="n", cex=1.1)
   mtext(main.text[2], cex=1.2, line=0.3)
   dev.off()
}

# -----------------------------------------------------------------------------
# Compare betweeen RT and LCL RT in sclc-wgs-rt.R
# Last Modified: 05/03/19
# -----------------------------------------------------------------------------
plotRTvsRTALL <- function(cors, file.name, main.text, ylab.text, xlab.text, ymin, ymax, col, c=NA, pos) {
   #png(paste0(file.name, ".png"), height=5, width=5, units="in", res=300)
   pdf(paste0(file.name, ".pdf"), height=5, width=5)
   plot(cors$cor ~ cors$chr, ylim=c(ymin, ymax), ylab=ylab.text, xlab=xlab.text, main=main.text, col=col, xaxt="n", yaxt="n", pch=19)
   lines(cors$cor, y=NULL, type="l", lwd=3, col=col)
   abline(h=0, lty=5)
 
   if (!is.na(c))
      text(cors[c,]$chr, cors[c,]$cor, round0(cors[c,]$cor, digits=2), cex=1.1, col=col, pos=pos)
   axis(side=1, at=seq(2, 22, by=2))
   axis(side=2, at=seq(-1, 1, by=0.2), labels=c(-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1))
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
   boxplot(cor ~ chr, data=cors.samples.plot, ylim=c(ymin, ymax), ylab="Spearman's rho", xlab="Chromosome", outline=T, xaxt="n", main=main.text[1], cex.lab=1.5, cex.axis=1.3, cex.main=1.6)#, medcol="red")
   axis(side=1, at=seq(2, 22, by=2), cex.axis=1.2)
   abline(h=0, lty=5)
   mtext(main.text[2], cex=1.1, line=0.3) 
   dev.off()
}

# -----------------------------------------------------------------------------
# Mean replication timing ratio (MRTR)
# Last Modified: 06/10/19
# -----------------------------------------------------------------------------
getYlim <- function(means, unit) {
   unit <- max(means)/unit
   ymax <- max(means) + unit
   
   ymin <- ((ymax - means[2]) - means[2]) * -1
   
   return(c(ymin, ymax))
}

plotSPR <- function(cors, file.name, main.text, cs=NULL, digits, unit, ylab.text) {
   xlab.text <- "Chromosome"
   cols <- c("red", "blue", "black")
   #ylim <- getYlim(cors$spr, unit)
   ylim <- c(-1, 1)
    
   pdf(paste0(file.name, ".pdf"), height=5, width=5)
   #plot(cors$skew ~ cors$chr, ylim=c(ymin, ymax), ylab=ylab.text, xlab=xlab.text, main=main.text, col=cols[3], xaxt="n", pch=19)   ## yaxt="n",
   plot(NULL, xlim=c(1, 22), ylim=ylim, xlab=xlab.text, ylab=ylab.text, main=main.text[1], col=cols[3], xaxt="n", pch=19)
   
   abline(h=cors$spr[2], lty=5)
   lines(cors$spr, y=NULL, type="l", lwd=3, col=cols[3])
   
   idx <- which(cors$spr > cors$spr[2])
   points(cors$spr[idx] ~ cors$chr[idx], col=cols[1], pch=19)
   idx <- which(cors$spr < cors$spr[2])
   points(cors$spr[idx] ~ cors$chr[idx], col=cols[2], pch=19)
   points(cors$spr[2] ~ cors$chr[2], col=cols[3], pch=19)
 
   text(cors$chr[2]+1.8, cors$spr[2], paste0("Chr2 (", round0(cors$spr[2], digits=digits), ")"), cex=1.1, col=cols[3], pos=3)
   if (!is.null(cs))
      for (c in 1:length(cs)) {
         c <- cs[c]
         if (cors$spr[c] > cors$spr[2])
            text(cors$chr[c]+1.8, cors$spr[c], paste0("Chr", c, " (", round0(cors$spr[c], digits=digits), ")"), cex=1.1, col=cols[1], pos=3)
         else
            text(cors$chr[c]+1.8, cors$spr[c], paste0("Chr", c, " (", round0(cors$spr[c], digits=digits), ")"), cex=1.1, col=cols[2], pos=1)
      }
   legend("topleft", "Earlier than chr2", text.col=cols[1], pch=16, col=cols[1])   
   legend("bottomleft", "Later than chr2", text.col=cols[2], pch=16, col=cols[2])

   axis(side=1, at=seq(2, 22, by=2))
   mtext(main.text[2], cex=1.2, line=0.3)
   dev.off()
}

plotSPRRDC <- function(means, cors, file.name, main.text, cs, xlab.text, unit, ylab.text) {
   cols <- c("red", "blue", "black", "purple")
   #ylim <- getYlim(cors$spr, unit)
   ylim <- c(-1, 1)
   
   unit <- (max(cors) - min(cors))/20
   xlim <- c(min(cors) - unit, max(cors) + unit)
   
   pdf(paste0(file.name, ".pdf"), height=5, width=5)
   if (is.null(xlim))
      plot(means ~ cors, ylim=ylim, ylab=ylab.text, xlab=xlab.text, main=main.text[1], col="black", pch=19)
   else
      plot(means ~ cors, ylim=ylim, xlim=xlim, ylab=ylab.text, xlab=xlab.text, main=main.text[1], col="black", pch=19)

   abline(h=means[2], lty=5)
   lm.fit <- lm(means ~ cors)
   abline(lm.fit, col=cols[4], lwd=3)

   idx <- which(means > means[2])
   points(means[idx] ~ cors[idx], col=cols[1], pch=19)
   idx <- which(means < means[2])
   points(means[idx] ~ cors[idx], col=cols[2], pch=19)
   points(means[2] ~ cors[2], col=cols[3], pch=19)
 
   text(cors[2], means[2], paste0("Chr", 2), cex=1.1, col=cols[3], pos=3)
   if (!is.null(cs))
      for (c in 1:length(cs)) {
         c <- cs[c]
         if (means[c] > means[2])
            text(cors[c], means[c], paste0("Chr", c), cex=1.1, col=cols[1], pos=3)
         else
            text(cors[c], means[c], paste0("Chr", c), cex=1.1, col=cols[2], pos=1)
      }
   
   cor <- cor.test(means, cors, method="spearman", exact=F)
   legends <- c("topright", "bottomleft")
   if (cor[[4]] > 0) legends[1] <- "topleft"
   legend(legends[1], "Earlier than chr2", text.col=cols[1], pch=16, col=cols[1])   ## bty="n"
   legend(legends[2], "Later than chr2", text.col=cols[2], pch=16, col=cols[2])
   
   legend("bottomright", c(paste0("rho = ", round0(cor[[4]], digits=2)), paste0("p-value = ", scientific(cor[[3]]))), text.col=cols[4], bty="n", cex=1.1)
   #legend("bottomright", c(paste0("R^2 = ", round0(summary(lm.fit)$r.squared, digits=2)), paste0("p-value = ", scientific(summary(lm.fit)$coefficients[2, 4]))), text.col=cols[4], bty="n", cex=1.1)
   #axis(side=1, at=seq(2, 22, by=2))
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
