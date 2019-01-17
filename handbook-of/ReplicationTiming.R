# =============================================================================
# Library      : DNA Replication Timing
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
   rt.chr <- toTable(NA, length(colnames), nrow(bed.gc.chr), colnames)
   rownames(rt.chr) <- rownames(bed.gc.chr)
 
   return(rt.chr)
}

# -----------------------------------------------------------------------------
# Methods: Bootstraps (in 2a_cmd-rt_rpkm.corr.gc.d_bstrp.R)
# Last Modified: 15/05/17
# -----------------------------------------------------------------------------
getDetectedRD <- function(rpkm) {   ## Find not dected (RPKM_CORR_GC = 0) windows in any of the samples 
   return(mapply(x = 1:nrow(rpkm), function(x) !any(as.numeric(rpkm[x,]) == 0)))
}

## Read in rpkm.corr.gc.txt.gz and getDetectedRD()
pipeGetDetectedRD <- function(wd.ngs.data, BASE, chr, PAIR, samples) {
   rpkms.chr <- readTable(file.path(wd.ngs.data, paste0(tolower(BASE), "_rpkm.corr.gc_", chr, "_", PAIR, ".txt.gz")), header=T, rownames=T, sep="")[, samples]   ## ADD samples 14/06/17
   rpkms.chr.d <- rpkms.chr[getDetectedRD(rpkms.chr),]   ## ADD getDetected() 13/06/17
 
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
   load(file.path(wd1.rt.data, b, paste0(type, "_ensGene.rt_bstrp.RData")))
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
# Determine replication timing from bootstrapped data (in 2c and 2d)
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
# Visualisation of bootstrapping data (Histogram, RO, and RT)
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
   ylab.text <- "Replication time (log2 FC)"
   if (!is.na(xmin) && !is.na(xmax)) file.name <- paste0(file.name, "_", xmin/1E6, "-", xmax/1E6, "Mb")
   if (is.na(xmin)) xmin <- 0
   if (is.na(xmax)) xmax <- subset(chromInfo, chrom == chr)$size
   ymin <- -ymax
 
   ## Initiation plot
   if (ext == "pdf") {
      pdf(paste0(file.name, ".pdf"), height=4, width=10)
   } else if (ext == "png")
      png(paste0(file.name, ".png"), height=4, width=10, units="in", res=300)   ## ADD 16/05/17: res=300
   plot(NULL, ylim=c(ymin, ymax), xlim=c(xmin/1E6, xmax/1E6), xlab=xlab.text, ylab=ylab.text, main=main.text)
   points(bed.gc.chr$START/1E6, rt.chr$RT, col="grey", cex=0.3)
   abline(h=0, lwd=0.5, col="grey")
 
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

# -----------------------------------------------------------------------------
# Visualisation of bootstrapping data (Expression level and gene length)
# Last Modified: 23/12/18
# -----------------------------------------------------------------------------
plotRD <- function(wd.rt.plots, sample, file, chr, PAIR, rpkm.chr, bed.gc.chr, isFlipped, ext) {
   main.text <- paste0(sample, " read depth (copy number-, GC-corrected RPKM)")
   xlab.text <- paste0("Chromosome ", gsub("chr", "", chr), " coordinate (Mb)")
   ylab.text <- "Read depth"
 
   xmin <- 0
   xmax <- subset(chromInfo, chrom == chr)$size
   ymin <- min(rpkm.chr)
   ymax <- max(rpkm.chr)
 
   if (ext == "pdf")
      pdf(paste0(wd.rt.plots, "/", chr, "/", sample, file, chr, "_", PAIR, ".pdf"), height=4, width=10)
   else if (ext == "png")
      png(paste0(wd.rt.plots, "/", chr, "/", sample, file, chr, "_", PAIR, ".png"), height=4, width=10, units="in", res=300)   ## ADD 16/05/17: res=300
 
   plot(NULL, ylim=c(ymin, ymax), xlim=c(xmin/1E6, xmax/1E6), xlab=xlab.text, ylab=ylab.text, main=main.text)
   points(bed.gc.chr$START/1E6, rpkm.chr, col="red", cex=0.3)
   lines(bed.gc.chr$START/1E6, smooth.spline(rpkm.chr)$y)
 
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
