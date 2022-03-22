# =============================================================================
# Library      : Replication Timing
# Name         : handbook-of/ReplicationTiming.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 29/11/20; 01/08/19; 11/07/19; 05/03/19; 15/11/18
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

getEnsGenesFromSegment <- function(seg) {
   ensGene.chr <- subset(ensGene, chromosome_name == seg$CHR)
   ensGene.chr.end <- ensGene.chr[which(ensGene.chr$end_position >= seg$START),]
   ensGene.chr.end.start <- ensGene.chr.end[which(ensGene.chr.end$start_position <= seg$END),]
 
   return(ensGene.chr.end.start)
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
   return(mapply(x = 1:nrow(rds), function(x) !any(as.numeric(rds[x,]) == 0)))   ## REMOVED rds[x, -1] 30/11/20; ADD 15/02/19; To skip the first column "BED"
}

## Read in rpkm.corr.gc.txt.gz and getDetectedRD()
pipeGetDetectedRD <- function(wd.ngs.data, BASE, chr, PAIR, method, samples) {
   nrds.chr <- readTable(file.path(wd.ngs.data, paste0(tolower(BASE), "_", method, ".gc.cn_", chr, "_", PAIR, ".txt.gz")), header=T, rownames=T, sep="")##[, samples]   ## REMOVED 15/02/19; if length(samples) == 1
   colnames(nrds.chr) <- gsub("\\.", "-", colnames(nrds.chr))    ## ADD 29/11/20; 05/10/19
   nrds.chr.d <- nrds.chr[getDetectedRD(nrds.chr[, samples]),]   ## ADD 13/06/17; getDetectedRD()
 
   return(nrds.chr.d)
}

## Read in rpkm.corr.gc.txt.gz and getDetectedRD()
#pipeGetDetectedRD <- function(wd.rt.data, BASE, chr, PAIR, n, method, samples) {
#   nrds.chr <- readTable(file.path(wd.rt.data, paste0(tolower(BASE), "_", method, ".gc.cn.m_", chr, "_", PAIR, "_n", n, ".txt.gz")), header=T, rownames=T, sep="")##[, samples]   ## REMOVED 15/02/19; if length(samples) == 1
#   colnames(nrds.chr) <- gsub("\\.", "-", colnames(nrds.chr))    ## ADD 29/11/20; 05/10/19
#   nrds.chr.d <- nrds.chr[getDetectedRD(nrds.chr[, samples]),]   ## ADD 13/06/17; getDetectedRD()
# 
#   return(nrds.chr.d)
#}

pipeMedianRD <- function(wd.ngs.data, BASE, chr, PAIR, method, samples=NULL) {
   nrds.chr <- readTable(file.path(wd.ngs.data, paste0(tolower(BASE), "_", method, ".gc.cn_", chr, "_", PAIR, ".txt.gz")), header=T, rownames=T, sep="")##[, samples]   ## REMOVED 15/02/19; if length(samples) == 1
   colnames(nrds.chr) <- gsub("\\.", "-", colnames(nrds.chr))   ## ADD 29/11/20; 05/10/19
   if (!is.null(samples))
      nrds.chr <- nrds.chr[,c("BED", samples)]
   
   nrds.chr$MEDIAN <- mapply(x = 1:nrow(nrds.chr), function(x) median(as.numeric(nrds.chr[x, -1])))   ## ADD 15/02/19; To skip the first column "BED"
   return(nrds.chr)
}

outputRT <- function(nrds.chr) {
   samples <- colnames(nrds.chr)
   nrds.chr$BED <- rownames(nrds.chr)
 
   return(nrds.chr[,c("BED", samples)])
}

# -----------------------------------------------------------------------------
# Get overall read depth
# Last Modified: 30/11/20
# -----------------------------------------------------------------------------
getOverallRD <- function(wd.rt.data, base, method, PAIR1, n1) {
   nrds.m <- readTable(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.m_", "chr1", "_", PAIR1, "_n", n1, ".txt.gz")), header=T, rownames=F, sep="\t")
   for (c in 2:22) {
      chr <- chrs[c]
      nrds.chr.m <- readTable(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.m_", chr, "_", PAIR1, "_n", n1, ".txt.gz")), header=T, rownames=F, sep="\t")
  
      nrds.m <- rbind(nrds.m, nrds.chr.m)
   }
 
   #colnames(nrds.m) <- gsub("\\.", "-", colnames(nrds.m))   ## ADD 29/11/20; 05/10/19
   #return(nrds.m[, -1])   ## BUG 2019/10/14: Remove column BED
   return(nrds.m)
}

getOverallDetectedRD <- function(wd.rt.data, base, method, PAIR1) {
   nrds.d <- pipeGetDetectedRD(wd.rt.data, base, "chr1", PAIR1, method)
   for (c in 2:22) {
      chr <- chrs[c]
      nrds.chr.d <- pipeGetDetectedRD(wd.rt.data, base, chr, PAIR1, method)
  
      nrds.d <- rbind(nrds.d, nrds.chr.d)
   }
 
   return(nrds.d)
}

# -----------------------------------------------------------------------------
# Plot RD and RT in sclc-wgs-rt.R (also refer to plotBootstrapsRT)
# Link(s): http://www.mun.ca/biology/scarr/2250_DNA_replication_&_transcription.html
#          https://stackoverflow.com/questions/43615469/how-to-calculate-the-slope-of-a-smoothed-curve-in-r
# Last Modified: 14/02/19
# -----------------------------------------------------------------------------
#getLog2ScaledDetectedRT <- function(wd.rt.data, base, method, PAIR1, PAIR0, n1, n0, chrs, bed.gc, isFlip=F) {
#   nrds <- toTable(NA, 4, 0, c("BED", "T", "N", "RT"))
#   for (c in 1:22) {
#      chr <- chrs[c]
#      bed.gc.chr <- subset(bed.gc, CHR == chr)
#      nrds.chr <- readTable(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
#  
#      nrds <- rbind(nrds, nrds.chr)
#   }
#   if (isFlip) {
#      nrds$RT <- nrds$N / nrds$T
#   }
#
#   nrds$RT <- log2(nrds$RT)   ## MUY MUY IMPORTANTE!! 2019/10/10
#   nrds$RT <- scale(nrds$RT)
#
#   return(nrds[,c("BED", "T", "N", "RT")])
#}

#getLog2ScaledRT <- function(wd.rt.data, base, method, PAIR1, PAIR0, n1, n0, chrs, bed.gc, isFlip=F) {
#   nrds <- toTable(NA, 4, 0, c("BED", "T", "N", "RATIO"))
#   for (c in 1:22) {
#      chr <- chrs[c]
#      bed.gc.chr <- subset(bed.gc, CHR == chr)
#      nrds.chr <- readTable(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.m.ratio_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
#  
#      nrds <- rbind(nrds, nrds.chr)
#   }
#   if (isFlip) {
#      nrds$RATIO <- nrds$N / nrds$T
#   }
#
#   noreads <- sort(unique(c(which(nrds$T == 0), which(nrds$N == 0), which(is.infinite(nrds$RATIO) == T), which(is.na(nrds$RATIO) == T))))
#   reads   <- setdiff(seq(1:nrow(nrds)), noreads)
#   
#   nrds$RT <- log2(nrds$RATIO)   ## MUY MUY IMPORTANTE!! 2019/10/10
#   nrds$RT[reads] <- scale(nrds$RT[reads])
#   nrds$RT[noreads] <- NA
#   
#   return(nrds[,c("BED", "T", "N", "RT")])
#}

getLog2ScaledRT <- function(wd.rt.data, base, method, PAIR1, PAIR0, n1, n0, chrs, bed.gc, isFlip=F) {
   nrds <- toTable(NA, 4, 0, c("BED", "T", "N", "RT"))
   for (c in 1:22) {
      chr <- chrs[c]
      bed.gc.chr <- subset(bed.gc, CHR == chr)
      nrds.chr <- readTable(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.m.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
  
      nrds <- rbind(nrds, nrds.chr)
   }
   if (isFlip) {
      nrds$RT <- nrds$N / nrds$T
   }
   nrds$RT <- log2(nrds$RT)   ## MUY MUY IMPORTANTE!! 2019/10/10
   nrds$RT <- scale(nrds$RT)
 
   return(nrds)
}

setSpline0 <- function(nrds.chr, bed.gc.chr, column) {
   overlaps <- intersect(rownames(bed.gc.chr), nrds.chr$BED)
   bed.gc.chr.o <- bed.gc.chr[overlaps,]
   nrds.chr.o   <- nrds.chr[overlaps,]
 
   spline <- smooth.spline(nrds.chr.o[, column])
   nrds.chr.o$SPLINE <- spline$y
 
   return(nrds.chr.o[, c("BED", column, "SPLINE")])
}

setSpline <- function(nrds.chr, bed.gc.chr, column) {
   overlaps <- intersect(rownames(bed.gc.chr), nrds.chr$BED)
   bed.gc.chr.o <- bed.gc.chr[overlaps,]
   nrds.chr.o   <- nrds.chr[overlaps,]
 
   spline <- smooth.spline(x=bed.gc.chr.o$START, y=nrds.chr.o[, column])
   if (nrow(nrds.chr.o) == length(spline$y)) {
      nrds.chr.o$SPLINE <- spline$y
   } else {
      bed.gc.chr.o.spline <- subset(bed.gc.chr.o, START %in% spline$x)
      nrds.chr.o$SPLINE <- NA
      nrds.chr.o[rownames(bed.gc.chr.o.spline),]$SPLINE <- spline$y
      
      diff <- setdiff(rownames(bed.gc.chr.o), rownames(bed.gc.chr.o.spline))
      idx <- which(rownames(nrds.chr.o) %in% diff)
      nrds.chr.o[idx,]$SPLINE <- nrds.chr.o[idx-1,]$SPLINE
   }
   
   return(nrds.chr.o[, c("BED", column, "SPLINE")])
}

setSpline2 <- function(nrds.chr, bed.gc.chr, column, kb=1, returnAll=F) {
   nrds.chr <- setSpline(nrds.chr, bed.gc.chr, column)
   overlaps <- intersect(rownames(bed.gc.chr), nrds.chr$BED)
   bed.gc.chr.o <- bed.gc.chr[overlaps,]
   nrds.chr.o   <- nrds.chr[overlaps,]
 
   ## https://stackoverflow.com/questions/43615469/how-to-calculate-the-slope-of-a-smoothed-curve-in-r
   slopes <- diff(nrds.chr.o$SPLINE)/diff(bed.gc.chr.o$START/1E6)   ## ADD 31/10/18
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

## http://genome.ucsc.edu/blog/the-ucsc-genome-browser-coordinate-counting-systems/
plotRT <- function(file.name, main.text, chr, xmin, xmax, nrds.chr, bed.gc.chr, colours, legends, colours2, legends2, ext, width, peaks, ylim=NULL, lcl.rt.chr=NULL, nrds.lcl.chr=NULL, legend="topright") {
   ## Colours (was "lightcoral", "skyblue3")
   ## http://r.789695.n4.nabble.com/plot-function-color-transparency-td4682424.html
   adjustcolor.red  <- adjustcolor(colours2[1], alpha.f=0.01)
   adjustcolor.blue <- adjustcolor(colours2[2], alpha.f=0.01)

   ## Read in replicaiton time
   overlaps <- intersect(rownames(bed.gc.chr), rownames(nrds.chr))   ## 29/11/19: Changed from intersect(rownames(nrds), rownames(bed.gc.chr))
   nrds.chr   <- nrds.chr[overlaps,]
   bed.gc.chr <- bed.gc.chr[overlaps,]
   
   nrds.chr.T  <- setSpline0(nrds.chr, bed.gc.chr, "T")
   nrds.chr.N  <- setSpline0(nrds.chr, bed.gc.chr, "N")
   nrds.chr.RT <- setSpline0(nrds.chr, bed.gc.chr, "RT")
   if (!is.null(nrds.lcl.chr))
      nrds.lcl.chr.RT <- setSpline(nrds.lcl.chr, bed.gc.chr, "RT")
    
   #nrds.chr.RT$SPLINE <- scale(nrds.chr.RT$SPLINE)
   #bed.gc.chr <- bed.gc.chr[rownames(nrds.chr.RT),]   ## NOT HERE?
 
   xlab.text <- paste0("Chromosome ", gsub("chr", "", chr), " position [Mb]")
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
      plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=c(min(rds), max(rds)), xlab="", ylab=ylab.text, main=main.text, xaxt="n", yaxt="n", cex.axis=1.2, cex.lab=1.3, cex.main=1.35)
   } else
      plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=ylim, xlab="", ylab=ylab.text, main=main.text, xaxt="n", yaxt="n", cex.axis=1.2, cex.lab=1.3, cex.main=1.35)
   #points(bed.gc.chr$START/1E6, nrds.chr, col=colours[1], cex=0.3)
   #abline(h=0, lty=5, lwd=1.5, col="lightgrey")
 
   ## Plot cytobands (before smoothing spline)
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.5, col="lightgrey") 

   ## Plot smoothing splines
   bed.gc.chr.lcl <- NULL
   if (!is.null(nrds.lcl.chr))
      bed.gc.chr.lcl <- bed.gc.chr[rownames(nrds.lcl.chr.RT),]
   bed.gc.chr <- bed.gc.chr[rownames(nrds.chr.RT),]   ## BUT HERE?
   points(bed.gc.chr$START/1E6, nrds.chr.N$SPLINE, col=colours[2], pch=16, cex=0.3)   ## G1 first
   points(bed.gc.chr$START/1E6, nrds.chr.T$SPLINE, col=colours[1], pch=16, cex=0.3)   ## S next (on top)!
   axis(side=2, at=seq(0.2, 0.6, by=0.2), labels=c(0.2, 0.4, 0.6), cex.axis=1.2)
   axis(side=2, at=seq(0.3, 0.5, by=0.2), labels=c("", ""), cex.axis=1.2)
   
   ## Plot legend and peaks
   legend(legend, legends, col=colours, lty=1, lwd=3, bty="n", horiz=T, cex=1.3)
   if (length(peaks) != 0)
      for (p in 1:length(peaks))
         abline(v=peaks[p]/1E6, lty=5, lwd=1, col="black")
   
   ### 
   ## Initiate RT plot
   par(mar=c(5.5,4,0,1))
   ylab.text <- "RT [log2]"
   
   plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=c(-2, 2), xlab=xlab.text, ylab=ylab.text, main="", yaxt="n", cex.axis=1.2, cex.lab=1.3, cex.main=1.35)
   idx <- which(nrds.chr.RT$RT == 0)
   points(bed.gc.chr[idx,]$START/1E6, nrds.chr.RT[idx,]$RT, col="lightgrey", cex=0.3)
   idx <- which(nrds.chr.RT$RT < 0)
   points(bed.gc.chr[idx,]$START/1E6, nrds.chr.RT[idx,]$RT, col=adjustcolor.blue, cex=0.3)
   idx <- which(nrds.chr.RT$RT > 0)
   points(bed.gc.chr[idx,]$START/1E6, nrds.chr.RT[idx,]$RT, col=adjustcolor.red, cex=0.3)
   
   ## Plot cytobands (before smoothing spline)
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.5, col="lightgrey") 
  
   ## Plot Koren 2012 (before smoothing spline)
   if (!is.null(lcl.rt.chr))
      points(lcl.rt.chr$POS/1E6, lcl.rt.chr$RT, col=colours[3], pch=16, cex=0.2)
   if (!is.null(nrds.lcl.chr)) {
      points(bed.gc.chr.lcl$START/1E6, nrds.lcl.chr.RT$SPLINE, col=colours[3], pch=16, cex=0.35)
   }
   ## Plot smoothing spline
   points(bed.gc.chr$START/1E6, nrds.chr.RT$SPLINE, col="black", pch=16, cex=0.3)
   abline(h=0, lty=5, lwd=1.5, col="black")
   axis(side=2, at=seq(-2, 2, by=4), labels=c("\u22122", 2), cex.axis=1.2)
   axis(side=2, at=seq(-1, 1, by=1), labels=c("\u22121", 0, 1), cex.axis=1.2)
   axis(side=2, at=0, labels=0, cex.axis=1.2)
   
   ## Plot legend and peaks
   legend("topleft", paste0("Early (", legends2[1], " > ", legends2[2], ")"), bty="n", text.col="black", pt.cex=0.9, pt.lwd=1.25, pch=1, col=colours2[1], cex=1.3)   
   legend("bottomleft", paste0("Late (", legends2[1], " < ", legends2[2], ")"), bty="n", text.col="black", pt.cex=0.9, pt.lwd=1.25, pch=1, col=colours2[2], cex=1.3)
   if (width != 5)
      legend("topright", paste0("", legends2[1], "/", legends2[2], " read depth ratio"), col="black", lty=1, lwd=3, bty="n", horiz=T, cex=1.3)
   else
      legend("topright", paste0(legends2[1], "/", legends2[2], " ratio"), col="black", lty=1, lwd=3, bty="n", horiz=T, cex=1.3)
   if (!is.null(lcl.rt.chr))
      legend("bottomright", "Koren et al. (~2 kb)", col=colours[3], lty=1, lwd=3, bty="n", horiz=T, cex=1.3)
   if (!is.null(nrds.lcl.chr))
      legend("bottomright", "LCL S/G1 ratio", col=colours[3], lty=1, lwd=3, bty="n", horiz=T, cex=1.3)
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
 
   jpeg(paste0(file.name, ".jpg"), height=5, width=5, units="in", res=300)
   plot(reads1 ~ timings, xlim=xlim, ylim=ylim, ylab=ylab.text, xlab=xlab.text, main=main.text[1], xaxt="n", col=colours[2], cex.axis=1.2, cex.lab=1.25, cex.main=1.35)
   #pdf(paste0(file.name, ".pdf"), height=5, width=5)
   #plot(NULL, xlim=xlim, ylim=ylim, ylab=ylab.text, xlab=xlab.text, main=main.text[1], xaxt="n", cex.axis=1.1, cex.lab=1.2, cex.main=1.25)
   abline(v=0, lty=5, lwd=2)
 
   lm.fit1 <- lm(reads1 ~ timings)
   abline(lm.fit1, col=colours[1], lwd=2)
   cor1 <- getCor(reads1, timings, method)

   RT <- "RT"   #paste0(legends[1], "/", legends[2])
   if (cor1 > 0) {
      legend("topright", paste0(legends[1], " vs. ", RT, " (", cor, " = ", round0(cor1, digits=2), ")"), text.col=colours[1], bty="n", cex=1.25)        
   } else {
      legend("bottomright", paste0(legends[2], " vs. ", RT, " (", cor, " = \u2212", round0(abs(cor1), digits=2), ")"), text.col=colours[1], bty="n", cex=1.25)
   }
   axis(side=1, at=seq(-3, 3, by=0.5), labels=c(-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3), cex.axis=1.2)
   mtext(main.text[2], line=0.3, cex=1.2)
   dev.off()
}

plotRD2vsRT <- function(reads1, reads2, timings, file.name, main.text, ylab.text, xlab.text, colours, legends, method) {
   xlim <- c(min(timings), max(timings))
   ylim <- c(min(c(reads1, reads2)), max(c(reads1, reads2)))
   cor <- "rho"
   main.text2 <- "Spearman's rho"
   if (method == "pearson") {
      cor <- "r" 
      main.text2 <- "Pearson correlation"
   }

   #png(paste0(file.name, ".png"), height=5, width=5, units="in", res=300)
   pdf(paste0(file.name, ".pdf"), height=5, width=5)
   plot(NULL, xlim=xlim, ylim=ylim, ylab=ylab.text, xlab=xlab.text, main=main.text[1], xaxt="n", cex.axis=1.25, cex.lab=1.3, cex.main=1.35)
   abline(v=0, lty=5, lwd=1.5)
   
   lm.fit2 <- lm(reads2 ~ timings)   ## G1 first
   abline(lm.fit2, col=colours[2], lwd=3)
   cor2 <- getCor(reads2, timings, method)

   lm.fit1 <- lm(reads1 ~ timings)   ## S next (on top)!
   abline(lm.fit1, col=colours[1], lwd=3)
   cor1 <- getCor(reads1, timings, method)
   
   RT <- paste0(legends[1], "/", legends[2])   # "RT"
   if (cor1 < 0 && cor2 < 0) {
      legend("bottomright", c(paste0(legends[1], " vs. ", RT, " (", cor, " = ", round0(cor1, digits=2), ")"), paste0(legends[2], " vs. ", RT, " (", cor, " = ", round0(cor2, digits=2), ")")), text.col=colours, text.font=2, bty="n", cex=1.3)
   } else if (as.numeric(cor1) > 0 && as.numeric(cor2) < 0) {
      legend("topright", paste0(legends[1], " vs. ", RT, " (", cor, " = ", round0(cor1, digits=2), ")"), text.col=colours[1], text.font=2, bty="n", cex=1.3)        
      legend("bottomright", paste0(legends[2], " vs. ", RT, " (", cor, " = ", round0(cor2, digits=2), ")"), text.col=colours[2], text.font=2, bty="n", cex=1.3)
   }
   axis(side=1, at=seq(-3, 3, by=0.5), labels=c(-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3), cex.axis=1.25)
   mtext(main.text[2], line=0.3, cex=1.3)
   dev.off()
}

plotRD2vsRTALL <- function(cors, file.name, main.text, ymin, ymax, cols, cols2, legends, c=NA) {
   ylab.text <- "Correlation [rho]"
   xlab.text <- "Chromosome"
 
   #png(paste0(file.name, ".png"), height=5, width=5, units="in", res=300)
   pdf(paste0(file.name, ".pdf"), height=6, width=6)
   plot(cors$cor1 ~ cors$chr, ylim=c(ymin, ymax), ylab="", xlab=xlab.text, main=main.text, col=cols[1], xaxt="n", yaxt="n", pch=19, cex=1.5, cex.lab=1.8, cex.main=2)
   lines(cors$cor1, y=NULL, type="l", lwd=3, col=cols[1])
   points(cors$chr, cors$cor2, col=cols[2], pch=19, cex=1.5)
   lines(cors$cor2, y=NULL, type="l", lwd=3, col=cols[2])
   abline(h=0, lty=5, lwd=2)
   
   RT <- "RT"   #paste0(legends[1], "/", legends[2])
   legend("topright", paste0(legends[1], " vs. ", RT), text.col=cols[1], bty="n", cex=1.9, text.font=2)        
   legend("bottomright", paste0(legends[2], " vs. ", RT), text.col=cols[2], bty="n", cex=1.9, text.font=2)
   
   if (!is.na(c)) {
      text(c, cors$cor1[c], round0(cors$cor1[c], digits=2), col=cols2[1], pos=3, cex=1.8)   ##, offset=1.3)
      text(c, cors$cor2[c], round0(cors$cor2[c], digits=2), col=cols2[2], pos=3, cex=1.8)
   }
   axis(side=1, at=seq(2, 22, by=4), cex.axis=1.7)
   axis(side=1, at=seq(4, 20, by=4), cex.axis=1.7)
   axis(side=2, at=seq(-0.8, 0.8, by=0.4), labels=c(-0.8, -0.4, 0, 0.4, 0.8), cex.axis=1.75)
   mtext(ylab.text, side=2, line=2.75, cex=1.8)
   dev.off()
}

plotRD3vsRTALL <- function(cors, file.name, main.text, ymin, ymax, cols, legends, c=NA, isRT=F) {
   ylab.text <- "Correlation [rho]"
   xlab.text <- "Chromosome"
 
   #png(paste0(file.name, ".png"), height=5, width=5, units="in", res=300)
   pdf(paste0(file.name, ".pdf"), height=5, width=5)
   plot(cors$cor ~ cors$chr, ylim=c(ymin, ymax), ylab=ylab.text, xlab=xlab.text, main=main.text, col="black", xaxt="n", yaxt="n", pch=19, cex.axis=1.2, cex.lab=1.25, cex.main=1.35)
   lines(cors$cor, y=NULL, type="l", lwd=2, col="black")
   points(cors$chr, cors$cor2, col=cols[2], pch=19)
   lines(cors$cor2, y=NULL, type="l", lwd=2, col=cols[2])
   points(cors$chr, cors$cor1, col=cols[1], pch=19)
   lines(cors$cor1, y=NULL, type="l", lwd=2, col=cols[1])
   
   abline(h=0, lty=5, lwd=1.5)
 
   RT <- "S/G1"
   #legend("topright", c(paste0(legends[3], " vs. ", RT), paste0("      ", legends[1], " vs. ", RT)), text.col=c(cols[3], cols[1]), bty="n", cex=1.25, text.font=c(1,2)) 
   #legend("bottomright", paste0(legends[2], " vs. ", RT), text.col=cols[2], bty="n", cex=1.25)
   legend("topright", paste0(legends[3], " vs. ", RT), text.col=cols[3], bty="n", cex=1.25) 
   legend("bottomright", c(paste0(legends[1], " vs. ", RT), paste0(legends[2], " vs. ", RT)), text.col=c(cols[1], cols[2]), bty="n", cex=1.25, text.font=c(2,1))        
   ##legend("bottomright", paste0(legends[2], " vs. ", RT), text.col=cols[2], bty="n", cex=1.1)
 
   if (!is.na(c)) {
      text(c, cors$cor1[c], round0(cors$cor1[c], digits=2), col=cols[1], pos=3, cex=1.2)   ##, offset=1.3)
      text(c, cors$cor2[c], round0(cors$cor2[c], digits=2), col=cols[2], pos=3, cex=1.2)
   }
   axis(side=1, at=seq(2, 22, by=4), cex.axis=1.2)
   axis(side=1, at=seq(4, 20, by=4), cex.axis=1.2)
   axis(side=2, at=seq(-0.8, 0.8, by=1.6), labels=c(-0.8, 0.8), cex.axis=1.2)
   axis(side=2, at=seq(-0.4, 0.4, by=0.4), labels=c(-0.4, 0, 0.4), cex.axis=1.2)
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

#plotRTvsRTALL <- function(cors, file.name, main.text, ylab.text, xlab.text, ymin, ymax, col, c=NA, pos) {
#   #png(paste0(file.name, ".png"), height=5, width=5, units="in", res=300)
#   pdf(paste0(file.name, ".pdf"), height=5, width=5)
#   plot(cors$cor ~ cors$chr, ylim=c(ymin, ymax), ylab=ylab.text, xlab=xlab.text, main=main.text, col=col, xaxt="n", yaxt="n", pch=19, cex.lab=1.2)
#   lines(cors$cor, y=NULL, type="l", lwd=3, col=col)
#   abline(h=0, lty=5)
# 
#   if (!is.na(c))
#      text(cors[c,]$chr, cors[c,]$cor, round0(cors[c,]$cor, digits=2), col=col, pos=pos, cex=1.2)
#   axis(side=1, at=seq(2, 22, by=2))
#   axis(side=2, at=seq(-1, 1, by=0.2), labels=c(-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1))
#   dev.off()
#}

# -----------------------------------------------------------------------------
# Compare betweeen M2/M1, Q4/Q4, LCL and SCLC-NL RTs in *-wgs-rt.R
# Last Modified: 05/02/20
# -----------------------------------------------------------------------------
plotRTvsRT2 <- function(cors, file.name, main.text, ymin, ymax, cols, legends) {
   ylab.text <- "Spearman's rho"
   xlab.text <- "Chromosome"
 
   #png(paste0(file.name, ".png"), height=5, width=5, units="in", res=300)
   pdf(paste0(file.name, ".pdf"), height=6, width=6)
   plot(cors$cor ~ cors$chr, ylim=c(ymin, ymax), ylab="", xlab=xlab.text, main=main.text, col=cols[1], xaxt="n", yaxt="n", pch=19, cex=1.5, cex.axis=1.8, cex.lab=1.9, cex.main=2)
   lines(cors$cor, y=NULL, type="l", lwd=3, col=cols[1])
   points(cors$chr, cors$cor1, col=cols[2], pch=19, cex=1.5)
   lines(cors$cor1, y=NULL, type="l", lwd=3, col=cols[2])

   legend("bottomright", c(legends[1], legends[2]), text.col=c(cols[1], cols[2]), bty="n", cex=1.9, text.font=c(2,2))       
   #legend("bottomright", c("", legends[2]), text.col=c(cols[1], cols[2]), bty="n", cex=1.9)
   
   axis(side=1, at=seq(2, 22, by=4), cex.axis=1.8)
   axis(side=1, at=seq(4, 20, by=4), cex.axis=1.8)

   mtext(ylab.text, side=2, line=2.75, cex=1.8)
   axis(side=2, at=seq(0.8, 1, by=0.05), labels=c(0.8, 0.85, 0.9, 0.95, 1), cex.axis=1.8)
   #axis(side=2, at=seq(0.2, 1, by=0.2), labels=c(0.2, 0.4, 0.6, 0.8, 1), cex.axis=1.7)   ## NBL-CL 
   #axis(side=2, at=seq(0.6, 0.9, by=0.1), labels=c(0.6, 0.7, 0.8, 0.9), cex.axis=1.8)   ## LUAD 0.6 ~ 0.9
   dev.off()
}

plotRTvsRT3 <- function(cors, file.name, main.text, ymin, ymax, cols, legends) {
   ylab.text <- "Spearman's rho"
   xlab.text <- "Chromosome"
 
   #png(paste0(file.name, ".png"), height=5, width=5, units="in", res=300)
   pdf(paste0(file.name, ".pdf"), height=6, width=6)
   plot(cors$cor ~ cors$chr, ylim=c(ymin, ymax), ylab="", xlab=xlab.text, main=main.text, col=cols[1], xaxt="n", yaxt="n", pch=19, cex=1.5, cex.axis=1.7, cex.lab=1.8, cex.main=2)
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

###
##
plotRD2 <- function(cors, file.name, main.text, ymin, ymax) {
   ylab.text <- "LCL S vs. RT"
   xlab.text <- "LCL G1 vs. RT"
 
   #png(paste0(file.name, ".png"), height=5, width=5, units="in", res=300)
   pdf(paste0(file.name, ".pdf"), height=5, width=5)
   plot(cors$cor1 ~ abs(cors$cor2), ylim=c(ymin, ymax), xlim=c(ymin, ymax), ylab="", xlab=xlab.text, yaxt="n", xaxt="n", main=main.text[1], col=NULL, pch=19, cex.axis=1.25, cex.lab=1.3, cex.main=1.35)
   lines(c(0, 1.1), c(0, 1.1), pch=NULL, col="black", type="l", lty=5, lwd=1.5)
   #lines(abs(cors$cor2), abs(cors$cor2), pch=NULL, col="purple", type="l", lty=2, lwd=3)
   #lm.fit <- lm(abs(cors$cor2) ~ abs(cors$cor2))
   #abline(lm.fit, col="purple", lty=2, lwd=2)
 
   for (c in 1:22)
      if (cors$cor1[c] > abs(cors$cor2[c])) {
         text(abs(cors$cor2[c]), cors$cor1[c], paste0("Chr", c), col="red", pos=3, cex=1.3)
         points(abs(cors$cor2[c]), cors$cor1[c], col=red, pch=19, cex=1.3)
      } else
         points(abs(cors$cor2[c]), cors$cor1[c], col=blue, pch=19, cex=1.3)
   text(abs(cors$cor2[2]), cors$cor1[2], paste0("Chr", 2), col="blue", pos=3, cex=1.3)
 
   axis(side=2, at=seq(0.1, 0.7, by=0.1), labels=c(0.1, "", 0.3, "", 0.5, "", 0.7), cex.axis=1.25)
   axis(side=1, at=seq(0.1, 0.7, by=0.1), labels=c(-0.1, "", -0.3, "", -0.5, "", -0.7), cex.axis=1.25)
   mtext(ylab.text, side=2, line=2.85, cex=1.3)
   mtext(main.text[2], line=0.3, cex=1.3)
   dev.off()
}

plotRD2 <- function(cors, file.name, main.text, ymin, ymax) {
   ylab.text <- "NBL Q4 vs. RT"
   xlab.text <- "NBL Q1 vs. RT"
 
   #png(paste0(file.name, ".png"), height=5, width=5, units="in", res=300)
   pdf(paste0(file.name, ".pdf"), height=5, width=5)
   plot(abs(cors$cor1) ~ abs(cors$cor2), ylim=c(ymin, ymax), xlim=c(ymin, ymax), ylab="", xlab=xlab.text, yaxt="n", xaxt="n", main=main.text[1], col=NULL, pch=19, cex.axis=1.25, cex.lab=1.3, cex.main=1.35)
   lines(c(0, 1.1), c(0, 1.1), pch=NULL, col="black", type="l", lty=5, lwd=1.5)
   #lines(abs(cors$cor2), abs(cors$cor2), pch=NULL, col="purple", type="l", lty=2, lwd=3)
   #lm.fit <- lm(abs(cors$cor2) ~ abs(cors$cor2))
   #abline(lm.fit, col="purple", lty=2, lwd=2)
 
   for (c in 1:22)
      if (abs(cors$cor1[c]) > abs(cors$cor2[c])) {
         text(abs(cors$cor2[c]), abs(cors$cor1[c]), paste0("Chr", c), col="red", pos=3, cex=1.3)
         points(abs(cors$cor2[c]), abs(cors$cor1[c]), col=red, pch=19, cex=1.3)
      } else
   points(abs(cors$cor2[c]), abs(cors$cor1[c]), col=blue, pch=19, cex=1.3)
   text(abs(cors$cor2[2]), abs(cors$cor1[2]), paste0("Chr", 2), col="blue", pos=3, cex=1.3)
 
   axis(side=2, at=seq(0.4, 1, by=0.1), labels=c(0.4, "", 0.6, "", 0.8, "", 1), cex.axis=1.25)
   axis(side=1, at=seq(0.4, 1, by=0.1), labels=c(-0.4, "", -0.6, "", -0.8, "", -1), cex.axis=1.25)   
   mtext(ylab.text, side=2, line=2.85, cex=1.3)
   mtext(main.text[2], line=0.3, cex=1.3)
   dev.off()
}

plotRD2 <- function(cors, file.name, main.text, ymin, ymax) {
   ylab.text <- "NBL M2 vs. RT"
   xlab.text <- "NBL M1 vs. RT"
 
   #png(paste0(file.name, ".png"), height=5, width=5, units="in", res=300)
   pdf(paste0(file.name, ".pdf"), height=5, width=5)
   plot(abs(cors$cor1) ~ abs(cors$cor2), ylim=c(ymin, ymax), xlim=c(ymin, ymax), ylab="", xlab=xlab.text, yaxt="n", xaxt="n", main=main.text[1], col=NULL, pch=19, cex.axis=1.25, cex.lab=1.3, cex.main=1.35)
   lines(c(0, 1.1), c(0, 1.1), pch=NULL, col="black", type="l", lty=5, lwd=1.5)
   #lines(abs(cors$cor2), abs(cors$cor2), pch=NULL, col="purple", type="l", lty=2, lwd=3)
   #lm.fit <- lm(abs(cors$cor2) ~ abs(cors$cor2))
   #abline(lm.fit, col="purple", lty=2, lwd=2)
 
   for (c in 1:22)
      if (abs(cors$cor1[c]) > abs(cors$cor2[c])) {
         text(abs(cors$cor2[c]), abs(cors$cor1[c]), paste0("Chr", c), col="red", pos=3, cex=1.3)
         points(abs(cors$cor2[c]), abs(cors$cor1[c]), col=red, pch=19, cex=1.3)
      } else
   points(abs(cors$cor2[c]), abs(cors$cor1[c]), col=blue, pch=19, cex=1.3)
   text(abs(cors$cor2[2]), abs(cors$cor1[2]), paste0("Chr", 2), col="blue", pos=3, cex=1.3)
 
   axis(side=2, at=seq(0.4, 1, by=0.1), labels=c(-0.4, "", -0.6, "", -0.8, "", -1), cex.axis=1.25)
   axis(side=1, at=seq(0.4, 1, by=0.1), labels=c(-0.4, "", -0.6, "", -0.8, "", -1), cex.axis=1.25)   
   mtext(ylab.text, side=2, line=2.85, cex=1.3)
   mtext(main.text[2], line=0.3, cex=1.3)
   dev.off()
}

# -----------------------------------------------------------------------------
# Comparison betweeen tumour sample's RD and LCL RT in sclc-wgs-rt.R
# Last Modified: 30/11/20; 03/06/19
# -----------------------------------------------------------------------------
getSAMPLEvsRT <- function(wd.rt.data, samples1) {
   cors.samples <- toTable(0, length(samples1)+4, 22, c("chr", "mean", "var", "cv2", samples1))
   cors.samples$chr <- 1:22
   for (s in 1:length(samples1)) {
      sample <- samples1[s]
      load(file.path(wd.rt.data, "samples", paste0("rd-vs-rt_", sample, "-vs-lcl_spline_spearman.RData")))
  
      cors.samples[, sample] <- cors$cor
   }
 
   for (c in 1:22) {
      cors.samples$mean[c] <- mean(as.numeric(cors.samples[c, samples1]))
      cors.samples$var[c]  <- var(as.numeric(cors.samples[c, samples1]))
      cors.samples$cv2[c]  <- cors.samples$var[c]/cors.samples$mean[c]^2
   }
   
   return(cors.samples)
}

plotSAMPLEvsRT <- function(cors.samples, samples, file.name, main.text=NA, ymin=NA, ymax=NA) {
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
   boxplot(cor ~ chr, data=cors.samples.plot, ylim=c(ymin, ymax), ylab="", xlab="Chromosome", outline=T, xaxt="n", yaxt="n", main=main.text[1], cex.lab=1.8, cex.main=2)#, medcol="red")
   
   axis(side=1, at=seq(2, 22, by=4), cex.axis=1.7)
   axis(side=1, at=seq(4, 20, by=4), cex.axis=1.7)
   axis(side=2, at=seq(-0.8, 0.8, by=0.4), labels=c(-0.8, -0.4, 0, 0.4, 0.8), cex.axis=1.7)
   abline(h=0, lty=5, lwd=2)
   mtext("Correlation [rho]", side=2, line=2.75, cex=1.8)
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

plotRTS <- function(sprs, file.name, main.text, cs=NULL, digits, unit, ylab.text, cex) {
   xlab.text <- "Chromosome"
   cols <- c(red, blue, "black")
   #ylim <- getYlim(sprs$spr, unit)
   ylim <- c(-1.1, 1.1)
    
   pdf(paste0(file.name, ".pdf"), height=5, width=5)
   #plot(sprs$skew ~ sprs$chr, ylim=c(ymin, ymax), ylab=ylab.text, xlab=xlab.text, main=main.text, col=cols[3], xaxt="n", pch=19)   ## yaxt="n",
   plot(NULL, xlim=c(1, 22), ylim=ylim, xlab=xlab.text, ylab=ylab.text, main=main.text[1], col=cols[3], xaxt="n", yaxt="n", pch=19, cex.axis=1.2, cex.lab=1.3, cex.main=1.35)
   
   abline(h=sprs$spr[2], lty=5, lwd=1.3)
   lines(sprs$spr, y=NULL, type="l", lwd=2, col=cols[3])
   
   idx <- which(sprs$spr > sprs$spr[2])
   points(sprs$spr[idx] ~ sprs$chr[idx], col=cols[1], pch=19, cex=cex)
   idx <- which(sprs$spr < sprs$spr[2])
   points(sprs$spr[idx] ~ sprs$chr[idx], col=cols[2], pch=19, cex=cex)
   points(sprs$spr[2] ~ sprs$chr[2], col=cols[3], pch=19, cex=cex)
 
   text(sprs$chr[2], sprs$spr[2], "Chr2", col=cols[3], pos=3, cex=1.3)
   #text(sprs$chr[2], sprs$spr[2], paste0(offset, "Chr2 (", chr2, ")"), col=cols[3], pos=3, cex=1.3)
   #text(sprs$chr[2]+1.8, sprs$spr[2], paste0("Chr2 (", round0(sprs$spr[2], digits=digits), ")"), cex=1.1, col=cols[3], pos=3)
   if (!is.null(cs))
      for (c in 1:length(cs)) {
         c <- cs[c]
         if (sprs$spr[c] > sprs$spr[2])
            text(sprs$chr[c], sprs$spr[c], paste0("Chr", c), col="red", pos=3, cex=1.3)
            #text(sprs$chr[c]+1.8, sprs$spr[c], paste0("Chr", c, " (", round0(sprs$spr[c], digits=digits), ")"), cex=1.1, col=cols[1], pos=3)
         else
            text(sprs$chr[c], sprs$spr[c], paste0("Chr", c), col="blue", pos=1, cex=1.3)
            #text(sprs$chr[c]+1.8, sprs$spr[c], paste0("Chr", c, " (", round0(sprs$spr[c], digits=digits), ")"), cex=1.1, col=cols[2], pos=1)
      }
   legend("topleft", "Earlier than chr2", text.col=cols[1], pch=19, pt.cex=1.5, col=cols[1], cex=1.3)   
   legend("bottomleft", "Later than chr2", text.col=cols[2], pch=19, pt.cex=1.5, col=cols[2], cex=1.3)

   axis(side=1, at=seq(2, 22, by=4), cex.axis=1.2)
   axis(side=1, at=seq(4, 20, by=4), cex.axis=1.2)
   axis(side=2, at=seq(-1, 1, by=0.5), labels=c(-1, -0.5, 0, 0.5, 1), cex.axis=1.2)
   mtext(main.text[2], line=0.3, cex=1.3)
   dev.off()
}

plotRTS2 <- function(sprs, means, file.name, main.text, cs, xlab.text, unit, ylab.text, cex) {
   cols <- c(red, blue, "black", green)
   ylim <- c(-1.1, 1.1)
   
   unit <- (max(means) - min(means))/15
   xlim <- c(min(means) - unit, max(means) + unit)
   
   pdf(paste0(file.name, ".pdf"), height=5, width=5)
   if (is.null(xlim))
      plot(sprs ~ means, ylim=ylim, ylab=ylab.text, xlab=xlab.text, main=main.text[1], col="white", yaxt="n", xaxt="n", pch=19, cex.axis=1.2, cex.lab=1.3, cex.main=1.35)
   else
      plot(sprs ~ means, ylim=ylim, xlim=xlim, ylab=ylab.text, xlab=xlab.text, main=main.text[1], col="white", yaxt="n", xaxt="n", pch=19, cex.axis=1.2, cex.lab=1.3, cex.main=1.35)
   
   abline(h=sprs[2], lty=5, lwd=1.3)

   idx <- which(sprs > sprs[2])
   points(sprs[idx] ~ means[idx], col=cols[1], pch=19, cex=cex)
   idx <- which(sprs < sprs[2])
   points(sprs[idx] ~ means[idx], col=cols[2], pch=19, cex=cex)
   points(sprs[2] ~ means[2], col=cols[3], pch=19, cex=cex)
 
   lm.fit <- lm(sprs ~ means)
   abline(lm.fit, col=cols[4], lwd=3)
   
   text(means[2], sprs[2], paste0("Chr", 2), col=cols[3], pos=3, cex=1.3)
   if (!is.null(cs))
      for (c in 1:length(cs)) {
         c <- cs[c]
         if (sprs[c] > sprs[2])
            #if (sprs[c] == sprs[22])
            #   text(means[c], sprs[c], paste0("Chr", c), col=cols[1], pos=1, cex=1.2)
            #else
               text(means[c], sprs[c], paste0("Chr", c), col="red", pos=3, cex=1.3)
         else
            text(means[c], sprs[c], paste0("Chr", c), col="blue", pos=1, cex=1.3)
      }
   
   cor <- cor.test(sprs, means, method="spearman", exact=F)
   legends <- c("topright", "bottomleft")
   if (cor[[4]] > 0) legends[1] <- "topleft"
   legend(legends[1], "Earlier than chr2", text.col=cols[1], pch=19, pt.cex=1.5, col=cols[1], cex=1.3)   ## bty="n"
   legend(legends[2], "Later than chr2", text.col=cols[2], pch=19, pt.cex=1.5, col=cols[2], cex=1.3)
   
   pvalue <- scientific(cor[[3]])
   #legend("bottomright", paste0("rho = ", round0(cor[[4]], digits=2)), text.col=cols[4], bty="n", cex=1.2)
   legend("bottomright", c(paste0("rho = ", round0(cor[[4]], digits=2)), expression(italic('P')~'= 5.68E-09')), text.col=cols[4], text.font=2, bty="n", cex=1.3)
   #legend("bottomright", c(paste0("R^2 = ", round0(summary(lm.fit)$r.squared, digits=2)), paste0("p-value = ", scientific(summary(lm.fit)$coefficients[2, 4]))), text.col=cols[4], bty="n", cex=1.1)
   #axis(side=1, at=seq(2, 22, by=2))
   axis(side=1, at=seq(-0.4, 0.8, by=0.4), labels=c(-0.4, 0, 0.4, 0.8), cex.axis=1.2)
   axis(side=1, at=seq(-0.2, 0.6, by=0.4), labels=c(-0.2, 0.2, 0.6), cex.axis=1.2)
   axis(side=2, at=seq(-1, 1, by=0.5), labels=c(-1, -0.5, 0, 0.5, 1), cex.axis=1.2)
   mtext(main.text[2], cex=1.3, line=0.3)
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
 
   cors.all$M2[which(cors.all$Q4 %in% c(3, 4))] <- 2
   cors.all$M2[which(cors.all$Q4 %in% c(1, 2))] <- 1

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

read.peiflyne.icgc.cn.txt <- function(cn.file) {
   cn <- readTable(cn.file, header=F, rownames=T, sep="")
   colnames(cn) <- c("BED", "CHR", "START", "END", "V5", "V6", "V7", "RD_T", "RD_N", "V10")
 
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

#setSamplesSG1 <- function(wd.rt.data, samples1, cors.samples) {
#   cors.all <- toTable(0, 4, length(samples1), c("SAMPLE_ID", "COR", "SG1", "M2"))
#   rownames(cors.all) <- samples1
#   cors.all$SAMPLE_ID <- samples1
#   for (s in 1:length(samples1)) {
#      sample <- samples1[s]
#      load(file.path(wd.rt.data, "samples", paste0("rd-vs-rt_", sample, "-vs-lcl_spline_spearman.RData")))
#  
#      cors.all$COR[s] <- cor
#   }
#   #print(quantile(as.numeric(cors.all$COR)))
 
   ##
#   s_likes <- c()
#   g1_likes <- c()
#   for (c in 1:22) {
#      cors <- cors.samples[c, -c(1:4)]
#      median <- median(as.numeric(cors))
      #median <- 0
#      
#      if (c == 1) {
#         s_likes  <- colnames(cors)[which(as.numeric(cors) > median)]
#         g1_likes <- colnames(cors)[which(as.numeric(cors) <= median)]
#      } else {
#         s_likes  <- intersect(s_likes, colnames(cors)[which(as.numeric(cors) > median)])
#         g1_likes <- intersect(g1_likes, colnames(cors)[which(as.numeric(cors) <= median)])
#      }
#   }
#   print(length(s_likes))
#   print(length(g1_likes))
#   print(s_likes)
#   print(g1_likes)
#   
#   cors.all$SG1 <- NA
#   cors.all[s_likes, ]$SG1 <- "SL"
#   cors.all[g1_likes,]$SG1 <- "G1L"
#   
#   cors.all$M2 <- NA
#   cors.all$M2[which(cors.all$SG1 == "SL") ] <- 1
#   cors.all$M2[which(cors.all$SG1 == "G1L")] <- 0
#   return(cors.all)
#}

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

#plotRTvsRT3 <- function(cors, file.name, main.text, ymin, ymax, cols, legends) {
#   ylab.text <- "Spearman's rho"
#   xlab.text <- "Chromosome"
# 
#   #png(paste0(file.name, ".png"), height=5, width=5, units="in", res=300)
#   pdf(paste0(file.name, ".pdf"), height=6, width=6)
#   plot(cors$cor ~ cors$chr, ylim=c(ymin, ymax), ylab="", xlab=xlab.text, main=main.text, col=cols[1], xaxt="n", yaxt="n", pch=19, cex.axis=1.7, cex.lab=1.7, cex.main=1.8)
#   lines(cors$cor, y=NULL, type="l", lwd=3, col=cols[1])
#   points(cors$chr, cors$cor1, col=cols[2], pch=19)
#   lines(cors$cor1, y=NULL, type="l", lwd=3, col=cols[2])
#   points(cors$chr, cors$cor2, col=cols[3], pch=19)
#   lines(cors$cor2, y=NULL, type="l", lwd=3, col=cols[3])
#
#   legend("bottomright", c(legends[1], legends[2], legends[3]), text.col=c(cols[1], cols[2], cols[3]), bty="n", cex=1.75)
#
#   axis(side=1, at=seq(2, 22, by=4), cex.axis=1.6)
#   axis(side=1, at=seq(4, 20, by=4), cex.axis=1.6)
#   mtext(ylab.text, side=2, line=2.85, cex=1.7)
#   axis(side=2, at=seq(0.5, 1, by=0.1), labels=c(0.5, 0.6, 0.7, 0.8, 0.9, 1), cex.axis=1.7)
#   dev.off()
#}
