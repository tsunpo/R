# =============================================================================
# Library      : DNA Replication Timing
# Name         : handbook-of/ReplicationTiming.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 13/06/17
# =============================================================================

# -----------------------------------------------------------------------------
# Methods: Calculate ratio of raw sequencing reads between tumour and matched normal
# Last Modified: 17/05/17
# -----------------------------------------------------------------------------
initRPKM <- function(bam, bed, pair) {
   bam.bed <- bam[intersect(rownames(bam), rownames(bed)),]   ## ADD 23/06/17: Autosomal-only RPKMs
 
   colnames <- c("BED", "BAM", "RPKM")
   rpkm <- toTable(NA, length(colnames), nrow(bam.bed), colnames)
   rpkm[,c("BED", "BAM")] <- bam.bed[,c("BED", paste0("BAM_", pair))]
 
   TOTAL <- sum(as.numeric(rpkm$BAM))   ## BUG FIX 17/05/17: Warning message: In sum(rpkm$BAM) : integer overflow - use sum(as.numeric(.))
   rpkm$RPKM <- rpkm$BAM / TOTAL * 1E6 * 1000
 
   return(rpkm)
}

# -----------------------------------------------------------------------------
# Methods: Calculate copy number-, GC-corrected read counts
# Last Modified: 15/05/17
# -----------------------------------------------------------------------------
## cmd-rt_bam.rpkm.corr.gc.R (command line version)
getBEDFromSegment <- function(bed.gc, seg) {
   #bed.gc.chr <- subset(bed.gc, CHR == seg$CHR)                      ## Using subset() takes ~15 MINS on 1,040 segments over 2,861,558 entries (S00022)
   #bed.gc.chr.end <- subset(bed.gc.chr, END > seg$START)             ## Important! Not >= seg$START)
   #bed.gc.chr.end.start <- subset(bed.gc.chr.end, START <= seg$END)
   bed.gc.chr <- bed.gc[which(bed.gc$CHR == seg$CHR),]                ## Using which() takes ~5 MINS on 1,040 segments over 2,861,558 entries (S00022)
   bed.gc.chr.end <- bed.gc.chr[which(bed.gc.chr$END > seg$START),]   ## Important! Not >= seg$START)
   bed.gc.chr.end.start <- bed.gc.chr.end[which(bed.gc.chr.end$START <= seg$END),]
 
   return(bed.gc.chr.end.start)
}

initRPKMCORRGC <- function(rpkm) {
   colnames <- c("BED", "BAM", "RPKM", "Ratio", "RPKM_CORR", "RPKM_CORR_GC")
   rpkm.corr.gc <- toTable(NA, length(colnames), nrow(rpkm), colnames)
   rpkm.corr.gc[,1:3] <- rpkm
   
   return(rpkm.corr.gc)
}

setRPKMCORR <- function(rpkm.corr.gc, segs, bed.gc) {
   #rpkm.corr.gc$Ratio <- mapply(r = 1:nrow(rpkm.corr.gc), function(r) getRatioFromSegment(bed.gc[r,], segs))
   #rpkm.corr.gc.ratio <- rpkm.corr.gc[!is.na(rpkm.corr.gc$Ratio),]    ## ~29 MINS on 2,861,558 entries; see deprecated method getRatioFromSegment()
   #bed.gc.ratio <- bed.gc[rownames(rpkm.corr.gc.ratio),]
   for (s in 1:nrow(segs)) {                                           ## ~5 MINS on 1,040 segments over 2,861,558 entries; see getBEDFromSegment()
      seg <- segs[s,]
      bed.gc.seg <- getBEDFromSegment(bed.gc, seg)
      bed.gc.seg.idx <- which(rpkm.corr.gc$BED %in% bed.gc.seg$BED)
  
      rpkm.corr.gc$Ratio[bed.gc.seg.idx] <- seg$Ratio
   }

   return(rpkm.corr.gc[which(rpkm.corr.gc$Ratio != 0),])   ## ADD 15/05/17: No information from chrY in female (e.g. S00035)
}                                                          ## 2852532 P2868190    0    0.0000000     0       NaN          NaN
                                                           ## 2852533 P2868191    1    0.9901062     0       Inf          Inf
getRPKMCORRGC <- function(rpkm, segs, bed, PAIR, CORR) {
   bed$BED <- rownames(bed)   ## To seed up line 56 and line 66
   overlaps <- intersect(rownames(rpkm), bed$BED)           ## ADD 24/06/17: Double-check if they have the same rows
   rpkm.gc <- rpkm[overlaps,]
   bed.gc <- bed[overlaps,]
   
   gc.mean <- mean(bed.gc$GC)
   rpkm.corr.gc <- initRPKMCORRGC(rpkm.gc)
 
   if (CORR) {   ## ADD 02/07/17: For LCL normals
      segs.gc <- subset(segs, CHR %in% unique(bed.gc$CHR))  ## ADD 24/06/17: If bed.gc.au, remove chrXY from segs accordingly
    
      rpkm.corr.gc <- setRPKMCORR(rpkm.corr.gc, segs.gc, bed.gc)
      rpkm.corr.gc$RPKM_CORR    <- rpkm.corr.gc$RPKM / (rpkm.corr.gc$Ratio / 2)  
      rpkm.corr.gc$RPKM_CORR_GC <- rpkm.corr.gc$RPKM_CORR / (bed.gc$GC / gc.mean)
   } else {
      rpkm.corr.gc$RPKM_CORR_GC <- rpkm.corr.gc$RPKM / (bed.gc$GC / gc.mean)
   }
    
   return(rpkm.corr.gc)
}

## For cmd-rt_rpkm.corr.gc_chr.R (command line version)
initReadDepthPerChromosome <- function(samples, bed.gc.chr) {
   colnames <- samples
   rt.chr <- toTable(NA, length(colnames), nrow(bed.gc.chr), colnames)
   rownames(rt.chr) <- rownames(bed.gc.chr)
 
   return(rt.chr)
}

getDetectedRD <- function(rpkm) {   ## Not dected (RPKM_CORR_GC = 0) windows in any of the samples 
   return(mapply(x = 1:nrow(rpkm), function(x) !any(as.numeric(rpkm[x,]) == 0)))
}

## pipeline read.rpkm.corr.gc.txt.gz() and getDetectedRD()
pipeGetDetectedRD <- function(BASE, chr, PAIR, samples) {
   rpkms.chr <- read.rpkm.corr.gc.txt.gz(paste0(wd.ngs.data, tolower(BASE), "_rpkm.corr.gc_", chr, "_", PAIR, ".txt.gz"))[,samples]   ## ADD samples 14/06/17
   rpkms.chr.d <- rpkms.chr[getDetectedRD(rpkms.chr),]   ## ADD getDetected() 13/06/17
 
   return(rpkms.chr.d)
}

outputRT <- function(rpkms.chr) {
   samples <- colnames(rpkms.chr)
   rpkms.chr$BED <- rownames(rpkms.chr)
   
   return(rpkms.chr[,c("BED", samples)])
}

plotRD <- function(wd.rt.plots, sample, file, chr, PAIR, rpkm.chr, bed.gc.chr, ext) {
   main.text <- paste0(sample, " read depth (copy number-, GC-corrected RPKM)")
   xlab.text <- paste0("Chromosome ", gsub("chr", "", chr), " coordinate (Mb)")
   ylab.text <- "Read depth"
   
   xmin <- 0
   xmax <- subset(chromInfo, chrom == chr)$size
   ymin <- min(rpkm.chr)
   ymax <- max(rpkm.chr)
   
   if (ext == "pdf")
      pdf(paste0(wd.rt.plots, chr, "/", sample, file, chr, "_", PAIR, ".pdf"), height=4, width=10)
   else if (ext == "png")
      png(paste0(wd.rt.plots, chr, "/", sample, file, chr, "_", PAIR, ".png"), height=4, width=10, units="in", res=300)   ## ADD 16/05/17: res=300
 
   plot(NULL, ylim=c(ymin, ymax), xlim=c(xmin/1E6, xmax/1E6), xlab=xlab.text, ylab=ylab.text, main=main.text)
   points(bed.gc.chr$START/1E6, rpkm.chr, col="red", cex=0.3)
   lines(bed.gc.chr$START/1E6, smooth.spline(rpkm.chr)$y)
 
   dev.off()
}

# -----------------------------------------------------------------------------
# Methods: RT
# Last Modified: 13/06/17
# -----------------------------------------------------------------------------
getLog2RDRatio <- function(rpkms.T.chr, rpkms.N.chr) {
   return(log2(as.numeric(rpkms.T.chr) + 0.01) - log2(as.numeric(rpkms.N.chr) + 0.01))
}

plotRT <- function(filename, title, chr, samples, xmin, xmax, rpkms.chr, bed.gc.chr, pair1, pair0, ext) {
   main.text <- paste0(title, " read depth (CN-, GC-corrected RPKM) ratio (", pair1, "/", pair0, ")")
   xlab.text <- paste0("Chromosome ", gsub("chr", "", chr), " coordinate (Mb)")
   ylab.text <- "Replication time (log2 FC)"
   
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   ymin <- min(rpkms.chr)
   ymax <- max(rpkms.chr)
   if (is.na(xmin)) xmin <- 0
   if (is.na(xmax)) xmax <- subset(chromInfo, chrom == chr)$size
 
   filename <- paste0(filename, "_", chr, "_n", samples)
   if (ext == "pdf")
      pdf(paste0(filename, ".pdf"), height=4, width=10)
   else if (ext == "png")
      png(paste0(filename, ".png"), height=4, width=10, units="in", res=300)   ## ADD 16/05/17: res=300
 
   plot(NULL, ylim=c(ymin, ymax), xlim=c(xmin/1E6, xmax/1E6), xlab=xlab.text, ylab=ylab.text, main=main.text)
   points(bed.gc.chr$START/1E6, rpkms.chr, col="red", cex=0.3)
   abline(h=0, lwd=0.5, col="grey")
   lines(bed.gc.chr$START/1E6, smooth.spline(rpkms.chr)$y)
   
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.4, col="lightgrey")
 
   dev.off()
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
# Inner Class  : In-house FileReader
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 08/05/17
# =============================================================================
read.bam.rpkm.txt.gz <- function(rpkm.file) {
   return(readTable(rpkm.file, header=T, rownames=T, sep=""))
}

read.rpkm.corr.gc.txt.gz <- function(rpkm.file) {
   return(readTable(rpkm.file, header=T, rownames=T, sep="")[,-1])
}

# =============================================================================
# Inner Class  : PeifLyne FileReader
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

read.peiflyne.mutcall.filtered.vcf <- function(vcf.file, pass, rs) {
   vcf <- readTable(vcf.file, header=F, rownames=F, sep="")
   colnames(vcf) <- c("CHROM",	"POS",	"ID",	"REF",	"ALT",	"QUAL",	"FILTER",	"INFO")
 
   if (pass)
      vcf <- subset(vcf, FILTER == "PASS")
   if (!rs)
      vcf <- subset(vcf, ID == ".")
   return(vcf)
}

# =============================================================================
# Inner Class  : Bedtools FileReader
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
# Inner Class  : Collections of test/obsolete/deprecated methods
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified:
# =============================================================================
getRatioFromSegment <- function(bed, segs) {
   #segs.chr <- subset(segs, CHR == bed$CHR)                          ## Using subset() takes ~48 MINS on 2,861,558 entries (S00022)
   #segs.chr.start <- subset(segs.chr, START <= bed$START)
   #segs.chr.start.end <- subset(segs.chr.start, END >= bed$END)
   segs.chr <- segs[which(segs$CHR == bed$CHR),]                    ## Using which() takes ~41 MINS on 2,861,558 entries (S00022)
   segs.chr.start <- segs.chr[which(segs.chr$START <= bed$START),]
   segs.chr.start.end <- segs.chr.start[which(segs.chr.start$END >= bed$END),]
 
   if (nrow(segs.chr.start.end) == 1)
      return(segs.chr.start.end$Ratio)
   else
      return(NA)
}

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
