#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
CHR  <- as.numeric(args[1])

# =============================================================================
# Manuscript   : 
# Chapter      : 
# Name         : manuscripts/replication/icgc-wgs-rt/icgc-wgs.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 19/02/22
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "ReplicationTiming.R", "Transcription.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.bed.gc.icgc.RData"))

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 18/02/22
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
#wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE  <- "ICGC"
base  <- tolower(BASE)

wd.ngs <- file.path(wd, BASE, "ngs/WGS")
wd.meta  <- file.path(wd, BASE, "metadata")

# -----------------------------------------------------------------------------
# Convert reference RT to ICGC partitioning
# Last Modified: 24/02/22
# -----------------------------------------------------------------------------
load(file.path(wd, "LCL/analysis/replication/lcl-wgs-rt/data/lcl_rpkm.gc.cn.d.rt.log2s_s-g1.RData"))
nrds.lcl <- nrds
# > nrow(nrds.lcl)
# [1] 2582940
load(file.path(wd.src.ref, "hg19.bed.gc.1kb.RData"))
# > nrow(bed.gc)
# [1] 2861558
bed.gc.o <- cbind(nrds.lcl, bed.gc[rownames(nrds.lcl),])

###
## Same as getNRDGCCN() in ReplicationTiming.R to remove chrXY
load(file.path(wd.src.ref, "hg19.bed.gc.icgc.RData"))
# > nrow(bed.gc)
# [1] 4062897
segs.gc  <- subset(bed.gc, CHR %in% chrs[1:22])
# > nrow(bed.gc.o)
# [1] 3810092

###
## Same as setNRDCN() in ReplicationTiming.R
nrds <- toTable(NA, 4, 0, c("BED", "T", "N", "RT"))
chr <- chrs[CHR]

segs.gc  <- subset(segs.gc, CHR == chr)
bed.gc.o <- subset(bed.gc.o, CHR == chr)

for (s in 1:nrow(segs.gc)) {
   seg.gc <- segs.gc[s,]
   bed.gc.seg <- getBEDFromSegment(bed.gc.o, seg.gc)
  
   if (nrow(bed.gc.seg) != 0) {
      nrds.s <- toTable(NA, 4, 1, c("BED", "T", "N", "RT"))
    
      if (nrow(bed.gc.seg) == 1) {
         nrds.s[1,] <- c(rownames(seg.gc), bed.gc.seg[, c("T", "N", "RT")])
      } else {
         nrds.s$BED[1] <- rownames(seg.gc)
         nrds.s$T[1]  <- sum(bed.gc.seg$T)/nrow(bed.gc.seg)  
         nrds.s$N[1]  <- sum(bed.gc.seg$N)/nrow(bed.gc.seg)  
         nrds.s$RT[1] <- sum(bed.gc.seg$RT)/nrow(bed.gc.seg)  
      }
      nrds<- rbind(nrds, nrds.s)
   }
}
nrds.chr <- nrds
save(nrds.chr, file=file.path(wd, paste0("LCL/analysis/replication/lcl-wgs-rt/data/lcl_rpkm.gc.cn.d.rt.log2s_s-g1_icgc_", chr, ".RData")))

###
##
load(file=file.path(wd, paste0("LCL/analysis/replication/lcl-wgs-rt/data/lcl_rpkm.gc.cn.d.rt.log2s_s-g1_icgc_chr1.RData")))
nrds <- nrds.chr
for (c in 2:22) {
   chr <- chrs[c]
   load(file=file.path(wd, paste0("LCL/analysis/replication/lcl-wgs-rt/data/lcl_rpkm.gc.cn.d.rt.log2s_s-g1_icgc_", chr, ".RData")))
   
   nrds <- rbind(nrds, nrds.chr)
}
save(nrds, file=file.path(wd, paste0("LCL/analysis/replication/lcl-wgs-rt/data/lcl_rpkm.gc.cn.d.rt.log2s_s-g1_icgc.RData")))
