#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
BASE   <- args[1]   ## Cancer type
SAMPLE <- args[2]   ## Sample ID
SAMPLE2 <- args[3]   ## Sample ID
PAIR   <- args[4]   ## T(umour) or N(ormal) pair
CN     <- as.logical(args[5])   ## T(rue) or F(alse) to correct CN
METHOD <- args[6]
base   <- tolower(BASE)
method <- tolower(METHOD)

# =============================================================================
# Name: 
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 04/12/23; 11/06/17
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"         ## tyang2@cheops
#wd.src <- "/re/home/tyang2/dev/R"                ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"             ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "ReplicationTiming.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.bed.gc.icgc.RData"))

# -----------------------------------------------------------------------------
# 1a: To calculate normalised read counts (NRD)
# Last Modified: 02/08/19; 01/05/17
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"   ## tyang2@cheops
wd.ngs <- file.path(wd, "ICGC", BASE, "ngs/WGS")
wd.ngs.data <- file.path(wd.ngs, "data")

#samples <- readTable(file.path(wd.ngs, LIST), header=T, rownames=2, sep="")

## Read coverage from PeifLyne
#for (s in 1:length(samples$wgs_id)) {
   wgs_id      <- SAMPLE2
   specimen_id <- SAMPLE
 
   ## INPUT: *_cn.txt.gz (read counts)
   rd <- read.peiflyne.icgc.cn.txt(file.path("/projects/cangen/PCAWG-Repository/PCAWG.raw.data/copy_number/converted_data", paste0(wgs_id, "_CONVERTED"), paste0(wgs_id, "_cn.txt")))   ## See ReplicationTiming.R / Inner Class / PeifLyne File Reader
   nrd <- initNRD(rd, bed.gc, pair="T", method=METHOD)   ## See ReplicationTiming.R
   #nrd.N <- initNRD(rd, bed.gc, pair="N", method=METHOD)
   rm(rd)
 
   #writeTable(nrd.T, gzfile(file.path(wd.ngs.data, specimen_id, paste0(specimen_id, "_T.", method, ".txt.gz"))), colnames=T, rownames=F, sep="\t")
   #writeTable(nrd.N, gzfile(file.path(wd.ngs.data, sample, paste0(sample, "_N.", method, ".txt.gz"))), colnames=T, rownames=F, sep="\t")
   #rm(nrd.T)
#}

# -----------------------------------------------------------------------------
# 1b: Calculate absolute RPKM
# Last Modified: 14/05/17
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"   ## tyang2@cheops
wd.ngs      <- file.path(wd, "ICGC", BASE, "ngs/WGS")
wd.ngs.data <- file.path(wd.ngs, "data")
wd.anlys   <- file.path(wd, "ICGC", BASE, "analysis")
wd.rt      <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data <- file.path(wd.rt, "data")
   
## Read copy number information from PeifLyne
#for (s in 1:length(samples)) {
   #sample <- samples[s]
   sample  <- SAMPLE 
   sample2 <- SAMPLE2
   
   ## INPUT: RPKM (from 1a) & SEQ
   #nrd <- readTable(file.path(wd.ngs.data, sample, paste0(sample, "_", PAIR, ".", method, ".txt.gz")), header=T, rownames=T, sep="")
   segs <- NA
   if (CN)
      segs <- read.peiflyne.cn.seg(file.path("/projects/cangen/PCAWG-Repository/PCAWG.raw.data/copy_number/sclust_final_copy_number_analysis_files", sample2, paste0(sample2, "_cn.seg")))   ## See ReplicationTiming.R / Inner Class / PeifLyne File Reader
   
   ## Correcte read counts fpr copy number
   nrd.cn <- getNRDGCCN(nrd, segs, bed.gc, PAIR, F, CN)   ## See ReplicationTiming.R
   #writeTable(nrd.cn, gzfile(file.path(wd.ngs.data, sample, paste0(sample, "_", PAIR, ".", method, ".cn.txt.gz"))), colnames=T, rownames=F, sep="\t")
#}
#nrd.cn <- readTable(file.path(wd.ngs.data, sample, paste0(sample, "_", PAIR, ".", method, ".gc.cn.txt.gz")), header=T, rownames=F, sep="\t")
rm(segs)

# -----------------------------------------------------------------------------
# Replication timing
# Last Modified: 21/04/19; 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"   ## tyang2@cheops
wd.ngs    <- file.path(wd, "ICGC", BASE, "ngs/WGS")
wd.ngs.data <- file.path(wd.ngs, "data") 
   
wd.anlys <- file.path(wd, "ICGC", BASE, "analysis")
wd.rt    <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data <- file.path(wd.rt, "data")
   
#samples <- readTable(file.path(wd1.ngs, LIST), header=F, rownames=F, sep="")
#samples <- toupper(gsub("-", "", samples))   ## ADD 03/07/17 for LCL (e.g. NA19240-2 to NA19240.2)   ## ADDED for NBL CCL
#n <- length(samples)
   
## LCL S/G1
load(file.path(wd, "LCL/analysis/replication/lcl-wgs-rt/data/rd-vs-rt_lcl-s-g1_spearman.RData"))
sprs.order <- sprs[order(abs(sprs$cor1), decreasing=T),]
rownames(sprs.order) <- 1:22

load(file.path(wd, "LCL/analysis/replication/lcl-wgs-rt/data/lcl_rpkm.cn.m.rt.log2s_s-g1.RData"))
nrds <- nrds[which(is.na(nrds$RT) == F),]   ## ADD 29/11/20

nrds.T.chr.m.sample <- nrd.cn[,c("BED", "NRD_GC_CN")]
rm(nrd.cn)
colnames(nrds.T.chr.m.sample) <- c("BED", "T")
rownames(nrds.T.chr.m.sample) <- nrds.T.chr.m.sample$BED

nrds.T.chr.m.sample.T.all <- NULL
nrds.chr.RT.all <- NULL
cors <- toTable(0, 3, 22, c("chr", "cor0", "cor"))
cors$chr <- sprs.order$chr
for (c in 1:nrow(sprs.order)) {
	  chr <- chrs[sprs.order$chr[c]]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
    
   ## Read depth
   #nrds.T.chr.m <- readTable(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.m_", chr, "_", PAIR, "_n", N, ".txt.gz")), header=T, rownames=T, sep="")
   #colnames(nrds.T.chr.m) <- gsub("\\.", "-", colnames(nrds.T.chr.m))   ## ADD 29/11/20; 05/10/19
   #nrds.T.chr.m <- nrds.T.chr.m[-which(nrds.T.chr.m$MEDIAN == 0),]      ## ADD 29/11/20
    
   ## Individual read depth
   #SAMPLE <- gsub("-", "", SAMPLE)     ## ADD 05/10/19
   #SAMPLE <- gsub("\\.", "", SAMPLE)   ## ADD 05/10/19
   #SAMPLE <- toupper(SAMPLE)           ## ADD 05/10/19
    
   #nrds.T.chr.d.sample$T <- log2(nrds.T.chr.d.sample$T + 0.01)
   overlaps <- intersect(rownames(bed.gc.chr), nrds.T.chr.m.sample$BED)
   nrds.T.chr.m.sample.T <- setSpline(nrds.T.chr.m.sample[overlaps,], bed.gc.chr[overlaps,], "T")
    
   ## Replication timing
   #load(file.path(wd, paste0("LCL/analysis/replication/lcl-wgs-rt/data/lcl_rpkm.cn.m.rt.log2s_s-g1_icgc_", chr, ".RData")))   ## 26/02/22: ADD
   #rownames(nrds.chr) <- nrds.chr$BED
   #overlaps <- intersect(rownames(bed.gc.chr), nrds.chr$BED)
   #nrds.chr.RT <- setSpline(nrds.chr[overlaps,], bed.gc.chr[overlaps,], "RT")   ## Reference LCL S/G1 ratio
   nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]   ## Reference LCL S/G1 ratio
   nrds.chr.RT <- setSpline(nrds.chr, bed.gc.chr, "RT")
   
   ## Keep 1kb slopes based on overlapping windows
   overlaps <- intersect(nrds.T.chr.m.sample.T$BED, nrds.chr.RT$BED)
   cors$cor0[c] <- getCor(nrds.chr.RT[overlaps,]$SPLINE, nrds.T.chr.m.sample.T[overlaps,]$SPLINE, method="spearman")
    
   ##
   if (is.null(nrds.T.chr.m.sample.T.all)) {
      nrds.T.chr.m.sample.T.all <- nrds.T.chr.m.sample.T[overlaps,]
   } else {
      nrds.T.chr.m.sample.T.all <- rbind(nrds.T.chr.m.sample.T.all, nrds.T.chr.m.sample.T[overlaps,])
   }
    
   if (is.null(nrds.chr.RT.all)) {
      nrds.chr.RT.all <- nrds.chr.RT[overlaps,]
   } else {
      nrds.chr.RT.all <- rbind(nrds.chr.RT.all, nrds.chr.RT[overlaps,])
   }
   
   cors$cor[c] <- getCor(nrds.chr.RT.all$RT, nrds.T.chr.m.sample.T.all$T, method="spearman")
   rm(nrds.T.chr.m.sample.T)
   rm(nrds.chr.RT)
}
cor.all <- getCor(nrds.chr.RT.all$SPLINE, nrds.T.chr.m.sample.T.all$SPLINE, method="spearman")
save(cor.all, cors, file=file.path(wd.rt.data, "samples", paste0("rd-vs-rt_", SAMPLE, "-vs-lcl_spearman.RData")))
