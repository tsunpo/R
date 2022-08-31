#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
BASE   <- args[1]   ## Cancer type
SAMPLE <- args[2]   ## Sample ID
PAIR   <- args[3]   ## T(umour) or N(ormal) pair
CN     <- as.logical(args[4])   ## T(rue) or F(alse) to correct CN
METHOD <- args[5]
base   <- tolower(BASE)
method <- tolower(METHOD)

# =============================================================================
# Name: 
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 11/06/17
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"         ## tyang2@cheops
#wd.src <- "/re/home/tyang2/dev/R"                ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"             ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "ReplicationTiming.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.bed.gc.peiflyne.wes.new.RData"))

# -----------------------------------------------------------------------------
# 1a: To calculate normalised read counts (NRD)
# Last Modified: 02/08/19; 01/05/17
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"   ## tyang2@cheops
wd.ngs <- file.path(wd, BASE, "ngs/WES_NEW")
wd.ngs.data <- file.path(wd.ngs, "data")

#samples <- readTable(file.path(wd.ngs, LIST), header=T, rownames=2, sep="")

## Read coverage from PeifLyne
#for (s in 1:length(samples$wgs_id)) {
   sample <- SAMPLE
 
   ## INPUT: *_cn.txt.gz (read counts)
   #rd <- read.peiflyne.wes.cn.txt(file.path("/projects/cangen/data/NB/WES/Cologne", SAMPLE, paste0(SAMPLE, "_ANALYSIS"), paste0(SAMPLE, "_cn.txt")))   ## See ReplicationTiming.R / Inner Class / PeifLyne File Reader
   #rd <- read.peiflyne.cn.txt(file.path("/projects/cangen/tyang2/SCLC/ngs/WES_1K", SAMPLE, paste0(SAMPLE, "_ANALYSIS"), paste0(SAMPLE, "_cn.txt")))   ## See ReplicationTiming.R / Inner Class / PeifLyne File Reader
   rd <- read.peiflyne.wes.cn.txt(file.path("/projects/cangen/tyang2/SCLC/ngs/WES", SAMPLE, paste0(SAMPLE, "_ANALYSIS"), paste0(SAMPLE, "_cn.txt")))   ## See ReplicationTiming.R / Inner Class / PeifLyne File Reader
   
   nrd <- initNRD(rd, bed.gc, pair=PAIR, method=METHOD)   ## See ReplicationTiming.R
   #nrd <- subset(nrd, RD > 200)   ## ADD 09/06/22
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
#wd <- "/projects/cangen/tyang2"   ## tyang2@cheops
#wd.ngs      <- file.path(wd, BASE, "ngs/WES")
#wd.ngs.data <- file.path(wd.ngs, "data")

wd.anlys   <- file.path(wd, BASE, "analysis")
wd.rt      <- file.path(wd.anlys, "replication", paste0(base, "-wes-rt"))
wd.rt.data <- file.path(wd.rt, "data_new")
   
## Read copy number information from PeifLyne
#for (s in 1:length(samples)) {
   #sample <- samples[s]
   sample  <- SAMPLE 
   
   ## INPUT: RPKM (from 1a) & SEQ
   #nrd <- readTable(file.path(wd.ngs.data, sample, paste0(sample, "_", PAIR, ".", method, ".txt.gz")), header=T, rownames=T, sep="")
   segs <- NA
   if (CN)
      #segs <- read.peiflyne.cn.seg(file.path("/projects/cangen/data/NB/WES/Cologne", SAMPLE, paste0(SAMPLE, "_ANALYSIS"), paste0(SAMPLE, "_cn.seg")))   ## See ReplicationTiming.R / Inner Class / PeifLyne File Reader
      #segs <- read.peiflyne.cn.seg(file.path("/projects/cangen/tyang2/SCLC/ngs/WES_1K", SAMPLE, paste0(SAMPLE, "_ANALYSIS"), paste0(SAMPLE, "_uncorr_cn.seg")))   ## See ReplicationTiming.R / Inner Class / PeifLyne File Reader
      segs <- read.peiflyne.cn.seg(file.path("/projects/cangen/tyang2/SCLC/ngs/WES", SAMPLE, paste0(SAMPLE, "_ANALYSIS"), paste0(SAMPLE, "_uncorr_cn.seg")))   ## See ReplicationTiming.R / Inner Class / PeifLyne File Reader
   
   ## Correcte read counts fpr copy number
   nrd.gc.cn <- getNRDGCCN(nrd, segs, bed.gc, PAIR, CN)   ## See ReplicationTiming.R
   writeTable(nrd.gc.cn, gzfile(file.path(wd.ngs.data, SAMPLE, paste0(SAMPLE, "_", PAIR, ".", method, ".gc.cn.txt.gz"))), colnames=T, rownames=F, sep="\t")
#}
#nrd.gc.cn <- readTable(file.path(wd.ngs.data, sample, paste0(sample, "_", PAIR, ".", method, ".gc.cn.txt.gz")), header=T, rownames=F, sep="\t")
rm(segs)

# -----------------------------------------------------------------------------
# Replication timing
# Last Modified: 21/04/19; 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
#wd <- "/projects/cangen/tyang2"   ## tyang2@cheops
#wd.ngs    <- file.path(wd, BASE, "ngs/WES")
#wd.ngs.data <- file.path(wd.ngs, "data") 
   
#wd.anlys <- file.path(wd, BASE, "analysis")
#wd.rt    <- file.path(wd.anlys, "replication", paste0(base, "-wes-rt"))
#wd.rt.data <- file.path(wd.rt, "data")
   
#samples <- readTable(file.path(wd1.ngs, LIST), header=F, rownames=F, sep="")
#samples <- toupper(gsub("-", "", samples))   ## ADD 03/07/17 for LCL (e.g. NA19240-2 to NA19240.2)   ## ADDED for NBL CCL
#n <- length(samples)
   
## LCL S/G1
#load(file.path(wd, "LCL/analysis/replication/lcl-wgs-rt/data/lcl_rpkm.gc.cn.d.rt.log2s_s-g1.RData"))
#nrds <- nrds[which(is.na(nrds$RT) == F),]   ## ADD 29/11/20

nrds.T.chr.m.sample <- nrd.gc.cn[,c("BED", "NRD_GC_CN")]
rm(nrd.gc.cn)
colnames(nrds.T.chr.m.sample) <- c("BED", "T")
rownames(nrds.T.chr.m.sample) <- nrds.T.chr.m.sample$BED

nrds.T.chr.m.sample.T.all <- NULL
nrds.chr.RT.all <- NULL
cors <- toTable(0, 2, 22, c("chr", "cor"))
cors$chr <- 1:22
for (c in 1:22) {
   chr <- chrs[c]
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
   nrds.T.chr.m.sample.T <- setSpline0(nrds.T.chr.m.sample[overlaps,], bed.gc.chr[overlaps,], "T")
    
   ## Replication timing
   load(file.path(wd, paste0("LCL/analysis/replication/lcl-wgs-rt/data/lcl_rpkm.gc.cn.d.rt.log2s_s-g1_peiflyne_wes_new", chr, ".RData")))   ## 26/02/22: ADD
   rownames(nrds.chr) <- nrds.chr$BED
   #nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]   ## Reference LCL S/G1 ratio
   overlaps <- intersect(rownames(bed.gc.chr), nrds.chr$BED)
   nrds.chr.RT <- setSpline0(nrds.chr[overlaps,], bed.gc.chr[overlaps,], "RT")   ## Reference LCL S/G1 ratio
    
   ## Keep 1kb slopes based on overlapping windows
   overlaps <- intersect(nrds.T.chr.m.sample.T$BED, nrds.chr.RT$BED)
   cors$cor[c] <- getCor(nrds.chr.RT[overlaps,]$SPLINE, nrds.T.chr.m.sample.T[overlaps,]$T, method="spearman")
    
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
   
   rm(nrds.T.chr.m.sample.T)
   rm(nrds.chr.RT)
}
cor <- getCor(nrds.chr.RT.all$SPLINE, nrds.T.chr.m.sample.T.all$T, method="spearman")
save(cor, cors, file=file.path(wd.rt.data, "samples", paste0("rd-vs-rt_", SAMPLE, "-vs-lcl_spline_spearman.RData")))
   