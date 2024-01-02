#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
BASE <- args[1]   ## Cancer type
PAIR <- args[2]   ## T(umour) or S (phase)
SAMPLE <- args[3]
#SAMPLE <- toupper(gsub("-", "", SAMPLE))   ## ADDED for NBL CCL
#N      <- args[4]
METHOD <- "RPKM"
base   <- tolower(BASE)
method <- "rpkm"

# =============================================================================
# Name: cmd-rt_2a_nrd.gc.cn.d_sample.R (commandline mode)
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 29/11/20; 22/10/18
# =============================================================================
wd.src <- "/nfs/users/nfs_t/ty2/dev/R"            ## ty2@farm5
#wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/re/home/tyang2/dev/R"                ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"             ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "ReplicationTiming.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.bed.gc.1kb.RData"))

# -----------------------------------------------------------------------------
# Replication timing
# Last Modified: 21/04/19; 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
wd <- "/lustre/scratch127/casm/team294rr/ty2"   ## ty2@farm5
#wd <- "/projects/cangen/tyang2"   ## tyang2@cheops
wd.ngs    <- file.path(wd, BASE, "ngs/WGS")
wd.ngs.data <- file.path(wd.ngs, "data") 

wd.anlys <- file.path(wd, BASE, "analysis")
wd.rt    <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data <- file.path(wd.rt, "data")

#samples <- readTable(file.path(wd1.ngs, LIST), header=F, rownames=F, sep="")
#samples <- toupper(gsub("-", "", samples))   ## ADD 03/07/17 for LCL (e.g. NA19240-2 to NA19240.2)   ## ADDED for NBL CCL
#n <- length(samples)

## LCL S/G1
load(file.path(wd, "LCL/analysis/replication/lcl-wgs-rt/data/rd-vs-rt_lcl-s-g1_spline_spearman.RData"))
sprs.order <- sprs[order(sprs$cor1, decreasing=T),]

load(file.path(wd, "LCL/analysis/replication/lcl-wgs-rt/data/lcl_rpkm.cn.m.rt.log2s_s-g1.RData"))
nrds <- nrds[which(is.na(nrds$RT) == F),]   ## ADD 29/11/20

nrds.T.chr.m.sample.T.all <- NULL
nrds.chr.RT.all <- NULL
cors <- toTable(0, 3, 22, c("chr", "cor0", "cor"))
cors$chr <- 1:22
for (c in 1:nrow(sprs.order)) {
   chr <- chrs[sprs.order$chr[c]]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   
   ## Read depth
   nrds.T.chr.m <- readTable(file.path(wd.ngs.data, paste0(base, "_", method, ".cn_", chr, "_", PAIR, ".txt.gz")), header=T, rownames=T, sep="")
   colnames(nrds.T.chr.m) <- gsub("\\.", "-", colnames(nrds.T.chr.m))   ## ADD 29/11/20; 05/10/19
   #nrds.T.chr.m <- nrds.T.chr.m[-which(nrds.T.chr.m$MEDIAN == 0),]      ## ADD 29/11/20

   ## Individual read depth
   #SAMPLE <- gsub("-", "", SAMPLE)     ## ADD 05/10/19
   #SAMPLE <- gsub("\\.", "", SAMPLE)   ## ADD 05/10/19
   #SAMPLE <- toupper(SAMPLE)           ## ADD 05/10/19
   
   nrds.T.chr.m.sample <- nrds.T.chr.m[,c("BED", SAMPLE)]
   colnames(nrds.T.chr.m.sample) <- c("BED", "T")
   #nrds.T.chr.d.sample$T <- log2(nrds.T.chr.d.sample$T + 0.01)
   nrds.T.chr.m.sample.T <- setSpline(nrds.T.chr.m.sample, bed.gc.chr, "T")

   ## Replication timing
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
   
   cors$cor[c] <- getCor(nrds.chr.RT.all$SPLINE, nrds.T.chr.m.sample.T.all$SPLINE, method="spearman")
}
cor <- getCor(nrds.chr.RT.all$SPLINE, nrds.T.chr.m.sample.T.all$SPLINE, method="spearman")
save(cor, cors, file=file.path(wd.rt.data, "samples", paste0("rd-vs-rt_", SAMPLE, "-vs-lcl_spline_spearman.RData")))
