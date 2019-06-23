#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
BASE1 <- args[1]   ## Cancer type
PAIR1 <- args[2]   ## T(umour) or S (phase)
NUM1  <- as.numeric(args[3])
BASE0 <- args[4]   ## Normal type
PAIR0 <- args[5]   ## N(ormal) or B(lood) or G1 (phase)
NUM0  <- as.numeric(args[6])
BSTRP <- as.numeric(args[7])
#CUTOFF <- as.numeric(args[8])
base1 <- tolower(BASE1)
base0 <- tolower(BASE0)

# =============================================================================
# Name: 2b_cmd-rt_rpkm.corr.gc.d.rt_bstrap.R (commandline mode)
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 22/10/18
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"         ## tyang2@cheops
#wd.src <- "/re/home/tyang2/dev/R"                ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"             ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "DifferentialExpression.R", "ReplicationForkDirectionality.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.1kb.gc.RData"))

# -----------------------------------------------------------------------------
# Replication timing
# Last Modified: 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
#wd <- "/projects/cangen/tyang2"   ## tyang2@cheops
wd <- "/scratch/tyang2"           ## tyang2@scratch
wd1.anlys <- file.path(wd, BASE1, "analysis")
wd1.rt    <- file.path(wd1.anlys, "replication", paste0(base1, "-wgs-rt"))
wd1.rt.data  <- file.path(wd1.rt, "data/bstrps", BSTRP)

bed.gc.rt <- bed.gc[1,]
bed.gc.rt$SPLINE <- 0
bed.gc.rt$SLOPE  <- 0
bed.gc.rt <- bed.gc.rt[-1,]

ensGene.rt <- ensGene[1,]
ensGene.rt$SLOPE_START <- 0
ensGene.rt$SLOPE_END   <- 0
ensGene.rt <- ensGene.rt[-1,]
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)

   ## Read in replicaiton time
   rpkms.chr.rt <-readTable(file.path(wd1.rt.data, paste0(base1, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", NUM1, "-", NUM0, ".txt.gz")), header=T, rownames=T, sep="\t")
   #rpkms.chr.rt <- rpkms.chr[which(rpkms.chr$RT > -CUTOFF),]
   #rpkms.chr.rt <- rpkms.chr.rt[which(rpkms.chr.rt$RT < CUTOFF),]
   
   overlaps <- intersect(rownames(rpkms.chr.rt), rownames(bed.gc.chr))
   bed.gc.chr.rt <- bed.gc.chr[overlaps,]
   bed.gc.chr.rt$SPLINE <- 0
   bed.gc.chr.rt$SLOPE  <- NA
   
   ## Determine replication direction for each expressed gene
   ## https://stackoverflow.com/questions/43615469/how-to-calculate-the-slope-of-a-smoothed-curve-in-r
   #slopes <- diff(smooth.spline(rpkms.chr.rt$MEDIAN)$y)/diff((bed.gc.chr.rt$START)/1E7)   ## BIG BUG FIX!!! 31/10/18: Warning message for using bed.gc.chr
   spline <- smooth.spline(x=bed.gc.chr.rt$START, y=rpkms.chr.rt$RT)
   slopes <- diff(spline$y)/diff(bed.gc.chr.rt$START/1E6)   ## ADD 31/10/18
   
   bed.gc.chr.rt$SPLINE <- spline$y
   bed.gc.chr.rt$SLOPE[1:length(slopes)] <- slopes   ## length(slopes) is 1 less than nrow(bed.gc.chr.rt), as no slope for the last 1KB window 
   #plotRT0(wd.rt.plots, BASE, chr, 101, 92, NA, NA, rpkms.chr.rt, bed.gc.chr.rt, PAIR1, PAIR0, "png", spline)
   bed.gc.rt <- rbind(bed.gc.rt, bed.gc.chr.rt)
   
   ensGene.rt.chr <- subset(ensGene, chromosome_name == chr)
   ensGene.rt.chr$SLOPE_START <- NA
   ensGene.rt.chr$SLOPE_END   <- NA   ## BUG 29/10/18
   for (g in 1:nrow(ensGene.rt.chr)) {
      gene <- ensGene.rt.chr[g,]
      bed.s <- getEnsGeneBED(gene$start_position, bed.gc.chr)   ## Not bed.gc.chr.rt
      bed.e <- getEnsGeneBED(gene$end_position, bed.gc.chr)
    
      if (length(bed.s) != 0) ensGene.rt.chr$SLOPE_START[g] <- bed.gc.chr.rt[bed.s[1],]$SLOPE
      if (length(bed.e) != 0) ensGene.rt.chr$SLOPE_END[g]   <- bed.gc.chr.rt[bed.e[1],]$SLOPE
   }
   ensGene.rt <- rbind(ensGene.rt, ensGene.rt.chr)
}
bed.gc.rt <- bed.gc.rt[, c("SPLINE", "SLOPE")]
ensGene.rt <- ensGene.rt[, c("SLOPE_START", "SLOPE_END")]
save(bed.gc.rt, ensGene.rt, file=file.path(wd1.rt.data, paste0(base1, "_ensGene.rt.RData")))
