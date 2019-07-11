#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
BASE1 <- args[1]   ## Cancer type
PAIR1 <- args[2]   ## T(umour) or S (phase)
NUM1  <- as.numeric(args[3])
BASE0 <- args[4]   ## Normal type
PAIR0 <- args[5]   ## N(ormal) or B(lood) or G1 (phase)
NUM0  <- as.numeric(args[6])
CHR  <- as.numeric(args[7])
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
handbooks  <- c("Commons.R", "DifferentialExpression.R", "ReplicationTiming.R", "ReplicationForkDirectionality.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.1kb.gc.RData"))

# -----------------------------------------------------------------------------
# Replication timing
# Last Modified: 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"   ## tyang2@cheops
wd1.anlys <- file.path(wd, BASE1, "analysis")
wd1.rt    <- file.path(wd1.anlys, "replication", paste0(base1, "-wgs-rt"))
wd1.rt.data  <- file.path(wd1.rt, "data")

# -----------------------------------------------------------------------------
# TO-DO (require "ReplicationForkDirectionality.R")
# Last Modified: 24/06/19
# -----------------------------------------------------------------------------
#ensGene.rt <- ensGene[1,]
#ensGene.rt$SPLINE_START <- NA
#ensGene.rt$SPLINE_END   <- NA
#ensGene.rt <- ensGene.rt[-1,]
#for (c in 1:22) {
   #chr <- chrs[c]
   chr <- chrs[CHR]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   
   ## SCLC CM2/CM1 RT
   rpkms.chr.rt <- readTable(file.path(wd1.rt.data, paste0(base1, "_rpkm.corr.gc.d.rt.cm2_", chr, "_", PAIR1, "-", PAIR0, "_n", NUM1, "-", NUM0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, scaledRT=T) 
   rpkms.chr.rt.RT <- setSpline(rpkms.chr.rt, bed.gc.chr, "RT")
 
   ##
   ensGene.rt.chr <- subset(ensGene, chromosome_name == chr)
   ensGene.rt.chr$SPLINE_START <- NA
   ensGene.rt.chr$SPLINE_END   <- NA   ## BUG 29/10/18
   for (g in 1:nrow(ensGene.rt.chr)) {
      gene <- ensGene.rt.chr[g,]
      bed.s <- getEnsGeneBED(gene$start_position, bed.gc.chr)   ## Not bed.gc.chr.rt
      bed.e <- getEnsGeneBED(gene$end_position, bed.gc.chr)
  
     if (length(bed.s) != 0) ensGene.rt.chr$SPLINE_START[g] <- rpkms.chr.rt.RT[bed.s[1],]$SPLINE
     if (length(bed.e) != 0) ensGene.rt.chr$SPLINE_END[g]   <- rpkms.chr.rt.RT[bed.e[1],]$SPLINE
   }
   #ensGene.rt <- rbind(ensGene.rt, ensGene.rt.chr)
   save(ensGene.rt.chr, file=file.path(wd1.rt.data, "chrs", paste0(base1, "_ensGene.rt.chr", CHR, ".RData")))
#}
