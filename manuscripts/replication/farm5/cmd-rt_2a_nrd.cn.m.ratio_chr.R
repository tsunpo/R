#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
BASE  <- args[1]   ## Cancer type
PAIR1 <- args[2]   ## T(umour) or S (phase)
LIST1 <- args[3]
PAIR0 <- args[4]   ## N(ormal) or B(lood) or G1 (phase)
LIST0 <- args[5]
CHR   <- as.numeric(args[6])
METHOD <- args[7]
RATIO  <- args[8]  ## NA or RATIO
base   <- tolower(BASE)
method <- tolower(METHOD)

# =============================================================================
# Name: cmd-rt_2a_nrd.gc.cn.d_chr.R (commandline mode)
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 29/11/20; 23/04/19; 22/10/18
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
# Last Modified: 23/04/19; 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
wd <- "/lustre/scratch127/casm/team294rr/ty2"   ## ty2@farm5
#wd <- "/projects/cangen/tyang2"   ## tyang2@cheops
wd.ngs    <- file.path(wd, BASE, "ngs/WGS")
wd.ngs.data <- file.path(wd.ngs, "data")

wd.anlys <- file.path(wd, BASE, "analysis")
wd.rt    <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data <- file.path(wd.rt, "data")

samples1 <- readTable(file.path(wd.ngs, LIST1), header=F, rownames=F, sep="")
samples0 <- readTable(file.path(wd.ngs, LIST0), header=F, rownames=F, sep="")
samples1 <- gsub("\\.", "-", samples1)
samples0 <- gsub("\\.", "-", samples0)
n1 <- length(samples1)
n0 <- length(samples0)

#for (c in 1:22) {
   chr <- chrs[CHR]

   ## Read depth
   nrds.T.chr.m <- pipeMedianRD(wd.ngs.data, BASE, chr, PAIR1, method, samples1)   ## ADD BACK 09/04/19; REMOVED 15/02/19; if length(samples) == 1
   nrds.N.chr.m <- pipeMedianRD(wd.ngs.data, BASE, chr, PAIR0, method, samples0)
   writeTable(nrds.T.chr.m[, c("BED", samples1, "MEDIAN")], gzfile(file.path(wd.rt.data, paste0(base, "_", method, ".cn.m_", chr, "_", PAIR1, "_n", n1, ".txt.gz"))), colnames=T, rownames=F, sep="\t")
   writeTable(nrds.N.chr.m[, c("BED", samples0, "MEDIAN")], gzfile(file.path(wd.rt.data, paste0(base, "_", method, ".cn.m_", chr, "_", PAIR0, "_n", n0, ".txt.gz"))), colnames=T, rownames=F, sep="\t")
   
   ## Replication timing
   if (RATIO != "NA") {
      nrds.T.chr.m <- nrds.T.chr.m[-which(nrds.T.chr.m$MEDIAN == 0),]      ## ADD 29/11/20
      nrds.N.chr.m <- nrds.N.chr.m[-which(nrds.N.chr.m$MEDIAN == 0),]      ## ADD 29/11/20
    
      overlaps <- intersect(rownames(nrds.T.chr.m), rownames(nrds.N.chr.m))
      nrds.T.chr.m.o <- nrds.T.chr.m[overlaps,]
      nrds.N.chr.m.o <- nrds.N.chr.m[overlaps,]
   
      nrds.chr <- toTable(NA, 3, length(overlaps), c("T", "N", "RT"))
      rownames(nrds.chr) <- overlaps
      nrds.chr$T  <- nrds.T.chr.m.o$MEDIAN
      nrds.chr$N  <- nrds.N.chr.m.o$MEDIAN
      nrds.chr$RT <- mapply(x = 1:nrow(nrds.chr), function(x) (nrds.chr[x,]$T / nrds.chr[x,]$N))

      writeTable(outputRT(nrds.chr), gzfile(file.path(wd.rt.data, paste0(base, "_", method, ".cn.m.ratio_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ".txt.gz"))), colnames=T, rownames=F, sep="\t")
   }
#}
