#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
BASE <- args[1]   ## Cancer type
PAIR <- args[2]   ## T(umour) or S (phase)
LIST <- args[3]
CHR  <- as.numeric(args[4])
N    <- as.numeric(args[5])
base <- tolower(BASE)

# =============================================================================
# Name: 2a_cmd-rt_rpkm.corr.gc.d.ensGene_chr.R (commandline mode)
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 22/05/19
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"         ## tyang2@cheops
#wd.src <- "/re/home/tyang2/dev/R"                ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"             ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "ReplicationTiming.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.1kb.gc.RData"))

# -----------------------------------------------------------------------------
# Replication timing
# Last Modified: 22/05/19
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"   ## tyang2@cheops
wd.ngs    <- file.path(wd, BASE, "ngs/WGS")
wd.ngs.data <- file.path(wd.ngs, "data")

wd.anlys   <- file.path(wd, BASE, "analysis")
wd.rt      <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data <- file.path(wd.rt, "data")

samples <- readTable(file.path(wd.ngs, LIST), header=T, rownames=T, sep="")

#for (c in 1:22) {
   chr <- chrs[CHR]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   ensGene.chr <- subset(ensGene, chromosome_name == chr)
   
   ## Replication timing
   rpkms.chr.d <- readTable(file.path(wd.rt.data, paste0(base, "_rpkm.corr.gc.d_", chr, "_", PAIR, "_n", N, ".txt.gz")), header=T, rownames=T, sep="\t")

   overlaps <- intersect(rownames(bed.gc.chr), rownames(rpkms.chr.d))
   bed.gc.chr <- bed.gc.chr[overlaps,]
   rpkms.chr.d <- rpkms.chr.d[overlaps, rownames(samples)]
   
   rpkms.chr.d.genes <- rpkms.chr.d[1,][-1,]
   for (g in 1:nrow(ensGene.chr)) {
      ensGene.chr.gene <- ensGene.chr[g,]
      
      bed.gc.chr.start <- subset(bed.gc.chr, END >= ensGene.chr.gene$start_position)
      bed.gc.chr.start.end <- subset(bed.gc.chr.start, START <= ensGene.chr.gene$end_position)
      
      if (nrow(bed.gc.chr.start.end) != 0) {
         rpkms.chr.d.1kb  <- rpkms.chr.d[rownames(bed.gc.chr.start.end),]
         rpkms.chr.d.gene <- rpkms.chr.d.1kb[1,]
         rownames(rpkms.chr.d.gene) <- rownames(ensGene.chr.gene)
      
         rpkms.chr.d.gene[1,] <- 0
         rpkms.chr.d.gene[1,] <- mapply(x = 1:ncol(rpkms.chr.d.1kb), function(x) sum(rpkms.chr.d.1kb[,x]))
      
         rpkms.chr.d.genes <- rbind(rpkms.chr.d.genes, rpkms.chr.d.gene)
      }
   }
   
   writeTable(rpkms.chr.d.genes, gzfile(file.path(wd.rt.data, paste0(base, "_rpkm.corr.gc.d.ensGene_", chr, ".txt.gz"))), colnames=T, rownames=T, sep="\t")
#}
