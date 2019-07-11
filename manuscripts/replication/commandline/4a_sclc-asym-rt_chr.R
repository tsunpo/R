#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
CHR  <- as.numeric(args[1])

# =============================================================================
# Manuscript   : 
# Chapter      : Chromosome replication timing of the human genome
# Name         : manuscripts/asymmetries/cll-asym-tx.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 20/06/18
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "Asymmetry.R", "DifferentialExpression.R", "Mutation.R", "ReplicationTiming.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.1kb.gc.RData"))

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 29/01/18
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"                   ## tyang2@gauss
#wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE <- "SCLC"
PAIR1 <- "T"
PAIR0 <- "T"
base <- tolower(BASE)

wd.ngs   <- file.path(wd, BASE, "ngs/WGS")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.asym       <- file.path(wd.anlys, "asymmetries", paste0(base, "-asym-rt"))
wd.asym.files <- file.path(wd.asym, "files")
wd.asym.data  <- file.path(wd.asym, "data")
wd.asym.plots <- file.path(wd.asym, "plots")

wd.rt       <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

samples <- readTable(file.path(wd.ngs, "sclc_wgs_n101.list"), header=F, rownames=F, sep="")
n1 <- 50
n0 <- 51

# -----------------------------------------------------------------------------
# Step 1: Load all mutations
# Last Modified: 07/07/19
# -----------------------------------------------------------------------------
#colnames <- c("CHROM",	"POS",	"ID",	"REF",	"ALT",	"QUAL",	"FILTER",	"INFO", "SAMPLE")
#snvs <- toTable("", length(colnames), 0, colnames)
#for (s in 1:length(samples)) {
#   sample <- samples[s]
#   vcf <- read.peiflyne.mutcall.filtered.vcf(file.path(wd.ngs, sample, paste0(sample, "_ANALYSIS"), paste0(sample, "_mutcall_filtered.vcf")), pass=T, rs=F)
#   vcf$SAMPLE <- sample
#   
   ## Seperate SNVs
#   vcf.snv <- subset(vcf,     REF %in% c("A", "T", "C", "G"))
#   vcf.snv <- subset(vcf.snv, ALT %in% c("A", "T", "C", "G"))   
 
   #writeTable(vcf.snv, gzfile(file.path(wd.asym.files, paste0(sample, "_mut.snv.gz"))), colnames=T, rownames=F, sep="\t")
#   snvs <- rbind(snvs, vcf.snv)
#}
#colnames <- c("SAMPLE", "CHROM",	"POS",	"ID",	"REF",	"ALT",	"QUAL",	"FILTER",	"INFO")
#snvs <- snvs[,colnames]
#save(snvs, file=file.path(wd.asym.data, paste0(base, "_mut_snvs.RData")))

load(file=file.path(wd.asym.data, paste0(base, "_mut_snvs.RData")))

# -----------------------------------------------------------------------------
# Step 2: Build S6 table
# Last Modified: 07/07/19
# -----------------------------------------------------------------------------
# > nrow(snvs)
# [1] 5268240

isWithinWindow <- function(bed.gc.chr.rt, POS) {
   bed.gc.chr.rt.start <- subset(bed.gc.chr.rt, START <= POS)
   bed.gc.chr.rt.start.end <- subset(bed.gc.chr.rt.start, END >= POS)
   
   if (nrow(bed.gc.chr.rt.start.end) == 1)
      return(rownames(bed.gc.chr.rt.start.end))
   else
      return("")
}

muts.chr <- toTable(0, 9, 1, c("CHR", "LENGTH", "TOTAL", "C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))
#for (c in 1:22) {
   chr <- chrs[CHR]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
 
   rpkms.chr.rt <- readTable(file.path(wd.rt.data, paste0(base, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, scaledRT=T) 
   rpkms.chr.rt.RT <- setSpline(rpkms.chr.rt, bed.gc.chr, "RT")
   rpkms.chr.rt <- cbind(rpkms.chr.rt[rownames(rpkms.chr.rt.RT),], rpkms.chr.rt.RT[, "SPLINE"])
   colnames(rpkms.chr.rt) <- c("BED", "T", "N", "RT", "SPLINE")
   
   bed.gc.chr.rt <- bed.gc.chr[rownames(rpkms.chr.rt),]
   snvs.chr <- subset(snvs, CHROM == chr)
   muts.chr$LENGTH[1] <- nrow(rpkms.chr.rt) * 1000/1000000
   muts.chr$CHR[1] <- CHR
   
   ## Keep only SNVs on overlapping 1kb windows
   snvs.chr$BED <- mapply(x = 1:nrow(snvs.chr), function(x) isWithinWindow(bed.gc.chr.rt, snvs.chr[x,]$POS))
   
   muts.chr$TOTAL[1]  <- nrow(subset(snvs.chr, BED != ""))
   snvs.chr <- subset(snvs.chr, BED != "")
   snvs.chr <- cbind(snvs.chr, rpkms.chr.rt[snvs.chr$BED, 2:5])
   save(snvs.chr, file=file.path(wd.asym.data, paste0(base, "_mut_snvs_chr", CHR,".RData")))
   
   snvs.chr.s6 <- getTableS6SNV(snvs.chr[, c("SAMPLE", "CHROM", "POS", "REF", "ALT", "QUAL")])
   for (s in 1:6) {
      muts.chr[1, 3+s] <- nrow(snvs.chr.s6[[s]][[1]]) + nrow(snvs.chr.s6[[s]][[2]])
   }
#}
save(muts.chr, file=file.path(wd.asym.data, paste0(base, "_mut_snvs_s6_chr", CHR,".RData")))
