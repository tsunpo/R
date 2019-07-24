#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
ICGC <- args[1]   ## Cancer type
NUM  <- as.numeric(args[2])

# =============================================================================
# Manuscript   : 
# Chapter      : Chromosome replication timing of the human genome
# Name         : manuscripts/asymmetries/sclc-asym-rt.R
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
# Last Modified: 10/07/18
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"                   ## tyang2@gauss
#wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE <- "SCLC"
PAIR1 <- "T"
PAIR0 <- "T"
base <- tolower(BASE)
icgc <- tolower(ICGC)

wd.anlys <- file.path(wd, "ICGC", "analysis")
wd.asym       <- file.path(wd.anlys, "asymmetries", paste0(icgc, "-asym-rt"))
wd.asym.data  <- file.path(wd.asym, "data")
wd.asym.plots <- file.path(wd.asym, "plots")

## Use NBL SPR
wd.anlys <- file.path(wd, BASE, "analysis")
wd.rt       <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt-m2"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

wd.ngs   <- file.path(wd, BASE, "ngs/WGS")
samples <- readTable(file.path(wd.ngs, "sclc_wgs_n101.list"), header=F, rownames=F, sep="")
n1 <- 50
n0 <- 51

# -----------------------------------------------------------------------------
# Step 2: Build S6 table
# Last Modified: 07/07/19
# -----------------------------------------------------------------------------
###
## After copy files back from cheops
muts <- toTable(0, 9, 0, c("CHR", "LENGTH", "TOTAL", "C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))
for (c in 1:22) {
   load(file=file.path(wd.asym.data, paste0(icgc, "_mut_snvs_s6_chr", c,".RData")))
   muts <- rbind(muts, muts.chr)
}
save(muts, file=file.path(wd.asym.data, paste0(icgc, "_mut_snvs_s6.RData")))

# -----------------------------------------------------------------------------
# SCLC RD vs RD
# Last Modified: 07/07/19
# -----------------------------------------------------------------------------
load(file=file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-m2-m1_spline_spearman.RData")))   ## Load skews
skews <- cors
#load(file=file.path(wd.asym.data, paste0(base, "_mut_snvs_s6.RData")))
#load(file=file.path(wd.asym.data, paste0(base, "_mut_snvs_s6_early.RData")))
#load(file=file.path(wd.asym.data, paste0(base, "_mut_snvs_s6_late.RData")))
time <- ""

s6 <- c("C>A/G>T", "C>G/G>C", "C>T/G>A", "T>A/A>T", "T>C/A>G", "T>G/A>C")
rhos <- toTable(0, 2, 6, c("rho", "P"))
rownames(rhos) <- s6
for (s in 1:6) {
   snvs <- toTable(0, 2, 22, c("chr", "cor"))
   snvs$chr <- 1:22
   idx <- idxs[s]
 
   for (c in 1:22)
      snvs$cor[c] <- muts[c, 3+s] / muts[c, "LENGTH"]
 
   file.name <- file.path(wd.asym.plots, paste0(time, "RDS-vs-MUT-", gsub("/", "_", s6[s]), "_", ICGC))
   main.text <- c(paste0(s6[s], ""), paste0("ICGC ", ICGC, " (n=", NUM, ")"))
   xlab.text <- "SNVs/Mb"
   plotSPRRDCSNV(snvs, cors, file.name, main.text, cs=c(4, 13, 17, 19, 21, 22), xlab.text, unit=4.5, ylab.text="SCLC SPR")
                                                                                                5.5 for SCLC?
   ## TO-DO
   ins <- cbind(snvs, cors[, "skew"])
   colnames(ins) <- c("chr", "cor", "skew")
   
   cor <- cor.test(ins$skew, ins$cor, method="spearman", exact=F)
   rhos$rho[s] <- cor[[4]]
   rhos$P[s]   <- cor[[3]]
}
save(rhos, file=file.path(wd.asym.data, paste0(icgc, "_mut_snvs_rho_s6.RData")))

# -----------------------------------------------------------------------------
# Step 2.1: Build S6 table for EARLY replicated region
# Last Modified: 08/07/19
# -----------------------------------------------------------------------------
muts <- toTable(0, 9, 22, c("CHR", "LENGTH", "TOTAL", "C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))
muts$CHR <- 1:22
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
 
   rpkms.chr.rt <- readTable(file.path(wd.rt.data, paste0(base, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, scaledRT=T) 
   rpkms.chr.rt.RT <- setSpline(rpkms.chr.rt, bed.gc.chr, "RT")
 
   rpkms.chr.rt.RT <- subset(rpkms.chr.rt.RT, SPLINE > 0)
   muts$LENGTH[c] <- nrow(rpkms.chr.rt.RT) * 1000/1000000
 
   ##
   load(file=file.path(wd.asym.data, paste0(icgc, "_mut_snvs_chr", c,".RData")))   ## Load snvs.chr
   snvs.chr <- subset(snvs.chr, SPLINE > 0)
   muts$TOTAL[c] <- nrow(snvs.chr)
 
   snvs.chr.s6 <- getTableS6SNV(snvs.chr[, c("SAMPLE", "CHROM", "POS", "REF", "ALT", "QUAL")])
   for (s in 1:6) {
      muts[c, 3+s] <- nrow(snvs.chr.s6[[s]][[1]]) + nrow(snvs.chr.s6[[s]][[2]])
   }
}
save(muts, file=file.path(wd.asym.data, paste0(icgc, "_mut_snvs_s6_early.RData")))

time <- "EARLY_"
s6 <- c("C>A/G>T", "C>G/G>C", "C>T/G>A", "T>A/A>T", "T>C/A>G", "T>G/A>C")
for (s in 1:6) {
   snvs <- toTable(0, 2, 22, c("chr", "cor"))
   snvs$chr <- 1:22
   idx <- idxs[s]
 
   for (c in 1:22)
      snvs$cor[c] <- muts[c, 3+s] / muts[c, "LENGTH"]
 
   file.name <- file.path(wd.asym.plots, paste0(time, "RDS-vs-MUT-", gsub("/", "_", s6[s]), "_", ICGC))
   main.text <- c(paste0(s6[s], ""), paste0("ICGC ", ICGC, " (n=", NUM, ")"))
   xlab.text <- "SNVs/Mb"
   plotSPRRDCSNV(snvs, cors, file.name, main.text, cs=c(4, 13, 17, 19, 21, 22), xlab.text, unit=4.5, ylab.text="SCLC SPR")
}

# -----------------------------------------------------------------------------
# Step 2.1: Build S6 table for LATE replicated region
# Last Modified: 08/07/19
# -----------------------------------------------------------------------------
muts <- toTable(0, 9, 22, c("CHR", "LENGTH", "TOTAL", "C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))
muts$CHR <- 1:22
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
 
   rpkms.chr.rt <- readTable(file.path(wd.rt.data, paste0(base, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, scaledRT=T) 
   rpkms.chr.rt.RT <- setSpline(rpkms.chr.rt, bed.gc.chr, "RT")
 
   rpkms.chr.rt.RT <- subset(rpkms.chr.rt.RT, SPLINE < 0)
   muts$LENGTH[c] <- nrow(rpkms.chr.rt.RT) * 1000/1000000
 
   ##
   load(file=file.path(wd.asym.data, paste0(icgc, "_mut_snvs_chr", c,".RData")))   ## Load snvs.chr
   snvs.chr <- subset(snvs.chr, SPLINE < 0)
   muts$TOTAL[c] <- nrow(snvs.chr)
 
   snvs.chr.s6 <- getTableS6SNV(snvs.chr[, c("SAMPLE", "CHROM", "POS", "REF", "ALT", "QUAL")])
   for (s in 1:6) {
      muts[c, 3+s] <- nrow(snvs.chr.s6[[s]][[1]]) + nrow(snvs.chr.s6[[s]][[2]])
   }
}
save(muts, file=file.path(wd.asym.data, paste0(icgc, "_mut_snvs_s6_late.RData")))

time <- "LATE_"
s6 <- c("C>A/G>T", "C>G/G>C", "C>T/G>A", "T>A/A>T", "T>C/A>G", "T>G/A>C")
for (s in 1:6) {
   snvs <- toTable(0, 2, 22, c("chr", "cor"))
   snvs$chr <- 1:22
   idx <- idxs[s]
 
   for (c in 1:22)
      snvs$cor[c] <- muts[c, 3+s] / muts[c, "LENGTH"]
 
   file.name <- file.path(wd.asym.plots, paste0(time, "RDS-vs-MUT-", gsub("/", "_", s6[s]), "_", ICGC))
   main.text <- c(paste0(s6[s], ""), paste0("ICGC ", ICGC, " (n=", NUM, ")"))
   xlab.text <- "SNVs/Mb"
   plotSPRRDCSNV(snvs, cors, file.name, main.text, cs=c(4, 13, 17, 19, 21, 22), xlab.text, unit=4.5, ylab.text="SCLC SPR")
}
