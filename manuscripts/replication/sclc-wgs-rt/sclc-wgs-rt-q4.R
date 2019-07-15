# =============================================================================
# Manuscript   : 
# Chapter      : Chromosome replication timing of the human genome
# Name         : manuscripts/replication/sclc-wgs-rt-m2.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 24/04/19; 05/03/19; 25/02/19; 30/01/18
# =============================================================================
#wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "DifferentialExpression.R", "ReplicationTiming.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.1kb.gc.RData"))

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 30/01/18
# -----------------------------------------------------------------------------
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE  <- "SCLC"
PAIR1 <- "T"
PAIR0 <- "T"
base <- tolower(BASE)

wd.ngs   <- file.path(wd, BASE, "ngs/WGS")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.rt       <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt-q4"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

#samples1 <- readTable(file.path(wd.ngs, "sclc_wgs_n101.txt"), header=T, rownames=T, sep="")   ## M2/M1
samples1 <- readTable(file.path(wd.ngs, "sclc_wgs_n51.txt"), header=T, rownames=T, sep="")     ## Q4/Q1
samples1 <- subset(samples1, M2 == 1)[,1]
#samples0 <- readTable(file.path(wd.ngs, "sclc_wgs_n101.txt"), header=T, rownames=T, sep="")
samples0 <- readTable(file.path(wd.ngs, "sclc_wgs_n51.txt"), header=T, rownames=T, sep="")
samples0 <- subset(samples0, M2 == 0)[,1]
n1 <- length(samples1)
n0 <- length(samples0)

# -----------------------------------------------------------------------------
# Plot RD and RT (see ReplicationTiming.R)
# Last Modified: 14/02/19; 10/01/19; 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   
   rpkms.chr.rt <- readTable(file.path(wd.rt.data, paste0(base, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, scaledRT=T)
   
   ## Plot RT
   main.text <- paste0(BASE, " M2/M1 read depth ratio between tumour (n=", n1, ") and tumour (n=", n0, ") samples")
   file.name <- file.path(wd.rt.plots, paste0("RT_", base, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0))
   plotRT(file.name, main.text, chr, NA, NA, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("M2 tumour", "M1 tumour"), c(adjustcolor.red, adjustcolor.blue), c("M2", "M1"), "png", width=10, peaks=c(), 7.75, 9, 3, 3)
   #plotRT(file.name, paste0(BASE, " M2/M1 read depth ratio"), chr, 44000000, 45000000, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("M2 tumour", "M1 tumour"), c(adjustcolor.red, adjustcolor.blue), c("M2", "M1"), "png", width=5, peaks=c(74393001, 85508001), 7.5, 9.25, 3, 3)
}
# > 9 - 7.75
# [1] 1.25
# > 9.25 - 7.5
# [1] 1.75

# -----------------------------------------------------------------------------
# RD vs RT (RDS and SPR)
# Last Modified: 11/07/19; 31/05/19
# -----------------------------------------------------------------------------
cors <- toTable(0, 7, 22, c("chr", "length", "cor", "cor1", "cor2", "intercept1", "intercept2"))
cors$chr <- 1:22

for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
 
   rpkms.chr.rt <- readTable(file.path(wd.rt.data, paste0(base, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, scaledRT=T) 
   rpkms.chr.rt.T  <- setSpline(rpkms.chr.rt, bed.gc.chr, "T")
   rpkms.chr.rt.N  <- setSpline(rpkms.chr.rt, bed.gc.chr, "N")
   rpkms.chr.rt.RT <- setSpline(rpkms.chr.rt, bed.gc.chr, "RT")
   cors$length[c] <- nrow(rpkms.chr.rt.RT)
   
   cor <- getCor(rpkms.chr.rt.T$SPLINE, rpkms.chr.rt.N$SPLINE, method="spearman")
   cors$cor[c] <- cor
   
   main.text <- c(paste0("SCLC read depths vs. SCLC Q4/Q1 (", "Chr", c, ")"), paste0("rho = ", round0(cor, digits=2), " (Q4 vs. Q1)"))
   xlab.text <- "SCLC Q4/Q1"
   ylab.text <- "SCLC read depth [log2]"
   file.name <- file.path(wd.rt.plots, "chrs", paste0("RD-vs-RT_SCLC-Q4-Q1_chr", c, "_spline_spearman"))
   plotRD2vsRT(rpkms.chr.rt.T$SPLINE, rpkms.chr.rt.N$SPLINE, rpkms.chr.rt.RT$SPLINE, file.name, main.text, ylab.text, xlab.text, c("red", "blue"), c("Q4", "Q1"), method="spearman")
 
   cors$cor1[c] <- getCor(rpkms.chr.rt.T$SPLINE, rpkms.chr.rt.RT$SPLINE, method="spearman")
   cors$cor2[c] <- getCor(rpkms.chr.rt.N$SPLINE, rpkms.chr.rt.RT$SPLINE, method="spearman")
   cors$intercept1[c] <- lm(rpkms.chr.rt.T$SPLINE ~ rpkms.chr.rt.RT$SPLINE)[[1]][1]
   cors$intercept2[c] <- lm(rpkms.chr.rt.N$SPLINE ~ rpkms.chr.rt.RT$SPLINE)[[1]][1]
   
   ## Read depth skew (RDS)
   cors$skew <- (cors$intercept1 - cors$intercept2) / (cors$intercept1 + cors$intercept2)
}
save(cors, file=file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-q4-q1_spline_spearman.RData")))
writeTable(cors, file=file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-q4-q1_spline_spearman.txt")), colnames=T, rownames=F, sep="\t")

#load(file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-q4-q1_spline_spearman.RData")))
#file.name <- file.path(wd.rt.plots, "RD-vs-RT_SCLC_spline_spearman")
#main.text <- paste0("SCLC read depths vs. SCLC M2/M1")
#ymin <- 0.6
#ymax <- 1
#plotRD2vsRTALLREVERSED(cors, file.name, main.text, ymin, ymax, cols=c("red", "blue"), c("M2", "M1"), c=2)

## Read depth skew (RDS)
file.name <- file.path(wd.rt.plots, "RDS_SCLC-Q4-Q1_spline_spearman")
main.text <- c(paste0(BASE, " read depth imbalance"), "Y-axis intercept")
plotRDS(cors, file.name, main.text, ymin=8, ymax=9, cols=c("red", "blue"), c("Q4 tumour", "Q1 tumour"), c(2, 13, 17), digits=3)

## S-phase progression rate (SPR)
file.name <- file.path(wd.rt.plots, "RDS-SPR_SCLC-Q4-Q1_spline_spearman")
main.text <- c(paste0(BASE, " S-phase progression rate"), "SPR = (Q4-Q1)/(Q4+Q1)")
plotSPR(cors, file.name, main.text, c(13, 17), digits=3, unit=5.5)

## SPR vs Read depth correlation (RDC)
file.name <- file.path(wd.rt.plots, "RDS-SPR-vs-RDC_SCLC-Q4-Q1_spline_spearman")
main.text <- c(paste0(BASE, " SPR vs. Read depth correlation"), "SPR = (Q4-Q1)/(Q4+Q1)")
xlab.text <- "Read depth correlation [rho]"
plotSPRRDC(cors, file.name, main.text, c(4, 13, 17, 19, 21, 22), xlab.text, unit=5.5)





# -----------------------------------------------------------------------------
# SCLC M2/M1 vs SCLC Q4/Q1
# Last Modified: 11/07/19; 06/06/19
# -----------------------------------------------------------------------------
cors <- toTable(0, 3, 22, c("chr", "length", "cor"))
cors$chr <- 1:22

for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
 
   rpkms.chr.rt <- readTable(file.path(wd.rt.data, paste0(base, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, scaledRT=T) 
   rpkms.chr.rt.RT <- setSpline(rpkms.chr.rt, bed.gc.chr, "RT")
 
   #rpkms.chr.rt.lcl <-readTable(paste0("/Users/tpyang/Work/uni-koeln/tyang2/LCL/analysis/replication/lcl-wgs-rt/data/lcl_rpkm.corr.gc.d.rt_", chr, "_LCL-LCL_n7-7.txt.gz"), header=T, rownames=T, sep="\t")
   #rpkms.chr.rt.lcl <-readTable(paste0("/Users/tpyang/Work/uni-koeln/tyang2/SCLC/analysis/replication/sclc-wgs-rt/data/sclc_rpkm.corr.gc.d.rt_", chr, "_SCLC-SCLC_n101-92.txt.gz"), header=T, rownames=T, sep="\t")
   #rpkms.chr.rt.lcl <-readTable(paste0("/Users/tpyang/Work/uni-koeln/tyang2/SCLC/analysis/replication/sclc-wgs-rt-q4/data/sclc_rpkm.corr.gc.d.rt_", chr, "_T-T_n25-26.txt.gz"), header=T, rownames=T, sep="\t")
   rpkms.chr.rt.lcl <-readTable(paste0("/Users/tpyang/Work/uni-koeln/tyang2/NBL/analysis/replication/nbl-wgs-rt-m2/data/nbl_rpkm.corr.gc.d.rt_", chr, "_T-T_n28-28.txt.gz"), header=T, rownames=T, sep="\t")
   rpkms.chr.rt.lcl <- setScaledRT(rpkms.chr.rt.lcl, pseudocount=0.01, recaliRT=T, scaledRT=T) 
   rpkms.chr.rt.lcl.RT <- setSpline(rpkms.chr.rt.lcl, bed.gc.chr, "RT")
 
   ## Keep only overlapping 1kb windows
   overlaps <- intersect(rpkms.chr.rt.RT$BED, rpkms.chr.rt.lcl.RT$BED)
   cors$length[c] <- length(overlaps)
   cors$cor[c] <- getCor(rpkms.chr.rt.RT[overlaps,]$SPLINE, rpkms.chr.rt.lcl.RT[overlaps,]$SPLINE, method="spearman")
}
#save(cors, file=file.path(wd.rt.data, paste0("rt-vs-rt_", base, "-m2-m1-vs-sclc-t-n_spline_spearman.RData")))
#save(cors, file=file.path(wd.rt.data, paste0("rt-vs-rt_", base, "-m2-m1-vs-sclc-q4-q1_spline_spearman.RData")))
save(cors, file=file.path(wd.rt.data, paste0("rt-vs-rt_", base, "-m2-m1-vs-nbl-m2-m1_spline_spearman.RData")))

ylab.text <- "Spearman's rho"
xlab.text <- "Chromosome"
#file.name <- file.path(wd.rt.plots, "RT-vs-RT_SCLC-M2-M1-vs-SCLC-T-N_spline_spearman")
#file.name <- file.path(wd.rt.plots, "RT-vs-RT_SCLC-M2-M1-vs-SCLC-Q4-Q1_spline_spearman")
file.name <- file.path(wd.rt.plots, "RT-vs-RT_SCLC-M2-M1-vs-NBL-M2-M1_spline_spearman")
main.text <- paste0("SCLC M2/M1 vs. SCLC Q4/Q1")
ymin <- 0.25
ymax <- 1.05
plotRTvsRTALL(cors, file.name, main.text, ylab.text, xlab.text, ymin, ymax, col="black", c=2, pos=3)





# -----------------------------------------------------------------------------
# Plot RT for individual genes (see ReplicationTiming.R)
# Last Modified: 26/04/19; 14/02/19; 10/01/19; 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
genes <- c("XRCC6", "AHR", "IL33", "IL7", "CDK4", "ERCC6L2")
genes <- c("ERCC6L2")
genes <- c("GNAQ", "EIF2B1", "PTAR1", "FBXL17", "RP11-141C7.3", "PTPRG", "ZNF780B", "NDFIP1", "BPNT1", "ZDHHC2", "CCDC155", "PTPRD")
genes <- c("GNAQ", "VEGFA")
genes <- c("TBC1D32", "BTBD1")
genes <- c("IL6ST")
genes <- c("ERCC6L2", "GNAQ")
genes <- c("MYC", "MYCN", "MYCL", "EGFR", "TP53", "RB1", "IRS2", "CDKN2A", "FHIT", "FGFR1", "WHSC1L1")
genes <- c("COL22A1", "GNAQ", "COL4A3BP")

plotWholeChr <- T
ranges <- c(50000, 500000, 5000000)
for (g in 1:length(genes)) {
   gene  <- getGene(genes[g])
   chr   <- gene$chromosome_name
   start <- gene$start_position
   end   <- gene$end_position
   
   rpkms.chr.rt <- readTable(file.path(wd.rt.data, paste0(base1, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, flipRT=F, scaledRT=T) 
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   bed.gc.chr <- bed.gc.chr[rpkms.chr.rt$BED,]
   
   ## RD & RT 
   main.text <- paste0(BASE1, " M2/M1 read depth ratio between tumour (n=", n1, ") and tumour (n=", n0, ") cells")   

   file.name <- file.path(wd.rt.plots, "genes", paste0("RTD_", base0, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, "_", genes[g]))
   if (plotWholeChr)
      plotRD3(file.name, main.text, chr, NA, NA, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("M2 tumour", "M1 tumour"), c(adjustcolor.gray, adjustcolor.gray), c("M2", "M1"), "png", width=10, peaks=c(start, end), 7.75, 9, 3, 3)
   
   for (r in 1:length(ranges))
      plotRD3(file.name, paste0(BASE, " M2/M1 read depth ratio"), chr, start-ranges[r],	end+ranges[r], rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("M2", "M1"), c(adjustcolor.red, adjustcolor.blue), c("M2", "M1"), "png", width=5, peaks=c(start, end), 7.75, 9, 3, 3)
}
