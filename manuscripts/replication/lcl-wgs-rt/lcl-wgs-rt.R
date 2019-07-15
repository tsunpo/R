# =============================================================================
# Manuscript   : 
# Chapter      : Chromosome replication timing of the human genome
# Name         : manuscripts/replication/lcl-wgs-rt.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 16/06/19; 25/02/19; 16/05/18
# =============================================================================
#wd.src <- "/projects/cangen/tyang2/dev/R"            ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"                 ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"                  ## tpyang@localhost

wd.src.handbook <- file.path(wd.src, "handbook-of")   ## Required handbooks/libraries for the manuscript
handbooks <- c("Commons.R", "DifferentialExpression.R", "ReplicationTiming.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.handbook, x))))

wd.src.guide <- file.path(wd.src, "guide-to-the")     ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.guide, "hg19.RData"))
load(file.path(wd.src.guide, "hg19.1kb.gc.RData"))

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 30/01/18
# -----------------------------------------------------------------------------
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE  <- "LCL"
PAIR1 <- "S"
PAIR0 <- "G1"
base  <- tolower(BASE)

wd.ngs   <- file.path(wd, BASE, "ngs/WGS")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.rt       <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

samples1 <- readTable(file.path(wd.ngs, "lcl_wgs_n7.list"), header=F, rownames=F, sep="")
samples0 <- readTable(file.path(wd.ngs, "lcl_wgs_n7.list"), header=F, rownames=F, sep="")
n1 <- length(samples1)
n0 <- length(samples0)

# -----------------------------------------------------------------------------
# Plot RD and RT (see ReplicationTiming.R)
# Last Modified: 28/05/19; 14/02/19; 10/01/19; 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
 
   rpkms.chr.rt <- readTable(file.path(wd.rt.data, paste0(base, "_rpkm.corr.gc.d.rt_", chr, "_", BASE, "-", BASE, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, scaledRT=T)

   ## Plot RT
   main.text <- paste0(BASE, " S/G1 read depth ratio between S phase (n=", n1, ") and G1 phase (n=", n0, ") cells")  
   file.name <- file.path(wd.rt.plots, paste0("RT_", base, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, "_9.25"))   
   plotRT(file.name, main.text, chr, NA, NA, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("S phase", "G1 phase"), c(adjustcolor.red, adjustcolor.blue), c("S", "G1"), "png", width=10, peaks=c(), 7.25, 9.5, 3, 3)
}
# > 9.5 - 7.25
# [1] 2.25

#spline <- smooth.spline(x=bed.gc.chr$START, y=rpkms.chr.rt$RT)
#idx.1 <- which(spline$x > 73353001)
#idx.2 <- which(spline$x < 75353001)
#idx <- intersect(idx.1, idx.2)

#which(spline$y == max(spline$y[idx]))
#spline$x[72636]

# -----------------------------------------------------------------------------
# Plot RT for individual genes (see ReplicationTiming.R)
# Last Modified: 26/04/19; 14/02/19; 10/01/19; 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
genes <- c("MYC", "MYCN", "MYCL1")
genes <- c("NRG1", "DDX55", "KNTC1")
genes <- c("CCAR1", "SMARCD1")

plotWholeChr <- T
ranges <- c(50000, 500000, 5000000)
for (g in 1:length(genes)) {
   gene  <- getGene(genes[g])
   chr   <- gene$chromosome_name
   start <- gene$start_position
   end   <- gene$end_position
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   
   rpkms.chr.rt <- readTable(file.path(wd.rt.data, paste0(base, "_rpkm.corr.gc.d.rt_", chr, "_", BASE, "-", BASE, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, scaledRT=T) 

   ## RD & RT 
   main.text <- paste0(BASE, " S/G1 read depth ratio between S phase (n=", n1, ") and G1 phase (n=", n0, ") cells")   
   file.name <- file.path(wd.rt.plots, "genes", paste0("RT_", base, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, "_", genes[g]))
   if (plotWholeChr)
      plotRT(file.name, main.text, chr, NA, NA, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("S phase", "G1 phase"), c(adjustcolor.gray, adjustcolor.gray), c("S", "G1"), "png", width=10, peaks=c(start, end), 7.75, 9, 3, 3)

   for (r in 1:length(ranges))
      plotRT(file.name, paste0(BASE, " S/G1 read depth ratio"), chr, start-ranges[r],	end+ranges[r], rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("S", "G1"), c(adjustcolor.red, adjustcolor.blue), c("S", "G1"), "png", width=5, peaks=c(start, end), 7.75, 9, 3, 3)
}

# -----------------------------------------------------------------------------
# RD vs RT (RDS and SPR)
# Last Modified: 11/07/19; 27/05/19
# -----------------------------------------------------------------------------
cors <- toTable(0, 7, 22, c("chr", "length", "cor", "cor1", "cor2", "intercept1", "intercept2"))
cors$chr <- 1:22

for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
 
   rpkms.chr.rt <- readTable(file.path(wd.rt.data, paste0(base, "_rpkm.corr.gc.d.rt_", chr, "_", BASE, "-", BASE, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, scaledRT=T)
   rpkms.chr.rt.T  <- setSpline(rpkms.chr.rt, bed.gc.chr, "T")
   rpkms.chr.rt.N  <- setSpline(rpkms.chr.rt, bed.gc.chr, "N")
   rpkms.chr.rt.RT <- setSpline(rpkms.chr.rt, bed.gc.chr, "RT")
   cors$length[c] <- nrow(rpkms.chr.rt.RT)
   
   cor <- getCor(rpkms.chr.rt.T$SPLINE, rpkms.chr.rt.N$SPLINE, method="spearman")
   cors$cor[c] <- cor
   
   main.text <- c(paste0("LCL read depth correlation (", "Chr", c, ")"), paste0("rho = ", round0(cor, digits=2), " (S vs. G1)"))
   xlab.text <- "LCL S/G1"
   ylab.text <- "LCL read depth [log2]"
   file.name <- file.path(wd.rt.plots, "chrs", paste0("RD-vs-RT_LCL-S-G1_chr", c, "_spline_spearman"))
   plotRD2vsRT(rpkms.chr.rt.T$SPLINE, rpkms.chr.rt.N$SPLINE, rpkms.chr.rt.RT$SPLINE, file.name, main.text, ylab.text, xlab.text, c("red", "blue"), c("S", "G1"), method="spearman")
   
   cors$cor1[c] <- getCor(rpkms.chr.rt.T$SPLINE, rpkms.chr.rt.RT$SPLINE, method="spearman")
   cors$cor2[c] <- getCor(rpkms.chr.rt.N$SPLINE, rpkms.chr.rt.RT$SPLINE, method="spearman")
   cors$intercept1[c] <- lm(rpkms.chr.rt.T$SPLINE ~ rpkms.chr.rt.RT$SPLINE)[[1]][1]
   cors$intercept2[c] <- lm(rpkms.chr.rt.N$SPLINE ~ rpkms.chr.rt.RT$SPLINE)[[1]][1]
   
   ## Read depth skew (RDS)
   cors$skew <- (cors$intercept1 - cors$intercept2) / (cors$intercept1 + cors$intercept2)
}
save(cors, file=file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-s-g1_spline_spearman.RData")))
writeTable(cors, file=file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-s-g1_spline_spearman.txt")), colnames=T, rownames=F, sep="\t")

#load(file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-s-g1_spline_spearman.RData")))
file.name <- file.path(wd.rt.plots, "RD-vs-RT_LCL_spline_spearman")
main.text <- paste0("LCL read depth vs. LCL S/G1")
ymin <- -0.8
ymax <- 0.8
plotRD2vsRTALL(cors, file.name, main.text, ymin, ymax, cols=c("red", "blue"), c("S", "G1"), c=2)

## Read depth skew (RDS)
file.name <- file.path(wd.rt.plots, "RDS_LCL-S-G1_spline_spearman")
main.text <- c(paste0(BASE, " read depth imbalance"), "Y-axis intercept")
plotRDS(cors, file.name, main.text, ymin=8, ymax=9, cols=c("red", "blue"), c("S phase", "G1 phase"), c(2, 13, 17), digits=3)

## S-phase progression rate (SPR)
file.name <- file.path(wd.rt.plots, "RDS-SPR_LCL-S-G1_spline_spearman")
main.text <- c(paste0(BASE, " S-phase progression rate"), "SPR = (S-G1)/(S+G1)")
plotSPR(cors, file.name, main.text, c(13, 17), digits=3, unit=5)

## SPR vs Read depth correlation (RDC)
file.name <- file.path(wd.rt.plots, "RDS-SPR-RDC_LCL-S-G1_spline_spearman")
main.text <- c(paste0(BASE, " SPR vs. Read depth correlation"), "SPR = (S-G1)/(S+G1)")
xlab.text <- "Read depth correlation [rho]"
plotSPRRDC(cors, file.name, main.text, c(4, 13, 17, 19, 21, 22), xlab.text, unit=5)
