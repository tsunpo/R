# =============================================================================
# Manuscript   : 
# Chapter      : Reconstruction of replication timing profile in tumour cells
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
   file.name <- file.path(wd.rt.plots, paste0("RT_", base, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0))   
   plotRT(file.name, main.text, chr, NA, NA, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("S phase", "G1 phase"), c(adjustcolor.gray, adjustcolor.gray), c("S", "G1"), "png", width=10, peaks=c(), 7.5, 9.25, 3, 3)
   
   ## chr2
   #plotRT(file.name, main.text, chr, NA, NA, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("S phase", "G1 phase"), c(adjustcolor.gray, adjustcolor.gray), c("S", "G1"), "png", width=10, peaks=c(74393001, 85508001), 7.5, 9.25, 3, 3)
   #plotRT(file.name, paste0(BASE, " S/G1 read depth ratio"), chr, 71500000, 90500000, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("S phase", "G1 phase"), c(adjustcolor.red, adjustcolor.blue), c("S", "G1"), "png", width=5, peaks=c(74393001, 85508001), 7.5, 9.25, 3, 3)   ## T/N: 74353001, 85951001
   ## chr4
   #plotRT(file.name, main.text, chr, NA, NA, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("S phase", "G1 phase"), c(adjustcolor.gray, adjustcolor.gray), c("S", "G1"), "png", width=10, peaks=c(128986001, 140285001), 7.5, 9.25, 3, 3)
   #plotRT(file.name, paste0(BASE, " S/G1 read depth ratio"), chr, 125000000, 144000000, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("S phase", "G1 phase"), c(adjustcolor.red, adjustcolor.blue), c("S", "G1"), "png", width=5, peaks=c(128986001, 140285001), 7.5, 9.25, 3, 3)
   ## chr22
   #plotRT(file.name, main.text, chr, NA, NA, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("S phase", "G1 phase"), c(adjustcolor.gray, adjustcolor.gray), c("S", "G1"), "png", width=10, peaks=c(), 7.25, 9, 3, 3)
}

spline <- smooth.spline(x=bed.gc.chr$START, y=rpkms.chr.rt$RT)
idx.1 <- which(spline$x > 73353001)
idx.2 <- which(spline$x < 75353001)
idx <- intersect(idx.1, idx.2)

which(spline$y == max(spline$y[idx]))
spline$x[72636]

# -----------------------------------------------------------------------------
# Plot RT for individual genes (see ReplicationTiming.R)
# Last Modified: 26/04/19; 14/02/19; 10/01/19; 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
genes <- c("MYC", "MYCN", "MYCL1")
genes <- c("NRG1", "DDX55", "KNTC1")

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
# LCL RD vs RT
# Last Modified: 27/05/19
# -----------------------------------------------------------------------------
skews <- toTable(0, 6, 22, c("chr", "length", "cor1", "cor2", "intercept1", "intercept2"))
skews$chr <- 1:22

for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
 
   rpkms.chr.rt <- readTable(file.path(wd.rt.data, paste0(base, "_rpkm.corr.gc.d.rt_", chr, "_", BASE, "-", BASE, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, scaledRT=T)
   rpkms.chr.rt.T  <- setSpline(rpkms.chr.rt, bed.gc.chr, "T")
   rpkms.chr.rt.N  <- setSpline(rpkms.chr.rt, bed.gc.chr, "N")
   rpkms.chr.rt.RT <- setSpline(rpkms.chr.rt, bed.gc.chr, "RT")
   skews$length[c] <- nrow(rpkms.chr.rt.RT)
   
   main.text <- paste0("LCL read depths vs. LCL S/G1 (", "Chr", c, ")")
   xlab.text <- "LCL S/G1"
   ylab.text <- "LCL read depth [log2]"
   file.name <- file.path(wd.rt.plots, paste0("plot_RD-vs-RT_LCL-vs-LCL_chr", c, "_spline_spearman"))
   #xmin <- min(rpkms.chr.rt.RT$SPLINE)
   #xmax <- max(rpkms.chr.rt.RT$SPLINE)
   #ymin <- min(c(rpkms.chr.rt.T$SPLINE, rpkms.chr.rt.N$SPLINE))
   #ymax <- max(c(rpkms.chr.rt.T$SPLINE, rpkms.chr.rt.N$SPLINE))
   plotRD2vsRT(rpkms.chr.rt.T$SPLINE, rpkms.chr.rt.N$SPLINE, rpkms.chr.rt.RT$SPLINE, file.name, main.text, ylab.text, xlab.text, c("red", "blue"), c("S", "G1"), method="spearman")
   
   skews$cor1[c] <- getCor(rpkms.chr.rt.T$SPLINE, rpkms.chr.rt.RT$SPLINE, method="spearman")
   skews$cor2[c] <- getCor(rpkms.chr.rt.N$SPLINE, rpkms.chr.rt.RT$SPLINE, method="spearman")
   skews$intercept1[c] <- lm(rpkms.chr.rt.T$SPLINE ~ rpkms.chr.rt.RT$SPLINE)[[1]][1]
   skews$intercept2[c] <- lm(rpkms.chr.rt.N$SPLINE ~ rpkms.chr.rt.RT$SPLINE)[[1]][1]    
}
save(skews, file=file.path(wd.rt.data, paste0("rds-vs-rt_", base, "-vs-lcl_spline_spearman.RData")))

ylab.text <- "Spearman's rho"
xlab.text <- "Chromosome"
file.name <- file.path(wd.rt.plots, "plot_RD-vs-RT_LCL-vs-LCL_spline_spearman")
main.text <- paste0("LCL read depths vs. LCL S/G1")
ymin <- -0.8
ymax <- 0.8
plotRD2vsRTALL(skews, file.name, main.text, ylab.text, xlab.text, ymin, ymax, cols=c("red", "blue"), c=2)

file.name <- file.path(wd.rt.plots, "plot_RDS-vs-RT_LCL_spline_spearman")
main.text <- c("LCL read depths vs. LCL S/G1", "RDS = (S-G1)/(S+G1)")
ymax <- 0.015
plotRDS(skews, file.name, main.text, ymax, c(4, 17), digits=2)

# -----------------------------------------------------------------------------
# LCL RD vs RD
# Last Modified: 03/06/19
# -----------------------------------------------------------------------------
cors <- toTable(0, 2, 22, c("chr", "cor"))
cors$chr <- 1:22

for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
 
   rpkms.chr.rt <- readTable(file.path(wd.rt.data, paste0(base, "_rpkm.corr.gc.d.rt_", chr, "_", BASE, "-", BASE, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, scaledRT=T)
   rpkms.chr.rt.T  <- setSpline(rpkms.chr.rt, bed.gc.chr, "T")
   rpkms.chr.rt.N  <- setSpline(rpkms.chr.rt, bed.gc.chr, "N")
 
   cors$cor[c] <- getCor(rpkms.chr.rt.T$SPLINE, rpkms.chr.rt.N$SPLINE, method="spearman")
}
save(cors, file=file.path(wd.rt.data, paste0("rds-vs-rd_", base, "_spline_spearman.RData")))

#load(file.path(wd.rt.data, paste0("rds-vs-rt_", base, "-vs-lcl_spline_spearman.RData")))
file.name <- file.path(wd.rt.plots, "plot_RDS-vs-RD_LCL_spline_spearman")
main.text <- c(paste0("LCL S vs. LCL G1"), "Linear regression")
ymax <- 0.015
xlim <- c(0.15, 0.9425)
plotRDCorsRDS(cors, skews, file.name, main.text, ymax, xlim, c(1, 4, 13, 17, 19, 22))
