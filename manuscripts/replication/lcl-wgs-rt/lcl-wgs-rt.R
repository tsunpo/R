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

wd.src.ref <- file.path(wd.src, "guide-to-the")     ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.1kb.gc.RData"))
load(file.path(wd.src.ref, "hg19.rt.lcl.koren.woodfine.RData"))

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 30/01/18
# -----------------------------------------------------------------------------
#wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE  <- "LCL"
BASE1 <- "T"
BASE0 <- "N"
PAIR1 <- "S"
PAIR0 <- "G1"
base  <- tolower(BASE)
method <- "rpkm"

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
# Last Modified: 07/08/19; 28/05/19; 14/02/19; 10/01/19; 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
nrds <- toTable(NA, 4, 0, c("BED", "T", "N", "RT"))
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   nrds.chr <- readTable(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d.rt_", chr, "_", BASE1, "-", BASE0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   #nrds.chr <- setScaledRT(nrds.chr, pseudocount=0, recaliRT=T, scaledRT=F)
   
   nrds <- rbind(nrds, nrds.chr)
}
nrds$RT <- scale(nrds$RT)
save(nrds, file=file.path(wd.rt.data, paste0("nrds_", base, "-s-g1_", method, ".RData")))
# > nrow(nrds)
# [1] 2582940 - 22
nrds.lcl <- nrds

#load(file=file.path(wd.rt.data, paste0("nrds_", base, "-s-g1_", method, ".RData")))
ymax <- 0.6
ymin <- 0.14
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]
   lcl.rt.chr <- subset(lcl.rt, CHR == chr)   ## Koren 2012
   
   ## Plot RT
   main.text <- paste0(BASE, " S/G1 read depth ratio between S phase (n=", n1, ") and G1 phase (n=", n0, ") cells")  
   file.name <- file.path(wd.rt.plots, paste0("RT_", base, "_", method, ".d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ""))   
   plotRT(file.name, main.text, chr, NA, NA, nrds.chr, bed.gc.chr, c("red", "blue"), c("S phase", "G1 phase"), c("lightcoral", "lightskyblue3"), c("S", "G1"), "png", width=10, peaks=c(), ylim=c(ymin, ymax), lcl.rt.chr)
}

# -----------------------------------------------------------------------------
# RD vs RT (RDS and SPR)
# Last Modified: 11/07/19; 27/05/19
# -----------------------------------------------------------------------------
cors <- toTable(0, 8, 22, c("chr", "length", "cor", "cor1", "cor2", "mean", "intercept1", "intercept2"))
cors$chr <- 1:22
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]
   #nrds.chr <- readTable(file.path(wd.rt.data, paste0(base, "_rpkb.gc.cn.d.rt_", chr, "_", BASE1, "-", BASE0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   #nrds.chr <- setScaledRT(nrds.chr, pseudocount=0, recaliRT=T, scaledRT=T)
   
   nrds.chr.T  <- setSpline(nrds.chr, bed.gc.chr, "T")
   nrds.chr.N  <- setSpline(nrds.chr, bed.gc.chr, "N")
   nrds.chr.RT <- setSpline(nrds.chr, bed.gc.chr, "RT")
   #nrds.chr.RT$SPLINE <- scale(nrds.chr.RT$SPLINE)
   cors$length[c] <- nrow(nrds.chr.RT)
   cors$mean[c]   <- mean(nrds.chr.RT$SPLINE)
   
   cor <- getCor(nrds.chr.T$SPLINE, nrds.chr.N$SPLINE, method="spearman")
   cors$cor[c] <- cor
 
   main.text <- c(paste0("LCL read depths correlation (", "Chr", c, ")"), paste0("rho = ", round0(cor, digits=2), " (S vs. G1)"))
   xlab.text <- "LCL S/G1"
   ylab.text <- "LCL read depth [RPKM]"
   file.name <- file.path(wd.rt.plots, "chrs", paste0("RD-vs-RT_LCL-S-G1_chr", c, "_spline_spearman"))
   plotRD2vsRT(nrds.chr.T$SPLINE, nrds.chr.N$SPLINE, nrds.chr.RT$SPLINE, file.name, main.text, ylab.text, xlab.text, c("red", "blue"), c("S", "G1"), method="spearman", intercept=F)
 
   cors$cor1[c] <- getCor(nrds.chr.T$SPLINE, nrds.chr.RT$SPLINE, method="spearman")
   cors$cor2[c] <- getCor(nrds.chr.N$SPLINE, nrds.chr.RT$SPLINE, method="spearman")
   cors$intercept1[c] <- lm(nrds.chr.T$SPLINE ~ nrds.chr.RT$SPLINE)[[1]][1]
   cors$intercept2[c] <- lm(nrds.chr.N$SPLINE ~ nrds.chr.RT$SPLINE)[[1]][1]
 
   ## Read depth skew (RDS)
   cors$skew <- (cors$intercept1 - cors$intercept2) / (cors$intercept1 + cors$intercept2)
}
save(cors, file=file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-s-g1_spline_spearman.RData")))
writeTable(cors, file=file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-s-g1_spline_spearman.txt")), colnames=T, rownames=F, sep="\t")

#load(file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-s-g1_spline_spearman.RData")))
file.name <- file.path(wd.rt.plots, "RD-vs-RT_LCL_spline_spearman")
main.text <- paste0("LCL read depths vs. LCL S/G1")
ymin <- -0.8
ymax <- 0.8
plotRD2vsRTALL(cors, file.name, main.text, ymin, ymax, cols=c("red", "blue"), c("S", "G1"), c=2)

## S-phase progression rate (SPR)
file.name <- file.path(wd.rt.plots, "SPR_LCL-S-G1_spline_spearman")
main.text <- c(paste0(BASE, " S-phase progression rate (SPR)"), "")
plotSPR(cors, file.name, main.text, c(13, 17), digits=3, unit=5, ylab.text="Mean S/G1 ratio")

## SPR vs Read depth correlation (RDC)
file.name <- file.path(wd.rt.plots, "SPR-RDC_LCL-S-G1_spline_spearman")
main.text <- c(paste0(BASE, " SPR vs. Read depths correlation"), "")
xlab.text <- "S vs. G1 [rho]"
plotSPRRDC(cors$mean, cors$cor, file.name, main.text, c(4, 13, 17, 19, 22), xlab.text, unit=5, ylab.text="Mean S/G1 ratio")

## SPR vs Woodfine 2004
file.name <- file.path(wd.rt.plots, "SPR-Woodfine_LCL-S-G1_spline_spearman")
main.text <- c(paste0(BASE, " SPR vs. Woodfine 2004"), "")
xlab.text <- "Mean replication timing ratio"
plotSPRRDC(cors$mean, lcl.mean$Mean, file.name, main.text, c(4, 13, 17, 19, 22), xlab.text, unit=5, ylab.text="Mean S/G1 ratio")

## S read depth correlation (RDC)
file.name <- file.path(wd.rt.plots, "SPR-S_LCL-S-G1_spline_spearman")
main.text <- c("S-phase read correlation", "")
xlab.text <- "S vs. S/G1 [rho]"
plotSPRRDC(cors$mean, cors$cor1, file.name, main.text, c(4, 13, 17, 19, 22), xlab.text, unit=5, ylab.text="Mean S/G1 ratio")

file.name <- file.path(wd.rt.plots, "SPR-G1_LCL-S-G1_spline_spearman")
main.text <- c("G1-phase read correlation", "")
xlab.text <- "G1 vs. S/G1 [rho]"
plotSPRRDC(cors$mean, cors$cor2, file.name, main.text, c(4, 13, 17, 19, 22), xlab.text, unit=5, ylab.text="Mean S/G1 ratio")

###
## S read depth correlation (RDC)
file.name <- file.path(wd.rt.plots, "SG1-S_LCL-S-G1_spline_spearman")
main.text <- c("S-phase read correlation", "")
xlab.text <- "S vs. G1 [rho]"
plotSRDC(cors$cor1, cors$cor, file.name, main.text, c(4, 13, 17, 19, 22), xlab.text, unit=5, ylab.text="S vs. S/G1 [rho]", "bottomleft")

file.name <- file.path(wd.rt.plots, "SG1-G1_LCL-S-G1_spline_spearman")
main.text <- c("G1-phase read correlation", "")
xlab.text <- "S vs. G1 [rho]"
plotSRDC(cors$cor2, cors$cor, file.name, main.text, c(4, 13, 17, 19, 22), xlab.text, unit=5, ylab.text="G1 vs. S/G1 [rho]", "topleft")










# -----------------------------------------------------------------------------
# Plot RD and RT (see ReplicationTiming.R)
# Last Modified: 04/08/19; 28/05/19; 14/02/19; 10/01/19; 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
rts <- toTable(NA, 4, 0, c("BED", "T", "N", "RT"))
for (c in 1:22) {
   chr <- chrs[c]
 
   rts.chr <- readTable(file.path(wd.rt.data, paste0(base, "_rpkb.gc.cn.d.rt_", chr, "_", BASE1, "-", BASE0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rts.chr <- setScaledRT(rts.chr, pseudocount=0.01, recaliRT=T, scaledRT=F)
   rts <- rbind(rts, rts.chr)
}
rts$RT <- scale(rts$RT)

for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
 
   nrds.chr.rt <- readTable(file.path(wd.rt.data, paste0(base, "_rpkb.gc.cn.d.rt_", chr, "_", BASE1, "-", BASE0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   nrds.chr.rt <- setScaledRT(nrds.chr.rt, pseudocount=0, recaliRT=T, scaledRT=T)

   ## Plot RT
   main.text <- paste0(BASE, " S/G1 read depth ratio between S phase (n=", n1, ") and G1 phase (n=", n0, ") cells")  
   file.name <- file.path(wd.rt.plots, paste0("RT_", base, "_rpkb.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, "_scale"))   
   plotRT(file.name, main.text, chr, NA, NA, nrds.chr.rt, bed.gc.chr, c("red", "blue"), c("S phase", "G1 phase"), c(adjustcolor.red, adjustcolor.blue), c("S", "G1"), "png", width=10, peaks=c(), NA, NA, 3, 3, isKoren=T)
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
   
   rpkms.chr.rt <- readTable(file.path(wd.rt.data, paste0(base, "_rpkb.gc.cn.d.rt_", chr, "_", BASE, "-", BASE, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
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
 
   nrds.chr.rt <- readTable(file.path(wd.rt.data, paste0(base, "_rpkb.gc.cn.d.rt_", chr, "_", BASE1, "-", BASE0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   nrds.chr.rt <- setScaledRT(nrds.chr.rt, pseudocount=0, recaliRT=T, scaledRT=T)
   nrds.chr.rt.T  <- setSpline(nrds.chr.rt, bed.gc.chr, "T")
   nrds.chr.rt.N  <- setSpline(nrds.chr.rt, bed.gc.chr, "N")
   nrds.chr.rt.RT <- setSpline(nrds.chr.rt, bed.gc.chr, "RT")
   cors$length[c] <- nrow(nrds.chr.rt.RT)
   
   cor <- getCor(nrds.chr.rt.T$SPLINE, nrds.chr.rt.N$SPLINE, method="spearman")
   cors$cor[c] <- cor
   
   main.text <- c(paste0("LCL read depth correlation (", "Chr", c, ")"), paste0("rho = ", round0(cor, digits=2), " (S vs. G1)"))
   xlab.text <- "LCL S/G1"
   ylab.text <- "LCL read depth [log2]"
   file.name <- file.path(wd.rt.plots, "chrs", paste0("RD-vs-RT_LCL-S-G1_chr", c, "_spline_spearman"))
   plotRD2vsRT(nrds.chr.rt.T$SPLINE, nrds.chr.rt.N$SPLINE, nrds.chr.rt.RT$SPLINE, file.name, main.text, ylab.text, xlab.text, c("red", "blue"), c("S", "G1"), method="spearman")
   
   cors$cor1[c] <- getCor(nrds.chr.rt.T$SPLINE, nrds.chr.rt.RT$SPLINE, method="spearman")
   cors$cor2[c] <- getCor(nrds.chr.rt.N$SPLINE, nrds.chr.rt.RT$SPLINE, method="spearman")
   cors$intercept1[c] <- lm(nrds.chr.rt.T$SPLINE ~ nrds.chr.rt.RT$SPLINE)[[1]][1]
   cors$intercept2[c] <- lm(nrds.chr.rt.N$SPLINE ~ nrds.chr.rt.RT$SPLINE)[[1]][1]
   
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
plotRDS(cors, file.name, main.text, ymin=8, ymax=9, cols=c("red", "blue"), c("S-phase intercept (S)", "G1-phase intercept (G1)"), c(2, 13, 17), digits=3)

## S-phase progression rate (SPR)
file.name <- file.path(wd.rt.plots, "RDS-SPR_LCL-S-G1_spline_spearman")
main.text <- c(paste0(BASE, " S-phase progression rate"), "SPR = (S-G1)/(S+G1)")
plotSPR(cors, file.name, main.text, c(13, 17), digits=3, unit=5)

## SPR vs Read depth correlation (RDC)
file.name <- file.path(wd.rt.plots, "RDS-SPR-RDC_LCL-S-G1_spline_spearman")
main.text <- c(paste0(BASE, " SPR vs. Read depth correlation"), "SPR = (S-G1)/(S+G1)")
xlab.text <- "Read depth correlation [rho]"
plotSPRRDC(cors, file.name, main.text, c(4, 13, 17, 19, 21, 22), xlab.text, unit=5)

# -----------------------------------------------------------------------------
# Test
# Last Modified: 01/08/19
# -----------------------------------------------------------------------------
timings <- toTable(0, 2, 22, c("chr", "RT"))
timings$chr <- 1:22
timings.log2 <- toTable(0, 2, 22, c("chr", "RT"))
timings.log2$chr <- 1:22
for (c in 22:1) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
 
   rpkms.chr.rt <- readTable(file.path(wd.rt.data, paste0(base, "_nrd.cn.d.rt_", chr, "_", BASE1, "-", BASE0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   timings$RT[c] <- sum(rpkms.chr.rt$T)/sum(rpkms.chr.rt$N)
   
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, scaledRT=F)
   timings.log2$RT[c] <- sum(rpkms.chr.rt$T)/sum(rpkms.chr.rt$N)
}
timings <- timings[order(timings$RT, decreasing=T),]
timings.log2 <- timings.log2[order(timings.log2$RT, decreasing=T),]
