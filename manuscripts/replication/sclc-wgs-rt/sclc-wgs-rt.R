# =============================================================================
# Manuscript   : 
# Chapter      : Reconstruction of replication timing profile in tumour cells
# Name         : manuscripts/replication/sclc-wgs-rt.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 05/03/19; 25/02/19; 30/01/18
# =============================================================================
#wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "ReplicationTiming.R")
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
PAIR0 <- "N"
base  <- tolower(BASE)

wd.ngs   <- file.path(wd, BASE, "ngs/WGS")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.rt       <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

samples1 <- readTable(file.path(wd.ngs, "sclc_wgs_n101.list"), header=F, rownames=F, sep="")
samples0 <- readTable(file.path(wd.ngs, "sclc_wgs_n92_N.list"), header=F, rownames=F, sep="")
#samples0 <- readTable(file.path(wd.ngs, "sclc_wgs_n9_B.list"), header=F, rownames=F, sep="")
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
   main.text <- paste0(BASE, " T/N read depth ratio between tumour (n=", n1, ") and normal (n=", n0, ") samples")
   file.name <- file.path(wd.rt.plots, paste0("RT_", base0, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0))
   plotRT(file.name, main.text, chr, NA, NA, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("Tumour", "Normal"), c(adjustcolor.gray, adjustcolor.gray), c("T", "N"), "png", width=10, peaks=c(), 7.75, 9, 3, 3)
   
   ## chr2
   #plotRT(file.name, main.text, chr, NA, NA, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("Tumour", "Normal"), c(adjustcolor.gray, adjustcolor.gray), c("T", "N"), "png", width=10, peaks=c(74393001, 85508001), 7.75, 9, 3, 3)   ## T/N: 74353001, 85951001
   #plotRT(file.name, paste0(BASE, " T/N read depth ratio"), chr, 71500000, 90500000, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("Tumour", "Normal"), c(adjustcolor.red, adjustcolor.blue), c("T", "N"), "png", width=5, peaks=c(74393001, 85508001), 7.75, 9, 3, 3)
   ## chr4
   #plotRT(file.name, main.text, chr, NA, NA, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("Tumour", "Normal"), c(adjustcolor.gray, adjustcolor.gray), c("T", "N"), "png", width=10, peaks=c(128986001, 140285001), 7.75, 9, 3, 3)
   #plotRT(file.name, paste0(BASE, " T/N read depth ratio"), chr, 125000000, 144000000, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("Tumour", "Normal"), c(adjustcolor.red, adjustcolor.blue), c("T", "N"), "png", width=5, peaks=c(128986001, 140285001), 7.75, 9, 3, 3)
   ## chr22
   #plotRT(file.name, main.text, chr, NA, NA, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("Tumour", "Normal"), c(adjustcolor.gray, adjustcolor.gray), c("T", "N"), "png", width=10, peaks=c(), 7.75, 9, 3, 3)
}

# -----------------------------------------------------------------------------
# Plot RT for individual genes (see ReplicationTiming.R)
# Last Modified: 26/04/19; 14/02/19; 10/01/19; 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
genes <- c("ERCC6L2", "GNAQ")

plotWholeChr <- T
ranges <- c(50000, 500000, 5000000)
for (g in 1:length(genes)) {
   gene  <- getGene(genes[g])
   chr   <- gene$chromosome_name
   start <- gene$start_position
   end   <- gene$end_position
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   
   rpkms.chr.rt <- readTable(file.path(wd.rt.data, paste0(base1, "_rpkm.corr.gc.d.rt_", chr, "_", BASE, "-", BASE, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, scaledRT=T) 
   rpkms.chr.rt.T  <- setSpline(rpkms.chr.rt, bed.gc.chr, "T")
   rpkms.chr.rt.N  <- setSpline(rpkms.chr.rt, bed.gc.chr, "N")
   rpkms.chr.rt.RT <- setSpline(rpkms.chr.rt, bed.gc.chr, "RT")

   ## RD & RT 
   main.text <- paste0(BASE, "T/N read depth ratio between tumour (n=", n1, ") and normal (n=", n0, ") samples")
   file.name <- file.path(wd.rt.plots, "genes", paste0("RTD_", base, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, "_", genes[g]))
   if (plotWholeChr)
      plotRD3(file.name, main.text, chr, NA, NA, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("Tumour", "Normal"), c(adjustcolor.gray, adjustcolor.gray), c("T", "N"), "png", width=10, peaks=c(start, end), 7.75, 9, 3, 3)
 
   for (r in 1:length(ranges))
      plotRD3(file.name, paste0(BASE, " T/N read depth ratio"), chr, start-ranges[r],	end+ranges[r], rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("T", "N"), c(adjustcolor.red, adjustcolor.blue), c("T", "N"), "png", width=5, peaks=c(start, end), 7.75, 9, 3, 3)
}

# -----------------------------------------------------------------------------
# SCLC RD vs RT
# Last Modified: 18/02/19
# -----------------------------------------------------------------------------
cors <- toTable(0, 6, 22, c("chr", "length", "cor1", "cor2", "intercept1", "intercept2"))
cors$chr <- 1:22

for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   
   rpkms.chr.rt <- readTable(file.path(wd.rt.data, paste0(base, "_rpkm.corr.gc.d.rt_", chr, "_", BASE, "-", BASE, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, scaledRT=T) 
   rpkms.chr.rt.T  <- setSpline(rpkms.chr.rt, bed.gc.chr, "T")
   rpkms.chr.rt.N  <- setSpline(rpkms.chr.rt, bed.gc.chr, "N")
   rpkms.chr.rt.RT <- setSpline(rpkms.chr.rt, bed.gc.chr, "RT")

   rpkms.chr.rt.lcl <-readTable(paste0("/Users/tpyang/Work/uni-koeln/tyang2/LCL/analysis/replication/lcl-wgs-rt/data/lcl_rpkm.corr.gc.d.rt_", chr, "_LCL-LCL_n7-7.txt.gz"), header=T, rownames=T, sep="\t")
   rpkms.chr.rt.lcl <- setScaledRT(rpkms.chr.rt.lcl, pseudocount=0.01, recaliRT=T, scaledRT=T) 
   rpkms.chr.rt.lcl.RT <- setSpline(rpkms.chr.rt.lcl, bed.gc.chr, "RT")

   ## Keep only overlapping 1kb windows
   overlaps <- intersect(rpkms.chr.rt.T$BED, rpkms.chr.rt.lcl.RT$BED)
   cors$length[c] <- length(overlaps)
   
   main.text <- paste0("SCLC read depths vs. LCL RT (", "Chr", c, ")")
   xlab.text <- "LCL S/G1 RT"
   ylab.text <- "SCLC read depth [log2]"
   file.name <- file.path(wd.rt.plots, paste0("plot_RD-vs-RT_SCLC-vs-LCL_chr", c, "_spearman_spline"))
   #xmin <- min(rpkms.chr.rt.lcl.RT[overlaps,]$SPLINE)
   #xmax <- max(rpkms.chr.rt.lcl.RT[overlaps,]$SPLINE)
   #ymin <- min(c(rpkms.chr.rt.T[overlaps,]$SPLINE, rpkms.chr.rt.N[overlaps,]$SPLINE))
   #ymax <- max(c(rpkms.chr.rt.T[overlaps,]$SPLINE, rpkms.chr.rt.N[overlaps,]$SPLINE))
   plotRD2vsRT(rpkms.chr.rt.T[overlaps,]$SPLINE, rpkms.chr.rt.N[overlaps,]$SPLINE, rpkms.chr.rt.lcl.RT[overlaps,]$SPLINE, file.name, main.text, ylab.text, xlab.text, c("red", "blue"), c("T", "N"), method="spearman")

   cors$cor1[c] <- getCor(rpkms.chr.rt.T[overlaps,]$SPLINE, rpkms.chr.rt.lcl.RT[overlaps,]$SPLINE, method="spearman")
   cors$cor2[c] <- getCor(rpkms.chr.rt.N[overlaps,]$SPLINE, rpkms.chr.rt.lcl.RT[overlaps,]$SPLINE, method="spearman")
   cors$intercept1[c] <- lm(rpkms.chr.rt.T[overlaps,]$SPLINE ~ rpkms.chr.rt.lcl.RT[overlaps,]$SPLINE)[[1]][1]
   cors$intercept2[c] <- lm(rpkms.chr.rt.N[overlaps,]$SPLINE ~ rpkms.chr.rt.lcl.RT[overlaps,]$SPLINE)[[1]][1] 
}
save(cors, file=file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-vs-lcl_spearman_spline.RData")))

ylab.text <- "Spearman's rho"
xlab.text <- "Chromosome"
file.name <- file.path(wd.rt.plots, "plot_RD-vs-RT_SCLC-vs-LCL_spearman_spline")
main.text <- paste0("SCLC read depths vs. LCL RT")
ymin <- -0.75
ymax <- 0.75
plotRD2vsRTALL(cors, file.name, main.text, ylab.text, xlab.text, ymin, ymax, cols=c("red", "blue"), c=2)

ylab.text <- "Intercept difference"
xlab.text <- "Chromosome"
file.name <- file.path(wd.rt.plots, "plot_RD-vs-RT_SCLC-vs-LCL_spearman_spline_intercept_0.06")
main.text <- c(paste0("SCLC read depths vs. LCL RT"), "Linear regression")
ymin <- -0.06034766
ymax <- 0.06034766
plotInterceptALL(cors, file.name, main.text, ylab.text, xlab.text, ymin, ymax, c("red", "blue", "black"), c(13, 17))

# -----------------------------------------------------------------------------
# SCLC RT vs LCL RT
# Last Modified: 05/03/19; 18/02/19
# -----------------------------------------------------------------------------
cors <- toTable(0, 3, 22, c("chr", "length", "cor"))
cors$chr <- 1:22

for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   
   ## SCLC and LCL RTs
   rpkms.chr.rt <- readTable(file.path(wd.rt.data, paste0(base1, "_rpkm.corr.gc.d.rt_", chr, "_", BASE1, "-", BASE0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, scaledRT=T)
   #rpkms.chr.rt.T <- setSpline(rpkms.chr.rt, bed.gc.chr, "T")
   #rpkms.chr.rt.N <- setSpline(rpkms.chr.rt, bed.gc.chr, "N")
   rpkms.chr.rt.RT <- setSpline(rpkms.chr.rt, bed.gc.chr, "RT")
   
   rpkms.chr.rt.lcl <-readTable(paste0("/Users/tpyang/Work/uni-koeln/tyang2/LCL/analysis/replication/lcl-wgs-rt/data/lcl_rpkm.corr.gc.d.rt_", chr, "_LCL-LCL_n7-7.txt.gz"), header=T, rownames=T, sep="\t")
   rpkms.chr.rt.lcl <- setScaledRT(rpkms.chr.rt.lcl, pseudocount=0.01, recaliRT=T, scaledRT=T) 
   #rpkms.chr.rt.lcl.T  <- setSpline(rpkms.chr.rt.lcl, bed.gc.chr, "T")
   #rpkms.chr.rt.lcl.N  <- setSpline(rpkms.chr.rt.lcl, bed.gc.chr, "N")
   rpkms.chr.rt.lcl.RT <- setSpline(rpkms.chr.rt.lcl, bed.gc.chr, "RT")
   
   ## Keep only overlapping 1kb windows
   overlaps <- intersect(rpkms.chr.rt.T$BED, rpkms.chr.rt.lcl.RT$BED)
   #rpkms.chr.rt.RT     <- rpkms.chr.rt.RT[overlaps,]       ## Too slow
   #rpkms.chr.rt.lcl.RT <- rpkms.chr.rt.lcl.RT[overlaps,]   ## Too slow
   cors$length[c] <- length(overlaps)
   
   main.text <- paste0("SCLC T/N vs. LCL S/G1 (", "Chr", c, ")")
   xlab.text <- "LCL S/G1 RT"
   ylab.text <- "SCLC T/N RT"
   file.name <- file.path(wd.rt.plots, paste0("plot_RT-vs-RT_SCLC-vs-LCL_spearman_chr", c, "_spline"))
   plotRTvsRT(rpkms.chr.rt.RT[overlaps,]$SPLINE, rpkms.chr.rt.lcl.RT[overlaps,]$SPLINE, file.name, main.text, ylab.text, xlab.text, c(adjustcolor.gray, "black"), method="spearman")

   #cors$T[c]  <- getCor(rpkms.chr.rt.T[overlaps,]$SPLINE, rpkms.chr.rt.lcl.T[overlaps,]$SPLINE, method="spearman")
   #cors$N[c]  <- getCor(rpkms.chr.rt.T[overlaps,]$SPLINE, rpkms.chr.rt.lcl.N[overlaps,]$SPLINE, method="spearman")   
   #cors$RT[c] <- getCor(rpkms.chr.rt.T[overlaps,]$SPLINE, rpkms.chr.rt.lcl.RT[overlaps,]$SPLINE, method="spearman")
   cors$cor[c] <- getCor(rpkms.chr.rt.T[overlaps,]$SPLINE, rpkms.chr.rt.lcl.RT[overlaps,]$SPLINE, method="spearman")
}
save(cors, file=file.path(wd.rt.data, paste0("rt-vs-rt_", base, "-vs-lcl_spearman_spline.RData")))

ylab.text <- "Spearman's rho"
xlab.text <- "Chromosome"
file.name <- file.path(wd.rt.plots, "plot_RT-vs-RT_SCLC-vs-LCL_spearman_spline")
main.text <- paste0("SCLC T/N vs. LCL S/G1")
ymin <- 0.2
ymax <- 1
plotRTvsRTALL(cors, file.name, main.text, ylab.text, xlab.text, ymin, ymax, col="black", c=2)

# -----------------------------------------------------------------------------
# 
# Last Modified: 27/05/19
# -----------------------------------------------------------------------------
ylab.text <- "Pearson's r"
xlab.text <- "Chromosome"
file.name <- file.path(wd.rt.plots, "plot_rt-vs-rt_SCLC-T-N-vs-LCL-S-G1_cors-pearson")
main.text <- paste0("SCLC T/N vs. LCL S/G1")
ymin <- 0.25
ymax <- 1.05
plotRTvsRTALL(cors, file.name, main.text, ylab.text, xlab.text, ymin, ymax, line0=T, c=4)

# -----------------------------------------------------------------------------
# 
# http://pklab.med.harvard.edu/scw2014/subpop_tutorial.html
# Last Modified: 20/04/19; 06/03/19
# -----------------------------------------------------------------------------
cors.samples <- toTable(0, length(samples1)+4, 22, c("chr", "mean", "var", "cv2", samples1))
cors.samples$chr <- 1:22
for (s in 1:length(samples1)) {
   sample <- samples1[s]
   load(file.path(wd.rt.data, "samples", paste0("rd-vs-rt_", sample, "-vs-lcl_cors-pearson_spline.RData")))
 
   cors.samples[, sample] <- cors$cor
}

for (c in 1:22) {
   cors.samples$mean[c] <- mean(as.numeric(cors.samples[c, samples1]))
   cors.samples$var[c]  <- var(as.numeric(cors.samples[c, samples1]))
   cors.samples$cv2[c]  <- cors.samples$var[c]/cors.samples$mean[c]^2
}
save(cors.samples, file=file.path(wd.rt.data, paste0("rd-vs-rt_samples-vs-lcl_cors-pearson_spline.RData")))
# > min(cors.samples[,-c(1:4)])
# [1] -0.7974991
# > max(cors.samples[,-c(1:4)])
# [1] 0.8457572

file.name <- file.path(wd.rt.plots, "plot_RD-vs-RT_SAMPLES-vs-LCL_pearson")
main.text <- c("SCLC read depth profiles (n=101) vs. LCL RT", "SCLC T vs. LCL S/G1 RT")   ## TO-DO
ymin <- -0.8318379
ymax <- 0.8457572
plotSAMPLEvsRTALL(cors.samples, samples1, file.name, main.text, ymin, ymax)

# -----------------------------------------------------------------------------
# CM2 and CQ4
# Last Modified: 19/05/19
# -----------------------------------------------------------------------------
cors <- t(cors.samples[, -c(1:4)]) 
colnames(cors) <- paste0("chr", c(1:22))
cm2 <- cors
cq4 <- cors

for (c in 1:22) {
   cors.chr <- cors[,c]
   q <- quantile(as.numeric(cors.chr))
 
   cq4[which(cors[, c] > q[4]), c] <- 4
   cq4[intersect(as.vector(which(cors[, c] > q[3])), as.vector(which(cors[, c] <= q[4]))), c] <- 3
   cq4[intersect(as.vector(which(cors[, c] > q[2])), as.vector(which(cors[, c] <= q[3]))), c] <- 2
   cq4[which(cors[, c] <= q[2]), c] <- 1
   
   cm2[which(cors[, c] > q[3]), c] <- 1
   cm2[which(cors[, c] <= q[3]), c] <- 0
}
cq4 <- as.data.frame(cq4)
cq4$SAMPLE_ID <- ""
cq4$SAMPLE_ID <- rownames(cors)
cm2 <- as.data.frame(cm2)
cm2$SAMPLE_ID <- ""
cm2$SAMPLE_ID <- rownames(cors)
writeTable(cq4[, c("SAMPLE_ID", paste0("chr", c(1:22)))], file.path(wd.ngs, "sclc_wgs_n101.cq4"), colnames=T, rownames=F, sep="\t")
writeTable(cm2[, c("SAMPLE_ID", paste0("chr", c(1:22)))], file.path(wd.ngs, "sclc_wgs_n101.cm2"), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Divide tumour samples into Q4
# Last Modified: 06/03/19
# -----------------------------------------------------------------------------
cors.all <- toTable(0, 4, length(samples1), c("SAMPLE_ID", "COR", "Q4", "RT"))
rownames(cors.all) <- samples1
cors.all$SAMPLE_ID <- samples1
for (s in 1:length(samples1)) {
   sample <- samples1[s]
   load(file.path(wd.rt.data, "samples", paste0("rd-vs-rt_", sample, "-vs-lcl_cors-pearson.RData")))
 
   cors.all$COR[s] <- cor
}

q <- quantile(as.numeric(cors.all$COR))
# 0%        25%        50%        75%       100% 
# -0.2991067 -0.2308337 -0.1999576 -0.1358137  0.3605179

samples.q4 <- list()
samples.q4[[4]] <- rownames(subset(cors.all, COR > as.numeric(q[4])))
samples.q4[[3]] <- rownames(subset(subset(cors.all, COR > as.numeric(q[3])), COR <= as.numeric(q[4])))
samples.q4[[2]] <- rownames(subset(subset(cors.all, COR > as.numeric(q[2])), COR <= as.numeric(q[3])))
samples.q4[[1]] <- rownames(subset(cors.all, COR <= as.numeric(q[2])))

cors.all$Q4[which(cors.all$SAMPLE_ID %in% samples.q4[[4]])] <- 4
cors.all$Q4[which(cors.all$SAMPLE_ID %in% samples.q4[[3]])] <- 3
cors.all$Q4[which(cors.all$SAMPLE_ID %in% samples.q4[[2]])] <- 2
cors.all$Q4[which(cors.all$SAMPLE_ID %in% samples.q4[[1]])] <- 1

cors.all$RT[which(cors.all$Q4 %in% c(3, 4))] <- 1
cors.all$RT[which(cors.all$Q4 %in% c(1, 2))] <- 0
writeTable(cors.all, file.path(wd.ngs, "sclc_wgs_n101.txt"), colnames=T, rownames=F, sep="\t")

cors.all <- subset(cors.all, Q4 %in% c(4,1))
writeTable(cors.all, file.path(wd.ngs, "sclc_wgs_n51.txt"), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# PCA
# Last Modified: 21/04/19
# -----------------------------------------------------------------------------
samples  <- readTable(file.path(wd.ngs, "sclc_wgs_n101.txt"), header=T, rownames=T, sep="")

trait <- as.numeric(samples$Q4)
trait[which(trait == 4)] <- "Q4 (-0.14 < r < 0.36)"
trait[which(trait == 3)] <- "Q3 (-0.20 < r < -0.14)"
trait[which(trait == 2)] <- "Q2 (-0.23 < r < -0.20)"
trait[which(trait == 1)] <- "Q1 (-0.30 < r < -0.23)"

rpkms.T.chr.d.all <- NA
## Copy from 2a_cmd-rt_rpkm.corr.gc.d_sample.R (commandline mode)
for (c in 22:1) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
 
   ## Read depth
   rpkms.T.chr.d <- pipeGetDetectedRD(wd1.ngs.data, BASE1, chr, PAIR1)
   #rpkms.N.chr.d <- pipeGetDetectedRD(wd0.ngs.data, BASE0, chr, PAIR0) 
   #overlaps <- intersect(rpkms.T.chr.d$BED, rpkms.N.chr.d$BED)
 
   ##
   test <- rpkms.T.chr.d[, samples$SAMPLE_ID]
   pca.de <- getPCA(t(test))

   file.main <- c(paste0("SCLC tumours on chr", c), "")
   plotPCA(1, 2, pca.de, trait, wd.rt.plots, paste0("pca_sclc_T_chr", c), size=6.5, file.main, "topleft", c("blue", "skyblue3", "lightcoral", "red"), NULL, flip.x=1, flip.y=1)
   save(pca.de, file=file.path(wd.rt.data, paste0("pca_sclc_T_chr", c, ".RData")))
 
   if (is.na(rpkms.T.chr.d.all)) {
      rpkms.T.chr.d.all <- test
   } else
      rpkms.T.chr.d.all <- rbind(rpkms.T.chr.d.all, test)
}

##
test <- rpkms.T.chr.d.all
pca.de <- getPCA(t(test))
save(pca.de, file=file.path(wd.rt.data, paste0("pca_sclc_T_chrs.RData")))

file.main <- c("SCLC read depth profiles (n=101)", "")
plotPCA(1, 2, pca.de, trait, wd.rt.plots, "pca_sclc_T_chrs", size=6, file.main, "bottomright", c("red", "lightcoral", "skyblue3", "blue"), NULL, flip.x=1, flip.y=1, legend.title="Overall corr. with LCL RT", decreasing=T)

# -----------------------------------------------------------------------------
# Before PCA (Beeswarm plots)
# Last Modified: 21/04/19
# -----------------------------------------------------------------------------
samples.sclc <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/SCLC/ngs/WGS/sclc_wgs_n101.txt", header=T, rownames=T, sep="")
samples.nbl  <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/NBL/ngs/WGS/nbl_wgs_n57-1.txt", header=T, rownames=T, sep="")
samples.cll  <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/CLL/ngs/WGS/cll_wgs_n96.txt", header=T, rownames=T, sep="")
n.sclc <- nrow(samples.sclc)
n.nbl  <- nrow(samples.nbl)
n.cll  <- nrow(samples.cll)

samples <- toTable(0, 3, n.sclc+n.nbl+n.cll, c("CANCER", "COR", "Q4"))
samples$CANCER[1:n.sclc] <- 0
samples$CANCER[(1+n.sclc):(n.sclc+n.nbl)] <- 1
samples$CANCER[(1+n.sclc+n.nbl):(n.sclc+n.nbl+n.cll)] <- 2
samples$COR <- c(samples.sclc$COR, samples.nbl$COR, samples.cll$COR)
samples$Q4 <- c(samples.sclc$Q4, samples.nbl$Q4, samples.cll$Q4)

#install.packages('beeswarm')
library(beeswarm)

pdf(file.path(wd.rt.plots, "beeswarm_sclc+nbl+cll.pdf"), height=5.2, width=5.2)
ymax <- max(samples$COR)
ymin <- -ymax
boxplot(COR ~ CANCER, data=samples, outline=F, names=c("SCLC (n=101)", "NBL (n=56)", "CLL (n=96)"), ylim=c(ymin, ymax), ylab="Pearson's r", main="Overall correlation with LCL RT")
abline(h=0, lty=5)

beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 1), col="blue", pch=16, add=T)
beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 2), col="skyblue3", pch=16, add=T)
beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 3), col="lightcoral", pch=16, add=T)
beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 4), col="red", pch=16, add=T)

legend("topright", legend = c("Q4", "Q3", "Q2", "Q1"), pch=16, col=c("red", "lightcoral", "skyblue3", "blue"))
dev.off()

# -----------------------------------------------------------------------------
# Find S-like and G1-like tumour samples
# Last Modified: 22/04/19
# -----------------------------------------------------------------------------
s_likes <- c()
g1_likes <- c()
for (c in 1:22) {
   samples <- cors.samples[c, -c(1:4)]
   median <- median(as.numeric(samples))
 
   if (c == 1) {
      s_likes  <- colnames(samples)[which(as.numeric(samples) > median)]
      g1_likes <- colnames(samples)[which(as.numeric(samples) <= median)]
   } else {
      s_likes  <- intersect(s_likes, colnames(samples)[which(as.numeric(samples) > median)])
      g1_likes <- intersect(g1_likes, colnames(samples)[which(as.numeric(samples) <= median)])
   }
}
# > length(s_likes)
# [1] 20
# > length(g1_likes)
# [1] 12
# > s_likes
# [1] "S00827" "S00832" "S00833" "S00837" "S00838" "S00941" "S01020" "S01022" "S01023" "S01366" "S01524" "S01864" "S02120" "S02284" "S02285" "S02287" "S02289" "S02290" "S02292" "S02295"
# > g1_likes
# [1] "S00050" "S00947" "S01516" "S01578" "S01861" "S01873" "S02139" "S02237" "S02241" "S02248" "S02256" "S02378"

samples$G1S <- NA
samples[s,]$G1S <- "S-like"
samples[g1,]$G1S <- "G1-like"

writeTable(samples, file.path(wd.ngs, "sclc_wgs_n101.txt"), colnames=T, rownames=F, sep="\t")









# -----------------------------------------------------------------------------
# SCLC-CLL-N RT vs LCL RT
# Last Modified: 18/02/19
# -----------------------------------------------------------------------------
BASE0 <- "CLL"
base0 <- tolower(BASE0)
n0 <- 96

cors <- toTable(0, 2, 22, c("chr", "cor"))
cors$chr <- 1:22
for (c in 1:22) {
 chr <- chrs[c]
 bed.gc.chr <- subset(bed.gc, CHR == chr)
 
 rpkms.chr.rt <- readTable(file.path(wd.rt.data, "hybrid", paste0(base1, "_rpkm.corr.gc.d.rt_", chr, "_", BASE1, "-", BASE0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
 rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, flipRT=F, scaledRT=T) 
 rpkms.chr.rt.RT <- setSlope(rpkms.chr.rt, bed.gc.chr, "RT")
 
 rpkms.chr.rt.lcl <-readTable(paste0("/Users/tpyang/Work/uni-koeln/tyang2/LCL/analysis/replication/lcl-wgs-rt/data/lcl_rpkm.corr.gc.d.rt_", chr, "_LCL-LCL_n7-7.txt.gz"), header=T, rownames=T, sep="\t")
 rpkms.chr.rt.lcl <- setScaledRT(rpkms.chr.rt.lcl, pseudocount=0.01, recaliRT=T, flipRT=F, scaledRT=T) 
 rpkms.chr.rt.lcl.RT <- setSlope(rpkms.chr.rt.lcl, bed.gc.chr, "RT")
 
 overlaps <- intersect(rpkms.chr.rt.RT$BED, rpkms.chr.rt.lcl.RT$BED)
 #rpkms.chr.rt.RT     <- rpkms.chr.rt.RT[overlaps,]
 #rpkms.chr.rt.lcl.RT <- rpkms.chr.rt.lcl.RT[overlaps,]
 
 cors$cor[c] <- getCor(rpkms.chr.rt.lcl.RT[overlaps,]$SLOPE, rpkms.chr.rt.RT[overlaps,]$SLOPE)
}
save(cors, file=file.path(wd.rt.data, paste0("rt-vs-rt_", base1, "&", base0, "-N-vs-lcl_cors-pearson.RData")))

ylab.text <- "Pearson's r"
xlab.text <- "Chromosome"
file.name <- file.path(wd.rt.plots, "plot_RT-vs-RT_SCLC&CLL-N-vs-LCL_pearson")
main.text <- paste0("SCLC–CLL T/N vs. LCL RT")
ymin <- -0.2
ymax <- 0.85
plotRTvsRTALL(cors, file.name, main.text, ylab.text, xlab.text, ymin, ymax, line0=T)


# -----------------------------------------------------------------------------
# 
# Last Modified: 06/03/19
# -----------------------------------------------------------------------------
# > length(which(as.numeric(cors.samples[2, -c(1:4)]) > 0))
# [1] 67
# > length(which(as.numeric(cors.samples[2, -c(1:4)]) < 0))
# [1] 25

idx.pos.chrs <- which(as.numeric(cors.samples[1, -c(1:4)]) > 0) + 4
idx.neg.chrs <- which(as.numeric(cors.samples[1, -c(1:4)]) < 0) + 4
for (c in 2:22) {
 idx.pos.chr <- which(as.numeric(cors.samples[c, -c(1:4)]) > 0) + 4
 idx.neg.chr <- which(as.numeric(cors.samples[c, -c(1:4)]) < 0) + 4
 
 idx.pos.chrs <- intersect(idx.pos.chrs, idx.pos.chr)
 idx.neg.chrs <- intersect(idx.neg.chrs, idx.neg.chr)
}
# > length(idx.pos.chrs)
# [1] 56
# > length(idx.neg.chrs)
# [1] 21
idx.var.chrs <- c(5:96)
idx.var.chrs <- setdiff(idx.var.chrs, c(idx.pos.chrs, idx.neg.chrs))
# > length(idx.var.chrs)
# [1] 15
# > 56+21+15
# [1] 92

##
samples.pos <- colnames(cors.samples[, idx.pos.chrs])
samples.neg <- colnames(cors.samples[, idx.neg.chrs])
samples.var <- colnames(cors.samples[, idx.var.chrs])

writeTable(samples.pos, file.path(wd.ngs, "sclc_wgs_n92-rt56.list"), colnames=F, rownames=F, sep="")
writeTable(samples.neg, file.path(wd.ngs, "sclc_wgs_n92-wt21.list"), colnames=F, rownames=F, sep="")
writeTable(samples.var, file.path(wd.ngs, "sclc_wgs_n92-var15.list"), colnames=F, rownames=F, sep="")











# -----------------------------------------------------------------------------
# SCLC RT vs LCL RT
# Last Modified: 18/02/19
# -----------------------------------------------------------------------------
for (c in 1:22) {
   chr <- chrs[c]

   BASE0 <- "NBL"
   base0 <- tolower(BASE0)
   n0 <- 56
   rpkms.chr.rt <- readTable(file.path(wd.rt.data, "hybrid", paste0(base1, "_rpkm.corr.gc.d.rt_", chr, "_", BASE1, "-", BASE0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   getCor
   #rpkms.chr.rt <- readTable(file.path(wd.rt.data, paste0(base1, "_rpkm.corr.gc.d.rt_", chr, "_", BASE1, "-", BASE0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, scaledRT=T) 
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   bed.gc.chr <- bed.gc.chr[rpkms.chr.rt$BED,]
   rpkms.chr.rt.RT     <- setSlopes(rpkms.chr.rt, bed.gc.chr, "RT")
   
   rpkms.chr.rt.lcl <-readTable(paste0("/Users/tpyang/Work/uni-koeln/tyang2/LCL/analysis/replication/lcl-wgs-rt/data/lcl_rpkm.corr.gc.d.rt_", chr, "_LCL-LCL_n7-7.txt.gz"), header=T, rownames=T, sep="\t")
   rpkms.chr.rt.lcl <- setScaledRT(rpkms.chr.rt.lcl, pseudocount=0.01, recaliRT=T, scaledRT=T) 
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   bed.gc.chr <- bed.gc.chr[rpkms.chr.rt.lcl$BED,]
   rpkms.chr.rt.lcl.RT <- setSlopes(rpkms.chr.rt.lcl, bed.gc.chr, "RT")

   overlaps <- intersect(rpkms.chr.rt.RT$BED, rpkms.chr.rt.lcl.RT$BED)
   rpkms.chr.rt.RT     <- rpkms.chr.rt.RT[overlaps,]
   rpkms.chr.rt.lcl.RT <- rpkms.chr.rt.lcl.RT[overlaps,]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   bed.gc.chr <- bed.gc.chr[overlaps,]
   
   ylab.text <- "LCL RT (S/G1)"
   xlab.text <- "SCLC RT (T/N)"
   file.name <- file.path(wd.rt.plots, paste0("plot_RT-vs-RT_SCLC&SCLC-vs-LCL_pearson_chr", c))
   main.text <- paste0("SCLC RT vs. LCL RT (", "Chr", c, ")")
   cors$cor[c] <- plotRDvsRT(rpkms.chr.rt.lcl.RT$SLOPE, rpkms.chr.rt.RT$SLOPE, file.name, main.text, ylab.text, xlab.text, c(adjustcolor.gray, "black"))
}



ylab.text <- "LCL (S/G1)"
xlab.text <- "SCLC (T/N)"
file.name <- file.path(wd.rt.plots, paste0("plot_RT-LCL_vs_RT-SCLC_slope_1kb_chr2"))
main.text <- paste0("Correlation between RT profiles (", "Chr2", ")")
cor <- plotRDvsRT(rpkms.chr.rt.lcl.RT$SLOPE, rpkms.chr.rt.RT$SLOPE, file.name, main.text, ylab.text, xlab.text, c("grey", "black"))




   ## LCL RD vs RT
   rpkms.chr.rt.lcl <-readTable(paste0("/Users/tpyang/Work/uni-koeln/tyang2/LCL/analysis/replication/lcl-wgs-rt/data/lcl_rpkm.corr.gc.d.rt_", chr, "_LCL-LCL_n7-7.txt.gz"), header=T, rownames=T, sep="\t")
   rpkms.chr.rt.lcl$RT <- scale(rpkms.chr.rt.lcl$RT)   ## ADD 19/02/19

bed.gc.chr <- subset(bed.gc, CHR == chr)
bed.gc.chr <- bed.gc.chr[rpkms.chr.rt.lcl$BED,]
overlaps <- intersect(rpkms.chr.rt.lcl$BED, rownames(bed.gc.chr))
idxes <- seq(1, length(overlaps), 1)
bed.gc.chr <- bed.gc.chr[overlaps[idxes],]

rpkms.chr.rt.lcl.RT <- setSlopes(rpkms.chr.rt.lcl, bed.gc.chr)
rpkms.chr.rt.lcl.RD  <- setSlopes2(rpkms.chr.rt.lcl, bed.gc.chr)

ylab.text <- "Read depth (G1)"
xlab.text <- "RT (LCL)"
file.name <- file.path(wd.rt.plots, "plot_RD_vs_RT_LCL_G1-LCL_slope_rho_1kb_log2")
main.text <- paste0("Slopes of each 1kb window (Chr2)")
plotRDvsRT(rpkms.chr.rt.lcl.RD$SLOPE_N, rpkms.chr.rt.lcl.RT$SLOPE, file.name, main.text, ylab.text, xlab.text, c("lightskyblue1", "blue"))

ylab.text <- "Read depth (S)"
file.name <- file.path(wd.rt.plots, "plot_RD_vs_RT_LCL_S-LCL_slope_rho_1kb_log2")
plotRDvsRT(rpkms.chr.rt.lcl.RD$SLOPE_T, rpkms.chr.rt.lcl.RT$SLOPE, file.name, main.text, ylab.text, xlab.text, c("pink", "red"))










# -----------------------------------------------------------------------------
# Insert size (each sample)
# Last Modified: 20/09/18
# -----------------------------------------------------------------------------
wd.meta <- file.path(wd.ngs, "insert_size")

for (s in 1:length(samples)) {
   table <- readTable(file.path(wd.meta, paste0(samples[s], "_qc.txt")), header=F, rownames=F, sep="")   
}



# -----------------------------------------------------------------------------
# Insert size
# Last Modified: 10/09/18
# -----------------------------------------------------------------------------
wd.meta <- file.path(wd, BASE, "metadata/George 2015")
table <- readTable(file.path(wd.meta, "nature14664-s1_ST2.txt"), header=T, rownames=T, sep="\t")
table <- table[samples,]

table1 <- table[,c("Insert.Size", "Mean.Coverage")]
table1$Group <- 1
table0 <- table[,c("Insert.Size.1", "Mean.Coverage.1")]
colnames(table0) <- c("Insert.Size", "Mean.Coverage")
table0$Group <- 0
table <- rbind(table1, table0)
table$Group <- as.factor(table$Group)

file.name <- file.path(wd.rt.plots, paste0("boxplot_", base, "_wgs_insert-size.pdf"))
pdf(file.name, height=6, width=3)
boxplot(Insert.Size ~ Group, data=table, outline=F, names=c("Normal", "Tumour"), col=c("dodgerblue", "red"), ylab="Insert size", main=BASE)
dev.off()

median(table1$Insert.Size)
# [1] 312.312
median(table0$Insert.Size)
# [1] 308.56
median(table0[normals,]$Insert.Size)

sd(table1$Insert.Size)
# [1] 16.19834
sd(table0$Insert.Size)
# [1] 13.34431

## Tumour vs. Normal
testW(table1$Insert.Size, table0$Insert.Size)
# [1] 0.5489026


# -----------------------------------------------------------------------------
# Step 6.1: Define replicaiton timing direction for expressed genes (Following Step 4 in "asym-sclc-tx.R" and Step 5 from rt-sclc-wgs.R)
# Link(s):  http://www.mun.ca/biology/scarr/2250_DNA_replication_&_transcription.html
#           https://stackoverflow.com/questions/43615469/how-to-calculate-the-slope-of-a-smoothed-curve-in-r
# Last Modified: 29/01/18
# -----------------------------------------------------------------------------
plotRT0 <- function(wd.rt.plots, BASE, chr, n1, n0, xmin, xmax, rpkms.chr.rt, bed.gc.chr.rt, pair1, pair0, ext, spline) {
   file.name  <- file.path(wd.rt.plots, paste0(tolower(BASE), "_wgs_rt_bstrp1000_", chr, "_", pair1, "-", pair0, "_n", n1, "-", n0, "_spline"))
   main.text <- paste0("Bootstrapped read depth (CN-, GC-corrected) ratio (", pair1, "/", pair0, ") in ", BASE)
   xlab.text <- paste0("Chromosome ", gsub("chr", "", chr), " coordinate (Mb)")
   ylab.text <- "Replication time (log2 FC)"
 
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   ymin <- min(rpkms.chr.rt$MEDIAN)
   ymax <- max(rpkms.chr.rt$MEDIAN)
   if (!is.na(xmin) && !is.na(xmax)) file.name <- paste0(file.name, "_", xmin/1E6, "-", xmax/1E6, "Mb")
   if (is.na(xmin)) xmin <- 0
   if (is.na(xmax)) xmax <- subset(chromInfo, chrom == chr)$size

   if (ext == "pdf") {
      pdf(paste0(file.name, ".pdf"), height=4, width=10)
   } else if (ext == "png")
      png(paste0(file.name, ".png"), height=4, width=10, units="in", res=300)   ## ADD 16/05/17: res=300
   
   plot(NULL, ylim=c(ymin, ymax), xlim=c(xmin/1E6, xmax/1E6), xlab=xlab.text, ylab=ylab.text, main=main.text)
   points(bed.gc.chr.rt$START/1E6, rpkms.chr.rt$MEDIAN, col="red", cex=0.3)
   abline(h=0, lwd=0.5, col="grey")
   lines(bed.gc.chr.rt$START/1E6, smooth.spline(rpkms.chr.rt$MEDIAN)$y)
   lines(spline$x/1E6, spline$y, col="blue")
   
   #slopes <- diff(smooth.spline(rpkms.chr)$y)/diff((bed.gc.chr$START)/1E6)
   #slopes2 <- diff(smooth.spline(rpkms.chr)$y)/diff(smooth.spline(rpkms.chr)$x)
   
   #temp <- loess.smooth(bed.gc.chr$START, rpkms.chr)
   #slopes3 <- diff(temp$y)/diff(temp$x)
   
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.4, col="lightgrey")
 
   dev.off()
}

getEnsGeneBED <- function(pos, bed.gc.chr) {
   bed.gc.chr.start <- subset(bed.gc.chr, pos >= START)
   bed.gc.chr.start.end <- subset(bed.gc.chr.start, pos <= END)
 
   return(rownames(bed.gc.chr.start.end))
}

BASE  <- "SCLC"
PAIR1 <- "T"
PAIR0 <- "N"
PAIR  <- paste0(PAIR1, "-", PAIR0)
CHR   <- 6
CUTOFF <- 0.15

###
##
bed.gc <- bed[which(bed$GC > 0),]   ## Only keep partitions (in the BED file) with a GC content
ensGene.tx <- ensGene[rownames(tpm.gene.input),]

ensGene.tx.rt <- ensGene.tx[1,]
ensGene.tx.rt$SLOPE_START <- 0
ensGene.tx.rt$SLOPE_END <- 0
ensGene.tx.rt <- ensGene.tx.rt[-1,]
for (c in 1:22) {
   #chr <- chrs[CHR]
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)

   ## Replication timing
   #rpkms.chr <- readTable(file.path(wd.rt.data, paste0(base, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR, "_n", length(samples), "-", length(normals), ".txt.gz")), header=T, rownames=T, sep="\t") 
   rpkms.chr <- readTable(file.path(wd.rt.data, paste0(base, "_rpkm.corr.gc.d.rt_bstrp1000_", chr, "_T-N_n101-92.txt.gz")), header=T, rownames=T, sep="\t") 
   rpkms.chr$MEDIAN <- rpkms.chr$RT   ## REMOVED 07/10/18
   
   ##
   rpkms.chr.rt <- rpkms.chr[which(rpkms.chr$MEDIAN > -CUTOFF),]
   rpkms.chr.rt <- rpkms.chr.rt[which(rpkms.chr.rt$MEDIAN < CUTOFF),]
   overlaps <- intersect(rownames(rpkms.chr.rt), rownames(bed.gc.chr))
   bed.gc.chr.rt <- bed.gc.chr[overlaps,]
   
   #plotRT0(wd.rt.plots, BASE, chr, 101, 92, NA, NA, rpkms.chr.rt, bed.gc.chr.rt, PAIR1, PAIR0, "png")
   plotRT0(wd.rt.plots, BASE, chr, 101, 92, 87730060,98505459, rpkms.chr.rt, bed.gc.chr.rt, PAIR1, PAIR0, "png")   ## CCDC67
   #plotRT0(wd.rt.plots, BASE, chr, 101, 92, 96754840, 131794839, rpkms.chr.rt, bed.gc.chr.rt, PAIR1, PAIR0, "png")   ## HDAC2
   #plotRT0(wd.rt.plots, BASE, chr, length(samples), 140813453, 153118090, rpkms.chr.rt, bed.gc.chr, PAIR1, PAIR0, "png")   ## CNTNAP2
   #plotRT0(wd.rt.plots, BASE, chr, length(samples), 43877887, 54056122, rpkms.chr.rt, bed.gc.chr, PAIR1, PAIR0, "png")   ## RB1
   #plotRT0(wd.rt.plots, BASE, chr, length(samples), 147504475, 149581413, rpkms.chr.rt, bed.gc.chr, PAIR1, PAIR0, "png")   ## EZH2
   #plotRT0(wd.rt.plots, BASE, chr, length(samples), 7314246, 11612723, rpkms.chr.rt, bed.gc.chr, PAIR1, PAIR0, "png")   ## PTPRD
   #plotRT0(wd.rt.plots, BASE, chr, length(samples), 114521235, 118716095, rpkms.chr.rt, bed.gc.chr, PAIR1, PAIR0, "png")   ## LSAMP
   
   ## Determin replication direction for each expressed gene
   slopes <- diff(smooth.spline(rpkms.chr.rt$MEDIAN)$y)/diff((bed.gc.chr$START)/1E7)   ## WHY?
 
   ensGene.tx.chr <- subset(ensGene.tx, chromosome_name == chr)
   ensGene.tx.chr$SLOPE_START <- NA
   ensGene.tx.chr$SLOPE_START <- NA
   for (g in 1:nrow(ensGene.tx.chr)) {
      gene <- ensGene.tx.chr[g,]
      bed.s <- getEnsGeneBED(gene$start_position, bed.gc.chr)
      bed.e <- getEnsGeneBED(gene$end_position, bed.gc.chr)
      
      if (length(bed.s) != 0) ensGene.tx.chr$SLOPE_START[g] <- slopes[which(rownames(bed.gc.chr) == bed.s[1])]
      if (length(bed.e) != 0) ensGene.tx.chr$SLOPE_END[g] <- slopes[which(rownames(bed.gc.chr) == bed.e[1])]
   }
   ensGene.tx.rt <- rbind(ensGene.tx.rt, ensGene.tx.chr)
}
save(ensGene.tx.rt, file=file.path(wd.asym.data, paste0(base, "_asym_tx_rt_bstrp1000.RData")))
# > nrow(ensGene.tx.rt)   ## All genes
# [1] 18440
# > nrow(ensGene.tx.rt)   ## All pcg genes
# [1] 16410
# > nrow(ensGene.tx.rt)   ## All non-pcg genes
# [1] 10604

