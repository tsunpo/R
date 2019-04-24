# =============================================================================
# Manuscript   : 
# Chapter      : Reconstruction of replication timing profile in tumour cells
# Name         : manuscripts/replication/nbl-wgs-rt.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 25/02/19
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
BASE1 <- "NBL"
PAIR1 <- "T"
BASE0 <- "NBL"
PAIR0 <- "T"
base1 <- tolower(BASE1)
base0 <- tolower(BASE0)

wd1.ngs    <- file.path(wd, BASE1, "ngs/WGS")
wd1.ngs.data <- file.path(wd1.ngs, "data") 
wd0.ngs    <- file.path(wd, BASE0, "ngs/WGS")
wd0.ngs.data <- file.path(wd0.ngs, "data")

wd.anlys <- file.path(wd, BASE1, "analysis")
wd.rt       <- file.path(wd.anlys, "replication", paste0(base1, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

wd.ngs <- file.path(wd, BASE1, "ngs/WGS")
samples1 <- readTable(file.path(wd.ngs, "nbl_wgs_n57-1.txt"), header=T, rownames=T, sep="")
#samples1 <- readTable(file.path(wd.ngs, "nbl_wgs_n28.txt"), header=T, rownames=T, sep="")
samples1 <- subset(samples1, RT == 1)[,1]
samples0 <- readTable(file.path(wd.ngs, "nbl_wgs_n57-1.txt"), header=T, rownames=T, sep="")
#samples0 <- readTable(file.path(wd.ngs, "nbl_wgs_n28.txt"), header=T, rownames=T, sep="")
samples0 <- subset(samples0, RT == 0)[,1]
n1 <- length(samples1)
n0 <- length(samples0)

# -----------------------------------------------------------------------------
# Plot RD and RT (see ReplicationTiming.R)
# Last Modified: 14/02/19; 10/01/19; 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
for (c in 1:22) {
   chr <- chrs[c]
   
   rpkms.chr.rt <- readTable(file.path(wd.rt.data, paste0(base1, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, flipRT=F, scaledRT=T) 
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   bed.gc.chr <- bed.gc.chr[rpkms.chr.rt$BED,]

   ## Colours (was "lightcoral", "skyblue3")
   adjustcolor.red  <- adjustcolor("lightcoral", alpha.f=0.3)
   adjustcolor.blue <- adjustcolor("skyblue3", alpha.f=0.3)
   adjustcolor.gray <- adjustcolor("gray", alpha.f=0.3)
   
   ## RD 
   #ylab.text <- "Read depth"
   #file.name <- file.path(wd.rt.plots, paste0("RD_", base1, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "_n", n1))
   #main.text <- paste0("Read depth in ", BASE1, " tumour cells (n=", n1, ")")
   #plotRD(file.name, main.text, ylab.text, chr, NA, NA, rpkms.chr.rt$T, bed.gc.chr, c(adjustcolor.red, "red"), "png", 7.75, 9.25)   #(7.75, 9.25) for chr2
   
   #file.name <- file.path(wd.rt.plots, paste0("RD_", base0, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR0, "_n", n0))
   #main.text <- paste0("Read depth in ", BASE0, " normal cells (n=", n0, ")")
   #plotRD(file.name, main.text, ylab.text, chr, NA, NA, rpkms.chr.rt$N, bed.gc.chr, c(adjustcolor.blue, "blue"), "png", 7.75, 9.25)
   
   #file.name <- file.path(wd.rt.plots, paste0("RD_", base0, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "+", PAIR0, "_n", n1, "-", n0))
   #main.text <- paste0("Read depth in ", BASE1, " tumour (n=", n1, ") and normal (n=", n0, ") cells")
   #plotRD2(file.name, main.text, ylab.text, chr, NA, NA, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("Tumour", "Normal"), "png", 7.75, 9)

   ## RD & RT 
   main.text <- paste0(BASE1, " Q4&3/Q2&1 read depth ratio between tumour (n=", n1, ") and tumour (n=", n0, ") cells")   
   file.name <- file.path(wd.rt.plots, paste0("RTD_", base0, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0))   
   plotRD3(file.name, main.text, chr, NA, NA, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("Q4&3 tumour", "Q2&1 tumour"), c(adjustcolor.gray, adjustcolor.gray), c("Q4&3", "Q2&1"), "png", width=10, peaks=c(74353001, 85951001), 7.75, 9, 3, 3)
   plotRD3(file.name, paste0(BASE1, " Q4&3/Q2&1 read depth ratio"), chr, 71500000, 90500000, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("Q4&3", "Q2&1"), c(adjustcolor.red, adjustcolor.blue), c("Q4&3", "Q2&1"), "png", width=5, peaks=c(74353001, 85951001), 7.75, 9, 3, 3)
   
   ## RT
   #ylab.text <- "Replication timing"
   #file.name <- file.path(wd.rt.plots, paste0("RT_", base0, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0))
   #main.text <- paste0(BASE1, " T/N read depth ratio between tumour (n=", n1, ") and normal (n=", n0, ") cells")
   #plotRT(file.name, main.text, ylab.text, chr, NA, NA, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c(adjustcolor.gray, adjustcolor.gray), c("T", "N"), "png", 3, 3)
}

# -----------------------------------------------------------------------------
# NBL RD vs RT
# Last Modified: 18/02/19
# -----------------------------------------------------------------------------
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   
   #rpkms.chr.rt <-readTable(paste0("/Users/tpyang/Work/uni-koeln/tyang2/LCL/analysis/replication/lcl-wgs-rt/data/lcl_rpkm.corr.gc.d.rt_", chr, "_LCL-LCL_n7-7.txt.gz"), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- readTable(file.path(wd.rt.data, paste0(base1, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, flipRT=F, scaledRT=T) 

   rpkms.chr.rt.T  <- setSlope(rpkms.chr.rt, bed.gc.chr, "T")
   rpkms.chr.rt.N  <- setSlope(rpkms.chr.rt, bed.gc.chr, "N")
   rpkms.chr.rt.RT <- setSlope(rpkms.chr.rt, bed.gc.chr, "RT")
   
   main.text <- paste0("NBL RT (Q4&3/Q2&1) vs. Read depths (", "Chr2", ")")
   xlab.text <- "NBL RT (Q4&3/Q2&1)"
   ylab.text <- "NBL read depth"
   file.name <- file.path(wd.rt.plots, paste0("plot_RT-Q4&3-Q2&1-vs-RD_NBL_pearson_chr2"))
   xmin <- min(rpkms.chr.rt.RT$SLOPE)
   xmax <- max(rpkms.chr.rt.RT$SLOPE)
   ymin <- min(c(rpkms.chr.rt.T$SLOPE, rpkms.chr.rt.N$SLOPE))
   ymax <- max(c(rpkms.chr.rt.T$SLOPE, rpkms.chr.rt.N$SLOPE))
   plotRD2vsRT(rpkms.chr.rt.T$SLOPE, rpkms.chr.rt.N$SLOPE, rpkms.chr.rt.RT$SLOPE, file.name, main.text, ylab.text, xlab.text, c("red", "blue"), c("Q4&3", "Q2&1"), xmin, xmax, ymin, ymax)
}

# -----------------------------------------------------------------------------
# CLL RT vs SCLC/LCL RT
# Last Modified: 05/03/19; 18/02/19
# -----------------------------------------------------------------------------
cors <- toTable(0, 2, 22, c("chr", "cor"))
cors$chr <- 1:22
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
 
   ## CLL and SCLC/LCL RTs
   rpkms.chr.rt <- readTable(file.path(wd.rt.data, paste0(base1, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, flipRT=F, scaledRT=T)
   rpkms.chr.rt.RT <- setSlope(rpkms.chr.rt, bed.gc.chr, "RT")
 
   rpkms.chr.rt.lcl <-readTable(paste0("/Users/tpyang/Work/uni-koeln/tyang2/LCL/analysis/replication/lcl-wgs-rt/data/lcl_rpkm.corr.gc.d.rt_", chr, "_LCL-LCL_n7-7.txt.gz"), header=T, rownames=T, sep="\t")
   #rpkms.chr.rt.lcl <-readTable(paste0("/Users/tpyang/Work/uni-koeln/tyang2/SCLC/analysis/replication/sclc-wgs-rt/data/sclc_rpkm.corr.gc.d.rt_", chr, "_SCLC-SCLC_n101-92.txt.gz"), header=T, rownames=T, sep="\t")
   rpkms.chr.rt.lcl <- setScaledRT(rpkms.chr.rt.lcl, pseudocount=0.01, recaliRT=T, flipRT=F, scaledRT=T) 
   rpkms.chr.rt.lcl.RT <- setSlope(rpkms.chr.rt.lcl, bed.gc.chr, "RT")
 
   ## Keep 1kb slopes based on overlapping windows
   overlaps <- intersect(rpkms.chr.rt.RT$BED, rpkms.chr.rt.lcl.RT$BED)
   #rpkms.chr.rt.RT     <- rpkms.chr.rt.RT[overlaps,]       ## Too slow
   #rpkms.chr.rt.lcl.RT <- rpkms.chr.rt.lcl.RT[overlaps,]   ## Too slow
 
   ##
   #main.text <- paste0("NBL RT (Q4/Q1) vs. LCL RT (", "Chr", c, ")")
   #xlab.text <- "NBL RT (Q4/Q1)"
   #ylab.text <- "LCL RT (S/G1)"
   #file.name <- file.path(wd.rt.plots, paste0("plot_RT-vs-RT_CLL-N-T-vs-LCL_pearson_chr", c, "_n"))
   #plotRTvsRT(rpkms.chr.rt.lcl.RT[overlaps,]$SLOPE, rpkms.chr.rt.RT[overlaps,]$SLOPE, file.name, main.text, ylab.text, xlab.text, c(adjustcolor.gray, "black"))
 
   cors$cor[c] <- getCor(rpkms.chr.rt.lcl.RT[overlaps,]$SLOPE, rpkms.chr.rt.RT[overlaps,]$SLOPE)
}
save(cors, file=file.path(wd.rt.data, paste0("rt-vs-rt_", "nbl-Q4&3-Q2&1", "-vs-lcl_cors-pearson.RData")))

ylab.text <- "Pearson's r"
xlab.text <- "Chromosome"
file.name <- file.path(wd.rt.plots, "plot_RT-Q4&3-Q2&1-vs-RT_NBL-vs-LCL_pearson")
main.text <- paste0("NBL RT (Q4&3/Q2&1) vs. LCL RT")
ymin <- 0
ymax <- 0.8
plotRTvsRTALL(cors, file.name, main.text, ylab.text, xlab.text, ymin, ymax, line0=T)







cors <- toTable(0, 2, 22, c("chr", "cor"))
cors$chr <- 1:22
for (c in 1:22) {
   chr <- chrs[c]
 
   rpkms.chr.rt <- readTable(file.path(wd.rt.data, paste0(base1, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
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
 
   cors$cor[c] <- getCor(rpkms.chr.rt.lcl.RT$SLOPE, rpkms.chr.rt.RT$SLOPE)
}
ylab.text <- "Pearson's r"
xlab.text <- "Chromosome"
file.name <- file.path(wd.rt.plots, "plot_RT-vs-RT_NBL&SCLC-vs-LCL_pearson_All")
main.text <- paste0("NBL-T/SCLC-N RT vs. LCL RT (All)")
ymin <- 0
ymax <- 0.8
plotRTvsRTALL(cors, file.name, main.text, ylab.text, xlab.text, ymin, ymax)
save(cors, file=file.path(wd.rt.data, paste0("rt-vs-rt_", base1, "&", base0, "-vs-lcl_cors-pearson.RData")))

###
## NBL RT vs LCL RT   ## TO-DO
BASE0 <- "SCLC"
base0 <- tolower(BASE0)
n0 <- 9

cors <- toTable(0, 2, 22, c("chr", "cor"))
cors$chr <- 1:22
for (c in 1:22) {
   chr <- chrs[c]

   rpkms.chr.rt <- readTable(file.path(wd.rt.data, "hybrid", paste0(base1, "_rpkm.corr.gc.d.rt_", chr, "_", BASE1, "-", BASE0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, flipRT=F, scaledRT=T) 
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   bed.gc.chr <- bed.gc.chr[rpkms.chr.rt$BED,]
   rpkms.chr.rt.RT     <- setSlope(rpkms.chr.rt, bed.gc.chr, "RT")
   
   rpkms.chr.rt.lcl <-readTable(paste0("/Users/tpyang/Work/uni-koeln/tyang2/LCL/analysis/replication/lcl-wgs-rt/data/lcl_rpkm.corr.gc.d.rt_", chr, "_LCL-LCL_n7-7.txt.gz"), header=T, rownames=T, sep="\t")
   rpkms.chr.rt.lcl <- setScaledRT(rpkms.chr.rt.lcl, pseudocount=0.01, recaliRT=T, flipRT=F, scaledRT=T) 
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   bed.gc.chr <- bed.gc.chr[rpkms.chr.rt.lcl$BED,]
   rpkms.chr.rt.lcl.RT <- setSlope(rpkms.chr.rt.lcl, bed.gc.chr, "RT")

   overlaps <- intersect(rpkms.chr.rt.RT$BED, rpkms.chr.rt.lcl.RT$BED)
   rpkms.chr.rt.RT     <- rpkms.chr.rt.RT[overlaps,]
   rpkms.chr.rt.lcl.RT <- rpkms.chr.rt.lcl.RT[overlaps,]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   bed.gc.chr <- bed.gc.chr[overlaps,]
   
   #ylab.text <- "LCL RT (S/G1)"
   #xlab.text <- "NBL RT (T/N)"
   #file.name <- file.path(wd.rt.plots, paste0("plot_RT-vs-RT_NBL-vs-LCL_pearson_chr", c))
   #main.text <- paste0("NBL RT vs. LCL RT (", "Chr", c, ")")
   #cors$cor[c] <- plotRDvsRT(rpkms.chr.rt.lcl.RT$SLOPE, rpkms.chr.rt.RT$SLOPE, file.name, main.text, ylab.text, xlab.text, c(adjustcolor.gray, "black"))
   cors$cor[c] <- getCor(rpkms.chr.rt.lcl.RT$SLOPE, rpkms.chr.rt.RT$SLOPE)
}
save(cors, file=file.path(wd.rt.data, paste0("rt-vs-rt_", base1, "&", base0, "-B-vs-lcl_cors-pearson.RData")))

ylab.text <- "Pearson's r"
xlab.text <- "Chromosome"
file.name <- file.path(wd.rt.plots, "plot_RT-vs-RT_NBL&SCLC-B-vs-LCL_pearson_All")
main.text <- paste0("NBL–SCLC T/N vs. LCL RT")
ymin <- -0.2
ymax <- 0.85
plotRTvsRTALL(cors, file.name, main.text, ylab.text, xlab.text, ymin, ymax, line0=T)

# -----------------------------------------------------------------------------
# 
# http://pklab.med.harvard.edu/scw2014/subpop_tutorial.html
# Last Modified: 06/03/19
# -----------------------------------------------------------------------------
cors.samples <- toTable(0, length(samples1)+4, 22, c("chr", "mean", "var", "cv2", samples1))
cors.samples$chr <- 1:22
for (s in 1:length(samples1)) {
   sample <- samples1[s]
   load(file.path(wd.rt.data, "samples", paste0("rd-vs-rt_", sample, "-vs-lcl_cors-pearson.RData")))
 
   cors.samples[, sample] <- cors$cor
}

for (c in 1:22) {
   cors.samples$mean[c] <- mean(as.numeric(cors.samples[c, samples1]))
   cors.samples$var[c]  <- var(as.numeric(cors.samples[c, samples1]))
   cors.samples$cv2[c]  <- cors.samples$var[c]/cors.samples$mean[c]^2
}
save(cors.samples, file=file.path(wd.rt.data, paste0("rd-vs-rt_samples-vs-lcl_cors-pearson.RData")))

file.name <- file.path(wd.rt.plots, "plot_RD-vs-RT_SAMPLES-vs-LCL_pearson")
main.text <- c("NBL (n=56) RD (T) vs. LCL RT", "NBL RD (T) vs. LCL RT")   ## TO-DO
ymin <- -0.6920546
ymax <- 0.6621544
plotSAMPLEvsRTALL(cors.samples, samples1, file.name, main.text, ymin, ymax, line0=T)

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
# -0.2669331 -0.2217933 -0.2002513 -0.1150665  0.4184586 

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
writeTable(cors.all, file.path(wd.ngs, "nbl_wgs_n57-1.txt"), colnames=T, rownames=F, sep="\t")

cors.all <- subset(cors.all, Q4 %in% c(4,1))
writeTable(cors.all, file.path(wd.ngs, "nbl_wgs_n28.txt"), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# PCA
# Last Modified: 21/04/19
# -----------------------------------------------------------------------------
samples  <- readTable(file.path(wd.ngs, "nbl_wgs_n57-1.txt"), header=T, rownames=T, sep="")

trait <- as.numeric(samples$Q4)
trait[which(trait == 4)] <- "Q4 (-0.12 < r < 0.42)"
trait[which(trait == 3)] <- "Q3 (-0.20 < r < -0.12)"
trait[which(trait == 2)] <- "Q2 (-0.22 < r < -0.20)"
trait[which(trait == 1)] <- "Q1 (-0.27 < r < -0.22)"

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
 
   file.main <- c(paste0("NBL tumours on chr", c), "")
   plotPCA(1, 2, pca.de, trait, wd.rt.plots, paste0("pca_nbl_T_chr", c), size=6.5, file.main, "topleft", c("blue", "skyblue3", "lightcoral", "red"), NULL, flip.x=1, flip.y=1)
   save(pca.de, file=file.path(wd.rt.data, paste0("pca_nbl_T_chr", c, ".RData")))
 
   if (is.na(rpkms.T.chr.d.all)) {
      rpkms.T.chr.d.all <- test
   } else
      rpkms.T.chr.d.all <- rbind(rpkms.T.chr.d.all, test)
}

##
test <- rpkms.T.chr.d.all
pca.de <- getPCA(t(test))

file.main <- c("NBL tumours on all chromosomes", "")
plotPCA(1, 2, pca.de, trait, wd.rt.plots, "pca_nbl_T_chrs", size=6.5, file.main, "topleft", c("blue", "skyblue3", "lightcoral", "red"), NULL, flip.x=1, flip.y=1)
save(pca.de, file=file.path(wd.rt.data, paste0("pca_nbl_T_chrs.RData")))





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
# [1] 15
# > length(g1_likes)
# [1] 4
# > s_likes
# [1] "P21702" "P21924" "P22496" "P23103" "P23206" "P23229" "P23267" "P24478" "P24632" "P24679" "P24702" "P24885" "P25114" "P25262" "P25376" 
# > g1_likes
# [1] "P1695"  "P20471" "P21442" "P22388"











# -----------------------------------------------------------------------------
# 
# Last Modified: 06/03/19
# -----------------------------------------------------------------------------
# > length(which(as.numeric(cors.samples[2, -c(1:4)]) > 0))
# [1] 13
# > length(which(as.numeric(cors.samples[2, -c(1:4)]) < 0))
# [1] 43

idx.pos.chrs <- which(as.numeric(cors.samples[1, -c(1:4)]) > 0) + 4
idx.neg.chrs <- which(as.numeric(cors.samples[1, -c(1:4)]) < 0) + 4
for (c in 2:22) {
   idx.pos.chr <- which(as.numeric(cors.samples[c, -c(1:4)]) > 0) + 4
   idx.neg.chr <- which(as.numeric(cors.samples[c, -c(1:4)]) < 0) + 4
 
   idx.pos.chrs <- intersect(idx.pos.chrs, idx.pos.chr)
   idx.neg.chrs <- intersect(idx.neg.chrs, idx.neg.chr)
}
# > length(idx.pos.chrs)
# [1] 7
# > length(idx.neg.chrs)
# [1] 33
idx.var.chrs <- c(5:60)
idx.var.chrs <- setdiff(idx.var.chrs, c(idx.pos.chrs, idx.neg.chrs))
# > length(idx.var.chrs)
# [1] 16
# > 7+33+16
# [1] 56

samples.pos.chr5 <- colnames(cors.samples[, which(as.numeric(cors.samples[5, -c(1:4)]) > 0) + 4])
samples.pos <- colnames(cors.samples[, idx.pos.chrs])
samples.neg <- colnames(cors.samples[, idx.neg.chrs])
samples.var <- colnames(cors.samples[, idx.var.chrs])


