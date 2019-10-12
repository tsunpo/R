# =============================================================================
# Manuscript   : 
# Chapter      : Chromosome replication timing of the human genome
# Name         : manuscripts/replication/sclc-wgs-rt.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 05/03/19; 25/02/19; 30/01/18
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "ReplicationTiming.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.1kb.gc.RData"))
load(file.path(wd.src.ref, "hg19.rt.lcl.koren.woodfine.RData"))

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 30/01/18
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
#wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE  <- "SCLC"
PAIR1 <- "T"
PAIR0 <- "N"
base  <- tolower(BASE)
method <- "rpkm"

wd.ngs   <- file.path(wd, BASE, "ngs/WGS")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.rt       <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

wd.ngs.data <- file.path(wd.ngs, "data")
samples1 <- readTable(file.path(wd.ngs, "sclc_wgs_n101.list"), header=F, rownames=F, sep="")
samples0 <- readTable(file.path(wd.ngs, "sclc_wgs_n92_N.list"), header=F, rownames=F, sep="")
#samples0 <- readTable(file.path(wd.ngs, "sclc_wgs_n9_B.list"), header=F, rownames=F, sep="")
n1 <- length(samples1)
n0 <- length(samples0)

# -----------------------------------------------------------------------------
# Plot RD and RT (see ReplicationTiming.R)
# Last Modified: 28/05/19; 14/02/19; 10/01/19; 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
nrds <- getLog2ScaledRT(wd.rt.data, base, method, PAIR1, PAIR0, chrs, bed.gc)
save(nrds, file=file.path(wd.rt.data, paste0("nrds_", base, "-t-n_", method, ".log2s.RData")))
# > nrow(nrds)
# [1] 2656170 - 22

#load(file.path(wd.rt.data, paste0("nrds_", base, "-t-n_", method, ".log2s.RData")))
ymax <- 0.6
ymin <- 0.14
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]
   lcl.rt.chr <- subset(lcl.rt, CHR == chr)   ## Koren 2012
 
   ## Plot RT
   main.text <- paste0(BASE, " T/N read depth ratio between tumour (n=", n1, ") and normal (n=", n0, ") samples")  
   file.name <- file.path(wd.rt.plots, paste0("RT_", base, "_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ""))   
   plotRT(file.name, main.text, chr, NA, NA, nrds.chr, bed.gc.chr, c("red", "blue"), c("Tumour", "Normal"), c("lightcoral", "lightskyblue3"), c("T", "N"), "png", width=10, peaks=c(), ylim=c(ymin, ymax), lcl.rt.chr)
}

# -----------------------------------------------------------------------------
# RD vs RT (RDS and SPR)
# Last Modified: 11/07/19; 31/05/19
# -----------------------------------------------------------------------------
cors <- toTable(0, 7, 22, c("chr", "cor", "cor1", "cor2", "e", "l", "spr"))
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
   
   cors$cor[c]  <- getCor(nrds.chr.T$SPLINE, nrds.chr.N$SPLINE,  method="spearman")
   cors$cor1[c] <- getCor(nrds.chr.T$SPLINE, nrds.chr.RT$SPLINE, method="spearman")
   cors$cor2[c] <- getCor(nrds.chr.N$SPLINE, nrds.chr.RT$SPLINE, method="spearman")
   
   e <- nrow(subset(nrds.chr.RT, SPLINE > 0))
   l <- nrow(subset(nrds.chr.RT, SPLINE < 0))
   cors$e[c] <- e
   cors$l[c] <- l
   cors$spr[c] <- (e - l)/(e + l)
   
   main.text <- c(paste0("SCLC read depth correlation (", "Chr", c, ")"), paste0("rho = ", round0(cors$cor[c], digits=2), " (T vs. N)"))
   xlab.text <- "T/N read depth ratio [log2]"
   ylab.text <- "Read depth [RPKM]"
   file.name <- file.path(wd.rt.plots, "chrs", paste0("RD-vs-RT_SCLC-T-N_chr", c, "_spline_spearman"))
   plotRD2vsRT(nrds.chr.T$SPLINE, nrds.chr.N$SPLINE, nrds.chr.RT$SPLINE, file.name, main.text, ylab.text, xlab.text, c("red", "blue"), c("T", "N"), method="spearman")
}
save(cors, file=file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-t-n_spline_spearman.RData")))
writeTable(cors, file=file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-t-n_spline_spearman.txt")), colnames=T, rownames=F, sep="\t")

## S-phase progression rate (SPR)
ylab.text <- "SPR"
file.name <- file.path(wd.rt.plots, "SPR_SCLC-T-N_spline_spearman")
main.text <- c(paste0(BASE, " S-phase progression rate"), "SPR = (E-L)/(E+L)")
plotSPR(cors, file.name, main.text, c(13, 17), digits=3, unit=5, ylab.text)

## SPR vs Read depth correlation
file.name <- file.path(wd.rt.plots, "SPR-RDC_SCLC-T-N_spline_spearman")
main.text <- c(paste0(BASE, " SPR vs. Read depths correlation"), "")
xlab.text <- "T vs. N [rho]"
plotSPRRDC(cors$spr, cors$cor, file.name, main.text, c(4, 13, 17, 19, 22), xlab.text, unit=5, ylab.text)

## SPR vs Woodfine 2004
file.name <- file.path(wd.rt.plots, "SPR-Woodfine_SCLC-T-N_spline_spearman")
main.text <- c(paste0(BASE, " SPR vs. Woodfine et al. 2004"), "")
xlab.text <- "Mean replication timing ratio"
plotSPRRDC(cors$spr, lcl.mean$Mean, file.name, main.text, c(4, 13, 17, 19, 22), xlab.text, unit=5, ylab.text)






# -----------------------------------------------------------------------------
# SCLC T/N vs LCL S/G1
# Last Modified: 27/05/19
# -----------------------------------------------------------------------------
## LCL S/G1
nrds.tmp <- nrds
load(file.path(wd, "LCL/analysis/replication/lcl-wgs-rt/data", paste0("nrds_lcl-s-g1_", method, ".RData")))
nrds.lcl <- nrds
nrds <- nrds.tmp

cors <- toTable(0, 3, 22, c("chr", "length", "cor"))
cors$chr <- 1:22
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]
   nrds.chr.RT <- setSpline(nrds.chr, bed.gc.chr, "RT")
   
   nrds.lcl.chr <- nrds.lcl[intersect(nrds.lcl$BED, rownames(bed.gc.chr)),]  ## Reference LCL S/G1 ratio
   nrds.lcl.chr.RT <- setSpline(nrds.lcl.chr, bed.gc.chr, "RT")
   
   ## Keep only overlapping 1kb windows
   overlaps <- intersect(nrds.chr.RT$BED, nrds.lcl.chr.RT$BED)
   cors$length[c] <- length(overlaps)
   cors$cor[c] <- getCor(nrds.chr.RT[overlaps,]$SPLINE, nrds.lcl.chr.RT[overlaps,]$SPLINE, method="spearman")
}
save(cors, file=file.path(wd.rt.data, paste0("rt-vs-rt_", base, "-t-n-vs-lcl-s-g1_spline_spearman.RData")))

ylab.text <- "Spearman's rho"
xlab.text <- "Chromosome"
file.name <- file.path(wd.rt.plots, "RT-vs-RT_SCLC-T-N-vs-LCL-S-G1_spline_spearman")
main.text <- paste0("SCLC T/N vs. LCL S/G1")
ymin <- 0.25
ymax <- 1.05
plotRTvsRTALL(cors, file.name, main.text, ylab.text, xlab.text, ymin, ymax, col="black", c=2, pos=1)

# -----------------------------------------------------------------------------
# SCLC T vs LCL S/G1
# http://pklab.med.harvard.edu/scw2014/subpop_tutorial.html
# Last Modified: 05/06/19; 20/04/19; 06/03/19
# -----------------------------------------------------------------------------
cors.samples <- toTable(0, length(samples1)+4, 22, c("chr", "mean", "var", "cv2", samples1))
cors.samples$chr <- 1:22
for (s in 1:length(samples1)) {
   sample <- samples1[s]
   load(file.path(wd.rt.data, "samples", paste0("rd-vs-rt_", sample, "-vs-lcl_spline_spearman.RData")))
 
   cors.samples[, sample] <- cors$cor
}

for (c in 1:22) {
   cors.samples$mean[c] <- mean(as.numeric(cors.samples[c, samples1]))
   cors.samples$var[c]  <- var(as.numeric(cors.samples[c, samples1]))
   cors.samples$cv2[c]  <- cors.samples$var[c]/cors.samples$mean[c]^2
}
save(cors.samples, file=file.path(wd.rt.data, paste0("samples-vs-rt_sclc-vs-lcl_spline_spearman.RData")))
# > min(cors.samples[,-c(1:4)])
# [1] -0.8466193
# > max(cors.samples[,-c(1:4)])
# [1] 0.8112673

#load(file.path(wd.rt.data, paste0("samples-vs-rt_sclc-vs-lcl_spline_spearman.RData")))
file.name <- file.path(wd.rt.plots, "SAMPLES-vs-RT_SCLC-vs-LCL_spline_spearman")
main.text <- c("SCLC (n=101) read depth vs. LCL S/G1", "")
ymin <- -0.8789273
ymax <- 0.8433154
plotSAMPLEvsRTALL(cors.samples, samples1, file.name, main.text, ymin, ymax)

# -----------------------------------------------------------------------------
# Overall correlation with LCL S/G1
# Last Modified: 16/06/19; 04/06/19; 06/03/19
# -----------------------------------------------------------------------------
samples.sclc <- setSamplesQ4(wd.rt.data, samples1)
writeTable(samples.sclc, file.path(wd.ngs, "sclc_wgs_n101.txt"), colnames=T, rownames=F, sep="\t")
# 0%        25%        50%        75%       100% 
# -0.7504529 -0.6907507 -0.6546286 -0.5624756  0.6811825

writeTable(subset(samples.sclc, Q4 %in% c(4,1)), file.path(wd.ngs, "sclc_wgs_q4_n51.txt"), colnames=T, rownames=F, sep="\t")
writeTable(subset(samples.sclc, Q4 %in% c(3,1)), file.path(wd.ngs, "sclc_wgs_q3_n51.txt"), colnames=T, rownames=F, sep="\t")
writeTable(subset(samples.sclc, Q4 %in% c(2,1)), file.path(wd.ngs, "sclc_wgs_q2_n51.txt"), colnames=T, rownames=F, sep="\t")
#samples.sclc <- setSamplesQ4(wd.rt.data, overlaps)
#writeTable(samples.sclc, file.path(wd.ngs, "sclc_wgs_n70.txt"), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# PCA
# Last Modified: 04/06/19; 21/04/19
# -----------------------------------------------------------------------------
## Copy from 2a_cmd-rt_rpkm.corr.gc.d_sample.R (commandline mode)
nrds.T.chr.d.all <- NULL
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
 
   ## Read depth
   nrds.T.chr.d <- pipeGetDetectedRD(wd.ngs.data, BASE, chr, PAIR1, method)
   #nrds.N.chr.d <- pipeGetDetectedRD(wd0.ngs.data, BASE, chr, PAIR0) 
   #overlaps <- intersect(nrds.T.chr.d$BED, nrds.N.chr.d$BED)
 
   test <- nrds.T.chr.d[, samples.sclc$SAMPLE_ID]
   if (is.null(nrds.T.chr.d.all)) {
      nrds.T.chr.d.all <- test
   } else
      nrds.T.chr.d.all <- rbind(nrds.T.chr.d.all, test)
}

##
test <- nrds.T.chr.d.all
pca.de <- getPCA(t(test))
save(pca.de, file=file.path(wd.rt.data, paste0("pca_sclc_chrs_spline_spearman.RData")))

#load(file.path(wd.rt.data, paste0("pca_sclc_chrs_spline_spearman.RData")))
file.main <- c("SCLC (n=101) read depth profiles", "")
trait <- as.numeric(samples.sclc$Q4)
trait[which(trait == 4)] <- "Q4"   ##(-0.56 < rho < 0.68)"
trait[which(trait == 3)] <- "Q3"   ##(-0.65 < rho < -0.56)"
trait[which(trait == 2)] <- "Q2"   ##(-0.69 < rho < -0.65)"
trait[which(trait == 1)] <- "Q1"   ##(-0.75 < rho < -0.69)"
plotPCA(1, 2, pca.de, trait, wd.rt.plots, "PCA_SCLC_chrs_spline_spearman", size=6, file.main, "bottomright", c("red", "lightcoral", "lightskyblue3", "blue"), NULL, flip.x=1, flip.y=1, legend.title=NA)

## SG1
#trait <- samples.sclc.sg1$SG1
#plotPCA(1, 2, pca.de, trait, wd.rt.plots, "pca_sclc_T_chrs_spline_spearman_SG1", size=6, file.main, "bottomright", c("red", "lightgray", "blue"), NULL, flip.x=1, flip.y=1, legend.title="Consist. CM in all chrs")

# -----------------------------------------------------------------------------
# Beeswarm plots
# Last Modified: 21/04/19
# -----------------------------------------------------------------------------
samples.sclc <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/SCLC/ngs/WGS/sclc_wgs_n101.txt", header=T, rownames=T, sep="")
samples.nbl  <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/NBL/ngs/WGS/nbl_wgs_n57-1.txt", header=T, rownames=T, sep="")
samples.cll  <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/CLL/ngs/WGS/cll_wgs_n96.txt", header=T, rownames=T, sep="")
n.sclc <- nrow(samples.sclc)
n.nbl  <- nrow(samples.nbl)
n.cll  <- nrow(samples.cll)

samples <- toTable(0, 3, n.cll+n.sclc+n.nbl, c("CANCER", "COR", "Q4"))
samples$CANCER[1:n.sclc] <- 0
samples$CANCER[(1+n.sclc):(n.sclc+n.nbl)] <- 1
samples$CANCER[(1+n.sclc+n.nbl):(n.sclc+n.nbl+n.cll)] <- 2
samples$COR <- c(samples.sclc$COR, samples.nbl$COR, samples.cll$COR)
samples$Q4  <- c(samples.sclc$Q4, samples.nbl$Q4, samples.cll$Q4)

#install.packages('beeswarm')
library(beeswarm)

pdf(file.path(wd.rt.plots, "beeswarm_sclc+nbl+cll.pdf"), height=6, width=6)
ymax <- max(samples$COR)
ymin <- -ymax
boxplot(COR ~ CANCER, data=samples, outline=F, names=c("SCLC", "NBL", "CLL"), ylim=c(ymin, ymax), ylab="Spearman's rho", main="Overall correlation with LCL S/G1", cex.lab=1.5, cex.axis=1.4, cex.main=1.6)
abline(h=0, lty=5)

beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 1), col="blue", pch=16, add=T)
beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 2), col="skyblue3", pch=16, add=T)
beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 3), col="lightcoral", pch=16, add=T)
beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 4), col="red", pch=16, add=T)

legend("topright", legend = c("Q4", "Q3", "Q2", "Q1"), pch=16, col=c("red", "lightcoral", "skyblue3", "blue"), cex=1.5)
mtext("", cex=1.2, line=0.3) 
dev.off()










# -----------------------------------------------------------------------------
# CM2 and CQ4
# Last Modified: 19/05/19
# -----------------------------------------------------------------------------
cors <- t(cors.samples[, -c(1:4)])
#cors <- cors[overlaps,]   ## Overlapping between WGS and RNA (n=70); after sclc-tpm.de.R
colnames(cors) <- paste0("chr", c(1:22))
cq4 <- cors
cm2 <- cors

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
# Find S-like and G1-like tumour samples
# Last Modified: 04/06/19; 06/03/19
# -----------------------------------------------------------------------------
samples.sclc.sg1 <- setSamplesSG1(wd.rt.data, samples1, cors.samples)
writeTable(samples.sclc.sg1, file.path(wd.ngs, "sclc_wgs_n101.sg1"), colnames=T, rownames=F, sep="\t")
# > length(s_likes)
# [1] 24
# > length(g1_likes)
# [1] 14
# > s_likes
# [1] "S00339" "S00825" "S00827" "S00833" "S00837" "S00838" "S00938" "S00941" "S01020" "S01022" "S01366" "S01524" "S01864" "S02120" "S02219" "S02285" "S02286" "S02287" "S02289" "S02290" "S02291" "S02292"
# [23] "S02295" "S02296"
# > g1_likes
# [1] "S00050" "S00472" "S00831" "S00944" "S01861" "S02139" "S02237" "S02241" "S02248" "S02342" "S02344" "S02352" "S02353" "S02378"
