# =============================================================================
# Manuscript   : 
# Chapter      : Chromosome replication timing of the human genome
# Name         : manuscripts/replication/sclc-wgs-rt-normal.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 24/10/19; 14/10/19; 05/03/19; 25/02/19; 30/01/18
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "DifferentialExpression.R", "ReplicationTiming.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.bed.gc.1kb.RData"))
load(file.path(wd.src.ref, "hg19.lcl.koren.woodfine.RData"))

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 30/01/18
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
#wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE  <- "SCLC"
PAIR1 <- "N"
base  <- tolower(BASE)
method <- "rpkm"

wd.ngs   <- file.path(wd, BASE, "ngs/WGS")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.rt       <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt-normal"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

#wd.ngs.data <- file.path(wd.ngs, "data")
wd.ngs.data <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"), "data")
samples1 <- readTable(file.path(wd.ngs, "sclc_wgs_n92_NL.list"), header=F, rownames=F, sep="")
n1 <- length(samples1)
## samples1 are modified to the normals in this script

# -----------------------------------------------------------------------------
# SCLC N vs LCL S/G1
# http://pklab.med.harvard.edu/scw2014/subpop_tutorial.html
# Last Modified: 24/10/19; 05/06/19; 20/04/19; 06/03/19
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
save(cors.samples, file=file.path(wd.rt.data, paste0("samples-vs-rt_sclc-tn-vs-lcl_spline_spearman.RData")))
# > min(cors.samples[,-c(1:4)])
# [1] -0.8553739
# > max(cors.samples[,-c(1:4)])
# [1] 0.790692

#load(file.path(wd.rt.data, paste0("samples-vs-rt_sclc-vs-lcl_spline_spearman.RData")))
file.name <- file.path(wd.rt.plots, "SAMPLES-vs-RT_SCLC-TN-vs-LCL_spline_spearman")
main.text <- c("SCLC-TN read depth vs. LCL RT", "")
ymin <- -0.8773492
ymax <- 0.8392611
plotSAMPLEvsRTALL(cors.samples, samples1, file.name, main.text, ymin, ymax)

# -----------------------------------------------------------------------------
# Overall correlation with LCL S/G1
# Last Modified: 16/06/19; 04/06/19; 06/03/19
# -----------------------------------------------------------------------------
samples.sclc.nl <- setSamplesQ4(wd.rt.data, samples1)
writeTable(samples.sclc.nl, file.path(wd.ngs, "sclc_wgs_n92_NL.txt"), colnames=T, rownames=F, sep="\t")
#         0%        25%        50%        75%       100% 
# -0.7585423 -0.7439594 -0.7293444 -0.6681727  0.6753625 

writeTable(subset(samples.sclc.nl, Q4 %in% c(4,1)), file.path(wd.ngs, "sclc_wgs_n92_NL_q4.txt"), colnames=T, rownames=F, sep="\t")

## Random 23/23
m2.23 <- sort(rownames(subset(samples.sclc.tn, M2 == 1))[sample(1:46, round(46/2), replace=F)])
m1.23 <- sort(rownames(subset(samples.sclc.tn, M2 == 0))[sample(1:46, round(46/2), replace=F)])

random1 <- sort(c(m2.23, m1.23))
writeTable(samples.sclc.tn[random1,], file.path(wd.ngs, "sclc_wgs_n92_TN-random1.txt"), colnames=T, rownames=F, sep="\t")

random2 <- sort(setdiff(rownames(samples.sclc.tn), random1))
writeTable(samples.sclc.tn[random2,], file.path(wd.ngs, "sclc_wgs_n92_TN-random2.txt"), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# PCA
# Last Modified: 04/06/19; 21/04/19
# -----------------------------------------------------------------------------
## Refer to cmd-rt_2a_nrd.gc.cn.d_sample.R (commandline mode)
nrds.T.chr.d.all <- readTable(file.path(wd.ngs.data, paste0(base, "_", method, ".gc.cn.d_", "chr1", "_", PAIR1, "_n", n1, ".txt.gz")), header=T, rownames=F, sep="\t")
for (c in 2:22) {
   chr <- chrs[c]
 
   ## Read depth
   #nrds.T.chr.d <- pipeGetDetectedRD(wd.ngs.data, BASE, chr, PAIR1, method)
   nrds.T.chr.d <- readTable(file.path(wd.ngs.data, paste0(base, "_", method, ".gc.cn.d_", chr, "_", PAIR1, "_n", n1, ".txt.gz")), header=T, rownames=F, sep="\t")
 
   nrds.T.chr.d.all <- rbind(nrds.T.chr.d.all, nrds.T.chr.d)
}

##
test <- nrds.T.chr.d.all[, -1]   ## BUG 2019/10/14: Remove column BED
pca.de <- getPCA(t(test))
save(pca.de, file=file.path(wd.rt.data, paste0("pca_sclc-nl_chrs.RData")))

#load(file.path(wd.rt.data, paste0("pca_sclc-nl_chrs.RData")))
file.main <- c("SCLC-NL", "")
trait <- as.numeric(samples.sclc.nl$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
plotPCA(1, 2, pca.de, trait, wd.rt.plots, "PCA_SCLC-NL", size=6, file.main, "bottomright", c(red, red.lighter, blue.lighter, blue), NULL, flip.x=1, flip.y=1, legend.title=NA)

## SG1
#trait <- samples.sclc.sg1$SG1
#plotPCA(1, 2, pca.de, trait, wd.rt.plots, "pca_sclc_T_chrs_spline_spearman_SG1", size=6, file.main, "bottomright", c("red", "lightgray", "blue"), NULL, flip.x=1, flip.y=1, legend.title="Consist. CM in all chrs")

# -----------------------------------------------------------------------------
# Beeswarm plots
# Last Modified: 21/04/19
# -----------------------------------------------------------------------------
samples.sclc.nl <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/SCLC/ngs/WGS/sclc_wgs_n92_NL.txt", header=T, rownames=T, sep="")
samples.nbl.wb  <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/NBL/ngs/WGS/nbl_wgs_n57-1_N.txt", header=T, rownames=T, sep="")
samples.cll.wb  <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/CLL/ngs/WGS/cll_wgs_n96_N.txt", header=T, rownames=T, sep="")
n.sclc.nl <- nrow(samples.sclc.nl)
n.nbl.wb  <- nrow(samples.nbl.wb)
n.cll.wb  <- nrow(samples.cll.wb)

samples <- toTable(0, 3, n.sclc.nl+n.nbl.wb+n.cll.wb, c("CANCER", "COR", "Q4"))
samples$CANCER[1:n.sclc.nl] <- 0
samples$CANCER[(1+n.sclc.nl):(n.sclc.nl+n.nbl.wb)] <- 1
samples$CANCER[(1+n.sclc.nl+n.nbl.wb):(n.sclc.nl+n.nbl.wb+n.cll.wb)] <- 2
samples$COR <- c(samples.sclc.nl$COR, samples.nbl.wb$COR, samples.cll.wb$COR)
samples$Q4  <- c(samples.sclc.nl$Q4, samples.nbl.wb$Q4, samples.cll.wb$Q4)

#install.packages('beeswarm')
library(beeswarm)

pdf(file.path(wd.rt.plots, "beeswarm_sclc+nbl+cll_normal.pdf"), height=6, width=6)
ymax <- 0.8   #max(samples$COR)
ymin <- -ymax
boxplot(COR ~ CANCER, data=samples, outline=F, names=c("", "", ""), ylim=c(ymin, ymax), ylab="", main="Normal read depth vs. RT", yaxt="n", cex.axis=1.8, cex.lab=1.9, cex.main=2.1)
abline(h=0, lty=5, lwd=2)

legend("topright", legend=c("Q4", "Q3", "Q2", "Q1"), pch=16, pt.cex=2.5, col=c(red, red.lighter, blue.lighter, blue), cex=1.9)

beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 1), col=blue, pch=19, cex=1, add=T)
beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 2), col=blue.lighter, pch=19, cex=1, add=T)
beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 3), col=red.lighter, pch=19, cex=1, add=T)
beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 4), col=red, pch=19, cex=1, add=T)

axis(side=2, at=seq(-0.8, 0.8, by=0.4), labels=c(-0.8, -0.4, 0, 0.4, 0.8), cex.axis=1.8)
mtext("Spearman's rho", side=2, line=2.7, cex=1.9)
mtext("", cex=1.2, line=0.3)
mtext(text=c("SCLC-NL", "NBL-WB", "CLL-WB"), side=1, cex=1.9, line=1.3, at=c(1,2,3))
mtext(text=c("n=92", "n=56", "n=96"), side=1, cex=1.8, line=3, at=c(1,2,3))
dev.off()

## 
testU(samples.sclc.nl$COR, samples.sclc$COR)
# [1] 5.601622e-10
testU(samples.sclc.nl$COR, samples.nbl$COR)
# [1] 0.5027551
testU(samples.sclc.nl$COR, samples.cll$COR)
# [1] 2.72587e-15

## 
testU(samples.sclc$COR, samples.nbl$COR)
# [1] 4.788974e-05
testU(samples.sclc$COR, samples.cll$COR)
# [1] 9.387524e-31









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
