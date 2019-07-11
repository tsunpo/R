# =============================================================================
# Manuscript   : Cell cycle status in the neuroblastoma cell lines
# Chapter      : Chromosome replication timing of the human genome
# Name         : manuscripts/replication/nbl-cl-rt.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 04/07/19
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
# Last Modified: 04/07/19
# -----------------------------------------------------------------------------
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE  <- "NBL"
PAIR1 <- "T"
base  <- tolower(BASE)

wd.ngs   <- file.path(wd, BASE, "ngs/CL")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.rt       <- file.path(wd.anlys, "replication", paste0(base, "-cl-rt"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

wd.ngs.data <- file.path(wd.ngs, "data")
samples1 <- readTable(file.path(wd.ngs, "nbl_cl_n8.list"), header=F, rownames=F, sep="")
samples1 <- toupper(gsub("-", "", samples1)) 
n1 <- length(samples1)

# -----------------------------------------------------------------------------
# NBL CL vs LCL S/G1
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
save(cors.samples, file=file.path(wd.rt.data, paste0("samples-vs-rt_nbl-cl-vs-lcl_spline_spearman.RData")))

# > min(cors.samples[,-c(1:4)])   ## RT
# [1] -0.8373368
# > max(cors.samples[,-c(1:4)])
# [1] 0.8643419
file.name <- file.path(wd.rt.plots, "boxplot_SAMPLES-vs-RT_NBL-CL-vs-LCL_spline_spearman")
main.text <- c("NBL CL (n=8) vs. LCL S/G1", "Chromosomal median (CM2/CM1)")
ymin <- -0.8732989
ymax <- 0.8643419
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
writeTable(cq4[, c("SAMPLE_ID", paste0("chr", c(1:22)))], file.path(wd.ngs, "nbl_cl_n8.cq4"), colnames=T, rownames=F, sep="\t")
writeTable(cm2[, c("SAMPLE_ID", paste0("chr", c(1:22)))], file.path(wd.ngs, "nbl_cl_n8.cm2"), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Find S-like and G1-like tumour samples
# Last Modified: 04/06/19; 06/03/19
# -----------------------------------------------------------------------------
samples.nbl.sg1 <- setSamplesSG1(wd.rt.data, samples1, cors.samples)
writeTable(samples.nbl.sg1, file.path(wd.ngs, "nbl_wgs_n57-1.sg1"), colnames=T, rownames=F, sep="\t")
# [1] 3
# [1] 3
# [1] "CLBGA" "NGP"   "SKNAS"
# [1] "GIMEN" "LAN6"  "SKNFI"

# -----------------------------------------------------------------------------
# Overall correlation with LCL S/G1
# Last Modified: 16/06/19; 04/06/19; 06/03/19
# -----------------------------------------------------------------------------
samples.nbl <- setSamplesQ4(wd.rt.data, samples1)
writeTable(samples.nbl, file.path(wd.ngs, "nbl_cl_n8.txt"), colnames=T, rownames=F, sep="\t")
# 0%        25%        50%        75%       100% 
# -0.31416476 -0.18939137  0.01659257  0.33883294  0.38311250

writeTable(subset(samples.nbl, Q4 %in% c(4,1)), file.path(wd.ngs, "nbl_cl_n4.txt"), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# PCA
# Last Modified: 04/06/19; 21/04/19
# -----------------------------------------------------------------------------
## Copy from 2a_cmd-rt_rpkm.corr.gc.d_sample.R (commandline mode)
rpkms.T.chr.d.all <- NULL
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
 
   ## Read depth
   rpkms.T.chr.d <- pipeGetDetectedRD(wd.ngs.data, BASE, chr, PAIR1)
   colnames(rpkms.T.chr.d) <- toupper(gsub("\\.", "", colnames(rpkms.T.chr.d)))
   #rpkms.N.chr.d <- pipeGetDetectedRD(wd0.ngs.data, BASE, chr, PAIR0) 
   #overlaps <- intersect(rpkms.T.chr.d$BED, rpkms.N.chr.d$BED)
 
   test <- rpkms.T.chr.d[, samples.nbl$SAMPLE_ID]
   if (is.null(rpkms.T.chr.d.all)) {
      rpkms.T.chr.d.all <- test
   } else
      rpkms.T.chr.d.all <- rbind(rpkms.T.chr.d.all, test)
}

##
test <- rpkms.T.chr.d.all
pca.de <- getPCA(t(test))
save(pca.de, file=file.path(wd.rt.data, paste0("pca_nbl_CL_chrs_spline_spearman.RData")))

file.main <- c("NBL CL (n=8) read depth profiles", "")
trait <- as.numeric(samples.nbl$Q4)
trait[which(trait == 4)] <- "Q4 (0.34 < r < 0.38)"
trait[which(trait == 3)] <- "Q3 (0.02 < r < 0.34)"
trait[which(trait == 2)] <- "Q2 (-0.19 < r < 0.02)"
trait[which(trait == 1)] <- "Q1 (-0.31 < r < -0.19)"
plotPCA(1, 2, pca.de, trait, wd.rt.plots, "pca_nbl_CL_chrs_spline_spearman", size=6, file.main, "topleft", c("red", "lightcoral", "skyblue3", "blue"), samples1, flip.x=-1, flip.y=1, legend.title="Overall corr. with LCL S/G1")

## SG1
trait <- samples.nbl.sg1$SG1
plotPCA(1, 2, pca.de, trait, wd.rt.plots, "pca_nbl_CL_chrs_spline_spearman_SG1", size=6, file.main, "topright", c("red", "lightgray", "blue"), samples1, flip.x=-1, flip.y=1, legend.title="Consist. CM in all chrs")
