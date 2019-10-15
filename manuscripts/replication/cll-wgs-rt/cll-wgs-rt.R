# =============================================================================
# Manuscript   : 
# Chapter      : Chromosome replication timing of the human genome
# Name         : manuscripts/replication/cll-wgs-rt.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 14/10/19; 26/02/19
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "DifferentialExpression.R", "ReplicationTiming.R")
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
BASE  <- "CLL"
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
samples1 <- readTable(file.path(wd.ngs, "cll_wgs_n96.list"), header=F, rownames=F, sep="")
samples0 <- readTable(file.path(wd.ngs, "cll_wgs_n96.list"), header=F, rownames=F, sep="")
n1 <- length(samples1)
n0 <- length(samples0)

# -----------------------------------------------------------------------------
# CLL T vs LCL S/G1
# http://pklab.med.harvard.edu/scw2014/subpop_tutorial.html
# Last Modified: 06/03/19
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
save(cors.samples, file=file.path(wd.rt.data, paste0("samples-vs-rt_cll-vs-lcl_spline_spearman.RData")))
# > min(cors.samples[,-c(1:4)])
# [1] -0.8773492
# > max(cors.samples[,-c(1:4)])
# [1] 0.6832058

#load(file.path(wd.rt.data, paste0("samples-vs-rt_cll-vs-lcl_spline_spearman.RData")))
file.name <- file.path(wd.rt.plots, "SAMPLES-vs-RT_CLL-vs-LCL_spline_spearman")
main.text <- c("CLL read depth vs. LCL S/G1", "")
ymin <- -0.8773492
ymax <- 0.8392611
plotSAMPLEvsRTALL(cors.samples, samples1, file.name, main.text, ymin, ymax)

# -----------------------------------------------------------------------------
# Overall correlation with LCL S/G1
# Last Modified: 16/06/19; 04/06/19; 06/03/19
# -----------------------------------------------------------------------------
samples.cll <- setSamplesQ4(wd.rt.data, samples1)
writeTable(samples.cll, file.path(wd.ngs, "cll_wgs_n96.txt"), colnames=T, rownames=F, sep="\t")
#         0%        25%        50%        75%       100% 
# -0.7952126 -0.7638682 -0.7514207 -0.7399593  0.5381444 

writeTable(subset(samples.cll, Q4 %in% c(4,1)), file.path(wd.ngs, "cll_wgs_q4_n48.txt"), colnames=T, rownames=F, sep="\t")
writeTable(subset(samples.cll, Q4 %in% c(3,1)), file.path(wd.ngs, "cll_wgs_q3_n48.txt"), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# PCA
# Last Modified: 04/06/19; 21/04/19
# -----------------------------------------------------------------------------
## Refer to cmd-rt_2a_nrd.gc.cn.d_sample.R (commandline mode)
nrds.T.chr.d.all <- readTable(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d_", "chr1", "_", PAIR1, "_n", n1, ".txt.gz")), header=T, rownames=F, sep="\t")
for (c in 2:22) {
   chr <- chrs[c]

   ## Read depth
   #nrds.T.chr.d <- pipeGetDetectedRD(wd.ngs.data, BASE, chr, PAIR1, method)
   nrds.T.chr.d <- readTable(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d_", chr, "_", PAIR1, "_n", n1, ".txt.gz")), header=T, rownames=F, sep="\t")
   
   nrds.T.chr.d.all <- rbind(nrds.T.chr.d.all, nrds.T.chr.d)
}

## Q4
test <- nrds.T.chr.d.all[, -1]   ## BUG 2019/10/14: Remove column BED
pca.de <- getPCA(t(test))
save(pca.de, file=file.path(wd.rt.data, paste0("pca_cll_chrs.RData")))

#load(file.path(wd.rt.data, paste0("pca_cll_chrs.RData")))
file.main <- c("CLL (n=96) read depth profiles", "")
trait <- as.numeric(samples.cll$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
plotPCA(1, 2, pca.de, trait, wd.rt.plots, "PCA_CLL_chrs", size=6, file.main, "bottomright", c("red", "lightcoral", "skyblue3", "blue"), NULL, flip.x=-1, flip.y=1, legend.title=NA)











## SG1
#trait <- samples.cll.sg1$SG1
#plotPCA(1, 2, pca.de, trait, wd.rt.plots, "pca_cll_T_chrs_spline_spearman_SG1", size=6, file.main, "bottomright", c("red", "lightgray", "blue"), NULL, flip.x=-1, flip.y=1, legend.title="Consist. CM in all chrs")

# -----------------------------------------------------------------------------
# Find S-like and G1-like tumour samples
# Last Modified: 04/06/19; 06/03/19
# -----------------------------------------------------------------------------
samples.cll.sg1 <- setSamplesSG1(wd.rt.data, samples1, cors.samples)
writeTable(samples.cll.sg1, file.path(wd.ngs, "cll_wgs_n96.sg1"), colnames=T, rownames=F, sep="\t")
# > length(s_likes)
# [1] 7
# > length(g1_likes)
# [1] 10
# > s_likes
# [1] "SP115064" "SP116746" "SP116766" "SP13291"  "SP13323"  "SP13365"  "SP13490"
# > g1_likes
# [1] "SP115060" "SP116744" "SP116754" "SP116764" "SP116776" "SP116780" "SP13442"  "SP13478"  "SP13624"  "SP96882"
