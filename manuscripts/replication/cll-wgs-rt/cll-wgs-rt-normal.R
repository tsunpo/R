# =============================================================================
# Manuscript   : 
# Chapter      : Chromosome replication timing of the human genome
# Name         : manuscripts/replication/cll-wgs-rt.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 24/10/19; 14/10/19; 26/02/19
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
BASE  <- "CLL"
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
samples1 <- readTable(file.path(wd.ngs, "cll_wgs_n96.list"), header=F, rownames=F, sep="")
n1 <- length(samples1)
## samples1 are modified to the normals in this script

# -----------------------------------------------------------------------------
# CLL N vs LCL S/G1
# http://pklab.med.harvard.edu/scw2014/subpop_tutorial.html
# Last Modified: 24/10/19; 06/03/19
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
# [1] -0.869272
# > max(cors.samples[,-c(1:4)])
# [1] 0.4091748

#load(file.path(wd.rt.data, paste0("samples-vs-rt_cll-vs-lcl_spline_spearman.RData")))
file.name <- file.path(wd.rt.plots, "SAMPLES-vs-RT_CLL-WB-vs-LCL_spline_spearman")
main.text <- c("CLL-WB read depth vs. LCL RT", "")
ymin <- -0.8773492
ymax <- 0.8392611
plotSAMPLEvsRTALL(cors.samples, samples1, file.name, main.text, ymin, ymax)

# -----------------------------------------------------------------------------
# Overall correlation with LCL S/G1
# Last Modified: 16/06/19; 04/06/19; 06/03/19
# -----------------------------------------------------------------------------
samples.cll.wb <- setSamplesQ4(wd.rt.data, samples1)
writeTable(samples.cll.wb, file.path(wd.ngs, "cll_wgs_n96_WB.txt"), colnames=T, rownames=F, sep="\t")
#         0%        25%        50%        75%       100% 
# -0.7829053 -0.7551694 -0.7481736 -0.7389845  0.2639248 

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

## Q4
test <- nrds.T.chr.d.all[, -1]   ## BUG 2019/10/14: Remove column BED
pca.de <- getPCA(t(test))
save(pca.de, file=file.path(wd.rt.data, paste0("pca_cll-wb_chrs.RData")))

#load(file.path(wd.rt.data, paste0("pca_cll-wb_chrs.RData")))
file.main <- c("CLL-WB overall read depth", "")
trait <- as.numeric(samples.cll.wb$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
plotPCA(1, 2, pca.de, trait, wd.rt.plots, "PCA_CLL-WB_chrs", size=6, file.main, "bottomright", c("red", "lightpink1", "lightskyblue2", "blue"), NULL, flip.x=-1, flip.y=1, legend.title=NA)











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
