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
handbooks  <- c("Commons.R", "ReplicationTiming.R", "Transcription.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.bed.gc.1kb.RData"))
#load(file.path(wd.src.ref, "hg19.lcl.koren.woodfine.RData"))

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
cors.samples <- getSAMPLEvsRT(wd.rt.data, samples1)
save(cors.samples, file=file.path(wd.rt.data, paste0("samples-vs-rt_", base, "-vs-lcl_spline_spearman.RData")))
# > min(cors.samples[,-c(1:4)])
# [1] -0.8773492
# > max(cors.samples[,-c(1:4)])
# [1] 0.6832058

#load(file.path(wd.rt.data, paste0("samples-vs-rt_cll-vs-lcl_spline_spearman.RData")))
file.name <- file.path(wd.rt.plots, "SAMPLES-vs-RT_CLL-vs-LCL_spline_spearman")
main.text <- c("CLL read depth vs. RT", "")
ymin <- -0.8773492
ymax <- 0.8392611
plotSAMPLEvsRT(cors.samples, samples1, file.name, main.text, ymin, ymax)

# -----------------------------------------------------------------------------
# Last Modified: 30/11/20; 04/06/19; 21/04/19
# -----------------------------------------------------------------------------
nrds <- getOverallReadDepth(wd.rt.data, base, method, PAIR1, n1)
# > nrow(nrds)
# [1] 2684771
nrds.m <- subset(nrds, MEDIAN != 0)   ## ADD 30/11/20
# > nrow(nrds.m)
# [1] 2674348
nrds.d <- getOverallDetectedRD(wd.rt.data, base, method, PAIR1, n1, samples1)
# > nrow(nrds.d)
# [1] 

# -----------------------------------------------------------------------------
# Overall correlation with LCL S/G1
# Last Modified: 19/11/19; 16/06/19; 04/06/19; 06/03/19
# -----------------------------------------------------------------------------
samples.cll <- setSamplesQ4(wd.rt.data, samples1)
writeTable(samples.cll, file.path(wd.ngs, "cll_wgs_n96.txt"), colnames=T, rownames=F, sep="\t")
#         0%        25%        50%        75%       100% 
# -0.7952126 -0.7638682 -0.7514207 -0.7399593  0.5381444 

writeTable(subset(samples.cll, Q4 %in% c(4,1)), file.path(wd.ngs, "cll_wgs_q4_n48.txt"), colnames=T, rownames=F, sep="\t")
writeTable(subset(samples.cll, Q4 %in% c(3,1)), file.path(wd.ngs, "cll_wgs_q3_n48.txt"), colnames=T, rownames=F, sep="\t")

## Random 50/50
m2.24 <- sort(rownames(subset(samples.cll, M2 == 1))[sample(1:48, round(48/2), replace=F)])
m1.24 <- sort(rownames(subset(samples.cll, M2 == 0))[sample(1:48, round(48/2), replace=F)])

random1 <- sort(c(m2.24, m1.24))
writeTable(samples.cll[random1,], file.path(wd.ngs, "cll_wgs_n96-random1.txt"), colnames=T, rownames=F, sep="\t")

random2 <- sort(setdiff(rownames(samples.cll), random1))
writeTable(samples.cll[random2,], file.path(wd.ngs, "cll_wgs_n96-random2.txt"), colnames=T, rownames=F, sep="\t")

## Random 25/25
m2 <- sort(rownames(subset(samples.cll, M2 == 1)))
m1 <- sort(rownames(subset(samples.cll, M2 == 0)))

m2.12 <- sort(m2[sample(1:48, round(48/4), replace=F)])
m1.12 <- sort(m1[sample(1:48, round(48/4), replace=F)])
random3 <- sort(c(m2.12, m1.12))
writeTable(samples.cll[random3,], file.path(wd.ngs, "cll_wgs_n96-random3.txt"), colnames=T, rownames=F, sep="\t")

m2 <- sort(setdiff(m2, m2.12))
m1 <- sort(setdiff(m1, m1.12))

m2.12 <- sort(m2[sample(1:36, round(48/4), replace=F)])
m1.12 <- sort(m1[sample(1:36, round(48/4), replace=F)])
random4 <- sort(c(m2.12, m1.12))
writeTable(samples.cll[random4,], file.path(wd.ngs, "cll_wgs_n96-random4.txt"), colnames=T, rownames=F, sep="\t")








## Random 14/14
m2 <- sort(rownames(subset(samples.nbl, M2 == 1)))
m1 <- sort(rownames(subset(samples.nbl, M2 == 0)))

m2.7 <- sort(m2[sample(1:28, round(28/4), replace=F)])
m1.7 <- sort(m1[sample(1:28, round(28/4), replace=F)])
random3 <- sort(c(m2.7, m1.7))
writeTable(samples.nbl[random3,], file.path(wd.ngs, "nbl_wgs_n56-random3.txt"), colnames=T, rownames=F, sep="\t")

m2 <- sort(setdiff(m2, m2.7))
m1 <- sort(setdiff(m1, m1.7))

m2.7 <- sort(m2[sample(1:21, round(28/4), replace=F)])
m1.7 <- sort(m1[sample(1:21, round(28/4), replace=F)])
random4 <- sort(c(m2.7, m1.7))
writeTable(samples.nbl[random4,], file.path(wd.ngs, "nbl_wgs_n56-random4.txt"), colnames=T, rownames=F, sep="\t")








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
file.main <- c("CLL", "")
trait <- as.numeric(samples.cll$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
plotPCA(1, 2, pca.de, trait, wd.rt.plots, "PCA_CLL", size=6, file.main, "bottomright", c(red, red.lighter, blue.lighter, blue), NULL, flip.x=-1, flip.y=1, legend.title=NA)











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
