# =============================================================================
# Manuscript   : 
# Chapter      : Chromosome replication timing of the human genome
# Name         : manuscripts/replication/nbl-wgs-rt.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 30/11/20; 14/10/19; 25/02/19
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
#wd <- "/Users/ty2/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE  <- "NBL"
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
samples1 <- readTable(file.path(wd.ngs, "nbl_wgs_n57-1.list"), header=F, rownames=F, sep="")
samples0 <- readTable(file.path(wd.ngs, "nbl_wgs_n57-1.list"), header=F, rownames=F, sep="")
n1 <- length(samples1)
n0 <- length(samples0)

# -----------------------------------------------------------------------------
# NBL T vs LCL S/G1
# http://pklab.med.harvard.edu/scw2014/subpop_tutorial.html
# Last Modified: 23/06/22; 05/06/19; 20/04/19; 06/03/19
# -----------------------------------------------------------------------------
cors.samples <- getSAMPLEvsRT(wd.rt.data, samples1)
save(cors.samples, file=file.path(wd.rt.data, paste0("samples-vs-rt_", base, "-vs-lcl_spline_spearman.RData")))
# > min(cors.samples[,-c(1:4)])
# [1] -0.8354356
# > max(cors.samples[,-c(1:4)])
# [1] 0.8392611

#load(file.path(wd.rt.data, paste0("samples-vs-rt_nbl-vs-lcl_spline_spearman.RData")))
file.name <- file.path(wd.rt.plots, "SAMPLES-vs-RT_NBL-vs-LCL_spline_spearman")
main.text <- c("NBL read depth vs. RT", "")
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
# [1] 2674140
nrds.d <- getOverallDetectedRD(wd.rt.data, base, method, PAIR1, n1, samples1)
nrow(nrds.d)
# > nrow(nrds.d)
# [1] 2659570

# -----------------------------------------------------------------------------
# Last Modified: 23/06/22
# -----------------------------------------------------------------------------
nrds <- getOverallReadDepth(wd.rt.data, base, method, PAIR1, n1)
nrow(nrds)
# [1] 
nrds.m <- subset(nrds, MEDIAN != 0)   ## ADD 30/11/20
nrow(nrds.m)
# [1] 
nrds.d <- getOverallDetectedRD(wd.rt.data, base, method, PAIR1, n1, samples1)
nrow(nrds.d)
# [1] 

# -----------------------------------------------------------------------------
# Overall correlation with LCL S/G1
# Last Modified: 19/11/19; 16/06/19; 04/06/19; 06/03/19
# -----------------------------------------------------------------------------
samples.nbl <- setSamplesM2(wd.rt.data, samples1, c=13)
writeTable(samples.nbl, file.path(wd.ngs, "nbl_wgs_n57-1_cor13.txt"), colnames=T, rownames=F, sep="\t")
writeTable(median(samples.nbl$COR), file.path(wd.ngs, "nbl_wgs_n57-1_COR.txt"), colnames=F, rownames=F, sep="\t")
writeTable(median(samples.sclc$COR), file.path(wd.ngs, "nbl_wgs_n57-1_COR_SCLC.txt"), colnames=F, rownames=F, sep="\t")
#         0%         25%         50%         75%        100% 
# -0.19619889 -0.06616305  0.01765322  0.19912003  0.35027779.  ## n=57-1 (NEW 02/06/24)
#         0%        25%        50%        75%       100%        ## n=56/57
# -0.7603965 -0.7321717 -0.7109754 -0.4246964  0.7624673
#         0%        25%        50%        75%       100%        ## n=57-1
# -0.7684906 -0.7397140 -0.7200663 -0.4456056  0.7665285

writeTable(subset(samples.nbl, Q4 %in% c(4,1)), file.path(wd.ngs, "nbl_wgs_q4_n28.txt"), colnames=T, rownames=F, sep="\t")
writeTable(subset(samples.nbl, Q4 %in% c(3,1)), file.path(wd.ngs, "nbl_wgs_q3_n28.txt"), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# 
# Last Modified: 03/06/22
# -----------------------------------------------------------------------------
y <- samples.nbl.old[rownames(samples.nbl),]$COR
x <- samples.nbl$COR
file.name <- file.path(wd.rt.plots, paste0("correlation_n56_vs_n57-1"))
plotCorrelation(file.name, expression(italic("In silico") ~ "sorting of NBL samples (n=55)"), "n=56", "n=57-1", x, y)

# -----------------------------------------------------------------------------
# 
# Last Modified: 
# -----------------------------------------------------------------------------
## Random 50/50
m2.14 <- sort(rownames(subset(samples.nbl, M2 == 1))[sample(1:28, round(28/2), replace=F)])
m1.14 <- sort(rownames(subset(samples.nbl, M2 == 0))[sample(1:28, round(28/2), replace=F)])

random1 <- sort(c(m2.14, m1.14))
writeTable(samples.nbl[random1,], file.path(wd.ngs, "nbl_wgs_n56-random1.txt"), colnames=T, rownames=F, sep="\t")

random2 <- sort(setdiff(rownames(samples.nbl), random1))
writeTable(samples.nbl[random2,], file.path(wd.ngs, "nbl_wgs_n56-random2.txt"), colnames=T, rownames=F, sep="\t")

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

##
test <- nrds.T.chr.d.all[, -1]   ## BUG 2019/10/14: Remove column BED
nrow(test)
# [1] 2684771
test.nona <- test[removeMissing(test),]
nrow(test.nona)
# [1] 2684771
test.no0 <- test[getExpressed(test),]
test.no0 <- test.no0[,-(ncol(test.no0))]
nrow(test.no0)
# [1] 2659482
pca.1kb <- getPCA(t(test.no0[, rownames(samples.nbl)]))
save(pca.1kb, file=file.path(wd.rt.data, paste0("PCA_1KB_NBL_n57.RData")))

#load(file.path(wd.rt.data, paste0("pca_nbl_chrs.RData")))
file.main <- c("NBL 1 kb read depth", "")
trait <- as.numeric(samples.nbl$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
#plotPCA(1, 2, pca.1kb.nbl, trait, wd.rt.plots, "PCA_NBL", size=6, file.main, "bottomright", c(red, red.lighter, blue.lighter, blue), samples.nbl, flip.x=1, flip.y=1, legend.title=NA)
plotPCA(1, 2, pca.1kb.nbl, trait, wd.rt.plots, "PCA_1KB_NBL_Q4", size=6, file.main, "bottomright", c("Q4", "Q3", "Q2", "Q1"), c(red, red.lighter, blue.lighter, blue), flip.x=1, flip.y=1)

trait <- as.numeric(samples.nbl$M2)
trait[which(trait == 2)] <- "Proliferative"
trait[which(trait == 1)] <- "Resting"
plotPCA(1, 2, pca.1kb.nbl, trait, wd.rt.plots, "PCA_1KB_NBL_M2", size=6, file.main, "bottomright", c("Proliferative", "Resting"), c(red, blue), flip.x=1, flip.y=1)

###
## WGD (28/06/22)
scores.1kb <- pcaScores(pca.1kb)
overlaps <- intersect(rownames(scores.1kb), rownames(wgds))

x <- scores.1kb[overlaps,]$PC1
y <- wgds[overlaps,]$decision_value
file.name <- file.path(wd.rt.plots, paste0("correlation_1KB-PC1-vs-WGD"))
plotCorrelation(file.name, "NB 1 kb read depth", paste0("PC", 1, " (", pcaProportionofVariance(pca.1kb, 1), "%)"), "WGD", x, y, size=6)

x <- scores.1kb[overlaps,]$PC1
y <- samples.nbl[overlaps,]$COR
file.name <- file.path(wd.rt.plots, paste0("correlation_1KB-PC1-vs-SORTING"))
plotCorrelation(file.name, "NB 1 kb read depth", paste0("PC", 1, " (", pcaProportionofVariance(pca.1kb, 1), "%)"), "Proliferation", x, y, size=6)

trait <- wgds[overlaps,]$wgd_predicted
file.main <- c("NB 1 kb read depth", "")
plotPCA(1, 2, pca.1kb, trait, wd.rt.plots, "PCA-1KB_NB_WGD", size=6, file.main, "bottomright", c("WGD", "non-WGD"), c("black", "lightgray"), flip.x=1, flip.y=1)

trait <- samples.nbl[overlaps,]$M2
trait[which(trait == 2)] <- "Proliferative"
trait[which(trait == 1)] <- "Resting"
file.main <- c("NB 1 kb read depth", "")
plotPCA(1, 2, pca.1kb, trait, wd.rt.plots, "PCA-1KB_NB_M2", size=6, file.main, "topleft", c("Proliferative", "Resting"), c(red, blue), flip.x=1, flip.y=1)

trait <- as.vector(samples.nbl[overlaps,]$SORTING)
trait[which(trait == "S")] <- "Proliferative"
trait[which(trait == "G1")] <- "Resting"
file.main <- c("NB 1 kb read depth", "")
plotPCA(1, 2, pca.1kb, trait, wd.rt.plots, "PCA-1KB_NB_SORTING", size=6, file.main, "topleft", c("Proliferative", "Resting"), c(red, blue), flip.x=1, flip.y=1)





###
## WGD (28/06/22)
scores.1kb.nbl <- pcaScores(pca.1kb.nbl)
wgds.nbl <- readTable(file.path(wd.ngs, "nbl_wgs_n57_WGD.txt"), header=T, rownames=T, sep="\t")
wgds.nbl[which(wgds.nbl$wgd_predicted == 0), ]$wgd_predicted <- "non-WGD"
wgds.nbl[which(wgds.nbl$wgd_predicted == 1), ]$wgd_predicted <- "WGD"
overlaps <- intersect(rownames(scores.1kb.nbl), rownames(wgds.nbl))

x <- scores.1kb.nbl[overlaps,]$PC1
y <- wgds.nbl[overlaps,]$decision_value
file.name <- file.path(wd.rt.plots, paste0("correlation_1KB-PC1-vs-WGD"))
plotCorrelation(file.name, "NB 1 kb read depth", paste0("PC", 1, " (", pcaProportionofVariance(pca.1kb.nbl, 1), "%)"), "WGD", x, y, size=6)

x <- scores.1kb.nbl[overlaps,]$PC1
y <- samples.nbl[overlaps,]$COR
file.name <- file.path(wd.rt.plots, paste0("correlation_NBL_PC1_1KB-vs-SORTING"))
plotCorrelation(file.name, "NBL 1 kb read depth", paste0("PC", 1, " (", pcaProportionofVariance(pca.1kb.nbl, 1), "%)"), expression(italic("In silico") ~ "sorting"), x, y, size=6)

trait <- as.numeric(samples.nbl$M2)
trait[which(trait == 2)] <- "Proliferative"
trait[which(trait == 1)] <- "Resting"
plotPCA(1, 2, pca.1kb.nbl, trait, wd.rt.plots, "PCA_1KB_NBL_M2", size=6, file.main, "bottomright", c("Proliferative", "Resting"), c(red, blue), flip.x=1, flip.y=1)

trait <- wgds.nbl$wgd_predicted
plotPCA(1, 2, pca.1kb.nbl, trait, wd.rt.plots, "PCA_1KB_NBL_WGD", size=6, file.main, "bottomright", c("WGD", "non-WGD"), c("black", "lightgray"), flip.x=1, flip.y=1)




y <- samples.sclc$COR
x <- scores$PC2
file.name <- file.path(wd.rt.plots, paste0("correlation_SORTING_vs_PC2"))
plotCorrelation(file.name, "SCLC 1 kb windows", paste0("PC", 2, " (", pcaProportionofVariance(pca.de, 2), "%)"), expression(italic("In silico") ~ "sorting"), x, y, size=6)










## SG1
#trait <- samples.nbl.sg1$SG1
#plotPCA(1, 2, pca.de, trait, wd.rt.plots, "pca_nbl_T_chrs_spline_spearman_SG1", size=6, file.main, "bottomright", c("red", "lightgray", "blue"), NULL, flip.x=1, flip.y=1, legend.title="Consist. CM in all chrs")

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
writeTable(cq4[, c("SAMPLE_ID", paste0("chr", c(1:22)))], file.path(wd.ngs, "nbl_wgs_n53.cq4"), colnames=T, rownames=F, sep="\t")
writeTable(cm2[, c("SAMPLE_ID", paste0("chr", c(1:22)))], file.path(wd.ngs, "nbl_wgs_n53.cm2"), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Find S-like and G1-like tumour samples
# Last Modified: 04/06/19; 06/03/19
# -----------------------------------------------------------------------------
samples.nbl.sg1 <- setSamplesSG1(wd.rt.data, samples1, cors.samples)
writeTable(samples.nbl.sg1, file.path(wd.ngs, "nbl_wgs_n57-1.sg1"), colnames=T, rownames=F, sep="\t")
# > length(s_likes)
# [1] 16
# > length(g1_likes)
# [1] 10
# > s_likes
# [1] "P15239" "P19537" "P21702" "P21776" "P21924" "P22496" "P23103" "P23267" "P24478" "P24632" "P24679" "P24702" "P24885" "P25114" "P25262" "P25376"
# > g1_likes
# [1] "P13967" "P16885" "P1695"  "P17344" "P17612" "P18478" "P18972" "P19743" "P20865" "P22283"
