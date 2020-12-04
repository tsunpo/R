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
#wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
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
# Last Modified: 05/06/19; 20/04/19; 06/03/19
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
# Overall correlation with LCL S/G1
# Last Modified: 19/11/19; 16/06/19; 04/06/19; 06/03/19
# -----------------------------------------------------------------------------
samples.nbl <- setSamplesQ4(wd.rt.data, samples1)
writeTable(samples.nbl, file.path(wd.ngs, "nbl_wgs_n57-1.txt"), colnames=T, rownames=F, sep="\t")
#         0%        25%        50%        75%       100% 
# -0.7684906 -0.7397140 -0.7200663 -0.4456056  0.7665285

writeTable(subset(samples.nbl, Q4 %in% c(4,1)), file.path(wd.ngs, "nbl_wgs_q4_n28.txt"), colnames=T, rownames=F, sep="\t")
writeTable(subset(samples.nbl, Q4 %in% c(3,1)), file.path(wd.ngs, "nbl_wgs_q3_n28.txt"), colnames=T, rownames=F, sep="\t")

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
pca.de <- getPCA(t(test))
save(pca.de, file=file.path(wd.rt.data, paste0("pca_nbl_chrs.RData")))

#load(file.path(wd.rt.data, paste0("pca_nbl_chrs.RData")))
file.main <- c("NBL overall read depth", "")
trait <- as.numeric(samples.nbl$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
plotPCA(1, 2, pca.de, trait, wd.rt.plots, "PCA_NBL_nasa.blue_google.red_lighter.50_pch=19_pt.cex=2.5_pch=19_2314", size=6, file.main, "bottomright", c(red, lighterred, lighterblue, blue), NULL, flip.x=1, flip.y=1, legend.title=NA)











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
