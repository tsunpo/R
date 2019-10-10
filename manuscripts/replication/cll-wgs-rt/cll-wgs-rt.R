# =============================================================================
# Manuscript   : 
# Chapter      : Chromosome replication timing of the human genome
# Name         : manuscripts/replication/cll-wgs-rt.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 26/02/19
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
load(file.path(wd.src.ref, "hg19.rt.lcl.koren.woodfine.RData"))

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 30/01/18
# -----------------------------------------------------------------------------
#wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE  <- "CLL"
PAIR1 <- "N"
PAIR0 <- "T"
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
# Plot RD and RT (see ReplicationTiming.R)
# Last Modified: 28/05/19; 14/02/19; 10/01/19; 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
nrds <- toTable(NA, 4, 0, c("BED", "T", "N", "RT"))
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   nrds.chr <- readTable(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   nrds.chr$RT <- nrds.chr$N / nrds.chr$T   ## ADD 12/08/19
   
   nrds <- rbind(nrds, nrds.chr)
}
nrds$RT <- scale(nrds$RT)
save(nrds, file=file.path(wd.rt.data, paste0("nrds_", base, "-n-t_", method, ".RData")))
# > nrow(nrds)
# [1] 2651906

#load(file.path(wd.rt.data, paste0("nrds_", base, "-n-t_", method, ".RData")))
ymax <- 0.6
ymin <- 0.14
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]
   lcl.rt.chr <- subset(lcl.rt, CHR == chr)   ## Koren 2012
 
   ## Plot RT
   main.text <- paste0(BASE, " N/T read depth ratio between normal (n=", n0, ") and tumour (n=", n1, ") samples")  
   file.name <- file.path(wd.rt.plots, paste0("RT_", base, "_", method, ".d.rt_", chr, "_", PAIR0, "-", PAIR1, "_n", n0, "-", n1, ""))   
   plotRT(file.name, main.text, chr, NA, NA, nrds.chr, bed.gc.chr, c("red", "blue"), c("Normal", "Tumour"), c("lightcoral", "lightskyblue3"), c("N", "T"), "png", width=10, peaks=c(), ylim=c(ymin, ymax), lcl.rt.chr)
}

# -----------------------------------------------------------------------------
# RD vs RT (RDS and SPR)
# Last Modified: 11/07/19; 31/05/19
# -----------------------------------------------------------------------------
cors <- toTable(0, 8, 22, c("chr", "length", "cor", "cor1", "cor2", "mean", "intercept1", "intercept2"))
cors$chr <- 1:22
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]
 
   nrds.chr.T  <- setSpline(nrds.chr, bed.gc.chr, "N")
   nrds.chr.N  <- setSpline(nrds.chr, bed.gc.chr, "T")
   nrds.chr.RT <- setSpline(nrds.chr, bed.gc.chr, "RT")
   #nrds.chr.RT$SPLINE <- scale(nrds.chr.RT$SPLINE)
   cors$length[c] <- nrow(nrds.chr.RT)
   cors$mean[c]   <- mean(nrds.chr.RT$SPLINE)
 
   cor <- getCor(nrds.chr.T$SPLINE, nrds.chr.N$SPLINE, method="spearman")
   cors$cor[c] <- cor
 
   main.text <- c(paste0("CLL read depth correlation (", "Chr", c, ")"), paste0("rho = ", round0(cor, digits=2), " (N vs. T)"))
   xlab.text <- "CLL N/T"
   ylab.text <- "CLL read depth [RPKM]"
   file.name <- file.path(wd.rt.plots, "chrs", paste0("RD-vs-RT_CLL-N-T_chr", c, "_spline_spearman"))
   plotRD2vsRT(nrds.chr.T$SPLINE, nrds.chr.N$SPLINE, nrds.chr.RT$SPLINE, file.name, main.text, ylab.text, xlab.text, c("red", "blue"), c("N", "T"), method="spearman")
 
   cors$cor1[c] <- getCor(nrds.chr.T$SPLINE, nrds.chr.RT$SPLINE, method="spearman")
   cors$cor2[c] <- getCor(nrds.chr.N$SPLINE, nrds.chr.RT$SPLINE, method="spearman")
   cors$intercept1[c] <- lm(nrds.chr.T$SPLINE ~ nrds.chr.RT$SPLINE)[[1]][1]
   cors$intercept2[c] <- lm(nrds.chr.N$SPLINE ~ nrds.chr.RT$SPLINE)[[1]][1]
 
   ## Read depth skew (RDS)
   cors$skew <- (cors$intercept1 - cors$intercept2) / (cors$intercept1 + cors$intercept2)   
}
save(cors, file=file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-n-t_spline_spearman.RData")))
writeTable(cors, file=file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-n-t_spline_spearman.txt")), colnames=T, rownames=F, sep="\t")

## S-phase progression rate (SPR)
file.name <- file.path(wd.rt.plots, "SPR_CLL-N-T_spline_spearman")
main.text <- c(paste0(BASE, " S-phase progression rate"), "SPR = Mean N/T ratio")
plotSPR(cors, file.name, main.text, c(13, 17), digits=3, unit=5, ylab.text=paste0(BASE, " SPR"))

## SPR vs Read depth correlation (RDC)
file.name <- file.path(wd.rt.plots, "SPR-RDC_CLL-N-T_spline_spearman")
main.text <- c(paste0(BASE, " SPR vs. Read depth correlation"), "")
xlab.text <- "CLL read depth correlation [rho]"
plotSPRRDC(cors, file.name, main.text, c(4, 13, 17, 19, 21, 22), xlab.text, unit=5, ylab.text=paste0(BASE, " SPR"), lcl.mean=NULL)

## SPR vs Woodfine 2004
file.name <- file.path(wd.rt.plots, "SPR-Woodfine_CLL-N-T_spline_spearman")
main.text <- c(paste0(BASE, " SPR vs. Woodfine 2004"), "Mean replication timing ratio ")
xlab.text <- "Woodfine et al 2004"
plotSPRRDC(cors, file.name, main.text, c(4, 13, 17, 19, 21, 22), xlab.text, unit=5, ylab.text=paste0(BASE, " SPR"), lcl.mean=lcl.mean)

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
# [1] -0.8789273
# > max(cors.samples[,-c(1:4)])
# [1] 0.6874687

#load(file.path(wd.rt.data, paste0("samples-vs-rt_cll-vs-lcl_spline_spearman.RData")))
file.name <- file.path(wd.rt.plots, "SAMPLES-vs-RT_CLL-vs-LCL_spline_spearman")
main.text <- c("CLL (n=96) read depth vs. LCL S/G1", "")
ymin <- -0.8789273
ymax <- 0.8433154
plotSAMPLEvsRTALL(cors.samples, samples1, file.name, main.text, ymin, ymax)

# -----------------------------------------------------------------------------
# Overall correlation with LCL S/G1
# Last Modified: 16/06/19; 04/06/19; 06/03/19
# -----------------------------------------------------------------------------
samples.cll <- setSamplesQ4(wd.rt.data, samples1)
writeTable(samples.cll, file.path(wd.ngs, "cll_wgs_n96.txt"), colnames=T, rownames=F, sep="\t")
# 0%        25%        50%        75%       100% 
# -0.8005515 -0.7701020 -0.7577076 -0.7466137  0.5436926 

writeTable(subset(samples.cll, Q4 %in% c(4,1)), file.path(wd.ngs, "cll_wgs_q4_n48.txt"), colnames=T, rownames=F, sep="\t")
writeTable(subset(samples.cll, Q4 %in% c(3,1)), file.path(wd.ngs, "cll_wgs_q3_n48.txt"), colnames=T, rownames=F, sep="\t")
writeTable(subset(samples.cll, Q4 %in% c(2,1)), file.path(wd.ngs, "cll_wgs_q2_n48.txt"), colnames=T, rownames=F, sep="\t")

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
 
   test <- nrds.T.chr.d[, samples.cll$SAMPLE_ID]
   if (is.null(nrds.T.chr.d.all)) {
      nrds.T.chr.d.all <- test
   } else
      nrds.T.chr.d.all <- rbind(nrds.T.chr.d.all, test)
}

## Q4
test <- nrds.T.chr.d.all
pca.de <- getPCA(t(test))
save(pca.de, file=file.path(wd.rt.data, paste0("pca_cll_chrs_spline_spearman.RData")))

#load(file.path(wd.rt.data, paste0("PCA_cll_chrs_spline_spearman.RData")))
file.main <- c("CLL (n=96) read depth profiles", "")
trait <- as.numeric(samples.cll$Q4)
trait[which(trait == 4)] <- "Q4"   ##(-0.75 < rho < 0.54)"
trait[which(trait == 3)] <- "Q3"   ##(-0.76 < rho < -0.75)"
trait[which(trait == 2)] <- "Q2"   ##(-0.77 < rho < -0.76)"
trait[which(trait == 1)] <- "Q1"   ##(-0.80 < rho < -0.77)"
plotPCA(1, 2, pca.de, trait, wd.rt.plots, "PCA_CLL_chrs_spline_spearman", size=6, file.main, "bottomright", c("red", "lightcoral", "skyblue3", "blue"), NULL, flip.x=-1, flip.y=1, legend.title=NA)

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
