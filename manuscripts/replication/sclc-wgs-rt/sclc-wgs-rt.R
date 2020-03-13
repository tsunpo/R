# =============================================================================
# Manuscript   : 
# Chapter      : Chromosome replication timing of the human genome
# Name         : manuscripts/replication/sclc-wgs-rt.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 14/10/19; 05/03/19; 25/02/19; 30/01/18
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
samples0 <- readTable(file.path(wd.ngs, "sclc_wgs_n92_NL.list"), header=F, rownames=F, sep="")
#samples0 <- readTable(file.path(wd.ngs, "sclc_wgs_n9_WB.list"), header=F, rownames=F, sep="")
n1 <- length(samples1)
n0 <- length(samples0)

# -----------------------------------------------------------------------------
# Plot RD and RT (see ReplicationTiming.R)
# Last Modified: 28/05/19; 14/02/19; 10/01/19; 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
nrds <- getLog2ScaledRT(wd.rt.data, base, method, PAIR1, PAIR0, n1, n0, chrs, bed.gc)
save(nrds, file=file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d.rt.log2s_", "t-n", ".RData")))
#load(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d.rt.log2s_", "t-n", ".RData")))
# > nrow(nrds)
# [1] 2656170 - 22
nrds.sclc <- nrds

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
sprs <- getSPR(nrds, bed.gc)
save(sprs, file=file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-t-n_spline_spearman.RData")))
writeTable(sprs, file=file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-t-n_spline_spearman.txt")), colnames=T, rownames=F, sep="\t")
#load(file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-t-n_spline_spearman.RData")))

for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
 
   nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]
   nrds.chr.T  <- setSpline(nrds.chr, bed.gc.chr, "T")
   nrds.chr.N  <- setSpline(nrds.chr, bed.gc.chr, "N")
   nrds.chr.RT <- setSpline(nrds.chr, bed.gc.chr, "RT")
 
   main.text <- c(paste0("SCLC read depth correlation (", "Chr", c, ")"), paste0("rho = ", round0(sprs$cor[c], digits=2), " (T vs. N)"))
   xlab.text <- "T/N read depth ratio [log2]"
   ylab.text <- "Read depth [RPKM]"
   file.name <- file.path(wd.rt.plots, "chrs", paste0("RD-vs-RT_SCLC-T-N_chr", c, "_spline_spearman"))
   plotRD2vsRT(nrds.chr.T$SPLINE, nrds.chr.N$SPLINE, nrds.chr.RT$SPLINE, file.name, main.text, ylab.text, xlab.text, c("red", "blue"), c("T", "N"), method="spearman")
}

## S-phase progression rate (SPR)
ylab.text <- "SPR"
file.name <- file.path(wd.rt.plots, "SPR_SCLC-T-N_spline_spearman")
main.text <- c(paste0(BASE, " S-phase progression rate"), "SPR = (E-L)/(E+L)")
plotSPR(sprs, file.name, main.text, c(13, 17), digits=3, unit=5, ylab.text)

## SPR vs Read depth correlation
file.name <- file.path(wd.rt.plots, "SPR-RDC_SCLC-T-N_spline_spearman")
main.text <- c(paste0(BASE, " SPR vs. Read depths correlation"), "")
xlab.text <- "T vs. N [rho]"
plotSPRRDC(sprs$spr, sprs$cor, file.name, main.text, c(4, 13, 17, 19, 22), xlab.text, unit=5, ylab.text)

## SPR vs Woodfine 2004
file.name <- file.path(wd.rt.plots, "SPR-Woodfine_SCLC-T-N_spline_spearman")
main.text <- c(paste0(BASE, " SPR vs. Woodfine et al. 2004"), "")
xlab.text <- "Mean replication timing ratio"
plotSPRRDC(sprs$spr, lcl.mean$Mean, file.name, main.text, c(4, 13, 17, 19, 22), xlab.text, unit=5, ylab.text)

# -----------------------------------------------------------------------------
# RT vs LCL S/G1
# Last Modified: 27/05/19
# -----------------------------------------------------------------------------
nrds.tmp <- nrds
load(file.path(wd, "LCL/analysis/replication/lcl-wgs-rt/data/lcl_rpkm.gc.cn.d.rt.log2s_s-g1.RData"))
nrds.lcl <- nrds
nrds <- nrds.tmp

cors <- getRTvsRT(nrds, nrds.lcl, bed.gc)
save(cors, file=file.path(wd.rt.data, paste0("rt-vs-rt_", base, "-t-n-vs-lcl-s-g1_spline_spearman.RData")))
#load(file.path(wd.rt.data, paste0("rt-vs-rt_", base, "-t-n-vs-lcl-s-g1_spline_spearman.RData")))

ylab.text <- "Spearman's rho"
xlab.text <- "Chromosome"
file.name <- file.path(wd.rt.plots, "RT-vs-RT_SCLC-T-N-vs-LCL-S-G1_spline_spearman")
main.text <- paste0("SCLC T/N vs. LCL S/G1")
ymin <- 0.25
ymax <- 1.05
plotRTvsRTALL(cors, file.name, main.text, ylab.text, xlab.text, ymin, ymax, col="black", c=2, pos=1)

##
file.name <- file.path(wd.rt.plots, "RTD-vs-RT_SCLC-T-N-vs-LCL-S-G1_spline_spearman")
ymin <- -1.1
ymax <- 1.1
plotRD3vsRTALL(cors, file.name, main.text, ymin, ymax, cols=c("red", "blue", "black"), c("T", "N", "T/N"), c=NA, isRT=T)

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
# [1] -0.8457817
# > max(cors.samples[,-c(1:4)])
# [1] 0.8088429

#load(file.path(wd.rt.data, paste0("samples-vs-rt_sclc-vs-lcl_spline_spearman.RData")))
file.name <- file.path(wd.rt.plots, "SAMPLES-vs-RT_SCLC-vs-LCL_spline_spearman")
main.text <- c("SCLC read depth vs. LCL RT", "")
ymin <- -0.8773492
ymax <- 0.8392611
plotSAMPLEvsRTALL(cors.samples, samples1, file.name, main.text, ymin, ymax)

# -----------------------------------------------------------------------------
# Overall correlation with LCL S/G1
# Last Modified: 18/11/19; 16/06/19; 04/06/19; 06/03/19
# -----------------------------------------------------------------------------
samples.sclc <- setSamplesQ4(wd.rt.data, samples1)
writeTable(samples.sclc, file.path(wd.ngs, "sclc_wgs_n101.txt"), colnames=T, rownames=F, sep="\t")
#         0%        25%        50%        75%       100% 
# -0.7438588 -0.6840244 -0.6474195 -0.5555826  0.6806205 

writeTable(subset(samples.sclc, Q4 %in% c(4,1)), file.path(wd.ngs, "sclc_wgs_q4_n51.txt"), colnames=T, rownames=F, sep="\t")
writeTable(subset(samples.sclc, Q4 %in% c(3,1)), file.path(wd.ngs, "sclc_wgs_q3_n51.txt"), colnames=T, rownames=F, sep="\t")

## Random 25/26
m2.25 <- sort(rownames(subset(samples.sclc, M2 == 1))[sample(1:50, round(50/2), replace=F)])
m1.26 <- sort(rownames(subset(samples.sclc, M2 == 0))[sample(1:51, round(51/2), replace=F)])

random51 <- sort(c(m2.25, m1.26))
writeTable(samples.sclc[random51,], file.path(wd.ngs, "sclc_wgs_n101-random51.txt"), colnames=T, rownames=F, sep="\t")

random50 <- sort(setdiff(rownames(samples.sclc), random51))
writeTable(samples.sclc[random50,], file.path(wd.ngs, "sclc_wgs_n101-random50.txt"), colnames=T, rownames=F, sep="\t")

## Random 12/13
m2 <- sort(rownames(subset(samples.sclc, M2 == 1)))
m1 <- sort(rownames(subset(samples.sclc, M2 == 0)))

m2.12 <- sort(m2[sample(1:50, round(50/4), replace=F)])
m1.13 <- sort(m1[sample(1:51, round(51/4), replace=F)])
random3 <- sort(c(m2.12, m1.13))
writeTable(samples.sclc[random3,], file.path(wd.ngs, "sclc_wgs_n101-random3.txt"), colnames=T, rownames=F, sep="\t")

m2 <- sort(setdiff(m2, m2.12))
m1 <- sort(setdiff(m1, m1.13))

m2.12 <- sort(m2[sample(1:38, round(50/4), replace=F)])
m1.13 <- sort(m1[sample(1:38, round(51/4), replace=F)])
random4 <- sort(c(m2.12, m1.13))
writeTable(samples.sclc[random4,], file.path(wd.ngs, "sclc_wgs_n101-random4.txt"), colnames=T, rownames=F, sep="\t")

## Random 6/7
m2 <- sort(rownames(subset(samples.sclc, M2 == 1)))
m1 <- sort(rownames(subset(samples.sclc, M2 == 0)))

m2.6 <- sort(m2[sample(1:50, round(50/8), replace=F)])
m1.6 <- sort(m1[sample(1:51, round(51/8), replace=F)])
random5 <- sort(c(m2.6, m1.6))
writeTable(samples.sclc[random5,], file.path(wd.ngs, "sclc_wgs_n101-random5.txt"), colnames=T, rownames=F, sep="\t")

m2 <- sort(setdiff(m2, m2.6))
m1 <- sort(setdiff(m1, m1.6))

m2.6 <- sort(m2[sample(1:44, round(50/8), replace=F)])
m1.6 <- sort(m1[sample(1:45, round(51/8), replace=F)])
random6 <- sort(c(m2.6, m1.6))
writeTable(samples.sclc[random6,], file.path(wd.ngs, "sclc_wgs_n101-random6.txt"), colnames=T, rownames=F, sep="\t")




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
save(pca.de, file=file.path(wd.rt.data, paste0("pca_sclc_chrs.RData")))

#load(file.path(wd.rt.data, paste0("pca_sclc_chrs.RData")))
file.main <- c("SCLC (n=101) read depth profiles", "")
trait <- as.numeric(samples.sclc$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
plotPCA(1, 2, pca.de, trait, wd.rt.plots, "PCA_SCLC_chrs_lightpink1_lightskyblue2_Q4", size=6, file.main, "bottomright", c("red", "lightpink1", "lightskyblue2", "blue"), NULL, flip.x=1, flip.y=1, legend.title=NA)

## SG1
#trait <- samples.sclc.sg1$SG1
#plotPCA(1, 2, pca.de, trait, wd.rt.plots, "pca_sclc_T_chrs_spline_spearman_SG1", size=6, file.main, "bottomright", c("red", "lightgray", "blue"), NULL, flip.x=1, flip.y=1, legend.title="Consist. CM in all chrs")

# -----------------------------------------------------------------------------
# Beeswarm plots
# Last Modified: 21/04/19
# -----------------------------------------------------------------------------
samples.sclc <- readTable("/projects/cangen/tyang2/SCLC/ngs/WGS/sclc_wgs_n101.txt", header=T, rownames=T, sep="")
samples.nbl  <- readTable("/projects/cangen/tyang2/NBL/ngs/WGS/nbl_wgs_n57-1.txt", header=T, rownames=T, sep="")
samples.cll  <- readTable("/projects/cangen/tyang2/CLL/ngs/WGS/cll_wgs_n96.txt", header=T, rownames=T, sep="")
n.sclc <- nrow(samples.sclc)
n.nbl  <- nrow(samples.nbl)
n.cll  <- nrow(samples.cll)

samples <- toTable(0, 3, n.sclc+n.nbl+n.cll, c("CANCER", "COR", "Q4"))
samples$CANCER[1:n.sclc] <- 0
samples$CANCER[(1+n.sclc):(n.sclc+n.nbl)] <- 1
samples$CANCER[(1+n.sclc+n.nbl):(n.sclc+n.nbl+n.cll)] <- 2
samples$COR <- c(samples.sclc$COR, samples.nbl$COR, samples.cll$COR)
samples$Q4  <- c(samples.sclc$Q4, samples.nbl$Q4, samples.cll$Q4)

#install.packages('beeswarm')
library(beeswarm)

pdf(file.path(wd.rt.plots, "beeswarm_sclc+nbl+cll_1.3.pdf"), height=6, width=6)
ymax <- 0.8   #max(samples$COR)
ymin <- -ymax
boxplot(COR ~ CANCER, data=samples, outline=F, names=c("", "", ""), ylim=c(ymin, ymax), ylab="", main="Overall correlation with LCL RT", yaxt="n", cex.axis=1.7, cex.lab=1.7, cex.main=1.8)
abline(h=0, lty=5)

beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 1), col="blue", pch=16, add=T)
beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 2), col="lightskyblue2", pch=16, add=T)
beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 3), col="lightpink1", pch=16, add=T)
beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 4), col="red", pch=16, add=T)

legend("topright", legend = c("Q4", "Q3", "Q2", "Q1"), pch=16, col=c("red", "lightpink1", "lightskyblue2", "blue"), cex=1.7)

axis(side=2, at=seq(-0.8, 0.8, by=0.4), labels=c(-0.8, -0.4, 0, 0.4, 0.8), cex.axis=1.6)
mtext("Spearman's rho", side=2, line=2.8, cex=1.7)
#mtext("", cex=1.2, line=0.3)
mtext(text=c("SCLC", "NBL", "CLL"), side=1, cex=1.7, line=1.3, at=c(1,2,3))
dev.off()

# -----------------------------------------------------------------------------
# PCA for all
# Last Modified: 11/03/20
# -----------------------------------------------------------------------------
# > dim(test.sclc)
# [1] 2657164     102
# > dim(test.nbl)
# [1] 2659570      57
# > dim(test.cll)
# [1] 2653842      97
# > dim(test.sclc.nl)
# [1] 2658679      93

overlaps <- intersect(intersect(intersect(rownames(test.sclc), rownames(test.nbl)), rownames(test.cll)), rownames(test.sclc.nl))
# > length(overlaps)
# [1] 2653842
test.sclc <- test.sclc[overlaps, -1]
test.nbl  <- test.nbl[overlaps, -1]
test.cll  <- test.cll[overlaps, -1]
test.sclc.nl <- test.sclc.nl[overlaps, -1]

test <- cbind(test.sclc, test.nbl, test.cll, test.sclc.nl)
pca.de <- getPCA(t(test))
save(pca.de, file=file.path(wd.rt.data, paste0("pca_sclc+nbl_cll+sclc.nl.RData")))

# -----------------------------------------------------------------------------
# PCA ALL+WB
# Last Modified: 12/03/20
# -----------------------------------------------------------------------------
#load(file.path(wd.rt.data, paste0("pca_sclc+nbl+cll+sclc.nl+nbl.wb+cll.wb_pca.de.RData")))

samples.sclc <- readTable("/projects/cangen/tyang2/SCLC/ngs/WGS/sclc_wgs_n101.txt", header=T, rownames=T, sep="")
samples.nbl  <- readTable("/projects/cangen/tyang2/NBL/ngs/WGS/nbl_wgs_n57-1.txt", header=T, rownames=T, sep="")
samples.cll  <- readTable("/projects/cangen/tyang2/CLL/ngs/WGS/cll_wgs_n96.txt", header=T, rownames=T, sep="")
samples.sclc.nl <- readTable("/projects/cangen/tyang2/SCLC/ngs/WGS/sclc_wgs_n92_NL.txt", header=T, rownames=T, sep="")
samples.nbl.wb  <- readTable("/projects/cangen/tyang2/NBL/ngs/WGS/nbl_wgs_n57-1_N.txt", header=T, rownames=T, sep="")
samples.cll.wb  <- readTable("/projects/cangen/tyang2/CLL/ngs/WGS/cll_wgs_n96_N.txt", header=T, rownames=T, sep="")
n.sclc <- nrow(samples.sclc)
n.nbl  <- nrow(samples.nbl)
n.cll  <- nrow(samples.cll)
n.sclc.nl <- nrow(samples.sclc.nl)
n.nbl.wb  <- nrow(samples.nbl.wb)
n.cll.wb  <- nrow(samples.cll.wb)

samples <- toTable(0, 3, n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb+n.cll.wb, c("CANCER", "COR", "Q4"))
samples$CANCER[1:n.sclc] <- 0
samples$CANCER[(1+n.sclc):(n.sclc+n.nbl)] <- 1
samples$CANCER[(1+n.sclc+n.nbl):(n.sclc+n.nbl+n.cll)] <- 2
samples$CANCER[(1+n.sclc+n.nbl+n.cll):(n.sclc+n.nbl+n.cll+n.sclc.nl)] <- 3
samples$CANCER[(1+n.sclc+n.nbl+n.cll+n.sclc.nl):(n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb)] <- 4
samples$CANCER[(1+n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb):(n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb+n.cll.wb)] <- 5
samples$COR <- c(samples.sclc$COR, samples.nbl$COR, samples.cll$COR, samples.sclc.nl$COR, samples.nbl.wb$COR, samples.cll.wb$COR)
samples$Q4  <- c(samples.sclc$Q4, samples.nbl$Q4, samples.cll$Q4, samples.sclc.nl$Q4, samples.nbl.wb$Q4, samples.cll.wb$Q4)
rownames(samples) <- c(rownames(samples.sclc), rownames(samples.nbl), rownames(samples.cll), paste0(rownames(samples.sclc.nl), "-NL"), paste0(rownames(samples.nbl.wb), "-WB"), paste0(rownames(samples.cll.wb), "-WB"))
samples$SAMPLE_ID <- rownames(samples)

## ALL
q <- quantile(as.numeric(samples$COR))
print(q)

samples.q4 <- list()
samples.q4[[4]] <- rownames(subset(samples, COR > as.numeric(q[4])))
samples.q4[[3]] <- rownames(subset(subset(samples, COR > as.numeric(q[3])), COR <= as.numeric(q[4])))
samples.q4[[2]] <- rownames(subset(subset(samples, COR > as.numeric(q[2])), COR <= as.numeric(q[3])))
samples.q4[[1]] <- rownames(subset(samples, COR <= as.numeric(q[2])))

samples$Q4[which(samples$SAMPLE_ID %in% samples.q4[[4]])] <- 4
samples$Q4[which(samples$SAMPLE_ID %in% samples.q4[[3]])] <- 3
samples$Q4[which(samples$SAMPLE_ID %in% samples.q4[[2]])] <- 2
samples$Q4[which(samples$SAMPLE_ID %in% samples.q4[[1]])] <- 1

file.main <- c("Total (n=497) tumour-normal profiles ", "")
trait <- as.numeric(samples$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
plotPCA(1, 3, pca.de, trait, wd.rt.plots, "PCA_ALL+WB", size=6, file.main, "topleft", c("red", "lightpink1", "lightskyblue2", "blue", "lightgray"), NULL, flip.x=-1, flip.y=1, legend.title=NA, "")

###
## SCLC
samples <- toTable(0, 3, n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb+n.cll.wb, c("CANCER", "COR", "Q4"))
samples$CANCER[1:n.sclc] <- 0
samples$CANCER[(1+n.sclc):(n.sclc+n.nbl)] <- 1
samples$CANCER[(1+n.sclc+n.nbl):(n.sclc+n.nbl+n.cll)] <- 2
samples$CANCER[(1+n.sclc+n.nbl+n.cll):(n.sclc+n.nbl+n.cll+n.sclc.nl)] <- 3
samples$CANCER[(1+n.sclc+n.nbl+n.cll+n.sclc.nl):(n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb)] <- 4
samples$CANCER[(1+n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb):(n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb+n.cll.wb)] <- 5
samples$COR <- c(samples.sclc$COR, samples.nbl$COR, samples.cll$COR, samples.sclc.nl$COR, samples.nbl.wb$COR, samples.cll.wb$COR)
samples$Q4  <- c(samples.sclc$Q4, samples.nbl$Q4, samples.cll$Q4, samples.sclc.nl$Q4, samples.nbl.wb$Q4, samples.cll.wb$Q4)

samples[which(samples$CANCER != 0),]$Q4 <- NA
file.main <- c("SCLC (n=101) read depth profiles", "")
trait <- as.numeric(samples$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
plotPCA(1, 3, pca.de, trait, wd.rt.plots, "PCA_ALL_SCLC", size=6, file.main, "topleft", c("red", "lightpink1", "lightskyblue2", "blue", "lightgray"), NULL, flip.x=-1, flip.y=1, legend.title=NA, "SCLC")

## NBL
samples <- toTable(0, 3, n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb+n.cll.wb, c("CANCER", "COR", "Q4"))
samples$CANCER[1:n.sclc] <- 0
samples$CANCER[(1+n.sclc):(n.sclc+n.nbl)] <- 1
samples$CANCER[(1+n.sclc+n.nbl):(n.sclc+n.nbl+n.cll)] <- 2
samples$CANCER[(1+n.sclc+n.nbl+n.cll):(n.sclc+n.nbl+n.cll+n.sclc.nl)] <- 3
samples$CANCER[(1+n.sclc+n.nbl+n.cll+n.sclc.nl):(n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb)] <- 4
samples$CANCER[(1+n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb):(n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb+n.cll.wb)] <- 5
samples$COR <- c(samples.sclc$COR, samples.nbl$COR, samples.cll$COR, samples.sclc.nl$COR, samples.nbl.wb$COR, samples.cll.wb$COR)
samples$Q4  <- c(samples.sclc$Q4, samples.nbl$Q4, samples.cll$Q4, samples.sclc.nl$Q4, samples.nbl.wb$Q4, samples.cll.wb$Q4)

samples[which(samples$CANCER != 1),]$Q4 <- NA
file.main <- c("NBL (n=56) read depth profiles", "")
trait <- as.numeric(samples$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
plotPCA(1, 3, pca.de, trait, wd.rt.plots, "PCA_ALL_NBL", size=6, file.main, "topleft", c("red", "lightpink1", "lightskyblue2", "blue", "lightgray"), NULL, flip.x=-1, flip.y=1, legend.title=NA, "NBL")

## CLL
samples <- toTable(0, 3, n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb+n.cll.wb, c("CANCER", "COR", "Q4"))
samples$CANCER[1:n.sclc] <- 0
samples$CANCER[(1+n.sclc):(n.sclc+n.nbl)] <- 1
samples$CANCER[(1+n.sclc+n.nbl):(n.sclc+n.nbl+n.cll)] <- 2
samples$CANCER[(1+n.sclc+n.nbl+n.cll):(n.sclc+n.nbl+n.cll+n.sclc.nl)] <- 3
samples$CANCER[(1+n.sclc+n.nbl+n.cll+n.sclc.nl):(n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb)] <- 4
samples$CANCER[(1+n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb):(n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb+n.cll.wb)] <- 5
samples$COR <- c(samples.sclc$COR, samples.nbl$COR, samples.cll$COR, samples.sclc.nl$COR, samples.nbl.wb$COR, samples.cll.wb$COR)
samples$Q4  <- c(samples.sclc$Q4, samples.nbl$Q4, samples.cll$Q4, samples.sclc.nl$Q4, samples.nbl.wb$Q4, samples.cll.wb$Q4)

samples[which(samples$CANCER != 2),]$Q4 <- NA
file.main <- c("CLL (n=96) read depth profiles", "")
trait <- as.numeric(samples$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
plotPCA(1, 3, pca.de, trait, wd.rt.plots, "PCA_ALL_CLL", size=6, file.main, "topleft", c("red", "lightpink1", "lightskyblue2", "blue", "lightgray"), NULL, flip.x=-1, flip.y=1, legend.title=NA, "CLL")

## SCLC-NL
samples <- toTable(0, 3, n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb+n.cll.wb, c("CANCER", "COR", "Q4"))
samples$CANCER[1:n.sclc] <- 0
samples$CANCER[(1+n.sclc):(n.sclc+n.nbl)] <- 1
samples$CANCER[(1+n.sclc+n.nbl):(n.sclc+n.nbl+n.cll)] <- 2
samples$CANCER[(1+n.sclc+n.nbl+n.cll):(n.sclc+n.nbl+n.cll+n.sclc.nl)] <- 3
samples$CANCER[(1+n.sclc+n.nbl+n.cll+n.sclc.nl):(n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb)] <- 4
samples$CANCER[(1+n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb):(n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb+n.cll.wb)] <- 5
samples$COR <- c(samples.sclc$COR, samples.nbl$COR, samples.cll$COR, samples.sclc.nl$COR, samples.nbl.wb$COR, samples.cll.wb$COR)
samples$Q4  <- c(samples.sclc$Q4, samples.nbl$Q4, samples.cll$Q4, samples.sclc.nl$Q4, samples.nbl.wb$Q4, samples.cll.wb$Q4)

samples[which(samples$CANCER != 3),]$Q4 <- NA
file.main <- c("SCLC-NL (n=92) read depth profiles", "")
trait <- as.numeric(samples$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
plotPCA(1, 3, pca.de, trait, wd.rt.plots, "PCA_ALL_SCLC-NL", size=6, file.main, "topleft", c("red", "lightpink1", "lightskyblue2", "blue", "lightgray"), NULL, flip.x=-1, flip.y=1, legend.title=NA, "SCLC-NL")

## NBL-WB
samples <- toTable(0, 3, n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb+n.cll.wb, c("CANCER", "COR", "Q4"))
samples$CANCER[1:n.sclc] <- 0
samples$CANCER[(1+n.sclc):(n.sclc+n.nbl)] <- 1
samples$CANCER[(1+n.sclc+n.nbl):(n.sclc+n.nbl+n.cll)] <- 2
samples$CANCER[(1+n.sclc+n.nbl+n.cll):(n.sclc+n.nbl+n.cll+n.sclc.nl)] <- 3
samples$CANCER[(1+n.sclc+n.nbl+n.cll+n.sclc.nl):(n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb)] <- 4
samples$CANCER[(1+n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb):(n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb+n.cll.wb)] <- 5
samples$COR <- c(samples.sclc$COR, samples.nbl$COR, samples.cll$COR, samples.sclc.nl$COR, samples.nbl.wb$COR, samples.cll.wb$COR)
samples$Q4  <- c(samples.sclc$Q4, samples.nbl$Q4, samples.cll$Q4, samples.sclc.nl$Q4, samples.nbl.wb$Q4, samples.cll.wb$Q4)

samples[which(samples$CANCER != 4),]$Q4 <- NA
file.main <- c("NBL-WB (n=56) read depth profiles", "")
trait <- as.numeric(samples$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
plotPCA(1, 3, pca.de, trait, wd.rt.plots, "PCA_ALL_NBL-WB", size=6, file.main, "topleft", c("red", "lightpink1", "lightskyblue2", "blue", "lightgray"), NULL, flip.x=-1, flip.y=1, legend.title=NA, "NBL-WB")

## CLL
samples <- toTable(0, 3, n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb+n.cll.wb, c("CANCER", "COR", "Q4"))
samples$CANCER[1:n.sclc] <- 0
samples$CANCER[(1+n.sclc):(n.sclc+n.nbl)] <- 1
samples$CANCER[(1+n.sclc+n.nbl):(n.sclc+n.nbl+n.cll)] <- 2
samples$CANCER[(1+n.sclc+n.nbl+n.cll):(n.sclc+n.nbl+n.cll+n.sclc.nl)] <- 3
samples$CANCER[(1+n.sclc+n.nbl+n.cll+n.sclc.nl):(n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb)] <- 4
samples$CANCER[(1+n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb):(n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb+n.cll.wb)] <- 5
samples$COR <- c(samples.sclc$COR, samples.nbl$COR, samples.cll$COR, samples.sclc.nl$COR, samples.nbl.wb$COR, samples.cll.wb$COR)
samples$Q4  <- c(samples.sclc$Q4, samples.nbl$Q4, samples.cll$Q4, samples.sclc.nl$Q4, samples.nbl.wb$Q4, samples.cll.wb$Q4)

samples[which(samples$CANCER != 5),]$Q4 <- NA
file.main <- c("CLL-WB (n=96) read depth profiles", "")
trait <- as.numeric(samples$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
plotPCA(1, 3, pca.de, trait, wd.rt.plots, "PCA_ALL_CLL-WB", size=6, file.main, "topleft", c("red", "lightpink1", "lightskyblue2", "blue", "lightgray"), NULL, flip.x=-1, flip.y=1, legend.title=NA, "CLL-WB")











# -----------------------------------------------------------------------------
# PCA ALL
# Last Modified: 11/03/20
# -----------------------------------------------------------------------------
#load(file.path(wd.rt.data, paste0("pca_sclc+nbl+cll+sclc.nl_pca.de.RData")))

samples.sclc <- readTable("/projects/cangen/tyang2/SCLC/ngs/WGS/sclc_wgs_n101.txt", header=T, rownames=T, sep="")
samples.nbl  <- readTable("/projects/cangen/tyang2/NBL/ngs/WGS/nbl_wgs_n57-1.txt", header=T, rownames=T, sep="")
samples.cll  <- readTable("/projects/cangen/tyang2/CLL/ngs/WGS/cll_wgs_n96.txt", header=T, rownames=T, sep="")
samples.sclc.nl <- readTable("/projects/cangen/tyang2/SCLC/ngs/WGS/sclc_wgs_n92_NL.txt", header=T, rownames=T, sep="")
n.sclc <- nrow(samples.sclc)
n.nbl  <- nrow(samples.nbl)
n.cll  <- nrow(samples.cll)
n.sclc.nl <- nrow(samples.sclc.nl)

samples <- toTable(0, 3, n.sclc+n.nbl+n.cll+n.sclc.nl, c("CANCER", "COR", "Q4"))
samples$CANCER[1:n.sclc] <- 0
samples$CANCER[(1+n.sclc):(n.sclc+n.nbl)] <- 1
samples$CANCER[(1+n.sclc+n.nbl):(n.sclc+n.nbl+n.cll)] <- 2
samples$CANCER[(1+n.sclc+n.nbl+n.cll):(n.sclc+n.nbl+n.cll+n.sclc.nl)] <- 3
samples$COR <- c(samples.sclc$COR, samples.nbl$COR, samples.cll$COR, samples.sclc.nl$COR)
samples$Q4  <- c(samples.sclc$Q4, samples.nbl$Q4, samples.cll$Q4, samples.sclc.nl$Q4)
rownames(samples) <- c(rownames(samples.sclc), rownames(samples.nbl), rownames(samples.cll), rownames(samples.sclc.nl))
samples$SAMPLE_ID <- rownames(samples)
 
## ALL
q <- quantile(as.numeric(samples$COR))
print(q)

samples.q4 <- list()
samples.q4[[4]] <- rownames(subset(samples, COR > as.numeric(q[4])))
samples.q4[[3]] <- rownames(subset(subset(samples, COR > as.numeric(q[3])), COR <= as.numeric(q[4])))
samples.q4[[2]] <- rownames(subset(subset(samples, COR > as.numeric(q[2])), COR <= as.numeric(q[3])))
samples.q4[[1]] <- rownames(subset(samples, COR <= as.numeric(q[2])))

samples$Q4[which(samples$SAMPLE_ID %in% samples.q4[[4]])] <- 4
samples$Q4[which(samples$SAMPLE_ID %in% samples.q4[[3]])] <- 3
samples$Q4[which(samples$SAMPLE_ID %in% samples.q4[[2]])] <- 2
samples$Q4[which(samples$SAMPLE_ID %in% samples.q4[[1]])] <- 1

file.main <- c("SCLC-NL, SCLC, NBL and CLL (n=345)", "")
trait <- as.numeric(samples$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
plotPCA(1, 2, pca.de, trait, wd.rt.plots, "PCA_ALL_Q4", size=6, file.main, "topleft", c("red", "lightpink1", "lightskyblue2", "blue", "lightgray"), NULL, flip.x=-1, flip.y=1, legend.title=NA)

###
## SCLC
samples <- toTable(0, 3, n.sclc+n.nbl+n.cll+n.sclc.nl, c("CANCER", "COR", "Q4"))
samples$CANCER[1:n.sclc] <- 0
samples$CANCER[(1+n.sclc):(n.sclc+n.nbl)] <- 1
samples$CANCER[(1+n.sclc+n.nbl):(n.sclc+n.nbl+n.cll)] <- 2
samples$CANCER[(1+n.sclc+n.nbl+n.cll):(n.sclc+n.nbl+n.cll+n.sclc.nl)] <- 3
samples$COR <- c(samples.sclc$COR, samples.nbl$COR, samples.cll$COR, samples.sclc.nl$COR)
samples$Q4  <- c(samples.sclc$Q4, samples.nbl$Q4, samples.cll$Q4, samples.sclc.nl$Q4)

samples[which(samples$CANCER != 0),]$Q4 <- NA
file.main <- c("SCLC (n=101) read depth profiles", "")
trait <- as.numeric(samples$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
plotPCA(1, 2, pca.de, trait, wd.rt.plots, "PCA_ALL_SCLC_Q4", size=6, file.main, "topleft", c("red", "lightpink1", "lightskyblue2", "blue", "lightgray"), NULL, flip.x=-1, flip.y=1, legend.title=NA, "SCLC")

## NBL
samples <- toTable(0, 3, n.sclc+n.nbl+n.cll+n.sclc.nl, c("CANCER", "COR", "Q4"))
samples$CANCER[1:n.sclc] <- 0
samples$CANCER[(1+n.sclc):(n.sclc+n.nbl)] <- 1
samples$CANCER[(1+n.sclc+n.nbl):(n.sclc+n.nbl+n.cll)] <- 2
samples$CANCER[(1+n.sclc+n.nbl+n.cll):(n.sclc+n.nbl+n.cll+n.sclc.nl)] <- 3
samples$COR <- c(samples.sclc$COR, samples.nbl$COR, samples.cll$COR, samples.sclc.nl$COR)
samples$Q4  <- c(samples.sclc$Q4, samples.nbl$Q4, samples.cll$Q4, samples.sclc.nl$Q4)

samples[which(samples$CANCER != 1),]$Q4 <- NA
file.main <- c("NBL (n=56) read depth profiles", "")
trait <- as.numeric(samples$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
plotPCA(1, 2, pca.de, trait, wd.rt.plots, "PCA_ALL_NBL_Q4", size=6, file.main, "topleft", c("red", "lightpink1", "lightskyblue2", "blue", "lightgray"), NULL, flip.x=-1, flip.y=1, legend.title=NA, "NBL")

## CLL
samples <- toTable(0, 3, n.sclc+n.nbl+n.cll+n.sclc.nl, c("CANCER", "COR", "Q4"))
samples$CANCER[1:n.sclc] <- 0
samples$CANCER[(1+n.sclc):(n.sclc+n.nbl)] <- 1
samples$CANCER[(1+n.sclc+n.nbl):(n.sclc+n.nbl+n.cll)] <- 2
samples$CANCER[(1+n.sclc+n.nbl+n.cll):(n.sclc+n.nbl+n.cll+n.sclc.nl)] <- 3
samples$COR <- c(samples.sclc$COR, samples.nbl$COR, samples.cll$COR, samples.sclc.nl$COR)
samples$Q4  <- c(samples.sclc$Q4, samples.nbl$Q4, samples.cll$Q4, samples.sclc.nl$Q4)

samples[which(samples$CANCER != 2),]$Q4 <- NA
file.main <- c("CLL (n=96) read depth profiles", "")
trait <- as.numeric(samples$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
plotPCA(1, 2, pca.de, trait, wd.rt.plots, "PCA_ALL_CLL_Q4", size=6, file.main, "topleft", c("red", "lightpink1", "lightskyblue2", "blue", "lightgray"), NULL, flip.x=-1, flip.y=1, legend.title=NA, "CLL")

## SCLC-NL
samples <- toTable(0, 3, n.sclc+n.nbl+n.cll+n.sclc.nl, c("CANCER", "COR", "Q4"))
samples$CANCER[1:n.sclc] <- 0
samples$CANCER[(1+n.sclc):(n.sclc+n.nbl)] <- 1
samples$CANCER[(1+n.sclc+n.nbl):(n.sclc+n.nbl+n.cll)] <- 2
samples$CANCER[(1+n.sclc+n.nbl+n.cll):(n.sclc+n.nbl+n.cll+n.sclc.nl)] <- 3
samples$COR <- c(samples.sclc$COR, samples.nbl$COR, samples.cll$COR, samples.sclc.nl$COR)
samples$Q4  <- c(samples.sclc$Q4, samples.nbl$Q4, samples.cll$Q4, samples.sclc.nl$Q4)

samples[which(samples$CANCER != 3),]$Q4 <- NA
file.main <- c("SCLC-NL (n=92) read depth profiles", "")
trait <- as.numeric(samples$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
plotPCA(1, 2, pca.de, trait, wd.rt.plots, "PCA_ALL_SCLC-NL_Q4", size=6, file.main, "topleft", c("red", "lightpink1", "lightskyblue2", "blue", "lightgray"), NULL, flip.x=-1, flip.y=1, legend.title=NA, "SCLC-NL")






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
