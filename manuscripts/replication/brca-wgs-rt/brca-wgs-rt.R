# =============================================================================
# Manuscript   : 
# Chapter      : 
# Name         : manuscripts/replication/brca-wgs-rt.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 29/01/21
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
# Last Modified: 04/07/19
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
#wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE  <- "BRCA"
PAIR1 <- "T"
base  <- tolower(BASE)
method <- "rpkm"

wd.ngs   <- file.path(wd, BASE, "ngs/WGS")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.rt       <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

wd.ngs.data <- file.path(wd.ngs, "data")
samples1 <- readTable(file.path(wd.ngs, "brca_wgs_n18.txt"), header=T, rownames=T, sep="\t")
n1 <- nrow(samples1)

# -----------------------------------------------------------------------------
# 
# Last Modified: 29/01/21
# -----------------------------------------------------------------------------
samples1$V5 <- paste0(samples1$V2, "_", samples1$V3, "_input.hg19.merged.nodup.", samples1$V4, ".bam")
writeTable(samples1, file.path(wd.ngs, "brca_wgs_n27.bam"), colnames=F, rownames=F, sep="\t")

samples3 <- readTable(file.path(wd.ngs, "brca_wgs_n27.list3"), header=F, rownames=F, sep="")
samples2 <- setdiff(samples1$V1, samples3)
writeTable(samples2, file.path(wd.ngs, "brca_wgs_n27.list2"), colnames=F, rownames=F, sep="\t")

samples4 <- readTable(file.path(wd.ngs, "brca_wgs_n18.list"), header=F, rownames=F, sep="")
samples1 <- samples1[samples4,]

# -----------------------------------------------------------------------------
# BRCA vs LCL S/G1
# http://pklab.med.harvard.edu/scw2014/subpop_tutorial.html
# Last Modified: 30/01/21
# -----------------------------------------------------------------------------
cors.samples <- getSAMPLEvsRT(wd.rt.data, samples1[,1])
save(cors.samples, file=file.path(wd.rt.data, paste0("samples-vs-rt_", base, "-vs-lcl_spline_spearman.RData")))
# > min(cors.samples[,-c(1:4)])
# [1] -0.835043
# > max(cors.samples[,-c(1:4)])
# [1] 0.7159929

#load(file.path(wd.rt.data, paste0("samples-vs-rt_brca-vs-lcl_spline_spearman.RData")))
file.name <- file.path(wd.rt.plots, "SAMPLES-vs-RT_BRCA-CL-vs-LCL_spline_spearman")
main.text <- c("BRCA read depth vs. LCL S/G1", "")
ymin <- -0.8773492
ymax <- 0.8392611
plotSAMPLEvsRT(cors.samples, samples1[,1], file.name, main.text, ymin, ymax)

# -----------------------------------------------------------------------------
# Overall correlation with LCL S/G1
# Last Modified: 30/01/21
# -----------------------------------------------------------------------------
samples.brca <- setSamplesQ4(wd.rt.data, samples1[,1])
samples1[,c("COR", "M2", "Q4")] <- samples.brca[,c("COR", "M2", "Q4")] 
writeTable(samples1, file.path(wd.ngs, "brca_wgs_n18.txt"), colnames=T, rownames=F, sep="\t")
#          0%         25%         50%         75%        100% 
# -0.28136620 -0.14326296 -0.07182227  0.08254313  0.19817708 

# -----------------------------------------------------------------------------
# Last Modified: 30/11/20; 04/06/19; 21/04/19
# -----------------------------------------------------------------------------
nrds <- getOverallRD(wd.rt.data, base, method, PAIR1, n1)
# > nrow(nrds)
# [1] 2684771
nrds.m <- subset(nrds, MEDIAN != 0)   ## ADD 30/11/20
# > nrow(nrds.m)
# [1] 2634215
nrds.d <- getOverallDetectedRD(wd.ngs.data, base, method, PAIR1)
nrow(nrds.d)
# > nrow(nrds.d)
# [1] 2404209

# -----------------------------------------------------------------------------
# PCA
# Last Modified: 15/11/19; 04/06/19; 21/04/19
# -----------------------------------------------------------------------------
## Refer to cmd-rt_2a_nrd.gc.cn.d_sample.R (commandline mode)
nrds.T.chr.d.all <- pipeGetDetectedRD(wd.ngs.data, BASE, "chr1", PAIR1, method, samples1[,1])
#nrds.T.chr.d.all <- readTable(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d_", "chr1", "_", PAIR1, "_n", n1, ".txt.gz")), header=T, rownames=F, sep="\t")
#colnames(nrds.T.chr.d.all) <- toupper(gsub("\\.", "", colnames(nrds.T.chr.d.all)))
for (c in 2:22) {
   chr <- chrs[c]
 
   ## Read depth
   nrds.T.chr.d <- pipeGetDetectedRD(wd.ngs.data, BASE, chr, PAIR1, method)
   #nrds.T.chr.d <- readTable(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d_", chr, "_", PAIR1, "_n", n1, ".txt.gz")), header=T, rownames=F, sep="\t")
   #colnames(nrds.T.chr.d) <- toupper(gsub("\\.", "", colnames(nrds.T.chr.d)))
   
   nrds.T.chr.d.all <- rbind(nrds.T.chr.d.all, nrds.T.chr.d)
}

##
test <- nrds.T.chr.d.all[, -1]   ## BUG 2019/10/14: Remove column BED
pca.de <- getPCA(t(test))
save(pca.de, file=file.path(wd.rt.data, paste0("pca_nbl-cl_chrs.RData")))

#load(file.path(wd.rt.data, paste0("pca_nbl-cl_chrs.RData")))
file.main <- c("NBL-CL (n=8) read depth profiles", "")
trait <- as.numeric(samples.nbl.cl$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
plotPCA(1, 2, pca.de, trait, wd.rt.plots, "PCA_NBL-CL", size=6, file.main, "top", c("red", "lightpink1", "lightskyblue2", "blue"), samples1, flip.x=-1, flip.y=1, legend.title=NA)

# -----------------------------------------------------------------------------
# Beeswarm plots
# Last Modified: 21/04/19
# -----------------------------------------------------------------------------
samples <- samples1
n.brca  <- nrow(samples)

#samples <- toTable(0, 4, n.nbl.cl, c("CANCER", "COR", "Q4", "SAMPLE_ID"))
samples$CANCER[1:n.brca] <- 1
#samples$COR <- samples$COR
#samples$Q4  <- samples$Q4
#samples$SAMPLE_ID <- samples$SAMPLE_ID
#rownames(samples) <- samples$SAMPLE_ID

pdf(file.path(wd.rt.plots, "boxplot_brca_n15.pdf"), height=6, width=4.2)
ymax <- 0.25
ymin <- -0.3
boxplot(COR ~ CANCER, data=samples, outline=F, names=c(""), ylim=c(ymin, ymax), ylab="", main="In silico prediction", yaxt="n", boxwex=0.75, cex.axis=1.5, cex.lab=1.6, cex.main=1.7)
abline(h=0, lty=5, lwd=2)

points(subset(samples, Q4 == 2)$CANCER, subset(samples, Q4 == 2)$COR, col=lighterblue, pch=19, cex=2)
points(subset(samples, Q4 == 1)$CANCER, subset(samples, Q4 == 1)$COR, col=blue, pch=19, cex=2)
points(subset(samples, Q4 == 3)$CANCER, subset(samples, Q4 == 3)$COR, col=lighterred, pch=19, cex=2)
points(subset(samples, Q4 == 4)$CANCER, subset(samples, Q4 == 4)$COR, col=red, pch=19, cex=2)
for (s in 1:nrow(samples)) {
   sample <- samples[s,]

   if (sample$SAMPLE_ID == "VHIO179-1" || sample$SAMPLE_ID == "AB521M" || sample$SAMPLE_ID == "PAR1006")
      text(sample$CANCER, sample$COR, sample$SAMPLE_ID2, col="black", adj=c(1.22, -0.55), cex=1.5)
   else if (sample$SAMPLE_ID == "AB577M")
      text(sample$CANCER, sample$COR, sample$SAMPLE_ID2, col="black", adj=c(1.23, -0.55), cex=1.5)
   else if (sample$SAMPLE_ID == "AB790")
      text(sample$CANCER, sample$COR, sample$SAMPLE_ID2, col="black", adj=c(1.3, 0.5), cex=1.5)
   else if (sample$SAMPLE_ID == "AB555M")
      text(sample$CANCER, sample$COR, sample$SAMPLE_ID2, col="black", adj=c(1.22, 1.17), cex=1.5)
   else if (sample$SAMPLE_ID == "STG139M-1" || sample$SAMPLE_ID == "STG282M")
      text(sample$CANCER, sample$COR, sample$SAMPLE_ID2, col="black", adj=c(1.19, 1.17), cex=1.5)
   else if (sample$SAMPLE_ID == "STG139M-2")
      text(sample$CANCER, sample$COR, sample$SAMPLE_ID2, col="black", adj=c(1.19, 0.5), cex=1.5)
   else
      text(sample$CANCER, sample$COR, sample$SAMPLE_ID2, col="black", adj=c(1.22, 0.5), cex=1.5)
}

legend("topright", legend = c("Q4", "Q3", "Q2", "Q1"), pch=19, pt.cex=2.5, col=c(red, lighterred, lighterblue, blue), cex=1.5)

axis(side=2, at=seq(-0.2, 0.2, by=0.2), labels=c(-0.2, 0, 0.2), cex.axis=1.5)
axis(side=2, at=seq(-0.1, 0.1, by=0.2), labels=c(-0.1, 0.1), cex.axis=1.5)
mtext("Overall read depth vs. RT [rho]", side=2, line=2.75, cex=1.6)
#mtext("", cex=1.2, line=0.3)
axis(side=1, at=1, labels="BRCA", cex.axis=1.6)
#mtext(text=c(), side=1, cex=1.4, line=0.9, at=c(1,2,3))
mtext(text=c("n=15"), side=1, cex=1.6, line=2.3, at=c(1,2,3))
dev.off()

# -----------------------------------------------------------------------------
# Plot signal-to-noice
# Last Modified: 10/12/19
# -----------------------------------------------------------------------------
plotFACS <- function(n1, snr1, n2, snr2, file.name, main.text, xlab.text, ylab.text, col, pos) {
   n <- c(n1, n2)
   snr <- c(snr1, snr2)
 
   unit <- (max(snr) - min(snr))/10
   #xlim <- c(min(snr) - unit, max(snr) + unit)
   xlim <- c(0, 100)
   unit <- (max(n) - min(n))/10
   #ylim <- c(min(n) - unit, max(n) + unit)
   ylim <- c(-0.367, 0.5)
 
   pdf(paste0(file.name, ".pdf"), height=6, width=6)
   plot(n1 ~ snr1, ylim=ylim, xlim=xlim, ylab="", xlab=xlab.text, main=main.text[1], col=col[1], pch=19, cex=2, cex.axis=1.5, cex.lab=1.6, cex.main=1.7)
   points(n2 ~ snr2, col=col[2], pch=19, cex=2)
 
   lm.fit <- lm(n1 ~ snr1)
   abline(lm.fit, col=col[1], lwd=3)
   lm.fit <- lm(n2 ~ snr2)
   abline(lm.fit, col=col[2], lwd=3)
 
   #samples <- c("SCLC-NL", "SCLC", "NBL", "CLL")
   #text(snr1, n1, samples, col=col[1], pos=c(1,2,2,2), cex=1.8)
   #text(snr2, n2, samples, col=col[2], pos=c(1,4,4,4), cex=1.8)
 
   cor <- cor.test(n1, snr1, method="spearman", exact=F)
   #legend(pos[1], c(paste0("rho = ", round0(cor[[4]], digits=1)), paste0("p-value = ", scientific(cor[[3]], digits=1))), text.col=col[1], bty="n", cex=1.75)
   legend(pos[1], paste0("G1 (rho = ", round0(cor[[4]], digits=1), ")      "), text.col=col[1], pch=c(NA), col=col[1], bty="n", cex=1.5)
 
   cor <- cor.test(n2, snr2, method="spearman", exact=F)
   #legend(pos[2], c(paste0("rho = ", round0(cor[[4]], digits=1)), paste0("p-value = ", scientific(cor[[3]], digits=1))), text.col=col[2], bty="n", cex=1.75)
   legend(pos[2], paste0("S (rho = ", round0(cor[[4]], digits=1), ")"), text.col=col[2], pch=c(NA), col=col[2], bty="n", cex=1.5)

   mtext(ylab.text, side=2, line=2.75, cex=1.6)
   #mtext(main.text[2], cex=1.2, line=0.3)
   dev.off()
}

# https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
plotFACS3 <- function(n3, snr3, file.name, main.text, xlab.text, ylab.text, col, col2, pos, xlim.max) {
   xlim <- c(0, xlim.max)
   ylim <- c(-0.3, 0.25)
 
   pdf(paste0(file.name, ".pdf"), height=6, width=6)
   plot(n3 ~ snr3, ylim=ylim, xlim=xlim, ylab="", xlab=xlab.text, main=main.text[1], yaxt="n", col=col2, pch=15, cex=2, lwd=0, cex.axis=1.5, cex.lab=1.6, cex.main=1.7)
   lm.fit <- lm(n3 ~ snr3)
   abline(lm.fit, col=col, lwd=4)
   
   cor3 <- cor.test(n3, snr3, method="spearman", exact=F)
   #cor3 <- round0(cor3[[4]], digits=2)
   #legend(pos, paste("rho = ", cor3), text.col=col, pch=15, col=col2, pt.cex=2.5, cex=1.5, pt.lwd=0, text.font=1)
   legend("bottomright", c(paste0("rho = ", round0(cor3[[4]], digits=2)), paste0("p-value = ", scientific(cor3[[3]]))), text.col=cols, text.font=2, bty="n", cex=1.5)
   
   axis(side=2, at=seq(-0.2, 0.2, by=0.2), labels=c(-0.2, 0, 0.2), cex.axis=1.5)
   axis(side=2, at=seq(-0.1, 0.1, by=0.2), labels=c(-0.1, 0.1), cex.axis=1.5)
   mtext(ylab.text, side=2, line=2.75, cex=1.6)
   dev.off()
}

## FACS
#facs <- readTable(file.path(wd.ngs, "nbl_cl_n8_FACS.txt"), header=T, rownames=T, sep="")
#samples <- samples[facs$SAMPLE_ID,]
#facs$SUM <- facs$G1 + facs$S + facs$G2
#facs$G1 <- (facs$G1/facs$SUM) * 100
#facs$S <- (facs$S/facs$SUM) * 100
#facs$G2 <- (facs$G2/facs$SUM) * 100

#file.name <- file.path(wd.rt.plots, "FACS_NBL-CL")
#main.text <- c("Flow cytometry validation", "")
#xlab.text <- "Proportion of cells"
#ylab.text <- "Proportion of S phase cells"
#plotFACS(samples$COR, facs$G1, samples$COR, facs$S, file.name, main.text, xlab.text, ylab.text, c("blue", "red"), c("right", "left"))

samples.tmp <- samples
samples <- samples[setdiff(rownames(samples), c("VHIO179-2", "STG139M-2", "STG201-2")),]
samples <- samples.tmp

file.name <- file.path(wd.rt.plots, "G4R_BRCA_3P_n15")
main.text <- c(paste("In silico vs. G4R"), "")
xlab.text <- expression(paste("Number of ", Delta, "G4R [#]"))
ylab.text <- "Overall read depth vs. RT [rho]"                                                                         ## "#619CFF", "#F8766D", "#00BA38"      "skyblue3", "lightcoral", "#59a523"
cols <- "black"
cols2 <- "darkgray"
plotFACS3(samples$COR, samples$G4R, file.name, main.text, xlab.text, ylab.text, cols, cols2, "topright", 16531)





###
## https://stackoverflow.com/questions/7588020/how-to-write-labels-in-barplot-on-x-axis-with-duplicated-names
file.name <- file.path(wd.rt.plots, "FACS_NBL-CL_barchart")
main.text <- c("Flow cytometry validation in vitro (Dean-Jett-Fox)", "")
xlab.text <- ""
ylab.text <- "Proportion of cells [%]"
#blue  <- "blue"   ## adjustcolor("#619CFF", alpha.f=0.9)
#red   <- "red"   ## adjustcolor("#F8766D", alpha.f=0.9)
#green <- "darkgray"   ## adjustcolor("#00BA38", alpha.f=0.9)
#cols <- c(blue, red, green)   ## #59a523 (Alcro wasabi)
cols <- c(blue, red, gray)
facs1 <- t(as.matrix(facs[,-1]))

pdf(paste0(file.name, ".pdf"), height=6, width=9.3)
par(mar=c(5.1, 4, 4.1, 3.7), xpd=TRUE)
barplot(facs1, col=cols, ylim=c(0, 100), main=main.text[1], cex.names=1.5, cex.axis=1.5, cex.lab=1.6, cex.main=1.7)

#mids <- barplot(facs1, xlab="")   ## To capture the midpoints
axis(1, at=1.9, labels="GIMEN", cex.axis=1.5, las=0, lwd.tick=0)

legend("right", rownames(facs1)[3:1], cex=1.6, fill=cols[3:1], horiz=F, bty="n", inset=c(-0.11, 0))
mtext(ylab.text, side=2, line=2.75, cex=1.6)
dev.off()





###
## FACS (14/11/20)
facs0 <- readTable(file.path(wd.ngs, "nbl_cl_n8_FACS_0.txt"), header=T, rownames=T, sep="")
facs0 <- facs0[samples$SAMPLE_ID,]
#ss <- facs$SAMPLE_ID
#facs0 <- facs
#for (s in 1:length(ss)) {
#   facs000 <- subset(facs00, model %in% ss[s])
#   
#   facs0$G1[s] <- sum(facs000$G1)
#   facs0$S[s] <- sum(facs000$S)
#   facs0$G2[s] <- sum(facs000$G2)
#}
#facs0$SUM <- facs0$G1 + facs0$S + facs0$G2
#facs0$G1 <- (facs0$G1/facs0$SUM) * 100
#facs0$S <- (facs0$S/facs0$SUM) * 100
#facs0$G2 <- (facs0$G2/facs0$SUM) * 100

file.name <- file.path(wd.rt.plots, "FACS_NBL-CL_2019_new")
main.text <- c("Flow cytometry validation", "")
xlab.text <- "Proportion of cells [%]"
ylab.text <- "Overall correlation with LCL RT"
plotFACS3(samples$COR, facs0$G1, samples$COR, facs0$S, samples$COR, facs0$G2, file.name, main.text, xlab.text, ylab.text, c("#619CFF", "#F8766D", "#00BA38"), "topright")










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
file.name <- file.path(wd.rt.plots, "SAMPLES-vs-RT_NBL-CL-vs-LCL_spline_spearman")
main.text <- c("NBL-CL (n=8) vs. LCL S/G1", "")
ymin <- -0.8732989
ymax <- 0.8643419
plotSAMPLEvsRTALL(cors.samples, samples1, file.name, main.text, ymin, ymax, size=6)

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
plotPCA(1, 2, pca.de, trait, wd.rt.plots, "PCA_NBL_CL_chrs_spline_spearman", size=6, file.main, "topleft", c("red", "lightpink1", "lightskyblue2", "blue"), samples1, flip.x=-1, flip.y=1, legend.title="Overall corr. with LCL S/G1")

## SG1
#trait <- samples.nbl.sg1$SG1
#plotPCA(1, 2, pca.de, trait, wd.rt.plots, "pca_nbl_CL_chrs_spline_spearman_SG1", size=6, file.main, "topright", c("red", "lightgray", "blue"), samples1, flip.x=-1, flip.y=1, legend.title="Consist. CM in all chrs")







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
