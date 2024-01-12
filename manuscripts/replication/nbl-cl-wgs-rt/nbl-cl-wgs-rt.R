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
#wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/ty2/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE  <- "NB-CL"
PAIR1 <- "T"
base  <- tolower(BASE)
method <- "rpkm"

wd.ngs   <- file.path(wd, BASE, "ngs/WGS")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.rt       <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

wd.ngs.data <- file.path(wd.ngs, "data")
samples1 <- readTable(file.path(wd.ngs, "nb-cl_wgs_n8.list"), header=F, rownames=F, sep="")
samples1 <- toupper(gsub("-", "", samples1)) 
n1 <- length(samples1)

# -----------------------------------------------------------------------------
# NBL-CL vs LCL S/G1
# http://pklab.med.harvard.edu/scw2014/subpop_tutorial.html
# Last Modified: 15/11/19; 05/06/19; 20/04/19; 06/03/19
# -----------------------------------------------------------------------------
cors.samples <- getSAMPLEvsRT(wd.rt.data, samples1)
colnames(cors.samples) <- toupper(gsub("-", "", colnames(cors.samples)))
save(cors.samples, file=file.path(wd.rt.data, paste0("samples-vs-rt_", base, "-vs-lcl_spline_spearman.RData")))

# -----------------------------------------------------------------------------
# Overall correlation with LCL S/G1
# Last Modified: 20/11/23
# -----------------------------------------------------------------------------
facs <- readTable(file.path(wd.ngs, "nbl_cl_n8_FACS.txt"), header=T, rownames=T, sep="")
for (s in 1:nrow(facs)) {
	  sum <- sum(facs[s, 2:4])
	  facs[s, 2:4] <- facs[s, 2:4] / sum * 100
}
facs <- facs[samples1,]

insilico <- toTable(0, 4, 22, c("CHR", "S", "G2", "G1"))
insilico$CHR <- sprs.order$chr
for (c in 1:nrow(sprs.order)) {
	  cors.all <- toTable(0, 4, length(samples1), c("SAMPLE_ID", "COR", "Q4", "M2"))
	  rownames(cors.all) <- samples1
	  cors.all$SAMPLE_ID <- samples1
	  for (s in 1:length(samples1)) {
		    sample <- samples1[s]
		    load(file.path(wd.rt.data, "samples", paste0("rd-vs-rt_", sample, "-vs-lcl_spline_spearman.RData")))
		
		    cors.all$COR[s] <- cors$cor[c]
	  }
	
	  cor <- cor.test(facs$S, cors.all$COR, method="spearman", exact=F)
	  insilico$S[c] <- cor[[4]]
	
	  cor <- cor.test(facs$G2, cors.all$COR, method="spearman", exact=F)
	  insilico$G2[c] <- cor[[4]]
	
  	cor <- cor.test(facs$G1, cors.all$COR, method="spearman", exact=F)
  	insilico$G1[c] <- cor[[4]]
}

##
ylim <- c(-0.4, 0.9)  #, dns.u$Freq)))
file.name <- file.path(wd.rt.plots, paste0("DNS_NB-CL_COR_-1_7_flowjo_2.5_red_FACS_22_3_3"))
#main.text <- c(expression(bold("NB-CL")~bolditalic('in silico')~bold("SCF estimation")), "")
main.text <- c("NB-CL SCF correlation", "")
ylab.text <- "Correlation to FACS cells"
xlab.text <- "Chromosome"
legends <- c("S", "G0/G1 (Inverted)", "G2/M")
#cols <- c(red, blue, "dimgray")
cols <- c(flowjo.red, flowjo.blue, flowjo.grey)

pdf(paste0(file.name, ".pdf"), height=4.5, width=9.8)
par(mar=c(5.1, 4.6, 4.1, 1.5))
plot(NULL, xlim=c(1, 22), ylim=ylim, xlab=xlab.text, ylab=ylab.text, main=main.text, yaxt="n", xaxt="n", pch=19, cex.axis=1.8, cex.lab=1.9, cex.main=2)

points(insilico$G2 ~ rownames(insilico), col=cols[3], pch=15, cex=2.5)
lines(rownames(insilico), y=insilico$G2, lty=5, lwd=2.5, col=cols[3])

points(insilico$G1*-1 ~ rownames(insilico), col=cols[2], pch=15, cex=2.5)
lines(rownames(insilico), y=insilico$G1*-1, lty=5, lwd=2.5, col=cols[2])

points(insilico$S ~ rownames(insilico), col=cols[1], pch=15, cex=2.5)
lines(rownames(insilico), y=insilico$S, lty=5, lwd=2.5, col=cols[1])

abline(v=7, lty=3, lwd=3, col=red)

axis(side=2, at=seq(-0.4, 0.8, by=0.2), labels=c(-0.4, "", 0, "", 0.4, "", 0.8), cex.axis=1.8)	  
axis(side=2, at=0.4, cex.axis=1.8)	
axis(side=2, at=0.8, cex.axis=1.8)	
axis(side=1, at=seq(1, 22, by=2), labels=insilico$CHR[seq(1, 22, by=2)], cex.axis=1.8)
axis(side=1, at=seq(2, 22, by=2), labels=insilico$CHR[seq(2, 22, by=2)], cex.axis=1.8)
legend("bottomright", legend=legends, col=cols, pch=15, lty=5, lwd=3, pt.cex=3, cex=1.9)
dev.off()


samples <- samples[order(samples$COR),]
#samples <- samples[facs$SAMPLE_ID,]




# -----------------------------------------------------------------------------
# Overall correlation with LCL S/G1
# Last Modified: 15/11/19; 16/06/19; 04/06/19; 06/03/19
# -----------------------------------------------------------------------------
samples.nbl.cl <- setSamplesM2(wd.rt.data, samples1, c=7)
samples.nbl.cl <- samples.nbl.cl[order(samples.nbl.cl$COR),]
samples.nbl.cl$SAMPLE_ID <- toupper(gsub("-", "", samples.nbl.cl$SAMPLE_ID))
rownames(samples.nbl.cl) <- samples.nbl.cl$SAMPLE_ID
writeTable(samples.nbl.cl, file.path(wd.ngs, "nb-cl_wgs_n8.txt"), colnames=T, rownames=F, sep="\t")
#          0%         25%         50%         75%        100% 
# -0.32832776 -0.23125788 -0.02602659  0.39228813  0.44813604

samples <- samples.nbl.cl[order(samples.nbl.cl$COR, decreasing=T),]

# -----------------------------------------------------------------------------
# Beeswarm plots
# Last Modified: 21/04/19
# -----------------------------------------------------------------------------
n.nbl.cl  <- nrow(samples.nbl.cl)

samples <- toTable(0, 5, n.nbl.cl, c("CANCER", "COR", "Q4", "M2", "SAMPLE_ID"))
samples$CANCER[1:n.nbl.cl] <- 1
samples$COR <- samples.nbl.cl$COR
samples$Q4  <- samples.nbl.cl$Q4
samples$M2  <- samples.nbl.cl$M2
samples$SAMPLE_ID <- samples.nbl.cl$SAMPLE_ID
rownames(samples) <- samples$SAMPLE_ID
#adjustcolor.gray <- adjustcolor("black", alpha.f=0.75)

pdf(file.path(wd.rt.plots, "boxplot_nb-cl_width=5_c=7_cex=3.3_5.5.pdf"), height=5.5, width=5)
par(mar=c(5.1, 4.6, 4.1, 1.5))
ymax <- 0.6
ymin <- 0.13
main.text <- c(expression(bolditalic('In silico')~bold("SCF estimation")), "")
boxplot(COR ~ CANCER, data=samples, outline=F, names=c(""), ylim=c(ymin, ymax), ylab="Correlation to RT", xlab="(n=8)", main=main.text, yaxt="n", col="white", boxwex=0.75, cex.axis=1.8, cex.lab=1.9, cex.main=2)
#abline(h=0, lty=5, lwd=2)

#adjustcolor.gray <- adjustcolor("black", alpha.f=0.625)
points(subset(samples, M2 == 2)$CANCER, subset(samples, M2 == 2)$COR, pch=19, col=red, cex=3, lwd=1)
points(subset(samples, M2 == 1)$CANCER, subset(samples, M2 == 1)$COR, pch=19, col=blue, cex=3, lwd=1)
#points(subset(samples, Q4 == 2)$CANCER, subset(samples, Q4 == 2)$COR, col=blue.lighter, pch=19, cex=2)
#points(subset(samples, Q4 == 1)$CANCER, subset(samples, Q4 == 1)$COR, col=blue, pch=19, cex=2)
#points(subset(samples, Q4 == 3)$CANCER, subset(samples, Q4 == 3)$COR, col=red.lighter, pch=19, cex=2)
#points(subset(samples, Q4 == 4)$CANCER, subset(samples, Q4 == 4)$COR, col=red, pch=19, cex=2)

#points(samples$CANCER, samples$COR, col=adjustcolor.gray, pch=19, cex=2.5)

for (s in 1:nrow(samples)) {
   sample <- samples[s,]

   if (sample$SAMPLE_ID == "NGP" || sample$SAMPLE_ID == "LAN6" || sample$SAMPLE_ID == "TR14")
   	  text(sample$CANCER, sample$COR, sample$SAMPLE_ID, col="black", adj=c(1.3, 0.5), cex=1.8)
   else if (sample$SAMPLE_ID == "CLBGA" || sample$SAMPLE_ID == "SKNAS")
   	  text(sample$CANCER, sample$COR, sample$SAMPLE_ID, col="black", adj=c(1.2, 0.5), cex=1.8)
   else if (sample$SAMPLE_ID == "SKNFI")
      text(sample$CANCER, sample$COR, sample$SAMPLE_ID, col="black", adj=c(1.3, -0.45), cex=1.8)
   else if (sample$SAMPLE_ID == "GIMEN")
   	  text(sample$CANCER, sample$COR, sample$SAMPLE_ID, col="black", adj=c(1.2, -0.2), cex=1.8)
   else if (sample$SAMPLE_ID == "LS")
      text(sample$CANCER, sample$COR, sample$SAMPLE_ID, col="black", adj=c(1.6, 0.5), cex=1.8)
   else if (sample$SAMPLE_ID == "")
      text(sample$CANCER, sample$COR, sample$SAMPLE_ID, col="black", adj=c(1.19, 1.1), cex=1.8)
   else
      text(sample$CANCER, sample$COR, sample$SAMPLE_ID, col="black", adj=c(1.15, 1.17), cex=1.8)
}

legend("topright", legend = c("M2", "M1"), pch=19, pt.cex=3.3, col=c(red, blue), cex=1.9)

axis(side=2, at=seq(-0.2, 0.6, by=0.1), labels=c(-0.2, "", 0, "", 0.2, "", 0.4, "", 0.6), cex.axis=1.8)
#mtext("Spearman's rho", side=2, line=2.73, cex=1.8)
#mtext("", cex=1.2, line=0.3)
axis(side=1, at=1, labels="NB-CL", cex.axis=1.9)
#mtext(text=c(), side=1, cex=1.4, line=0.9, at=c(1,2,3))
#mtext(text=c("(n=8)"), side=1, cex=1.9, line=2.8, at=c(1,2,3))
dev.off()

# -----------------------------------------------------------------------------
# 
# Last Modified: 10/10/23
# -----------------------------------------------------------------------------
facs <- readTable(file.path(wd.ngs, "nbl_cl_n8_FACS.txt"), header=T, rownames=T, sep="")
for (s in 1:nrow(facs)) {
	  sum <- sum(facs[s, 2:4])
	  facs[s, 2:4] <- facs[s, 2:4] / sum * 100
}
samples <- samples[order(samples$COR),]
#samples <- samples[facs$SAMPLE_ID,]
facs <- facs[samples$SAMPLE_ID,]

file.name <- file.path(wd.rt.plots, paste0("Cor_SCF-vs-S_c=7_NEW"))
main.text <- "S"
xlab.text <- expression(italic('In vitro')~"FACS")
ylab.text <- expression(italic('In silico')~"SCF")
plotCorrelation(file.name, main.text, xlab.text, ylab.text, x=facs$S, y=samples$COR, pos="bottomright", cols=c(flowjo.red, red), size=4.5, pch=15, cex=3)

file.name <- file.path(wd.rt.plots, paste0("Cor_SCF-vs-G1_c=7_NEW"))
main.text <- "G0/G1"
xlab.text <- expression(italic('In vitro')~"FACS")
ylab.text <- expression(italic('In silico')~"SCF")
plotCorrelation(file.name, main.text, xlab.text, ylab.text, x=facs$G1, y=samples$COR, pos="bottomleft", cols=c(flowjo.blue, blue), size=4.5, pch=15, cex=3)

file.name <- file.path(wd.rt.plots, paste0("Cor_SCF-vs-G2_c=7_NEW"))
main.text <- "G2/M"
xlab.text <- expression(italic('In vitro')~"FACS")
ylab.text <- expression(italic('In silico')~"SCF")
plotCorrelation(file.name, main.text, xlab.text, ylab.text, x=facs$G2, y=samples$COR, pos="bottomright", cols=c(flowjo.grey, "dimgray"), size=4.5, pch=15, cex=3)

# -----------------------------------------------------------------------------
# 
# Last Modified: 10/10/23
# -----------------------------------------------------------------------------
## https://stackoverflow.com/questions/7588020/how-to-write-labels-in-barplot-on-x-axis-with-duplicated-names
file.name <- file.path(wd.rt.plots, "FACS_NBL-CL_RDC_16_G0_G1")
main.text <- c(expression(bolditalic('In vitro')~bold("FACS validation")), "")
xlab.text <- ""
ylab.text <- "Fraction of cells"
#blue  <- "blue"   ## adjustcolor("#619CFF", alpha.f=0.9)
#red   <- "red"   ## adjustcolor("#F8766D", alpha.f=0.9)
#green <- "darkgray"   ## adjustcolor("#00BA38", alpha.f=0.9)
#cols <- c(blue, red, green)   ## #59a523 (Alcro wasabi)
cols <- c(blue, red, "darkgray")
cols2 <- c(flowjo.red, flowjo.grey, flowjo.blue)
facs1 <- t(as.matrix(facs[,-1]))
facs2 <- facs1[c(2,3,1),]

pdf(paste0(file.name, ".pdf"), height=5.5, width=12)
par(mar=c(5.1, 4.6, 4.1, 7.4), xpd=TRUE)
barplot(facs2, col=cols2, ylim=c(0, 100), ylab=ylab.text, xaxt="n", main=main.text[1], cex.names=1.7, cex.axis=1.8, cex.lab=1.9, cex.main=2)
text(labels=facs$SAMPLE_ID, x=c(0.8, 2, 3.2, 4.4, 5.6, 6.8, 8, 9.2), y=par("usr")[3] - 4, srt=45, adj=0.965, xpd=NA, cex=1.9)

#mids <- barplot(facs1, xlab="")   ## To capture the midpoints
#axis(1, at=1.9, labels="GIMEN", cex.axis=1.7, las=0, lwd.tick=0)
#axis(1, at=4.3, labels="SKNFI", cex.axis=1.7, las=0, lwd.tick=0)
#axis(1, at=7.9, labels="NGP", cex.axis=1.7, las=0, lwd.tick=0)

legend("right", c("G0/G1", "G2/M", "S"), text.col="black", pch=c(15, 15, 15), col=cols2[3:1], pt.cex=3, cex=1.9, horiz=F, bty="n", inset=c(-0.16, 0))
#mtext(ylab.text, side=2, line=2.75, cex=1.8)
dev.off()

























# -----------------------------------------------------------------------------
# NBL-CL vs LCL S/G1
# http://pklab.med.harvard.edu/scw2014/subpop_tutorial.html
# Last Modified: 15/11/19; 05/06/19; 20/04/19; 06/03/19
# -----------------------------------------------------------------------------
cors.samples <- getSAMPLEvsRT(wd.rt.data, samples1)
save(cors.samples, file=file.path(wd.rt.data, paste0("samples-vs-rt_", base, "-vs-lcl_spline_spearman.RData")))
# > min(cors.samples[,-c(1:4)])
# [1] -0.2733976
# > max(cors.samples[,-c(1:4)])
# [1] 0.8063618

#load(file.path(wd.rt.data, paste0("samples-vs-rt_nbl-cl-vs-lcl_spline_spearman.RData")))
#file.name <- file.path(wd.rt.plots, "SAMPLES-vs-RT_NBL-CL-vs-LCL_spline_spearman")
#main.text <- c("NBL-CL read depth vs. LCL S/G1", "")
#ymin <- -0.2733976
#ymax <- 0.8063618
#plotSAMPLEvsRT(cors.samples, samples1, file.name, main.text, ymin, ymax)

# -----------------------------------------------------------------------------
# Last Modified: 30/11/20; 04/06/19; 21/04/19
# -----------------------------------------------------------------------------
nrds <- getOverallReadDepth(wd.rt.data, base, method, PAIR1, n1)
# > nrow(nrds)
# [1] 2684771
nrds.m <- subset(nrds, MEDIAN != 0)   ## ADD 30/11/20
# > nrow(nrds.m)
# [1] 2673685
nrds.d <- getOverallDetectedRD(wd.rt.data, base, method, PAIR1, n1, samples1)
nrow(nrds.d)
# > nrow(nrds.d)
# [1] 2663308

# -----------------------------------------------------------------------------
# PCA
# Last Modified: 15/11/19; 04/06/19; 21/04/19
# -----------------------------------------------------------------------------
## Refer to cmd-rt_2a_nrd.gc.cn.d_sample.R (commandline mode)
nrds.T.chr.d.all <- readTable(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d_", "chr1", "_", PAIR1, "_n", n1, ".txt.gz")), header=T, rownames=F, sep="\t")
colnames(nrds.T.chr.d.all) <- toupper(gsub("\\.", "", colnames(nrds.T.chr.d.all)))
for (c in 2:22) {
	chr <- chrs[c]
	
	## Read depth
	#nrds.T.chr.d <- pipeGetDetectedRD(wd.ngs.data, BASE, chr, PAIR1, method)
	nrds.T.chr.d <- readTable(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d_", chr, "_", PAIR1, "_n", n1, ".txt.gz")), header=T, rownames=F, sep="\t")
	colnames(nrds.T.chr.d) <- toupper(gsub("\\.", "", colnames(nrds.T.chr.d)))
	
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
   par(mar=c(5.1, 4.6, 4.1, 1.5))
   plot(n1 ~ snr1, ylim=ylim, xlim=xlim, ylab=ylab.text, xlab=xlab.text, main=main.text[1], col=col[1], pch=19, cex=2.5, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
   points(n2 ~ snr2, col=col[2], pch=19, cex=2.5)

   lm.fit <- lm(n1 ~ snr1)
   abline(lm.fit, col=col[1], lwd=3)
   lm.fit <- lm(n2 ~ snr2)
   abline(lm.fit, col=col[2], lwd=3)
 
   #samples <- c("SCLC-NL", "SCLC", "NBL", "CLL")
   #text(snr1, n1, samples, col=col[1], pos=c(1,2,2,2), cex=1.8)
   #text(snr2, n2, samples, col=col[2], pos=c(1,4,4,4), cex=1.8)
 
   cor <- cor.test(n1, snr1, method="spearman", exact=F)
   #legend(pos[1], c(paste0("rho = ", round0(cor[[4]], digits=1)), paste0("p-value = ", scientific(cor[[3]], digits=1))), text.col=col[1], bty="n", cex=1.75)
   legend(pos[1], paste0("G1 (rho = ", round0(cor[[4]], digits=2), ")      "), text.col=col[1], pch=c(NA), col=col[1], bty="n", cex=1.5)
 
   cor <- cor.test(n2, snr2, method="spearman", exact=F)
   #legend(pos[2], c(paste0("rho = ", round0(cor[[4]], digits=1)), paste0("p-value = ", scientific(cor[[3]], digits=1))), text.col=col[2], bty="n", cex=1.75)
   legend(pos[2], paste0("S (rho = ", round0(cor[[4]], digits=2), ")"), text.col=col[2], pch=c(NA), col=col[2], bty="n", cex=1.5)
 
   #mtext(ylab.text, side=2, line=2.75, cex=1.7)
   #mtext(main.text[2], cex=1.2, line=0.3)
   dev.off()
}

# https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
plotFACS3 <- function(n1, snr1, n2, snr2, n3, snr3, file.name, main.text, xlab.text, ylab.text, col, col2, pos, xlim.max) {
   xlim <- c(0, xlim.max)
   ylim <- c(-0.367, 0.5)
 
   pdf(paste0(file.name, ".pdf"), height=6, width=6)
   par(mar=c(5.1, 4.6, 4.1, 1.5))
   plot(n3 ~ snr3, ylim=ylim, xlim=xlim, ylab=ylab.text, xlab=xlab.text, main=main.text[1], yaxt="n", col=col2[3], pch=15, cex=2.5, lwd=0, cex.axis=1.8, cex.lab=1.9, cex.main=2)
   lm.fit <- lm(n3 ~ snr3)
   abline(lm.fit, col=col[3], lwd=5)
   
   #legend(pos, c(expression(paste("S   (", rho, " = 0.8)")), expression(paste("G2 (", rho, " = 0.4)")), expression(paste("G1 (", rho, " = -0.8)"))), text.col=c(col[2], col[3], col[1]), pch=15, col="white", pt.cex=3, cex=1.7, pt.lwd=0, text.font=2, bty="n")
   #legend(pos, c(expression(paste(bold("S"), "   (", rho, " = 0.8)")), expression(paste(bold("G2"), " (", rho, " = 0.4)")), expression(paste(bold("G1"), " (", rho, " = -0.8)"))), text.col=c(col[2], col[3], col[1]), pch=15, col="white", pt.cex=3, cex=1.7, pt.lwd=0, text.font=2, bty="n")
   #legend(pos, c(expression(bold("S")), expression(bold("G2")), expression(bold("G1"))), text.col=c(col[2], col[3], col[1]), pch=c(15, 15, 15), col=c(col2[2], col2[3], col2[1]), pt.cex=3, cex=1.7, pt.lwd=0, text.font=c(2,2,2))
   legend(pos, c("S", "G2", "G1"), pch=c(15, 15, 15), col=c(col2[2], col2[3], col2[1]), pt.cex=3, cex=1.9, pt.lwd=0)
   #cor <- cor.test(n1, snr1, method="spearman", exact=F)
   #legend(pos[1], paste0("G1 (rho = ", round0(cor[[4]], digits=1), ")     "), text.col=col[1], pch=c(NA), col=col[1], bty="n", cex=1.5)
   #
   #cor <- cor.test(n2, snr2, method="spearman", exact=F)
   #legend(pos[2], paste0("S (rho = ", round0(cor[[4]], digits=1), ")"), text.col=col[2], pch=c(NA), col=col[2], bty="n", cex=1.5)
      
   points(n1 ~ snr1, col=col2[1], pch=15, cex=2.5, lwd=0)
   points(n2 ~ snr2, col=col2[2], pch=15, cex=2.5, lwd=0)
   lm.fit <- lm(n1 ~ snr1)
   abline(lm.fit, col=col[1], lwd=5)
   lm.fit <- lm(n2 ~ snr2)
   abline(lm.fit, col=col[2], lwd=5)

   par(xpd=T)
   cor3 <- cor.test(n3, snr3, method="spearman", exact=F)
   cor3 <- round0(cor3[[4]], digits=2)
   text(mean(snr3)-3, mean(n3)-0.025, label=paste("rho = ", cor3), col=col[3], cex=1.8, font=2)
   
   cor2 <- cor.test(n2, snr2, method="spearman", exact=F)
   cor2 <- round0(cor2[[4]], digits=2)
   text(mean(snr2), mean(n2)+0.05, label=paste("rho = ", cor2), col=col[2], cex=1.8, font=2)
   
   cor1 <- cor.test(n1, snr1, method="spearman", exact=F)
   cor1 <- round0(cor1[[4]], digits=2)
   text(mean(snr1), mean(n1)-0.05, label=paste("rho = ", cor1), col=col[1], cex=1.8, font=2)
   par(xpd=F)
   
   axis(side=2, at=seq(-0.4, 0.4, by=0.2), labels=c(-0.4, -0.2, 0, 0.2, 0.4), cex.axis=1.8)
   #mtext(ylab.text, side=2, line=2.75, cex=1.8)
   dev.off()
}

## FACS
## https://link.springer.com/protocol/10.1385/1-59259-227-9:129
facs <- readTable(file.path(wd.ngs, "nbl_cl_n8_FACS_1.8-1.9-2.txt"), header=T, rownames=T, sep="")
samples <- samples[facs$SAMPLE_ID,]
#facs$SUM <- facs$G1 + facs$S + facs$G2
#facs$G1 <- (facs$G1/facs$SUM) * 100
#facs$S <- (facs$S/facs$SUM) * 100
#facs$G2 <- (facs$G2/facs$SUM) * 100

#file.name <- file.path(wd.rt.plots, "FACS_NBL-CL_italic_Spearman's_count")
#main.text <- c(expression(bolditalic('In silico')~bold(" vs. ")~bolditalic('in vitro')), "")
#xlab.text <- "Cell count [%]"
#ylab.text <- "Spearmans's rho"
#plotFACS(samples$COR, facs$G1, samples$COR, facs$S, file.name, main.text, xlab.text, ylab.text, c("blue", "red"), c("right", "left"))

file.name <- file.path(wd.rt.plots, "FACS_NBL-CL_cex=1.8")
main.text <- c(expression(bolditalic('In silico')~bold("vs.")~bolditalic('in vitro')), "")
xlab.text <- "Fraction of cells"
ylab.text <- "Spearman's rho"                                                                         ## "#619CFF", "#F8766D", "#00BA38"      "skyblue3", "lightcoral", "#59a523"
cols <- c(blue, red, "darkgray")
flowjo.blue <- "#989aff"
flowjo.red  <- "#ff9899"
flowjo.grey <- "#b7b7b7"
#cols2 <- c(adjustcolor(blue, alpha.f=0.6), adjustcolor(red, alpha.f=0.6), adjustcolor(dimgray, alpha.f=0.6))
cols2 <- c(flowjo.blue, flowjo.red, flowjo.grey)
plotFACS3(samples$COR, facs$G1, samples$COR, facs$S, samples$COR, facs$G2, file.name, main.text, xlab.text, ylab.text, cols, cols2, "topright", 80)

















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
