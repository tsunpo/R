# =============================================================================
# Manuscript   : 
# Chapter      : 
# Name         : manuscripts/replication/brca-wgs-rt-g4s.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 02/02/21
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "GQuadruplex.R", "TranscriptionReplicationConflict.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.bed.gc.1kb.RData"))

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 04/07/19
# -----------------------------------------------------------------------------
#wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE  <- "BRCA"
PAIR1 <- "T"
base  <- tolower(BASE)
method <- "rpkm"

wd.ngs   <- file.path(wd, BASE, "ngs/WGS")
wd.anlys <- file.path(wd, BASE, "analysis")
wd.meta  <- file.path(wd, BASE, "metadata", "Hänsel-Hertsch 2020")   #, "DG4Rs")

wd.rt       <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

wd.ngs.data <- file.path(wd.ngs, "data")
samples <- readTable(file.path(wd.ngs, "brca_wgs_n22.txt"), header=T, rownames=T, sep="\t")
n <- nrow(samples)

helicases <- readTable(file.path(wd.meta, "Helicases.list"), header=F, rownames=F, sep="\t")
#helicases <- rownames(getGenes(helicases))
tfbs      <- readTable(file.path(wd.meta, "TFBS.list"), header=F, rownames=F, sep="\t")
#tfbs      <- rownames(getGenes(tfbs))

# -----------------------------------------------------------------------------
# G-quadruplex region skew (G4RS)
# Last Modified: 02/02/21
# -----------------------------------------------------------------------------
g4rs <- getG4RS(wd.meta, samples, bed.gc, nrds.lcl)
save(g4rs, file=file.path(wd.rt.data, paste0("brca_g4rs_n22.RData")))

samples$G4RS <- g4rs$G4RS

g4rs$rho <- 0
for (s in 1:nrow(g4rs)) {
   ## G4RS vs. LCL RTS
   g4s <- as.numeric(g4rs[s, paste0("chr", 1:22)])
   rts <- sprs.lcl$spr
   cor3 <- cor.test(rts, g4s, method="spearman", exact=F)
   g4rs$rho[s] <- cor3[[4]]

   sample <- samples$SAMPLE_ID[s]
   file.name <- file.path(wd.rt.plots, paste0("BRCA-G4RS_vs_LCL-RTS_", sample, "_NEW"))
   #main.text <- c(paste0(sample, " G-quadruplex region skew"), "G4RS = (E-L)/(E+L)")
   main.text <- c(paste0(sample, " G4RS correlation"), "")
   xlab.text <- "G-quadruplex region skew (G4RS)"
   ylab.text <- "RTS = (E-L)/(E+L)"
   cols <- c(red, blue, "black", green)
   plotG4RS(rts, g4s, file.name, main.text, xlab.text, ylab.text, cols, "topright")
}

## G4RS E vs. L
ylim <- c(min(c(g4rs$L, g4rs$E)), max(c(g4rs$L, g4rs$E)))
file.name <- paste0("boxplot_brca_g4rs_E+L_darkgray_NEW_mar=4.6")
main.txt <- c("BRCA G4R", "")
#plotBoxG4R(wd.rt.plots, file.name, g4rs$L, g4rs$E, main.txt, , ylim)
tpm.1 <- g4rs$L
tpm.2 <- g4rs$E
names <- c("Late", "Early")
cols <- c("darkgray", "darkgray")

trait <- rep(0, length(tpm.1))
trait <- c(trait, rep(1, length(tpm.2)))
trait <- as.factor(trait)
expr <- as.numeric(c(tpm.1, tpm.2))

pdf(file.path(wd.rt.plots, paste0(file.name, ".pdf")), height=6, width=4.2)
par(mar=c(5.1, 4.6, 4.1, 1.5))
boxplot(expr ~ trait, outline=T, xaxt="n", ylab="Number of G4R", xlab="", main=main.txt[1], col=cols, ylim=ylim, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)

p <- testU(tpm.1, tpm.2)
#text(1.5, ylim[2]-1250, paste0("p-value = ", scientific(wilcox.test(expr ~ trait, exact=F)$p.value)), col="black", cex=1.5)
#text(1.5, ylim[2], getPvalueSignificanceLevel(p), col="black", cex=3)
text(1.5, ylim[2]-1500, expression(italic('P')~'= 3.24E-08'), col="black", cex=1.8)
text(1.5, ylim[2], "*", col="black", cex=3)

##
axis(side=1, at=seq(1, 2, by=1), labels=names, font=2, cex.axis=1.8)
axis(side=1, at=1, labels=paste0("n=", format(sum(tpm.1), big.mark=",", scientific=F)), line=1.7, col=NA, cex.axis=1.7)
axis(side=1, at=2, labels=paste0("n=", format(sum(tpm.2), big.mark=",", scientific=F)), line=1.7, col=NA, cex.axis=1.7)

#mtext("Number of G4R", side=2, line=2.7, cex=1.8)   
#mtext("", line=0.3, cex=1.25)
dev.off()





## Replication timing skew (RTS)
file.name <- file.path(wd.rt.plots, "RTS_BRCA-M2-M1_spline_spearman_chr2")
main.text <- c("Replication timing skew", "RTS = (E-L)/(E+L)")
ylab.text <- "BRCA M2/M1"
plotRTS(sprs.brca, file.name, main.text, c(4, 13, 17, 19), digits=3, unit=5, ylab.text, cex=1.2, chr2="-0.13", offset="           ")

### Figure 4D
## SCLC vs. LCL
file.name <- file.path(wd.rt.plots, "RTS2_BRCA-M2-M1_vs_LCL_spline_spearman")
main.text <- c("Replication timing skew", "")
xlab.text <- "LCL S/G1"
ylab.text <- "BRCA M2/M1"
plotRTS2(sprs.brca$spr, sprs.lcl$spr, file.name, main.text, c(4, 13, 17, 19, 22), xlab.text, unit=5, ylab.text, cex=1.2)


# -----------------------------------------------------------------------------
# 
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
   ylim <- c(-0.34, 0.25)
 
   pdf(paste0(file.name, ".pdf"), height=6, width=6)

   plot(n3 ~ snr3, ylim=ylim, xlim=xlim, ylab=ylab.text, xlab=xlab.text, main=main.text[1], yaxt="n", col=col2, pch=15, cex=2, lwd=0, cex.axis=1.5, cex.lab=1.6, cex.main=1.7)
   lm.fit <- lm(n3 ~ snr3)
   abline(lm.fit, col=col, lwd=4)
   
   cor3 <- cor.test(n3, snr3, method="spearman", exact=F)
   #cor3 <- round0(cor3[[4]], digits=2)
   #legend(pos, paste("rho = ", cor3), text.col=col, pch=15, col=col2, pt.cex=2.5, cex=1.5, pt.lwd=0, text.font=1)
   legend("bottomright", c(paste0("rho = ", round0(cor3[[4]], digits=2)), paste0("p-value = ", scientific(cor3[[3]]))), text.col=cols, text.font=2, bty="n", cex=1.5)
   
   axis(side=2, at=seq(-0.2, 0.2, by=0.2), labels=c(-0.2, 0, 0.2), cex.axis=1.5)
   axis(side=2, at=seq(-0.3, 0.1, by=0.2), labels=c(-0.3, -0.1, 0.1), cex.axis=1.5)
   #mtext(ylab.text, side=2, line=2.75, cex=1.6)
   mtext(main.text[2], line=0.3, cex=1.6)
   dev.off()
}

plotFACS30 <- function(n3, snr3, file.name, main.text, xlab.text, ylab.text, col, col2, pos, xlim.max) {
   xlim <- c(0, xlim.max)
   #ylim <- c(-0.35, 0.25)
 
   pdf(paste0(file.name, ".pdf"), height=6, width=6)
   plot(n3 ~ snr3, xlim=xlim, ylab="", xlab=xlab.text, main=main.text[1], col=col2, pch=15, cex=2, lwd=0, cex.axis=1.5, cex.lab=1.6, cex.main=1.7)
   lm.fit <- lm(n3 ~ snr3)
   abline(lm.fit, col=col, lwd=4)
 
   cor3 <- cor.test(n3, snr3, method="spearman", exact=F)
   #cor3 <- round0(cor3[[4]], digits=2)
   #legend(pos, paste("rho = ", cor3), text.col=col, pch=15, col=col2, pt.cex=2.5, cex=1.5, pt.lwd=0, text.font=1)
   legend("bottomright", c(paste0("rho = ", round0(cor3[[4]], digits=2)), paste0("p-value = ", scientific(cor3[[3]]))), text.col=cols, text.font=2, bty="n", cex=1.5)
 
   axis(side=2, at=seq(-0.2, 0.2, by=0.2), labels=c(-0.2, 0, 0.2), cex.axis=1.5)
   axis(side=2, at=seq(-0.3, 0.1, by=0.2), labels=c(-0.3, -0.1, 0.1), cex.axis=1.5)
   mtext(ylab.text, side=2, line=2.75, cex=1.6)
   mtext(main.text[2], line=0.3, cex=1.6)
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
samples <- samples[setdiff(rownames(samples), c("STG139M-2", "STG143-2", "STG201-2", "VHIO098-2", "VHIO179-2")),]
samples <- samples.tmp

file.name <- file.path(wd.rt.plots, "G4R_BRCA_3P_n22_mar=4.6")
main.text <- c(paste("In silico vs. G4R"), "n=22")
xlab.text <- expression(paste("Number of ", Delta, "G4R [#]"))
ylab.text <- "Overall read depth vs. RT [rho]"                                                                         ## "#619CFF", "#F8766D", "#00BA38"      "skyblue3", "lightcoral", "#59a523"
cols <- "black"
cols2 <- "darkgray"
plotFACS3(samples$COR, samples$G4R, file.name, main.text, xlab.text, ylab.text, cols, cols2, "topright", 16531)

###
##
file.name <- file.path(wd.rt.plots, "G4R_READ")
main.text <- c(paste("Reads vs. G4R"), "n=27")
xlab.text <- expression(paste("Number of ", Delta, "G4R [#]"))
ylab.text <- "Number of reads [log10]"
cols <- "black"
cols2 <- "darkgray"
plotFACS30(log10(samples1$BAM_C), samples1$G4R, file.name, main.text, xlab.text, ylab.text, cols, cols2, "topright", 16531)
