# =============================================================================
# Library      : Transcription-Replication Conflict
# Name         : handbook-of/TranscriptionReplicationConflict.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 27/05/20
# =============================================================================

# -----------------------------------------------------------------------------
# Get NRFD information for TSS, TTS and Gene
# Last Modified: 27/05/20
# -----------------------------------------------------------------------------
getGeneNRFD <- function(tpm.gene.log2.m.rfd) {
   length <- abs(tpm.gene.log2.m.rfd$end_position - tpm.gene.log2.m.rfd$start_position) / 1000
   nrfd <- 0
 
   if (tpm.gene.log2.m.rfd$TSS_RFD < 0) {
      if (tpm.gene.log2.m.rfd$strand < 0) {   ## e.g. PIF1 (CD)
         nrfd <- tpm.gene.log2.m.rfd$TSS_RFD - tpm.gene.log2.m.rfd$TTS_RFD 
      } else {                         ## e.g. TOR1AIP1 (HO)
         nrfd <- tpm.gene.log2.m.rfd$TTS_RFD - tpm.gene.log2.m.rfd$TSS_RFD 
      }
   } else {
      if (tpm.gene.log2.m.rfd$strand > 0) {   ## e.g. BRCA2 (CD)
         nrfd <- tpm.gene.log2.m.rfd$TTS_RFD - tpm.gene.log2.m.rfd$TSS_RFD 
      } else {                         ## e.g. ? (HO)
         nrfd <- tpm.gene.log2.m.rfd$TSS_RFD - tpm.gene.log2.m.rfd$TTS_RFD 
      }
   }
 
   return(nrfd/length)
}

getTRC <- function(tpm.gene.log2.m, nrds.RT.NRFD) {
   ## Ensembl gene annotations
   annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
   tpm.gene.log2.m.rfd <- cbind(annot[rownames(tpm.gene.log2.m),], tpm.gene.log2.m[, "MEDIAN"], ensGene.bed[rownames(tpm.gene.log2.m),]) 
   colnames(tpm.gene.log2.m.rfd)[8] <- "MEDIAN"
 
   tpm.gene.log2.m.rfd$TSS_RFD <- nrds.RT.NRFD[tpm.gene.log2.m.rfd$TSS,]$RFD
   tpm.gene.log2.m.rfd$TTS_RFD <- nrds.RT.NRFD[tpm.gene.log2.m.rfd$TTS,]$RFD
   tpm.gene.log2.m.rfd$TSS_NRFD <- nrds.RT.NRFD[tpm.gene.log2.m.rfd$TSS,]$NRFD
   tpm.gene.log2.m.rfd$TTS_NRFD <- nrds.RT.NRFD[tpm.gene.log2.m.rfd$TTS,]$NRFD
   #tpm.gene.log2.m.rfd$TSS_SPLINE <- nrds.RT.NRFD[tpm.gene.log2.m.rfd$TSS,]$SPLINE
   #tpm.gene.log2.m.rfd$TTS_SPLINE <- nrds.RT.NRFD[tpm.gene.log2.m.rfd$TTS,]$SPLINE
   #tpm.gene.log2.m.rfd$RT          <- nrds.RT.NRFD[tpm.gene.log2.m.rfd$TSS,]$RT
   tpm.gene.log2.m.rfd$RT <- nrds.RT.NRFD[tpm.gene.log2.m.rfd$TSS,]$SPLINE
   
   tpm.gene.log2.m.rfd <- tpm.gene.log2.m.rfd[!is.na(tpm.gene.log2.m.rfd$TSS_RFD),]
   tpm.gene.log2.m.rfd <- tpm.gene.log2.m.rfd[!is.na(tpm.gene.log2.m.rfd$TTS_RFD),]
   tpm.gene.log2.m.rfd <- tpm.gene.log2.m.rfd[!is.na(tpm.gene.log2.m.rfd$TSS_NRFD),]
   tpm.gene.log2.m.rfd <- tpm.gene.log2.m.rfd[!is.na(tpm.gene.log2.m.rfd$TTS_NRFD),]
   
   ##
   tpm.gene.log2.m.rfd$GENE_NRFD <- mapply(x = 1:nrow(tpm.gene.log2.m.rfd), function(x) getGeneNRFD(tpm.gene.log2.m.rfd[x,]))
   tpm.gene.log2.m.rfd$TRC <- tpm.gene.log2.m.rfd$strand * tpm.gene.log2.m.rfd$TSS_RFD
   tpm.gene.log2.m.rfd$TRC[which(tpm.gene.log2.m.rfd$TRC > 0)] <- 1
   tpm.gene.log2.m.rfd$TRC[which(tpm.gene.log2.m.rfd$TRC < 0)] <- -1
   tpm.gene.log2.m.rfd$TRC[which(tpm.gene.log2.m.rfd$TRC == 0)] <- 0   ## ADD 27/05/20
   
   return(tpm.gene.log2.m.rfd)
}

setTRC <- function(tpm.gene.log2.m.rfd, rfd, file.name) {
   ## CTR 
   tpm.gene.log2.m.rfd.ctr <- subset(subset(tpm.gene.log2.m.rfd, TSS_RFD < rfd), TSS_RFD > -rfd)
   
   tpm.gene.log2.m.rfd.ctr.iz <- subset(tpm.gene.log2.m.rfd.ctr, TSS_NRFD > 0)
   tpm.gene.log2.m.rfd.ctr.tz <- subset(tpm.gene.log2.m.rfd.ctr, TSS_NRFD < 0)

   tpm.gene.log2.m.rfd.ctr.iz.cd <- subset(tpm.gene.log2.m.rfd.ctr.iz, TRC > 0)
   tpm.gene.log2.m.rfd.ctr.iz.ho <- subset(tpm.gene.log2.m.rfd.ctr.iz, TRC < 0)
   tpm.gene.log2.m.rfd.ctr.tz.cd <- subset(tpm.gene.log2.m.rfd.ctr.tz, TRC > 0)
   tpm.gene.log2.m.rfd.ctr.tz.ho <- subset(tpm.gene.log2.m.rfd.ctr.tz, TRC < 0)
 
   ##
   tpm.gene.log2.m.rfd.ctr.iz.e <- subset(tpm.gene.log2.m.rfd.ctr.iz, RT > 0)
   tpm.gene.log2.m.rfd.ctr.iz.l <- subset(tpm.gene.log2.m.rfd.ctr.iz, RT < 0)
   tpm.gene.log2.m.rfd.ctr.tz.e <- subset(tpm.gene.log2.m.rfd.ctr.tz, RT > 0)
   tpm.gene.log2.m.rfd.ctr.tz.l <- subset(tpm.gene.log2.m.rfd.ctr.tz, RT < 0)
   
   tpm.gene.log2.m.rfd.ctr.iz.e.cd <- subset(tpm.gene.log2.m.rfd.ctr.iz.e, TRC > 0)
   tpm.gene.log2.m.rfd.ctr.iz.e.ho <- subset(tpm.gene.log2.m.rfd.ctr.iz.e, TRC < 0)
   tpm.gene.log2.m.rfd.ctr.iz.l.cd <- subset(tpm.gene.log2.m.rfd.ctr.iz.l, TRC > 0)
   tpm.gene.log2.m.rfd.ctr.iz.l.ho <- subset(tpm.gene.log2.m.rfd.ctr.iz.l, TRC < 0)
   tpm.gene.log2.m.rfd.ctr.tz.e.cd <- subset(tpm.gene.log2.m.rfd.ctr.tz.e, TRC > 0)
   tpm.gene.log2.m.rfd.ctr.tz.e.ho <- subset(tpm.gene.log2.m.rfd.ctr.tz.e, TRC < 0)
   tpm.gene.log2.m.rfd.ctr.tz.l.cd <- subset(tpm.gene.log2.m.rfd.ctr.tz.l, TRC > 0)
   tpm.gene.log2.m.rfd.ctr.tz.l.ho <- subset(tpm.gene.log2.m.rfd.ctr.tz.l, TRC < 0)
   
   ## TTR
   tpm.gene.log2.m.rfd.ttr <- rbind(subset(tpm.gene.log2.m.rfd, TSS_RFD >= rfd), subset(tpm.gene.log2.m.rfd, TSS_RFD <= -rfd))
   
   tpm.gene.log2.m.rfd.ttr.e <- subset(tpm.gene.log2.m.rfd.ttr, RT > 0)
   tpm.gene.log2.m.rfd.ttr.l <- subset(tpm.gene.log2.m.rfd.ttr, RT < 0)

   tpm.gene.log2.m.rfd.ttr.e.cd <- subset(tpm.gene.log2.m.rfd.ttr.e, TRC > 0)
   tpm.gene.log2.m.rfd.ttr.e.ho <- subset(tpm.gene.log2.m.rfd.ttr.e, TRC < 0)
   tpm.gene.log2.m.rfd.ttr.l.cd <- subset(tpm.gene.log2.m.rfd.ttr.l, TRC > 0)
   tpm.gene.log2.m.rfd.ttr.l.ho <- subset(tpm.gene.log2.m.rfd.ttr.l, TRC < 0)
   
   save(file=file.name, tpm.gene.log2.m.rfd, tpm.gene.log2.m.rfd.ttr, tpm.gene.log2.m.rfd.ttr.e, tpm.gene.log2.m.rfd.ttr.l, tpm.gene.log2.m.rfd.ttr.e.cd, tpm.gene.log2.m.rfd.ttr.e.ho, tpm.gene.log2.m.rfd.ttr.l.cd, tpm.gene.log2.m.rfd.ttr.l.ho, tpm.gene.log2.m.rfd.ctr, tpm.gene.log2.m.rfd.ctr.iz, tpm.gene.log2.m.rfd.ctr.iz.cd, tpm.gene.log2.m.rfd.ctr.iz.ho, tpm.gene.log2.m.rfd.ctr.tz, tpm.gene.log2.m.rfd.ctr.tz.cd, tpm.gene.log2.m.rfd.ctr.tz.ho, tpm.gene.log2.m.rfd.ctr.iz.e, tpm.gene.log2.m.rfd.ctr.iz.l, tpm.gene.log2.m.rfd.ctr.tz.e, tpm.gene.log2.m.rfd.ctr.tz.l, tpm.gene.log2.m.rfd.ctr.iz.e.cd, tpm.gene.log2.m.rfd.ctr.iz.e.ho, tpm.gene.log2.m.rfd.ctr.iz.l.cd, tpm.gene.log2.m.rfd.ctr.iz.l.ho, tpm.gene.log2.m.rfd.ctr.tz.e.cd, tpm.gene.log2.m.rfd.ctr.tz.e.ho, tpm.gene.log2.m.rfd.ctr.tz.l.cd, tpm.gene.log2.m.rfd.ctr.tz.l.ho)
}

plotBox0 <- function(wd.de.plots, file.name, tpm.1, tpm.2, main, names, cols, ylim) {
   plotBox(wd.de.plots, file.name, tpm.1, tpm.2, main, names, cols, ylim, height=5, width=3.2)
}

plotBox <- function(wd.de.plots, file.name, tpm.1, tpm.2, main, names, cols, ylim, height=6, width=3) {
   trait <- rep(0, nrow(tpm.1))
   trait <- c(trait, rep(1, nrow(tpm.2)))
   trait <- as.factor(trait)
 
   expr <- as.numeric(c(tpm.1$MEDIAN, tpm.2$MEDIAN))
 
   pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=height, width=width)
   boxplot(expr ~ trait, outline=T, xaxt="n", ylab=paste0("log2(TPM + 0.01)"), main=main, boxcol=cols, whiskcol=cols, outcol=cols, medcol=cols, staplecol=cols, ylim=ylim, cex=1.3, cex.axis=1.1, cex.lab=1.2, cex.main=1.25)
   
   p <- testU(tpm.1$MEDIAN, tpm.2$MEDIAN)
   text <- ""
   if (p < 1E-3)  text <- "*"
   if (p < 1E-9)  text <- "**"
   if (p < 1E-15) text <- "***"
   text(1.5, ylim[2], text, col="black", cex=2.5)
 
   ##
   trait <- rep(0, nrow(tpm.1))
   trait <- c(trait, rep(1, nrow(tpm.2)))
   trait <- as.factor(trait)
   expr <- as.numeric(c(tpm.1$MEDIAN, tpm.2$MEDIAN))
 
   axis(side=1, at=seq(1, 2, by=1), labels=names, font=2, cex.axis=1.2)
   #axis(side=1, at=2, labels="Total n=30,978", line=1.3, col=NA, cex.axis=1.2)
   axis(side=1, at=1, labels=paste0("n=", format(nrow(tpm.1), big.mark=",", scientific=F)), line=1.3, col=NA, cex.axis=1.2)
   axis(side=1, at=2, labels=paste0("n=", format(nrow(tpm.2), big.mark=",", scientific=F)), line=1.3, col=NA, cex.axis=1.2)
   
   mtext(paste0("p-value = ", scientific(wilcox.test(expr ~ trait, exact=F)$p.value)), cex=1.25, line=0.3)
   dev.off()
}

plotBox2 <- function(wd.de.plots, file.name, tpm.1, tpm.2, main, names) {
   trait <- rep(0, length(tpm.1))
   trait <- c(trait, rep(1, length(tpm.2)))
   trait <- as.factor(trait)
 
   expr <- as.numeric(c(tpm.1, tpm.2))
   ylim <- c(min(expr), max(expr))
 
   pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=6, width=3)
   boxplot(expr ~ trait, outline=T, names=names, ylab=paste0("log2(TPM + 0.01)"), main=main, xaxt="n", ylim=ylim, cex.axis=1.1, cex.lab=1.2, cex.main=1.25)
 
   p <- wilcox.test(expr ~ trait, exact=F)$p.value
   text <- ""
   if (p < 1E-3)  text <- "*"
   if (p < 1E-9)  text <- "**"
   if (p < 1E-15) text <- "***"
   text(1.5, ylim[2], text, col="black", cex=2.5)
   
   axis(side=1, at=seq(1, 2, by=1), labels=names, font=2, cex.axis=1.2)
   #axis(side=2, at=seq(5, 8, by=1), cex.axis=1.2)
   mtext(paste0("p-value = ", scientific(p)), cex=1.25, line=0.3)
   dev.off()
}

plotBox20 <- function(wd.de.plots, file.name, tpm.1, tpm.2, main, names) {
   trait <- rep(0, length(tpm.1))
   trait <- c(trait, rep(1, length(tpm.2)))
   trait <- as.factor(trait)
 
   expr <- as.numeric(c(tpm.1, tpm.2))
   ylim <- c(min(expr), max(expr))
 
   pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=6, width=3)
   boxplot(expr ~ trait, outline=T, names=names, ylab=paste0("log2(TPM + 0.01)"), main=main, xaxt="n", yaxt="n", ylim=ylim, cex.axis=1.1, cex.lab=1.2, cex.main=1.25)
 
   axis(side=1, at=seq(1, 2, by=1), labels=names, font=2, cex.axis=1.2)
   axis(side=2, at=seq(5, 8, by=1), labels=c(5, 6, 7, 8), cex.axis=1.1)
   mtext(paste0("p-value = ", scientific(wilcox.test(expr ~ trait, exact=F)$p.value)), cex=1.25, line=0.3)
   dev.off()
}

plotBox3 <- function(wd.de.plots, file.name, tpm.1, tpm.2, tpm.3, main, names, cols, ylab.txt) {
   trait <- rep(0, length(tpm.1))
   trait <- c(trait, rep(1, length(tpm.2)))
   trait <- c(trait, rep(2, length(tpm.3)))
   trait <- as.factor(trait)
 
   expr <- as.numeric(c(tpm.1, tpm.2, tpm.3))
 
   pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=6, width=3.9)
   boxplot(expr ~ trait, outline=T, xaxt="n", ylab=ylab.txt, main=main, boxcol=cols, whiskcol=cols, outcol=cols, medcol=cols, staplecol=cols, cex.axis=1.1, cex.lab=1.2, cex.main=1.25)
   #text(2, 15, "***", col="black", cex=2.5)

   axis(side=1, at=seq(1, 3, by=1), labels=names, font=2, cex.axis=1.2)
   axis(side=1, at=1, labels=paste0("n=", format(length(tpm.1), big.mark=",", scientific=F)), line=1.3, col=NA, cex.axis=1.2)
   axis(side=1, at=2, labels=paste0("n=", format(length(tpm.2), big.mark=",", scientific=F)), line=1.3, col=NA, cex.axis=1.2)
   axis(side=1, at=3, labels=paste0("n=", format(length(tpm.3), big.mark=",", scientific=F)), line=1.3, col=NA, cex.axis=1.2)
 
   #mtext(paste0("p-value = ", scientific(wilcox.test(expr ~ trait, exact=F)$p.value)), cex=1.25, line=0.3)
   dev.off()
}

# -----------------------------------------------------------------------------
# Box plots of TSS_NRFD and GENE_NRFD
# Last Modified: 15/08/20; 31/08/20
# -----------------------------------------------------------------------------
plotBoxTSSNRFD <- function(wd.de.plots, file.name, tpm.1, tpm.2, main.txt, names, ylab.txt, cols, ylim, height=6, width=3, isFlip=F, outline=F) {
   trait <- rep(0, nrow(tpm.1))
   trait <- c(trait, rep(1, nrow(tpm.2)))
   trait <- as.factor(trait)
 
   expr <- as.numeric(c(tpm.1$TSS_NRFD, tpm.2$TSS_NRFD))
   if (isFlip)
      expr <- as.numeric(c(tpm.1$TSS_NRFD * -1, tpm.2$TSS_NRFD * -1))
 
   pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=height, width=width)
   if (outline)
      boxplot(expr ~ trait, outline=outline, xaxt="n", yaxt="n", ylab=ylab.txt, main=main.txt, col=cols, cex.axis=1.1, cex.lab=1.2, cex.main=1.25)
   else
      boxplot(expr ~ trait, outline=outline, xaxt="n", yaxt="n", ylab=ylab.txt, main=main.txt, col=cols, ylim=ylim, cex.axis=1.1, cex.lab=1.2, cex.main=1.25)
 
   p <- testU(tpm.1$TSS_NRFD, tpm.2$TSS_NRFD)
   text <- ""
   if (p < 1E-3)  text <- "*"
   if (p < 1E-9)  text <- "**"
   if (p < 1E-15) text <- "***"
   text(1.5, ylim[2], text, col="black", cex=2.5)
 
   ##
   trait <- rep(0, nrow(tpm.1))
   trait <- c(trait, rep(1, nrow(tpm.2)))
   trait <- as.factor(trait)
   expr <- as.numeric(c(tpm.1$TSS_NRFD, tpm.2$TSS_NRFD))
 
   if (isFlip)
      axis(side=2, at=seq(-0.01, 0.03, by=0.01), labels=c(0.01, 0, -0.01, -0.02, -0.03), cex.axis=1.1)
   else
      axis(side=2, at=seq(-0.01, 0.03, by=0.01), labels=c(-0.01, 0, 0.01, 0.02, 0.03), cex.axis=1.1)
 
   axis(side=1, at=seq(1, 2, by=1), labels=names, font=2, cex.axis=1.2)
   #axis(side=1, at=2, labels="Total n=30,978", line=1.3, col=NA, cex.axis=1.2)
   axis(side=1, at=1, labels=paste0("n=", format(nrow(tpm.1), big.mark=",", scientific=F)), line=1.3, col=NA, cex.axis=1.2)
   axis(side=1, at=2, labels=paste0("n=", format(nrow(tpm.2), big.mark=",", scientific=F)), line=1.3, col=NA, cex.axis=1.2)
 
   mtext(paste0("p-value = ", scientific(wilcox.test(expr ~ trait, exact=F)$p.value)), cex=1.25, line=0.3)
   dev.off()
}

plotBoxGENENRFD <- function(wd.de.plots, file.name, tpm.1, tpm.2, main.txt, names, ylab.txt, cols, ylim, height=6, width=3, isFlip=F, outline=F) {
   trait <- rep(0, nrow(tpm.1))
   trait <- c(trait, rep(1, nrow(tpm.2)))
   trait <- as.factor(trait)
 
   expr <- as.numeric(c(tpm.1$GENE_NRFD, tpm.2$GENE_NRFD))
   if (isFlip)
      expr <- as.numeric(c(tpm.1$GENE_NRFD * -1, tpm.2$GENE_NRFD * -1))
   
   pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=height, width=width)
   if (outline)
      boxplot(expr ~ trait, outline=outline, xaxt="n", yaxt="n", ylab=ylab.txt, main=main.txt, col=cols, cex.axis=1.1, cex.lab=1.2, cex.main=1.25)
   else
      boxplot(expr ~ trait, outline=outline, xaxt="n", yaxt="n", ylab=ylab.txt, main=main.txt, col=cols, ylim=ylim, cex.axis=1.1, cex.lab=1.2, cex.main=1.25)
   
   p <- testU(tpm.1$GENE_NRFD, tpm.2$GENE_NRFD)
   text <- ""
   if (p < 1E-3)  text <- "*"
   if (p < 1E-9)  text <- "**"
   if (p < 1E-15) text <- "***"
   text(1.5, ylim[2], text, col="black", cex=2.5)
 
   ##
   trait <- rep(0, nrow(tpm.1))
   trait <- c(trait, rep(1, nrow(tpm.2)))
   trait <- as.factor(trait)
   expr <- as.numeric(c(tpm.1$GENE_NRFD, tpm.2$GENE_NRFD))
 
   if (isFlip)
      axis(side=2, at=seq(-0.01, 0.03, by=0.01), labels=c(0.01, 0, -0.01, -0.02, -0.03), cex.axis=1.1)
   else
      axis(side=2, at=seq(-0.01, 0.03, by=0.01), labels=c(-0.01, 0, 0.01, 0.02, 0.03), cex.axis=1.1)
   
   axis(side=1, at=seq(1, 2, by=1), labels=names, font=2, cex.axis=1.2)
   #axis(side=1, at=2, labels="Total n=30,978", line=1.3, col=NA, cex.axis=1.2)
   axis(side=1, at=1, labels=paste0("n=", format(nrow(tpm.1), big.mark=",", scientific=F)), line=1.3, col=NA, cex.axis=1.2)
   axis(side=1, at=2, labels=paste0("n=", format(nrow(tpm.2), big.mark=",", scientific=F)), line=1.3, col=NA, cex.axis=1.2)
 
   mtext(paste0("p-value = ", scientific(wilcox.test(expr ~ trait, exact=F)$p.value)), cex=1.25, line=0.3)
   dev.off()
}

plotBoxNRFD <- function(base, BASE, ylim, ylim2=NA, tpm.gene.log2.m.rfd.ctr.iz.e.cd, tpm.gene.log2.m.rfd.ctr.iz.e.ho, tpm.gene.log2.m.rfd.ctr.tz.e.cd, tpm.gene.log2.m.rfd.ctr.tz.e.ho) {
   names <- c("HO", "CD")
   if (is.na(ylim2)[1]) ylim2 <- ylim
   
   ## IZ
   main.txt <- paste0(BASE, " early IZ genes")
   cols=c("red", "red")

   file.name <- paste0("boxplot0_", base, "_tpm.gene_median0_TSS-NRFD_IZ_E_HO+CD")
   plotBoxTSSNRFD(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.e.ho, tpm.gene.log2.m.rfd.ctr.iz.e.cd, main=main.txt, names, "IZ efficiency [TSS]", cols, ylim, height=5, width=3.2)
   file.name <- paste0("boxplot2_", base, "_tpm.gene_median0_TSS-NRFD_IZ_E_HO+CD")
   plotBoxTSSNRFD(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.e.ho, tpm.gene.log2.m.rfd.ctr.iz.e.cd, main=main.txt, names, "IZ efficiency [TSS]", cols, ylim2, height=5, width=3.2, outline=T)
   file.name <- paste0("boxplot1_", base, "_tpm.gene_median0_TSS-NRFD_IZ_E_HO+CD")
   plotBoxTSSNRFD(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.e.ho, tpm.gene.log2.m.rfd.ctr.iz.e.cd, main="TSS", names, "IZ efficiency", cols, ylim, height=5, width=3.2)

   file.name <- paste0("boxplot_", base, "_tpm.gene_median0_GENE-NRFD_IZ_E_HO+CD")
   plotBoxGENENRFD(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.e.ho, tpm.gene.log2.m.rfd.ctr.iz.e.cd, main=main.txt, names, "IZ efficiency [Gene body]", cols, ylim)
   file.name <- paste0("boxplot0_", base, "_tpm.gene_median0_GENE-NRFD_IZ_E_HO+CD")
   plotBoxGENENRFD(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.e.ho, tpm.gene.log2.m.rfd.ctr.iz.e.cd, main=main.txt, names, "IZ efficiency [Gene body]", cols, ylim, height=5, width=3.2)
   file.name <- paste0("boxplot2_", base, "_tpm.gene_median0_GENE-NRFD_IZ_E_HO+CD")
   plotBoxGENENRFD(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.e.ho, tpm.gene.log2.m.rfd.ctr.iz.e.cd, main=main.txt, names, "IZ efficiency [Gene body]", cols, ylim2, height=5, width=3.2, outline=T)
   file.name <- paste0("boxplot1_", base, "_tpm.gene_median0_GENE-NRFD_IZ_E_HO+CD")
   plotBoxGENENRFD(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.e.ho, tpm.gene.log2.m.rfd.ctr.iz.e.cd, main="Gene body", names, "IZ efficiency", cols, ylim, height=5, width=3.2)

   #file.name <- file.path(wd.de.plots, paste0("GENE_NRFD-vs-LENGTH_median0_IZ-E-CD.pdf"))
   #plotCYS0(file.name, paste0(BASE, " early IZ, CD genes"), "IZ efficiency", "Gene length [log10]", tpm.gene.log2.m.rfd.ctr.iz.e.cd$GENE_NRFD, log10(tpm.gene.log2.m.rfd.ctr.iz.e.cd$length), "topright")
   #file.name <- file.path(wd.de.plots, paste0("GENE_NRFD-vs-LENGTH_median0_IZ-E-HO.pdf"))
   #plotCYS0(file.name, paste0(BASE, " early IZ, HO genes"), "IZ efficiency", "Gene length [log10]", tpm.gene.log2.m.rfd.ctr.iz.e.ho$GENE_NRFD, log10(tpm.gene.log2.m.rfd.ctr.iz.e.ho$length), "topright")

   ## TZ
   main.txt <- paste0(BASE, " early TZ genes")
   cols=c("blue", "blue")

   file.name <- paste0("boxplot0_", base, "_tpm.gene_median0_TSS-NRFD_TZ_E_HO+CD")
   plotBoxTSSNRFD(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.e.ho, tpm.gene.log2.m.rfd.ctr.tz.e.cd, main=main.txt, names, "TZ efficiency [TSS]", cols, ylim, isFlip=T, height=5, width=3.2)
   file.name <- paste0("boxplot1_", base, "_tpm.gene_median0_TSS-NRFD_TZ_E_HO+CD")
   plotBoxTSSNRFD(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.e.ho, tpm.gene.log2.m.rfd.ctr.tz.e.cd, main="TSS", names, "TZ efficiency", cols, ylim, isFlip=T, height=5, width=3.2)

   file.name <- paste0("boxplot_", base, "_tpm.gene_median0_GENE-NRFD_TZ_E_HO+CD")
   plotBoxGENENRFD(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.e.ho, tpm.gene.log2.m.rfd.ctr.tz.e.cd, main=main.txt, names, "TZ efficiency [Gene body]", cols, ylim, isFlip=T)
   file.name <- paste0("boxplot0_", base, "_tpm.gene_median0_GENE-NRFD_TZ_E_HO+CD")
   plotBoxGENENRFD(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.e.ho, tpm.gene.log2.m.rfd.ctr.tz.e.cd, main=main.txt, names, "TZ efficiency [Gene body]", cols, ylim, isFlip=T, height=5, width=3.2)
   file.name <- paste0("boxplot1_", base, "_tpm.gene_median0_GENE-NRFD_TZ_E_HO+CD")
   plotBoxGENENRFD(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.e.ho, tpm.gene.log2.m.rfd.ctr.tz.e.cd, main="Gene body", names, "TZ efficiency", cols, ylim, isFlip=T, height=5, width=3.2)
}

plotBoxLength <- function(wd.de.plots, file.name, tpm.1, tpm.2, main, names, cols, ylim, height=6, width=3) {
   trait <- rep(0, nrow(tpm.1))
   trait <- c(trait, rep(1, nrow(tpm.2)))
   trait <- as.factor(trait)
  
   expr <- as.numeric(c(log10(tpm.1$length), log10(tpm.2$length)))
 
   pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=height, width=width)
   boxplot(expr ~ trait, outline=T, xaxt="n", ylab=paste0("Gene length [log10]"), main=main, boxcol=cols, whiskcol=cols, outcol=cols, medcol=cols, staplecol=cols, ylim=ylim, cex.axis=1.1, cex.lab=1.2, cex.main=1.25)
 
   p <- testU(log10(tpm.1$length), log10(tpm.2$length))
   text <- ""
   if (p < 1E-3)  text <- "*"
   if (p < 1E-9)  text <- "**"
   if (p < 1E-15) text <- "***"
   text(1.5, ylim[2], text, col="black", cex=2.5)
 
   ##
   trait <- rep(0, nrow(tpm.1))
   trait <- c(trait, rep(1, nrow(tpm.2)))
   trait <- as.factor(trait)
   expr <- as.numeric(c(log10(tpm.1$length), log10(tpm.2$length)))
 
   axis(side=1, at=seq(1, 2, by=1), labels=names, font=2, cex.axis=1.2)
   #axis(side=1, at=2, labels="Total n=30,978", line=1.3, col=NA, cex.axis=1.2)
   axis(side=1, at=1, labels=paste0("n=", format(nrow(tpm.1), big.mark=",", scientific=F)), line=1.3, col=NA, cex.axis=1.2)
   axis(side=1, at=2, labels=paste0("n=", format(nrow(tpm.2), big.mark=",", scientific=F)), line=1.3, col=NA, cex.axis=1.2)
 
   mtext(paste0("p-value = ", scientific(wilcox.test(expr ~ trait, exact=F)$p.value)), cex=1.25, line=0.3)
   dev.off()
}

# -----------------------------------------------------------------------------
# 
# Last Modified: 29/05/20
# -----------------------------------------------------------------------------
plotCYS0 <- function(file.name, main.text, ylab.text, xlab.text, cn, snr, pos) {
   #unit <- (max(snr) - min(snr))/10
   #xlim <- c(min(snr) - unit, max(snr) + unit)
   #unit <- (max(cn) - min(cn))/10
   #ylim <- c(min(cn) - unit, max(cn) + unit)
 
   pdf(paste0(file.name, ".pdf"), height=6, width=6)
   plot(cn ~ snr, ylab="", xlab=xlab.text, main=main.text, pch=1, cex=2, cex.axis=1.7, cex.lab=1.9, cex.main=1.9)
 
   lm.fit <- lm(cn ~ snr)
   abline(lm.fit, lwd=5)
 
   cor <- cor.test(cn, snr, method="spearman", exact=F)
   legend(pos, c(paste0("rho = ", round0(cor[[4]], digits=2)), paste0("p-value = ", scientific(cor[[3]], digits=2))), bty="n", cex=1.9)
 
   #axis(side=1, at=seq(-0.5, 0.5, by=0.5), labels=c(-0.5, 0, 0.5), cex.axis=1.7)
   #axis(side=2, at=seq(6, 8, by=1), labels=c(6, 7, 8), cex.axis=1.7)   ## MARS
   #axis(side=2, at=seq(4, 6, by=1), labels=c(4, 5, 6), cex.axis=1.7)   ## GTPBP3
   mtext(ylab.text, side=2, line=2.74, cex=1.85)
   #mtext(main.text[2], cex=1.2, line=0.3)
   dev.off()
}

plotCYS <- function(gene, cn, snr, pch, col, pos) {
   unit <- (max(snr) - min(snr))/10
   xlim <- c(min(snr) - unit, max(snr) + unit)
   unit <- (max(cn) - min(cn))/10
   ylim <- c(min(cn) - unit, max(cn) + unit)
 
   xlab.text <- "In-silico sorting [rho]"
   ylab.text <- "log2(TPM + 0.01)"
   id <- subset(ensGene, external_gene_name == gene)$ensembl_gene_id
   file.name <- file.path(wd.de.plots, paste0("TPM-vs-SORTING_", genes[g], ""))
 
   pdf(paste0(file.name, ".pdf"), height=6, width=6)
   plot(cn ~ snr, ylim=ylim, xlim=xlim, ylab="", xaxt="n", xlab=xlab.text, main=paste0(gene, " (", id, ")"), col=col, pch=pch, cex=2, cex.axis=1.7, cex.lab=1.9, cex.main=1.9)
 
   lm.fit <- lm(cn ~ snr)
   abline(lm.fit, col="gold", lwd=5)
 
   cor <- cor.test(cn, snr, method="spearman", exact=F)
   legend(pos, c(paste0("rho = ", round0(cor[[4]], digits=2)), paste0("p-value = ", scientific(cor[[3]], digits=2))), text.col=col, bty="n", cex=1.9)
 
   axis(side=1, at=seq(-0.5, 0.5, by=0.5), labels=c(-0.5, 0, 0.5), cex.axis=1.7)
   #axis(side=2, at=seq(6, 8, by=1), labels=c(6, 7, 8), cex.axis=1.7)   ## MARS
   #axis(side=2, at=seq(4, 6, by=1), labels=c(4, 5, 6), cex.axis=1.7)   ## GTPBP3
   mtext(ylab.text, side=2, line=2.74, cex=1.85)
   #mtext(main.text[2], cex=1.2, line=0.3)
   dev.off()
}

# -----------------------------------------------------------------------------
# Gene length vs. RFD slopes
# Last Modified: 10/01/20
# -----------------------------------------------------------------------------
plotTRC <- function(rfd, tpm, main.text, file.name, xlim, ylim, col, pos, ylab.text, isFlip=F) {
   xlab.text <- "Gene length [log10]"
 
   pdf(paste0(file.name, ".pdf"), height=6, width=6)
   plot(rfd ~ tpm, ylim=ylim, xlim=xlim, ylab="", xlab=xlab.text, main="", col=col[1], pch=1, cex=1.2, cex.axis=1.5, cex.lab=1.55, yaxt="n")
 
   lm.fit <- lm(rfd ~ tpm)
   abline(lm.fit, col=col[2], lwd=3)
 
   cor <- cor.test(rfd, tpm, method="spearman", exact=F)
   legend(pos, c(paste0("rho = ", round0(cor[[4]], digits=2)), paste0("p-value = ", scientific(cor[[3]], digits=2))), text.col=col[2], bty="n", cex=1.55)
 
   if (isFlip)
      axis(side=2, at=seq(0, 0.2, by=0.1), labels=c(0, -0.1, -0.2), cex.axis=1.5)
   else
      axis(side=2, at=seq(0, 0.2, by=0.1), labels=c(0, 0.1, 0.2), cex.axis=1.5)
   mtext(ylab.text, side=2, line=2.85, cex=1.5)
   mtext(main.text[1], line=1.8, cex=1.55, font=2)
   mtext(main.text[2], line=0.3, cex=1.55)
   dev.off()
}

plotTRC0 <- function(rfd, tpm, main.text, file.name, xlim, ylim, col, pos, ylab.text, isFlip=F) {
   xlab.text <- "Gene length [log10]"
 
   pdf(paste0(file.name, ".pdf"), height=6, width=6)
   plot(rfd ~ tpm, ylim=ylim, xlim=xlim, ylab="", xlab=xlab.text, main="", col="white", pch=1, cex=1.2, cex.axis=1.5, cex.lab=1.55, yaxt="n")
 
   lm.fit <- lm(rfd ~ tpm)
   abline(lm.fit, col=col[1], lwd=3)
 
   cor <- cor.test(rfd, tpm, method="spearman", exact=F)
   legend(pos, c(paste0("rho = ", round0(cor[[4]], digits=2)), paste0("p-value = ", scientific(cor[[3]], digits=2))), text.col=col[2], bty="n", cex=1.55)
 
   if (isFlip)
      axis(side=2, at=seq(0, 0.2, by=0.1), labels=c(0, -0.1, -0.2), cex.axis=1.5)
   else
      axis(side=2, at=seq(0, 0.2, by=0.1), labels=c(0, 0.1, 0.2), cex.axis=1.5)
   mtext(ylab.text, side=2, line=2.85, cex=1.5)
   mtext(main.text[1], line=1.8, cex=1.55, font=2)
   mtext(main.text[2], line=0.3, cex=1.55)
   dev.off()
}