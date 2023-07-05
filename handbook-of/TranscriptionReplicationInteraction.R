# =============================================================================
# Library      : Transcription-Replication Interaction
# Name         : handbook-of/TranscriptionReplicationInteraction.R
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
   
   tpm.gene.log2.m.rfd.ctr.iz <- subset(tpm.gene.log2.m.rfd.ctr, TSS_NRFD >= 0)
   tpm.gene.log2.m.rfd.ctr.tz <- subset(tpm.gene.log2.m.rfd.ctr, TSS_NRFD < 0)

   tpm.gene.log2.m.rfd.ctr.iz.cd <- subset(tpm.gene.log2.m.rfd.ctr.iz, TRC >= 0)
   tpm.gene.log2.m.rfd.ctr.iz.ho <- subset(tpm.gene.log2.m.rfd.ctr.iz, TRC < 0)
   tpm.gene.log2.m.rfd.ctr.tz.cd <- subset(tpm.gene.log2.m.rfd.ctr.tz, TRC >= 0)
   tpm.gene.log2.m.rfd.ctr.tz.ho <- subset(tpm.gene.log2.m.rfd.ctr.tz, TRC < 0)
 
   ##
   tpm.gene.log2.m.rfd.ctr.iz.e <- subset(tpm.gene.log2.m.rfd.ctr.iz, RT >= 0)
   tpm.gene.log2.m.rfd.ctr.iz.l <- subset(tpm.gene.log2.m.rfd.ctr.iz, RT < 0)
   tpm.gene.log2.m.rfd.ctr.tz.e <- subset(tpm.gene.log2.m.rfd.ctr.tz, RT >= 0)
   tpm.gene.log2.m.rfd.ctr.tz.l <- subset(tpm.gene.log2.m.rfd.ctr.tz, RT < 0)
   
   tpm.gene.log2.m.rfd.ctr.iz.e.cd <- subset(tpm.gene.log2.m.rfd.ctr.iz.e, TRC >= 0)
   tpm.gene.log2.m.rfd.ctr.iz.e.ho <- subset(tpm.gene.log2.m.rfd.ctr.iz.e, TRC < 0)
   tpm.gene.log2.m.rfd.ctr.iz.l.cd <- subset(tpm.gene.log2.m.rfd.ctr.iz.l, TRC >= 0)
   tpm.gene.log2.m.rfd.ctr.iz.l.ho <- subset(tpm.gene.log2.m.rfd.ctr.iz.l, TRC < 0)
   tpm.gene.log2.m.rfd.ctr.tz.e.cd <- subset(tpm.gene.log2.m.rfd.ctr.tz.e, TRC >= 0)
   tpm.gene.log2.m.rfd.ctr.tz.e.ho <- subset(tpm.gene.log2.m.rfd.ctr.tz.e, TRC < 0)
   tpm.gene.log2.m.rfd.ctr.tz.l.cd <- subset(tpm.gene.log2.m.rfd.ctr.tz.l, TRC >= 0)
   tpm.gene.log2.m.rfd.ctr.tz.l.ho <- subset(tpm.gene.log2.m.rfd.ctr.tz.l, TRC < 0)
   
   ## TTR
   tpm.gene.log2.m.rfd.ttr <- rbind(subset(tpm.gene.log2.m.rfd, TSS_RFD >= rfd), subset(tpm.gene.log2.m.rfd, TSS_RFD <= -rfd))
   
   tpm.gene.log2.m.rfd.ttr.e <- subset(tpm.gene.log2.m.rfd.ttr, RT >= 0)
   tpm.gene.log2.m.rfd.ttr.l <- subset(tpm.gene.log2.m.rfd.ttr, RT < 0)

   tpm.gene.log2.m.rfd.ttr.e.cd <- subset(tpm.gene.log2.m.rfd.ttr.e, TRC >= 0)
   tpm.gene.log2.m.rfd.ttr.e.ho <- subset(tpm.gene.log2.m.rfd.ttr.e, TRC < 0)
   tpm.gene.log2.m.rfd.ttr.l.cd <- subset(tpm.gene.log2.m.rfd.ttr.l, TRC >= 0)
   tpm.gene.log2.m.rfd.ttr.l.ho <- subset(tpm.gene.log2.m.rfd.ttr.l, TRC < 0)
   
   save(file=file.name, tpm.gene.log2.m.rfd, tpm.gene.log2.m.rfd.ttr, tpm.gene.log2.m.rfd.ttr.e, tpm.gene.log2.m.rfd.ttr.l, tpm.gene.log2.m.rfd.ttr.e.cd, tpm.gene.log2.m.rfd.ttr.e.ho, tpm.gene.log2.m.rfd.ttr.l.cd, tpm.gene.log2.m.rfd.ttr.l.ho, tpm.gene.log2.m.rfd.ctr, tpm.gene.log2.m.rfd.ctr.iz, tpm.gene.log2.m.rfd.ctr.iz.cd, tpm.gene.log2.m.rfd.ctr.iz.ho, tpm.gene.log2.m.rfd.ctr.tz, tpm.gene.log2.m.rfd.ctr.tz.cd, tpm.gene.log2.m.rfd.ctr.tz.ho, tpm.gene.log2.m.rfd.ctr.iz.e, tpm.gene.log2.m.rfd.ctr.iz.l, tpm.gene.log2.m.rfd.ctr.tz.e, tpm.gene.log2.m.rfd.ctr.tz.l, tpm.gene.log2.m.rfd.ctr.iz.e.cd, tpm.gene.log2.m.rfd.ctr.iz.e.ho, tpm.gene.log2.m.rfd.ctr.iz.l.cd, tpm.gene.log2.m.rfd.ctr.iz.l.ho, tpm.gene.log2.m.rfd.ctr.tz.e.cd, tpm.gene.log2.m.rfd.ctr.tz.e.ho, tpm.gene.log2.m.rfd.ctr.tz.l.cd, tpm.gene.log2.m.rfd.ctr.tz.l.ho)
}

# -----------------------------------------------------------------------------
# 
# Last Modified: 27/11/20
# -----------------------------------------------------------------------------
getPvalueSignificanceLevel <- function(p) {
   text <- ""
   if (p < 1E-3) text <- "*"
   if (p < 1E-6) text <- "**"
   if (p < 1E-9) text <- "***"
   
   return(text)
}

plotBox <- function(wd.de.plots, file.name, tpm.1, tpm.2, main, names, cols, ylim, height=6, width=3) {
   trait <- rep(0, nrow(tpm.1))
   trait <- c(trait, rep(1, nrow(tpm.2)))
   trait <- as.factor(trait)
   expr <- as.numeric(c(tpm.1$MEDIAN, tpm.2$MEDIAN))
 
   pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=height, width=width)
   boxplot(expr ~ trait, outline=T, xaxt="n", ylab=paste0("log2(TPM + 1)"), main=main, boxcol=cols, whiskcol=cols, outcol=cols, medcol=cols, staplecol=cols, ylim=ylim, cex.axis=1.2, cex.lab=1.25, cex.main=1.3)
   
   p <- testU(tpm.1$MEDIAN, tpm.2$MEDIAN)
   text(1.5, ylim[2], getPvalueSignificanceLevel(p), col="black", cex=2.5)
 
   ##
   axis(side=1, at=seq(1, 2, by=1), labels=names, font=2, cex.axis=1.25)
   #axis(side=1, at=2, labels="Total n=30,978", line=1.3, col=NA, cex.axis=1.2)
   axis(side=1, at=1, labels=paste0("n=", format(nrow(tpm.1), big.mark=",", scientific=F)), line=1.2, col=NA, cex.axis=1.2)
   axis(side=1, at=2, labels=paste0("n=", format(nrow(tpm.2), big.mark=",", scientific=F)), line=1.2, col=NA, cex.axis=1.2)
   
   mtext(paste0("p-value = ", scientific(wilcox.test(expr ~ trait, exact=F)$p.value)), cex=1.25, line=0.3)
   dev.off()
}

plotBox0 <- function(wd.de.plots, file.name, tpm.1, tpm.2, main, names, cols, ylim) {
   plotBox(wd.de.plots, file.name, tpm.1, tpm.2, main, names, cols, ylim, height=5, width=3.2)
}

plotBox00 <- function(wd.de.plots, file.name, tpm.1, tpm.2, main, names, cols, ylim, height=6, width=3) {
   trait <- rep(0, nrow(tpm.1))
   trait <- c(trait, rep(1, nrow(tpm.2)))
   trait <- as.factor(trait)
   expr <- as.numeric(c(tpm.1$MEDIAN, tpm.2$MEDIAN))
 
   pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=height, width=width)
   boxplot(expr ~ trait, outline=T, xaxt="n", ylab=paste0("log2(TPM + 0.01)"), main=main, boxcol=cols, whiskcol=cols, outcol=cols, medcol=cols, staplecol=cols, ylim=ylim, cex.axis=1.2, cex.lab=1.25, cex.main=1.3)
 
   p <- testU(tpm.1$MEDIAN, tpm.2$MEDIAN)
   text(1.5, ylim[2], getPvalueSignificanceLevel(p), col="black", cex=2.5)
 
   ##
   axis(side=1, at=seq(1, 2, by=1), labels=names, font=2, cex.axis=1.25)
   #axis(side=1, at=2, labels="Total n=30,978", line=1.3, col=NA, cex.axis=1.2)
   axis(side=1, at=1, labels=paste0("n=", format(nrow(tpm.1), big.mark=",", scientific=F)), line=1.2, col=NA, cex.axis=1.2)
   axis(side=1, at=2, labels=paste0("n=", format(nrow(tpm.2), big.mark=",", scientific=F)), line=1.2, col=NA, cex.axis=1.2)
 
   mtext(paste0("p-value = ", scientific(wilcox.test(expr ~ trait, exact=F)$p.value)), cex=1.25, line=0.3)
   dev.off()
}

plotBox02 <- function(wd.de.plots, file.name, tpm.1, tpm.2, main, names, cols) {
   trait <- rep(0, length(tpm.1))
   trait <- c(trait, rep(1, length(tpm.2)))
   trait <- as.factor(trait)
 
   expr <- as.numeric(c(tpm.1, tpm.2))
   ylim <- c(min(expr), max(expr))
 
   pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=6, width=4.5)
   boxplot(expr ~ trait, outline=T, names=names, main="", col=cols, xaxt="n", xlab="", ylab="", ylim=ylim, cex=2, cex.axis=1.8, cex.lab=1.9, cex.main=2)
 
   p <- wilcox.test(expr ~ trait, exact=F)$p.value
   offset <- (ylim[2] - ylim[1])/35
   text(1.5, ylim[2] - offset, getPvalueSignificanceLevel(p), col="black", cex=4)
 
   axis(side=1, at=seq(1, 2, by=1), labels=names, font=2, cex.axis=1.9)
   axis(side=1, at=1, labels=paste0("n=", length(tpm.1)), line=1.8, col=NA, cex.axis=1.9)
   axis(side=1, at=2, labels=paste0("n=", length(tpm.2)), line=1.8, col=NA, cex.axis=1.9)
   #axis(side=2, at=seq(5, 8, by=1), cex.axis=1.2)
   mtext(main, font=2, cex=2, line=2)
   mtext(paste0("p-value = ", scientific(p)), cex=1.9, line=0.3)
   mtext("Age at diagnosis", side=2, line=2.7, cex=1.9)
   dev.off()
}

plotBox020 <- function(wd.de.plots, file.name, tpm.1, tpm.2, main, names, cols, ylim) {
   trait <- rep(0, length(tpm.1))
   trait <- c(trait, rep(1, length(tpm.2)))
   trait <- as.factor(trait)
 
   expr <- as.numeric(c(tpm.1, tpm.2))

   pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=6, width=4.5)
   boxplot(expr ~ trait, outline=T, names=names, main="", col=cols, xaxt="n", ylab="", ylim=ylim, cex=2, cex.axis=1.7, cex.lab=1.9, cex.main=2)
 
   p <- wilcox.test(expr ~ trait, exact=F)$p.value
   offset <- (ylim[2] - ylim[1])/35
   text(1.5, ylim[2] - offset, getPvalueSignificanceLevel(p), col="black", cex=4)
 
   axis(side=1, at=seq(1, 2, by=1), labels=names, font=2, cex.axis=1.9)
   axis(side=1, at=1, labels=paste0("n=", length(tpm.1)), line=1.8, col=NA, cex.axis=1.9)
   axis(side=1, at=2, labels=paste0("n=", length(tpm.2)), line=1.8, col=NA, cex.axis=1.9)
   #axis(side=2, at=seq(5, 8, by=1), cex.axis=1.2)
   mtext(main, font=2, cex=2, line=2)
   mtext(paste0("p-value = ", scientific(p)), cex=2, line=0.3)
   mtext(expression(italic('in silico')~"sorting [rho]"), side=2, line=2.48, cex=1.85)
   dev.off()
}

plotStripchart <- function(wd.de.plots, file.name, phenos, main, names, cols, height, width) {
   pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=height, width=width)
   boxplot(COR ~ group, data=phenos, xaxt="n", ylab="", main=main, col="white", outline=F, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
   
   #stripchart(COR ~ group, data=subset(phenos, M2 == 0), method="jitter", cex=2, pch=19, col=cols[1], vertical=T, add=T)
   #stripchart(COR ~ group, data=subset(phenos, M2 == 1), method="jitter", cex=2, pch=19, col=cols[2], vertical=T, add=T, at=2:4)
   stripchart(COR ~ group, data=subset(phenos, Q4 == 1), method="jitter", cex=1.5, pch=19, col=cols[1], vertical=T, add=T)
   stripchart(COR ~ group, data=subset(phenos, Q4 == 2), method="jitter", cex=1.5, pch=19, col=cols[2], vertical=T, add=T, at=c(1,3,4))
   stripchart(COR ~ group, data=subset(phenos, Q4 == 3), method="jitter", cex=1.5, pch=19, col=cols[3], vertical=T, add=T, at=2:4)
   stripchart(COR ~ group, data=subset(phenos, Q4 == 4), method="jitter", cex=1.5, pch=19, col=cols[4], vertical=T, add=T, at=2:4)
   
   #p <- testU(tpm.1$COR, tpm.3$COR)
   #text(2, ylim[2], getPvalueSignificanceLevel(p), cex=2.5)
   #lines(c(1, 3), y=c(ylim[2]-0.8, ylim[2]-0.8), type="l", lwd=2)
 
   #p <- testU(tpm.1$MEDIAN, tpm.2$MEDIAN)
   #text(1.5, ylim[2]-1.6, getPvalueSignificanceLevel(p), cex=2.5)
   #lines(c(1, 2), y=c(ylim[2]-2.4, ylim[2]-2.4), type="l", lwd=2)
 
   #p <- testU(tpm.2$MEDIAN, tpm.3$MEDIAN)
   #text(2.5, ylim[2]-3.2, getPvalueSignificanceLevel(p), cex=2.5)
   #lines(c(2, 3), y=c(ylim[2]-4, ylim[2]-4), type="l", lwd=2)
 
   axis(side=1, at=seq(1, 4, by=1), labels=names, font=2, cex.axis=1.7)
   axis(side=1, at=1, labels=paste0("n=", nrow(subset(phenos, group == 0))), line=1.35, col=NA, cex.axis=1.7)
   axis(side=1, at=2, labels=paste0("n=", nrow(subset(phenos, group == 1))), line=1.35, col=NA, cex.axis=1.7)
   axis(side=1, at=3, labels=paste0("n=", nrow(subset(phenos, group == 2))), line=1.35, col=NA, cex.axis=1.7)
   axis(side=1, at=4, labels=paste0("n=", nrow(subset(phenos, group == 3))), line=1.35, col=NA, cex.axis=1.7)
   
   #mtext(paste0(text), cex=1.25, line=0.3)
   mtext(expression(italic('in silico')~"sorting [rho]"), side=2, line=2.38, cex=1.85)
   legend("topleft", legend=c("Q4", "Q3", "Q2", "Q1"), pch=19, pt.cex=2.5, col=c(red, red.light, blue.light, blue), cex=1.8)
   dev.off()
}

plotStripchartATRX <- function(wd.de.plots, file.name, phenos, main, names, cols, height, width) {
   pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=height, width=width)
   boxplot(COR ~ group, data=phenos, xaxt="n", ylab="", main=main, col="white", outline=F, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
 
   #stripchart(COR ~ group, data=subset(phenos, M2 == 0), method="jitter", cex=2, pch=19, col=cols[1], vertical=T, add=T)
   #stripchart(COR ~ group, data=subset(phenos, M2 == 1), method="jitter", cex=2, pch=19, col=cols[2], vertical=T, add=T, at=2:4)
   stripchart(COR ~ group, data=subset(phenos, Q4 == 1), method="jitter", cex=1.5, pch=19, col=cols[1], vertical=T, add=T, at=c(1,3,4))
   stripchart(COR ~ group, data=subset(phenos, Q4 == 2), method="jitter", cex=1.5, pch=19, col=cols[2], vertical=T, add=T, at=c(1,2,3,5))
   stripchart(COR ~ group, data=subset(phenos, Q4 == 3), method="jitter", cex=1.5, pch=19, col=cols[3], vertical=T, add=T, at=2:5)
   stripchart(COR ~ group, data=subset(phenos, Q4 == 4), method="jitter", cex=1.5, pch=19, col=cols[4], vertical=T, add=T, at=2:5)
 
   #p <- testU(tpm.1$COR, tpm.3$COR)
   #text(2, ylim[2], getPvalueSignificanceLevel(p), cex=2.5)
   #lines(c(1, 3), y=c(ylim[2]-0.8, ylim[2]-0.8), type="l", lwd=2)
 
   #p <- testU(tpm.1$MEDIAN, tpm.2$MEDIAN)
   #text(1.5, ylim[2]-1.6, getPvalueSignificanceLevel(p), cex=2.5)
   #lines(c(1, 2), y=c(ylim[2]-2.4, ylim[2]-2.4), type="l", lwd=2)
 
   #p <- testU(tpm.2$MEDIAN, tpm.3$MEDIAN)
   #text(2.5, ylim[2]-3.2, getPvalueSignificanceLevel(p), cex=2.5)
   #lines(c(2, 3), y=c(ylim[2]-4, ylim[2]-4), type="l", lwd=2)
 
   axis(side=1, at=seq(1, 5, by=1), labels=names, font=2, cex.axis=1.7)
   axis(side=1, at=1, labels=paste0("n=", nrow(subset(phenos, group == 0))), line=1.35, col=NA, cex.axis=1.7)
   axis(side=1, at=2, labels=paste0("n=", nrow(subset(phenos, group == 1))), line=1.35, col=NA, cex.axis=1.7)
   axis(side=1, at=3, labels=paste0("n=", nrow(subset(phenos, group == 2))), line=1.35, col=NA, cex.axis=1.7)
   axis(side=1, at=4, labels=paste0("n=", nrow(subset(phenos, group == 3))), line=1.35, col=NA, cex.axis=1.7)
   axis(side=1, at=5, labels=paste0("n=", nrow(subset(phenos, group == 4))), line=1.35, col=NA, cex.axis=1.7)
   
   #mtext(paste0(text), cex=1.25, line=0.3)
   mtext(expression(italic('in silico')~"sorting [rho]"), side=2, line=2.38, cex=1.85)
   legend("topleft", legend=c("Q4", "Q3", "Q2", "Q1"), pch=19, pt.cex=2.5, col=c(red, red.light, blue.light, blue), cex=1.8)
   dev.off()
}

plotBox2 <- function(wd.de.plots, file.name, tpm.1, tpm.2, main, names, cols, ylim, ylab=text.Log2.TPM.1) {
   trait <- rep(0, nrow(tpm.1))
   trait <- c(trait, rep(1, nrow(tpm.2)))
   trait <- as.factor(trait)
   expr <- as.numeric(c(tpm.1$MEDIAN, tpm.2$MEDIAN))
   #ylim <- c(min(expr), max(expr))
 
   pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=6, width=3.8)
   par(mar=c(4.5, 4.7, 4.1, 1.4))
   boxplot(expr ~ trait, outline=T, names=names, xlab="", ylab=ylab, main=main, xaxt="n", col="white", boxcol=cols, whiskcol=cols, outcol=cols, medcol=cols, staplecol=cols, ylim=ylim, yaxt="n", cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
 
   p <- testU(tpm.1$MEDIAN, tpm.2$MEDIAN)
   text(1.5, ylim[2]-0.2, getPvalueSignificanceLevel(p), cex=3)
   lines(c(1, 2), y=c(ylim[2]-0.8, ylim[2]-0.8), type="l", lwd=2)
   
   axis(side=2, at=seq(0, 12, by=2), labels=c(0, "", 4, "", 8, "", 12), cex.axis=1.7)
   
   axis(side=1, at=seq(1, 2, by=1), labels=c("", ""), cex.axis=1.8)
   axis(side=1, at=seq(1, 2, by=1), labels=names, line=0.5, col=NA, cex.axis=1.8)
   axis(side=1, at=1, labels=paste0("n=", separator(nrow(tpm.1))), line=2.2, col=NA, cex.axis=1.7)
   axis(side=1, at=2, labels=paste0("n=", separator(nrow(tpm.2))), line=2.2, col=NA, cex.axis=1.7)
   
   dev.off()
}

plotBox2length <- function(wd.de.plots, file.name, tpm.1, tpm.2, main, names, cols, ylim, ylab=text.Log2.TPM.1) {
	  trait <- rep(0, length(tpm.1))
	  trait <- c(trait, rep(1, length(tpm.2)))
	  trait <- as.factor(trait)
	  expr <- as.numeric(c(tpm.1, tpm.2))
	  #ylim <- c(min(expr), max(expr))
	
	  pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=6, width=3.8)
	  par(mar=c(4.5, 4.7, 4.1, 1.4))
	  boxplot(expr ~ trait, outline=T, names=names, xlab="", ylab=ylab, main=main, xaxt="n", col="white", boxcol=cols, whiskcol=cols, outcol=cols, medcol=cols, staplecol=cols, ylim=ylim, yaxt="n", cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
	
	  p <- testU(tpm.1, tpm.2)
	  text(1.5, ylim[2]-0.2, getPvalueSignificanceLevel(p), cex=3)
	  lines(c(1, 2), y=c(ylim[2]-0.8, ylim[2]-0.8), type="l", lwd=2)
	
	  axis(side=2, at=seq(0, 12, by=2), labels=c(0, "", 4, "", 8, "", 12), cex.axis=1.7)
	
	  axis(side=1, at=seq(1, 2, by=1), labels=c("", ""), cex.axis=1.8)
	  axis(side=1, at=seq(1, 2, by=1), labels=names, line=0.5, col=NA, cex.axis=1.8)
	  axis(side=1, at=1, labels=paste0("n=", separator(nrow(tpm.1))), line=2.2, col=NA, cex.axis=1.7)
	  axis(side=1, at=2, labels=paste0("n=", separator(nrow(tpm.2))), line=2.2, col=NA, cex.axis=1.7)
	
	  dev.off()
}

# plotBox20 <- function(wd.de.plots, file.name, tpm.1, tpm.2, main, names) {
#    trait <- rep(0, length(tpm.1))
#    trait <- c(trait, rep(1, length(tpm.2)))
#    trait <- as.factor(trait)
#  
#    expr <- as.numeric(c(tpm.1, tpm.2))
#    ylim <- c(min(expr), max(expr))
#  
#    pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=6, width=3)
#    boxplot(expr ~ trait, outline=T, names=names, ylab=paste0("log2(TPM + 0.01)"), main=main, xaxt="n", yaxt="n", ylim=ylim, cex.axis=1.1, cex.lab=1.2, cex.main=1.25)
#  
#    axis(side=1, at=seq(1, 2, by=1), labels=names, font=2, cex.axis=1.2)
#    axis(side=2, at=seq(5, 8, by=1), labels=c(5, 6, 7, 8), cex.axis=1.1)
#    mtext(paste0("p-value = ", scientific(wilcox.test(expr ~ trait, exact=F)$p.value)), cex=1.25, line=0.3)
#    dev.off()
# }

plotBox3 <- function(wd.de.plots, file.name, tpm.1, tpm.2, tpm.3, main, names, cols, ylim, height=6, width=4, text="", ylab=text.Log2.TPM.1) {
   trait <- rep(0, nrow(tpm.1))
   trait <- c(trait, rep(1, nrow(tpm.2)))
   trait <- c(trait, rep(2, nrow(tpm.3)))
   trait <- as.factor(trait)
   expr <- as.numeric(c(tpm.1$MEDIAN, tpm.2$MEDIAN, tpm.3$MEDIAN))
 
   pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=height, width=width)
   par(mar=c(4.5, 4.7, 4.1, 1.4))
   boxplot(expr ~ trait, outline=F, xlab="", xaxt="n", ylab=ylab, main=main, col="white", boxcol=cols, whiskcol=cols, outcol=cols, medcol=cols, staplecol=cols, ylim=ylim, yaxt="n", cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
 
   p <- testU(tpm.1$MEDIAN, tpm.3$MEDIAN)
   text(2, ylim[2]-0.2, getPvalueSignificanceLevel(p), cex=3)
   lines(c(1, 3), y=c(ylim[2]-0.8, ylim[2]-0.8), type="l", lwd=2)
   
   p <- testU(tpm.1$MEDIAN, tpm.2$MEDIAN)
   text(1.5, ylim[2]-1.4, getPvalueSignificanceLevel(p), cex=3)
   lines(c(1, 2), y=c(ylim[2]-2, ylim[2]-2), type="l", lwd=2)
   
   p <- testU(tpm.2$MEDIAN, tpm.3$MEDIAN)
   text(2.5, ylim[2]-2.6, getPvalueSignificanceLevel(p), cex=3)
   lines(c(2, 3), y=c(ylim[2]-3.2, ylim[2]-3.2), type="l", lwd=2)
   
   axis(side=2, at=seq(0, 12, by=2), labels=c(0, "", 4, "", 8, "", 12), cex.axis=1.7)
   axis(side=1, at=seq(1, 3, by=1), labels=c("", "", ""), cex.axis=1.8)
   axis(side=1, at=seq(1, 3, by=1), labels=names, line=0.5, col=NA, cex.axis=1.8)

   dev.off()
}

plotBox30 <- function(wd.de.plots, file.name, tpm.1, tpm.2, tpm.3, main, names, cols, ylim, height=6, width=4, text="", ylab=text.Log2.TPM.1) {
	  trait <- rep(0, nrow(tpm.1))
	  trait <- c(trait, rep(1, nrow(tpm.2)))
	  trait <- c(trait, rep(2, nrow(tpm.3)))
	  trait <- as.factor(trait)
	  expr <- as.numeric(c(tpm.1$MEDIAN, tpm.2$MEDIAN, tpm.3$MEDIAN))
	
	  pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=height, width=width)
	  par(mar=c(4.5, 4.7, 4.1, 1.4))
	  boxplot(expr ~ trait, outline=F, xlab="", xaxt="n", ylab=ylab, main=main, col="white", boxcol=cols, whiskcol=cols, outcol=cols, medcol=cols, staplecol=cols, ylim=ylim, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
	
	  unit <- (ylim[2] - ylim[1])/10
	  p <- testU(tpm.1$MEDIAN, tpm.3$MEDIAN)
	  text(2, ylim[2]-unit+unit/2, getPvalueSignificanceLevel(p), cex=3)
	  lines(c(1, 3), y=c(ylim[2]-unit, ylim[2]-unit), type="l", lwd=2)
	
	  p <- testU(tpm.1$MEDIAN, tpm.2$MEDIAN)
	  text(1.5, ylim[2]-unit*2+unit/2, getPvalueSignificanceLevel(p), cex=3)
	  lines(c(1, 2), y=c(ylim[2]-unit*2, ylim[2]-unit*2), type="l", lwd=2)
	
	  p <- testU(tpm.2$MEDIAN, tpm.3$MEDIAN)
	  text(2.5, ylim[2]-unit*3+unit/2, getPvalueSignificanceLevel(p), cex=3)
	  lines(c(2, 3), y=c(ylim[2]-unit*3, ylim[2]-unit*3), type="l", lwd=2)
	
	  #axis(side=2, at=seq(0, 12, by=2), labels=c(0, "", 4, "", 8, "", 12), cex.axis=1.7)
	  axis(side=1, at=seq(1, 3, by=1), labels=c("", "", ""), cex.axis=1.8)
	  axis(side=1, at=seq(1, 3, by=1), labels=names, line=0.5, col=NA, cex.axis=1.8)
	
  	dev.off()
}

plotBox4 <- function(wd.de.plots, file.name, tpm.1, tpm.2, tpm.3, tpm.4, main, names, cols, ylim, height=5, width=3.2) {
   trait <- rep(0, nrow(tpm.1))
   trait <- c(trait, rep(1, nrow(tpm.2)))
   trait <- c(trait, rep(2, nrow(tpm.3)))
   trait <- c(trait, rep(3, nrow(tpm.4)))
   trait <- as.factor(trait)
   expr <- as.numeric(c(tpm.1$MEDIAN, tpm.2$MEDIAN, tpm.3$MEDIAN, tpm.4$MEDIAN))
 
   pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=height, width=width)
   boxplot(expr ~ trait, outline=F, xaxt="n", ylab=paste0("log2(TPM + 1)"), main=main, boxcol=cols, whiskcol=cols, outcol=cols, medcol=cols, staplecol=cols, ylim=ylim, yaxt="n", cex.axis=1.25, cex.lab=1.3, cex.main=1.35)
 
   p <- testU(tpm.1$MEDIAN, tpm.2$MEDIAN)
   text(1.5, ylim[2]-0.2, getPvalueSignificanceLevel(p), cex=2.5)
   lines(c(1, 2), y=c(ylim[2]-0.8, ylim[2]-0.8), type="l", lwd=2)
 
   p <- testU(tpm.3$MEDIAN, tpm.4$MEDIAN)
   text(3.5, ylim[2]-0.2, getPvalueSignificanceLevel(p), cex=2.5)
   lines(c(3, 4), y=c(ylim[2]-0.8, ylim[2]-0.8), type="l", lwd=2)
 
   axis(side=2, at=seq(0, 12, by=2), labels=c(0, "", 4, "", 8, "", 12), cex.axis=1.25)
   #axis(side=1, at=seq(1, 4, by=1), labels=c("L", "E", "L", "E"), cex.axis=1.2)
   #axis(side=1, at=1.5, labels="TZ", line=1.2, col=NA, font=2, cex.axis=1.25)
   #axis(side=1, at=3.5, labels="IZ", line=1.2, col=NA, font=2, cex.axis=1.25)
   axis(side=1, at=seq(1.5, 3.5, by=2), labels=c("TZ", "IZ"), font=2, cex.axis=1.3)
   axis(side=1, at=seq(1, 4, by=1), labels=c("(L)", "(E)", "(L)", "(E)"), line=1.2, col=NA, cex.axis=1.3)
 
   mtext(paste0("(L)ate vs. (E)arly"), cex=1.25, line=0.3)
   dev.off()
}

plotBox6 <- function(wd.de.plots, file.name, tpm.1, tpm.2, tpm.3, tpm.4, tpm.5, tpm.6, main, names, cols, ylim, height=5, width=3.5) {
   trait <- rep(0, nrow(tpm.1))
   trait <- c(trait, rep(1, nrow(tpm.2)))
   trait <- c(trait, rep(2, nrow(tpm.3)))
   trait <- c(trait, rep(3, nrow(tpm.4)))
   trait <- c(trait, rep(4, nrow(tpm.5)))
   trait <- c(trait, rep(5, nrow(tpm.6)))
   trait <- as.factor(trait)
   expr <- as.numeric(c(tpm.1$MEDIAN, tpm.2$MEDIAN, tpm.3$MEDIAN, tpm.4$MEDIAN, tpm.5$MEDIAN, tpm.6$MEDIAN))
 
   pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=height, width=width)
   boxplot(expr ~ trait, outline=F, xaxt="n", ylab=paste0("log2(TPM + 1)"), main=main, boxcol=cols, whiskcol=cols, outcol=cols, medcol=cols, staplecol=cols, ylim=ylim, yaxt="n", cex.axis=1.25, cex.lab=1.3, cex.main=1.35)
 
   p <- testU(tpm.1$MEDIAN, tpm.2$MEDIAN)
   text(1.5, ylim[2]-0.2, getPvalueSignificanceLevel(p), cex=2.5)
   lines(c(1, 2), y=c(ylim[2]-0.8, ylim[2]-0.8), type="l", lwd=2)
 
   p <- testU(tpm.3$MEDIAN, tpm.4$MEDIAN)
   text(3.5, ylim[2]-0.2, getPvalueSignificanceLevel(p), cex=2.5)
   lines(c(3, 4), y=c(ylim[2]-0.8, ylim[2]-0.8), type="l", lwd=2)
 
   p <- testU(tpm.5$MEDIAN, tpm.6$MEDIAN)
   text(5.5, ylim[2]-0.2, getPvalueSignificanceLevel(p), cex=2.5)
   lines(c(5, 6), y=c(ylim[2]-0.8, ylim[2]-0.8), type="l", lwd=2)
   
   axis(side=2, at=seq(0, 12, by=2), labels=c(0, "", 4, "", 8, "", 12), cex.axis=1.25)
   #axis(side=1, at=seq(1, 6, by=1), labels=c("L", "E", "L", "E", "L", "E"), cex.axis=1.2)
   #axis(side=1, at=1.5, labels="TTR", line=1.2, col=NA, font=2, cex.axis=1.25)
   #axis(side=1, at=3.5, labels="TZ", line=1.2, col=NA, font=2, cex.axis=1.25)
   #axis(side=1, at=5.5, labels="IZ", line=1.2, col=NA, font=2, cex.axis=1.25)
   axis(side=1, at=seq(1.5, 5.5, by=2), labels=c("TTR", "TZ", "IZ"), font=2, cex.axis=1.3)
   axis(side=1, at=seq(1, 5, by=2), labels=c("(L)", "(L)", "(L)"), line=1.2, col=NA, cex.axis=1.3)
   axis(side=1, at=seq(2, 6, by=2), labels=c("(E)", "(E)", "(E)"), line=1.2, col=NA, cex.axis=1.3)
   
   mtext(paste0("(L)ate vs. (E)arly"), cex=1.25, line=0.3)
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
      boxplot(expr ~ trait, outline=outline, xaxt="n", yaxt="n", ylab=ylab.txt, main=main.txt, col=cols, cex.axis=1.2, cex.lab=1.25, cex.main=1.3)
   else
      boxplot(expr ~ trait, outline=outline, xaxt="n", yaxt="n", ylab=ylab.txt, main=main.txt, col=cols, ylim=ylim, cex.axis=1.2, cex.lab=1.25, cex.main=1.3)
 
   p <- testU(tpm.1$TSS_NRFD, tpm.2$TSS_NRFD)
   text <- ""
   if (p < 1E-3)  text <- "*"
   if (p < 1E-6)  text <- "**"
   if (p < 1E-9) text <- "***"
   text(1.5, ylim[2], text, col="black", cex=2.5)
 
   ##
   trait <- rep(0, nrow(tpm.1))
   trait <- c(trait, rep(1, nrow(tpm.2)))
   trait <- as.factor(trait)
   expr <- as.numeric(c(tpm.1$TSS_NRFD, tpm.2$TSS_NRFD))
 
   if (isFlip)
      axis(side=2, at=seq(-0.01, 0.03, by=0.01), labels=c(0.01, 0, -0.01, -0.02, -0.03), cex.axis=1.2)
   else
      axis(side=2, at=seq(-0.01, 0.03, by=0.01), labels=c(-0.01, 0, 0.01, 0.02, 0.03), cex.axis=1.2)
 
   axis(side=1, at=seq(1, 2, by=1), labels=names, font=2, cex.axis=1.25)
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
      boxplot(expr ~ trait, outline=outline, xaxt="n", yaxt="n", ylab=ylab.txt, main=main.txt, col=cols, cex.axis=1.2, cex.lab=1.25, cex.main=1.3)
   else
      boxplot(expr ~ trait, outline=outline, xaxt="n", yaxt="n", ylab=ylab.txt, main=main.txt, col=cols, ylim=ylim, cex.axis=1.2, cex.lab=1.25, cex.main=1.3)
   
   p <- testU(tpm.1$GENE_NRFD, tpm.2$GENE_NRFD)
   text <- ""
   if (p < 1E-3)  text <- "*"
   if (p < 1E-6)  text <- "**"
   if (p < 1E-9) text <- "***"
   text(1.5, ylim[2], text, col="black", cex=2.5)
 
   ##
   trait <- rep(0, nrow(tpm.1))
   trait <- c(trait, rep(1, nrow(tpm.2)))
   trait <- as.factor(trait)
   expr <- as.numeric(c(tpm.1$GENE_NRFD, tpm.2$GENE_NRFD))
 
   if (isFlip)
      axis(side=2, at=seq(-0.01, 0.03, by=0.01), labels=c(0.01, 0, -0.01, -0.02, -0.03), cex.axis=1.2)
   else
      axis(side=2, at=seq(-0.01, 0.03, by=0.01), labels=c(-0.01, 0, 0.01, 0.02, 0.03), cex.axis=1.2)
   
   axis(side=1, at=seq(1, 2, by=1), labels=names, font=2, cex.axis=1.25)
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
   cols=c(red, red)

   #file.name <- paste0("boxplot0_", base, "_tpm.gene_median0_TSS-NRFD_IZ_E_HO+CD")
   #plotBoxTSSNRFD(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.e.ho, tpm.gene.log2.m.rfd.ctr.iz.e.cd, main=main.txt, names, "IZ efficiency [TSS]", cols, ylim, height=5, width=3.2)
   #file.name <- paste0("boxplot2_", base, "_tpm.gene_median0_TSS-NRFD_IZ_E_HO+CD")
   #plotBoxTSSNRFD(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.e.ho, tpm.gene.log2.m.rfd.ctr.iz.e.cd, main=main.txt, names, "IZ efficiency [TSS]", cols, ylim2, height=5, width=3.2, outline=T)
   file.name <- paste0("boxplot1_", base, "_tpm.gene_median0_TSS-NRFD_IZ_E_HO+CD")
   plotBoxTSSNRFD(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.e.ho, tpm.gene.log2.m.rfd.ctr.iz.e.cd, main="TSS", names, "IZ efficiency", cols, ylim, height=5, width=3.2)

   #file.name <- paste0("boxplot_", base, "_tpm.gene_median0_GENE-NRFD_IZ_E_HO+CD")
   #plotBoxGENENRFD(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.e.ho, tpm.gene.log2.m.rfd.ctr.iz.e.cd, main=main.txt, names, "IZ efficiency [Gene body]", cols, ylim)
   #file.name <- paste0("boxplot0_", base, "_tpm.gene_median0_GENE-NRFD_IZ_E_HO+CD")
   #plotBoxGENENRFD(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.e.ho, tpm.gene.log2.m.rfd.ctr.iz.e.cd, main=main.txt, names, "IZ efficiency [Gene body]", cols, ylim, height=5, width=3.2)
   #file.name <- paste0("boxplot2_", base, "_tpm.gene_median0_GENE-NRFD_IZ_E_HO+CD")
   #plotBoxGENENRFD(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.e.ho, tpm.gene.log2.m.rfd.ctr.iz.e.cd, main=main.txt, names, "IZ efficiency [Gene body]", cols, ylim2, height=5, width=3.2, outline=T)
   file.name <- paste0("boxplot1_", base, "_tpm.gene_median0_GENE-NRFD_IZ_E_HO+CD")
   plotBoxGENENRFD(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.e.ho, tpm.gene.log2.m.rfd.ctr.iz.e.cd, main="Gene body", names, "IZ efficiency", cols, ylim, height=5, width=3.2)

   #file.name <- file.path(wd.de.plots, paste0("GENE_NRFD-vs-LENGTH_median0_IZ-E-CD.pdf"))
   #plotCYS0(file.name, paste0(BASE, " early IZ, CD genes"), "IZ efficiency", "Gene length [log10]", tpm.gene.log2.m.rfd.ctr.iz.e.cd$GENE_NRFD, log10(tpm.gene.log2.m.rfd.ctr.iz.e.cd$length), "topright")
   #file.name <- file.path(wd.de.plots, paste0("GENE_NRFD-vs-LENGTH_median0_IZ-E-HO.pdf"))
   #plotCYS0(file.name, paste0(BASE, " early IZ, HO genes"), "IZ efficiency", "Gene length [log10]", tpm.gene.log2.m.rfd.ctr.iz.e.ho$GENE_NRFD, log10(tpm.gene.log2.m.rfd.ctr.iz.e.ho$length), "topright")

   ## TZ
   main.txt <- paste0(BASE, " early TZ genes")
   cols=c(blue, blue)

   #file.name <- paste0("boxplot0_", base, "_tpm.gene_median0_TSS-NRFD_TZ_E_HO+CD")
   #plotBoxTSSNRFD(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.e.ho, tpm.gene.log2.m.rfd.ctr.tz.e.cd, main=main.txt, names, "TZ efficiency [TSS]", cols, ylim, isFlip=T, height=5, width=3.2)
   file.name <- paste0("boxplot1_", base, "_tpm.gene_median0_TSS-NRFD_TZ_E_HO+CD")
   plotBoxTSSNRFD(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.e.ho, tpm.gene.log2.m.rfd.ctr.tz.e.cd, main="TSS", names, "TZ efficiency", cols, ylim, isFlip=T, height=5, width=3.2)

   #file.name <- paste0("boxplot_", base, "_tpm.gene_median0_GENE-NRFD_TZ_E_HO+CD")
   #plotBoxGENENRFD(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.e.ho, tpm.gene.log2.m.rfd.ctr.tz.e.cd, main=main.txt, names, "TZ efficiency [Gene body]", cols, ylim, isFlip=T)
   #file.name <- paste0("boxplot0_", base, "_tpm.gene_median0_GENE-NRFD_TZ_E_HO+CD")
   #plotBoxGENENRFD(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.e.ho, tpm.gene.log2.m.rfd.ctr.tz.e.cd, main=main.txt, names, "TZ efficiency [Gene body]", cols, ylim, isFlip=T, height=5, width=3.2)
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
   abline(lm.fit, lwd=2.5)
 
   cor <- cor.test(cn, snr, method="spearman", exact=F)
   legend(pos, c(paste0("rho = ", round0(cor[[4]], digits=2)), paste0("p-value = ", scientific(cor[[3]], digits=2))), bty="n", cex=1.9)
 
   #axis(side=1, at=seq(-0.5, 0.5, by=0.5), labels=c(-0.5, 0, 0.5), cex.axis=1.7)
   #axis(side=2, at=seq(6, 8, by=1), labels=c(6, 7, 8), cex.axis=1.7)   ## MARS
   #axis(side=2, at=seq(4, 6, by=1), labels=c(4, 5, 6), cex.axis=1.7)   ## GTPBP3
   mtext(ylab.text, side=2, line=2.74, cex=1.85)
   #mtext(main.text[2], cex=1.2, line=0.3)
   dev.off()
}

plotSRC <- function(gene, cn, snr, pch, col, pos, xlab.text="") {
   unit <- (max(snr) - min(snr))/10
   xlim <- c(min(snr) - unit, max(snr) + unit)
   unit <- (max(cn) - min(cn))/10
   ylim <- c(min(cn) - unit, max(cn) + unit)
 
   if (xlab.text == "")
      xlab.text <- expression(italic('In silico')~'SCF index')
   ylab.text <- expression("log" * ""[2] * "(TPM + 1)")
   id <- subset(ensGene, external_gene_name == gene)$ensembl_gene_id
   file.name <- file.path(wd.de.plots, paste0("TPM-vs-SORTING_", genes[g], ""))
 
   pdf(paste0(file.name, ".pdf"), height=5, width=5)
   par(mar=c(5.1, 4.7, 4.1, 1.4))
   plot(cn ~ snr, ylim=ylim, xlim=xlim, ylab=ylab.text, xaxt="n", xlab=xlab.text, main=gene, col="black", pch=pch, cex=2, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
   
   lm.fit <- lm(cn ~ snr)
   abline(lm.fit, col=col, lwd=7)
 
   cor <- cor.test(cn, snr, method="spearman", exact=F)
   legend("bottomright", c(paste0("rho = ", round0(cor[[4]], digits=2)), paste0("P = 1.00E-00")), text.col=c(col, "white"), text.font=2, bty="n", cex=1.8)
   legend("bottomright", expression(bolditalic('P')~"                   "), text.col=col, text.font=2, bty="n", cex=1.8)
   legend("bottomright", paste0("   = ", scientific(cor[[3]])), text.col=col, text.font=2, bty="n", cex=1.8)
   
   axis(side=1, at=seq(-0.5, 0.5, by=0.5), labels=c(-0.5, 0, 0.5), cex.axis=1.7)
   dev.off()
}

plotPuritySRC <- function(column, gene, cn, snr, pch, col, pos) {
   unit <- (max(snr) - min(snr))/10
   ylim <- c(min(snr) - unit, max(snr) + unit)
   unit <- (max(cn) - min(cn))/10
   xlim <- c(min(cn) - unit, max(cn) + unit)
 
   column <- sub("_", " ", column)
   ylab.text <- expression(italic('In silico')~'sorting [rho]')
   xlab.text <- column
   file.name <- file.path(wd.de.plots, paste0("SORTING_vs_", column, ""))
 
   pdf(paste0(file.name, ".pdf"), height=6, width=6)
   #plot(cn ~ snr, ylim=ylim, xlim=xlim, ylab="", xaxt="n", xlab=xlab.text, main=paste0(gene, " (", id, ")"), col="black", pch=pch, cex=2, cex.axis=1.7, cex.lab=1.9, cex.main=2)
   plot(snr ~ cn, ylim=ylim, xlim=xlim, xlab="", yaxt="n", ylab="", main=gene, col="black", pch=pch, cex=2, cex.axis=1.7, cex.lab=1.9, cex.main=2)
 
   lm.fit <- lm(snr ~ cn)
   abline(lm.fit, col=col, lwd=5)
 
   cor <- cor.test(snr, cn, method="spearman", exact=F)
   legend(pos, c(paste0("rho = ", round0(cor[[4]], digits=2)), paste0("p-value = ", scientific(cor[[3]], digits=2))), text.col=col, text.font=2, bty="n", cex=1.9)
 
   axis(side=2, at=seq(-0.5, 0.5, by=0.5), labels=c(-0.5, 0, 0.5), cex.axis=1.7)
   #axis(side=2, at=seq(6, 8, by=1), labels=c(6, 7, 8), cex.axis=1.7)   ## MARS
   #axis(side=2, at=seq(4, 6, by=1), labels=c(4, 5, 6), cex.axis=1.7)   ## GTPBP3
   mtext(ylab.text, side=2, line=2.35, cex=1.9)
   mtext(xlab.text, side=1, line=3, cex=1.9)
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