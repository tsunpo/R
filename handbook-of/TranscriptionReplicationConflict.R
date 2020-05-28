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
 
   ## TTR
   tpm.gene.log2.m.rfd.ttr <- rbind(subset(tpm.gene.log2.m.rfd, TSS_RFD >= rfd), subset(tpm.gene.log2.m.rfd, TSS_RFD <= -rfd))
 
   save(tpm.gene.log2.m.rfd, tpm.gene.log2.m.rfd.ctr, tpm.gene.log2.m.rfd.ctr.iz, tpm.gene.log2.m.rfd.ctr.iz.cd, tpm.gene.log2.m.rfd.ctr.iz.ho, tpm.gene.log2.m.rfd.ctr.tz, tpm.gene.log2.m.rfd.ctr.tz.cd, tpm.gene.log2.m.rfd.ctr.tz.ho, tpm.gene.log2.m.rfd.ttr, file=file.name)
}

plotBox0 <- function(wd.de.plots, file.name, tpm.1, tpm.2, main, names, cols, ylim) {
   trait <- rep(0, nrow(tpm.1))
   trait <- c(trait, rep(1, nrow(tpm.2)))
   trait <- as.factor(trait)
 
   expr <- as.numeric(c(tpm.1$MEDIAN, tpm.2$MEDIAN))
 
   pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=5, width=3.2)
   boxplot(expr ~ trait, outline=T, xaxt="n", ylab=paste0("log2(TPM + 0.01)"), main=main, boxcol=cols, whiskcol=cols, outcol=cols, medcol=cols, staplecol=cols, ylim=ylim, cex.axis=1.1, cex.lab=1.2, cex.main=1.25)
   
   p <- testU(tpm.1$MEDIAN, tpm.2$MEDIAN)
   text <- ""
   if (p < 1E-5)  text <- "*"
   if (p < 1E-9)  text <- "**"
   if (p < 1E-13) text <- "***"
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
   boxplot(expr ~ trait, outline=T, names=names, ylab=paste0("log2(TPM + 0.01)"), main=main, ylim=ylim, cex.axis=1.1, cex.lab=1.2, cex.main=1.25)
 
   mtext(paste0("p-value = ", scientific(wilcox.test(expr ~ trait, exact=F)$p.value)), cex=1.25, line=0.3)
   dev.off()
}

plotBox <- function(wd.de.plots, file.name, tpm.1, tpm.2, main, names, cols, ylim) {
   trait <- rep(0, nrow(tpm.1))
   trait <- c(trait, rep(1, nrow(tpm.2)))
   trait <- as.factor(trait)
 
   expr <- as.numeric(c(tpm.1$MEDIAN, tpm.2$MEDIAN))
 
   pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=6, width=3)
   boxplot(expr ~ trait, outline=T, xaxt="n", ylab=paste0("log2(TPM + 0.01)"), main=main,boxcol=cols, whiskcol=cols, outcol=cols, medcol=cols, staplecol=cols, ylim=ylim, cex.axis=1.1, cex.lab=1.2, cex.main=1.25)
 
   axis(side=1, at=seq(1, 2, by=1), labels=names, font=2, cex.axis=1.2)
   axis(side=1, at=1, labels="n=428", line=1.3, col=NA, cex.axis=1.2)
   axis(side=1, at=2, labels="n=446", line=1.3, col=NA, cex.axis=1.2)
 
   mtext(paste0("p-value = ", scientific(wilcox.test(expr ~ trait, exact=F)$p.value)), cex=1.25, line=0.3)
   dev.off()
}
