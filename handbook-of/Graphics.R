# =============================================================================
# Library      : Transcription-Replication Conflict
# Name         : handbook-of/TranscriptionReplicationConflict.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 27/05/20
# =============================================================================

# -----------------------------------------------------------------------------
# 
# Last Modified: 06/03/22; 29/05/20
# -----------------------------------------------------------------------------
plotCorrelation <- function(file.name, main.text, xlab.text, ylab.text, x, y, pos="bottomright", cols=c("dimgray", "black"), size, ylim=ylim) {
   pdf(paste0(file.name, ".pdf"), height=size, width=size)
   par(mar=c(5.1, 4.7, 4.1, 1.4))  
   plot(y ~ x, ylab=ylab.text, ylim=ylim, xlab=xlab.text, main=main.text, pch=1, cex=2, col=cols[1], cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
 
   lm.fit <- lm(y ~ x)
   abline(lm.fit, lwd=5, col=cols[2])
 
   cor <- cor.test(y, x, method="spearman", exact=F)
   legend(pos, c(paste0("rho = ", round0(cor[[4]], digits=2)), paste0("P = 1.00E-00")), text.col=c(cols[2], "white"), text.font=2, bty="n", cex=1.8)
   legend(pos, expression(bolditalic('P')~"                   "), text.col=cols[2], text.font=2, bty="n", cex=1.8)
   legend(pos, paste0("   = ", scientific(cor[[3]])), text.col=cols[2], text.font=2, bty="n", cex=1.8)
   
   dev.off()
}

plotCorrelation5 <- function(file.name, main.text, xlab.text, ylab.text, x, y, pos="topright", cols=c("dimgray", "black"), size=5) {
	  pdf(paste0(file.name, ".pdf"), height=size, width=size)
	  par(mar=c(5.1, 4.7, 4.1, 1.4))  
	  plot(y ~ x, ylab=ylab.text, xlab=xlab.text, main=main.text, pch=1, cex=2, col=cols[1], cex.axis=1.9, cex.lab=2, cex.main=2.1)
	
	  lm.fit <- lm(y ~ x)
	  abline(lm.fit, lwd=5, col=cols[2])
	
	  cor <- cor.test(y, x, method="spearman", exact=F)
	  legend(pos, c(paste0("rho = ", round0(cor[[4]], digits=2)), paste0("P = 1.00E-00")), text.col=c(cols[2], "white"), text.font=1, bty="n", cex=2)
	  legend(pos, c("", expression(italic('P')~"                   ")), text.col=cols[2], text.font=1, bty="n", cex=2)
	  legend(pos, c("", paste0("   = ", scientific(cor[[3]]))), text.col=cols[2], text.font=1, bty="n", cex=2)
	
  	dev.off()
}

plotBox <- function(file.name, tpm.1, tpm.2, main, names, cols, ylab.txt="", h=6, w=4.5) {
   trait <- rep(0, length(tpm.1))
   trait <- c(trait, rep(1, length(tpm.2)))
   trait <- as.factor(trait)
 
   expr <- as.numeric(c(tpm.1, tpm.2))
   ylim <- c(min(expr), max(expr))
   if (ylab.txt == "")
      ylab.txt <- expression("Log" * ""[2] * "(TPM + 1)")
   
   pdf(paste0(file.name, ".pdf"), height=h, width=w)
   par(mar=c(5.1, 4.7, 4.1, 1.4))
   boxplot(expr ~ trait, outline=T, main=main, col=cols, xlab="", xaxt="n", ylab=ylab.txt, ylim=ylim, cex=2, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
 
   p <- wilcox.test(expr ~ trait, exact=F)$p.value
   offset <- (ylim[2] - ylim[1])/35
   text(1.5, ylim[2] - offset, expression(italic('P')~"                   "), col="black", cex=1.8)
   text(1.5, ylim[2] - offset, paste0("   = ", scientific(p)), col="black", cex=1.8)
   
   #axis(side=1, at=seq(1, 2, by=1), labels=names, font=2, cex.axis=1.8)
   axis(side=1, at=1, labels=names[1], font=2, cex.axis=1.8)
   axis(side=1, at=2, labels=names[2], font=2, cex.axis=1.8)
   axis(side=1, at=1, labels=paste0("n=", separator(length(tpm.1))), line=1.8, col=NA, cex.axis=1.8)
   axis(side=1, at=2, labels=paste0("n=", separator(length(tpm.2))), line=1.8, col=NA, cex.axis=1.8)

   dev.off()
}

# plotBox3 <- function(wd.de.plots, file.name, tpm.1, tpm.2, tpm.3, main, names, cols, ylim, height=6, width=3.5, text="") {
#    trait <- rep(0, nrow(tpm.1))
#    trait <- c(trait, rep(1, nrow(tpm.2)))
#    trait <- c(trait, rep(2, nrow(tpm.3)))
#    trait <- as.factor(trait)
#    expr <- as.numeric(c(tpm.1$MEDIAN, tpm.2$MEDIAN, tpm.3$MEDIAN))
#  
#    pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=height, width=width)
#    par(mar=c(4.5, 4.7, 4.1, 1.4)) 
#    boxplot(expr ~ trait, outline=F, xlab="", xaxt="n", ylab=text.Log2.TPM.1, main=main, col="white", boxcol=cols, whiskcol=cols, outcol=cols, medcol=cols, staplecol=cols, ylim=ylim, yaxt="n", cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
# 
#    p <- testU(tpm.1$MEDIAN, tpm.3$MEDIAN)
#    text(2, ylim[2]-0.2, getPvalueSignificanceLevel(p), cex=2.5)
#    lines(c(1, 3), y=c(ylim[2]-0.8, ylim[2]-0.8), type="l", lwd=2)
#  
#    p <- testU(tpm.1$MEDIAN, tpm.2$MEDIAN)
#    text(1.5, ylim[2]-1.4, getPvalueSignificanceLevel(p), cex=2.5)
#    lines(c(1, 2), y=c(ylim[2]-2, ylim[2]-2), type="l", lwd=2)
#  
#    p <- testU(tpm.2$MEDIAN, tpm.3$MEDIAN)
#    text(2.5, ylim[2]-2.6, getPvalueSignificanceLevel(p), cex=2.5)
#    lines(c(2, 3), y=c(ylim[2]-3.2, ylim[2]-3.2), type="l", lwd=2)
#  
#    #axis(side=2, at=seq(0, 12, by=2), labels=c(0, "", 4, "", 8, "", 12), cex.axis=1.7)
#    axis(side=1, at=seq(1, 3, by=1), labels=names, font=2, cex.axis=1.7)
# 
#    #mtext(paste0(text), cex=1.8, line=0.3)
#    dev.off()
# }

# plotSRC <- function(wd.plots, gene, cn, snr, pch, pos="bottomright", cols, size=6, xlab.text="SCF index") {
#    unit <- (max(snr) - min(snr))/10
#    xlim <- c(min(snr) - unit, max(snr) + unit)
#    unit <- (max(cn) - min(cn))/10
#    ylim <- c(min(cn) - unit, max(cn) + unit)
#  
#    xlab.text2 <- expression(italic('In silico')~'SCF index')
#    #xlab.text <- 
#    ylab.text <- expression("log" * ""[2] * "(TPM + 1)")
#    #ylab.text <- "Expression"
#    id <- subset(ensGene, external_gene_name == gene)$ensembl_gene_id
#    file.name <- file.path(wd.plots, paste0("TPM-vs-", xlab.text, "_", genes[g], ""))
#  
#    pdf(paste0(file.name, ".pdf"), height=size, width=size)
#    par(mar=c(5.1, 4.7, 4.1, 1.4))
#    if (xlab.text == "SCF index")
#       plot(cn ~ snr, ylim=ylim, xlim=xlim, ylab=ylab.text, xlab=xlab.text2, xaxt="n", main=paste0(gene, ""), col="black", pch=pch, cex=2, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
#    else
#       plot(cn ~ snr, ylim=ylim, xlim=xlim, ylab=ylab.text, xlab=xlab.text,           main=paste0(gene, ""), col="black", pch=pch, cex=2, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
#    
#    cor <- cor.test(cn, snr, method="spearman", exact=F)
#    col <- cols[1]
#    if (cor[[4]] < 0)
#       col <- cols[2]
#    legend(pos, c(paste0("rho = ", round0(cor[[4]], digits=2)), paste0("P = 1.00E-00")), text.col=c(col, "white"), text.font=2, bty="n", cex=1.8)
#    legend(pos, expression(bolditalic('P')~"                   "), text.col=col, text.font=2, bty="n", cex=1.8)
#    legend(pos, paste0("   = ", scientific(cor[[3]])), text.col=col, text.font=2, bty="n", cex=1.8)
#  
#    lm.fit <- lm(cn ~ snr)
#    abline(lm.fit, col=col, lwd=7)
#    
#    if (xlab.text == "SCF index")
#       axis(side=1, at=seq(-0.5, 0.5, by=0.5), labels=c(-0.5, 0, 0.5), cex.axis=1.7)
#    dev.off()
# }

plotCNA <- function(wd.plots, gene, cn, snr, pch, pos="bottomright", col, size=6) {
   unit <- (max(snr) - min(snr))/10
   xlim <- c(min(snr) - unit, max(snr) + unit)
   unit <- (max(cn) - min(cn))/10
   ylim <- c(min(cn) - unit, max(cn) + unit)
 
   xlab.text <- "CNA"
   ylab.text <- expression("Log" * ""[2] * "(TPM + 1)")
   #ylab.text <- "Expression"
   id <- subset(ensGene, external_gene_name == gene)$ensembl_gene_id
   file.name <- file.path(wd.plots, paste0("TPM-vs-CNA_", genes[g], ""))
 
   pdf(paste0(file.name, ".pdf"), height=size, width=size)
   par(mar=c(5.1, 4.7, 4.1, 1.4))
   plot(cn ~ snr, ylim=ylim, xlim=xlim, ylab=ylab.text, xlab=xlab.text, main=paste0(gene, ""), col="black", pch=pch, cex=2, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
 
   lm.fit <- lm(cn ~ snr)
   abline(lm.fit, col=col, lwd=7)
 
   cor <- cor.test(cn, snr, method="spearman", exact=F)
   legend(pos, c(paste0("rho = ", round0(cor[[4]], digits=2)), paste0("P = 1.00E-00")), text.col=c(col, "white"), text.font=2, bty="n", cex=1.8)
   legend(pos, expression(bolditalic('P')~"                   "), text.col=col, text.font=2, bty="n", cex=1.8)
   legend(pos, paste0("   = ", scientific(cor[[3]])), text.col=col, text.font=2, bty="n", cex=1.8)
 
   dev.off()
}

plotMYCN <- function(wd.plots, gene, cn, snr, pch, pos="bottomright", col, size=6) {
   unit <- (max(snr) - min(snr))/10
   xlim <- c(min(snr) - unit, max(snr) + unit)
   unit <- (max(cn) - min(cn))/10
   ylim <- c(min(cn) - unit, max(cn) + unit)
 
   xlab.text <- "MYCN expression"
   ylab.text <- expression("Log" * ""[2] * "(TPM + 1)")
   #ylab.text <- "Expression"
   id <- subset(ensGene, external_gene_name == gene)$ensembl_gene_id
   file.name <- file.path(wd.plots, paste0("TPM-vs-MYCN_", genes[g], ""))
 
   pdf(paste0(file.name, ".pdf"), height=size, width=size)
   par(mar=c(5.1, 4.7, 4.1, 1.4))
   plot(cn ~ snr, ylim=ylim, xlim=xlim, ylab=ylab.text, xlab=xlab.text, main=paste0(gene, ""), col="black", pch=pch, cex=2, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
 
   lm.fit <- lm(cn ~ snr)
   abline(lm.fit, col=col, lwd=7)
 
   cor <- cor.test(cn, snr, method="spearman", exact=F)
   legend(pos, c(paste0("rho = ", round0(cor[[4]], digits=2)), paste0("P = 1.00E-00")), text.col=c(col, "white"), text.font=2, bty="n", cex=1.8)
   legend(pos, expression(bolditalic('P')~"                   "), text.col=col, text.font=2, bty="n", cex=1.8)
   legend(pos, paste0("   = ", scientific(cor[[3]])), text.col=col, text.font=2, bty="n", cex=1.8)
 
   dev.off()
}

plotCNASorting <- function(wd.plots, gene, cn, snr, pch, pos="bottomright", col, size=6) {
   unit <- (max(snr) - min(snr))/10
   xlim <- c(min(snr) - unit, max(snr) + unit)
   unit <- (max(cn) - min(cn))/10
   ylim <- c(min(cn) - unit, max(cn) + unit)
 
   xlab.text <- "CNA"
   ylab.text <- "Proliferation index"
   id <- subset(ensGene, external_gene_name == gene)$ensembl_gene_id
   file.name <- file.path(wd.plots, paste0("TPM-vs-CNA_", genes[g], "_S"))
 
   pdf(paste0(file.name, ".pdf"), height=size, width=size)
   par(mar=c(5.1, 4.7, 4.1, 1.4))
   plot(cn ~ snr, ylim=ylim, xlim=xlim, ylab=ylab.text, xlab=xlab.text, main=paste0(gene, " (n=25)"), col="black", pch=pch, cex=2, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
 
   lm.fit <- lm(cn ~ snr)
   abline(lm.fit, col=col, lwd=7)
 
   cor <- cor.test(cn, snr, method="spearman", exact=F)
   legend(pos, c(paste0("rho = ", round0(cor[[4]], digits=2)), paste0("P = 1.00E-00")), text.col=c(col, "white"), text.font=2, bty="n", cex=1.8)
   legend(pos, expression(bolditalic('P')~"                   "), text.col=col, text.font=2, bty="n", cex=1.8)
   legend(pos, paste0("   = ", scientific(cor[[3]])), text.col=col, text.font=2, bty="n", cex=1.8)
 
   dev.off()
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

plotBox <- function(wd.de.plots, file.name, tpm.1, tpm.2, main, names, cols, ylim, height=6, width=3, cols) {
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

plotBox02 <- function(wd.de.plots, file.name, tpm.1, tpm.2, main, names, cols, ylab.text) {
   trait <- rep(0, length(tpm.1))
   trait <- c(trait, rep(1, length(tpm.2)))
   trait <- as.factor(trait)
   
   expr <- as.numeric(c(tpm.1, tpm.2))
   ylim <- c(min(expr), max(expr))
 
   pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=6, width=4.5)
   par(mar=c(5.1, 4.6, 4.1, 1.5))
   boxplot(expr ~ trait, outline=T, names=names, main=main, col=cols, xaxt="n", ylab=ylab.text, xlab="", ylim=ylim, cex=2, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
 
   p <- wilcox.test(expr ~ trait, exact=F)$p.value
   offset <- (ylim[2] - ylim[1])/35
   text(1.5, ylim[2], expression(italic('P')~'= 5.75E-01'), col="black", cex=1.8)
   
   axis(side=1, at=seq(1, 2, by=1), labels=names, font=2, cex.axis=1.9)
   axis(side=1, at=1, labels=paste0("n=", length(tpm.1)), line=1.8, col=NA, cex.axis=1.9)
   axis(side=1, at=2, labels=paste0("n=", length(tpm.2)), line=1.8, col=NA, cex.axis=1.9)

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

plotBox2 <- function(wd.de.plots, file.name, tpm.1, tpm.2, main, names, cols, ylim) {
   trait <- rep(0, length(tpm.1))
   trait <- c(trait, rep(1, length(tpm.2)))
   trait <- as.factor(trait)
 
   expr <- as.numeric(c(tpm.1, tpm.2))
   ylim <- c(min(expr), max(expr))
 
   pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=6, width=3)
   boxplot(expr ~ trait, outline=T, names=names, ylab=ylab, main=main, xaxt="n", ylim=ylim, cex.axis=1.1, cex.lab=1.2, cex.main=1.25)
 
   p <- wilcox.test(expr ~ trait, exact=F)$p.value
   text(1.5, ylim[2], getPvalueSignificanceLevel(p), col="black", cex=2.5)
   
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

plotBox3 <- function(wd.de.plots, file.name, tpm.1, tpm.2, tpm.3, main, names, cols, ylim, height=5, width=3.5, text="") {
   trait <- rep(0, nrow(tpm.1))
   trait <- c(trait, rep(1, nrow(tpm.2)))
   trait <- c(trait, rep(2, nrow(tpm.3)))
   trait <- as.factor(trait)
   expr <- as.numeric(c(tpm.1$MEDIAN, tpm.2$MEDIAN, tpm.3$MEDIAN))
 
   pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=height, width=width)
   boxplot(expr ~ trait, outline=F, xaxt="n", ylab=paste0("log2(TPM + 1)"), main=main, boxcol=cols, whiskcol=cols, outcol=cols, medcol=cols, staplecol=cols, ylim=ylim, yaxt="n", cex.axis=1.25, cex.lab=1.3, cex.main=1.35)
 
   p <- testU(tpm.1$MEDIAN, tpm.3$MEDIAN)
   text(2, ylim[2]-0.2, getPvalueSignificanceLevel(p), cex=2.5)
   lines(c(1, 3), y=c(ylim[2]-0.8, ylim[2]-0.8), type="l", lwd=2)
   
   p <- testU(tpm.1$MEDIAN, tpm.2$MEDIAN)
   text(1.5, ylim[2]-1.4, getPvalueSignificanceLevel(p), cex=2.5)
   lines(c(1, 2), y=c(ylim[2]-2, ylim[2]-2), type="l", lwd=2)
   
   p <- testU(tpm.2$MEDIAN, tpm.3$MEDIAN)
   text(2.5, ylim[2]-2.6, getPvalueSignificanceLevel(p), cex=2.5)
   lines(c(2, 3), y=c(ylim[2]-3.2, ylim[2]-3.2), type="l", lwd=2)
   
   axis(side=2, at=seq(0, 12, by=2), labels=c(0, "", 4, "", 8, "", 12), cex.axis=1.25)
   axis(side=1, at=seq(1, 3, by=1), labels=names, font=2, cex.axis=1.3)
   #axis(side=1, at=2, labels="Total n=30,978", line=1.3, col=NA, cex.axis=1.2)
   #axis(side=1, at=1, labels=paste0("n=", format(nrow(tpm.1), big.mark=",", scientific=F)), line=1.2, col=NA, cex.axis=1.25)
   #axis(side=1, at=3, labels=paste0("n=", format(nrow(tpm.3), big.mark=",", scientific=F)), line=1.2, col=NA, cex.axis=1.25)
   
   #axis(side=1, at=1, labels=paste0("n=", format(nrow(tpm.1), big.mark=",", scientific=F)), line=1.2, col=NA, cex.axis=1.25)
   #axis(side=1, at=2, labels=paste0("n=", format(nrow(tpm.2), big.mark=",", scientific=F)), line=1.2, col=NA, cex.axis=1.25)
 
   #mtext(paste0("p-value = ", scientific(wilcox.test(expr ~ trait, exact=F)$p.value)), cex=1.2, line=0.3)
   mtext(paste0(text), cex=1.3, line=0.3)
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

plotSRC <- function(gene, cn, snr, pch, col, pos, xlab.text="") {
   unit <- (max(snr) - min(snr))/10
   xlim <- c(min(snr) - unit, max(snr) + unit)
   unit <- (max(cn) - min(cn))/10
   ylim <- c(min(cn) - unit, max(cn) + unit)
 
   if (xlab.text == "")
      xlab.text <- expression(italic('In silico')~'sorting [rho]')
   ylab.text <- "log2(TPM + 1)"
   id <- subset(ensGene, external_gene_name == gene)$ensembl_gene_id
   file.name <- file.path(wd.de.plots, paste0("TPM-vs-SORTING_", genes[g], ""))
 
   pdf(paste0(file.name, ".pdf"), height=6, width=6)
   #plot(cn ~ snr, ylim=ylim, xlim=xlim, ylab="", xaxt="n", xlab=xlab.text, main=paste0(gene, " (", id, ")"), col="black", pch=pch, cex=2, cex.axis=1.7, cex.lab=1.9, cex.main=2)
   plot(cn ~ snr, ylim=ylim, xlim=xlim, ylab="", xaxt="n", xlab="", main=gene, col="black", pch=pch, cex=2, cex.axis=1.7, cex.lab=1.9, cex.main=2)
   
   lm.fit <- lm(cn ~ snr)
   abline(lm.fit, col=col, lwd=7)
 
   cor <- cor.test(cn, snr, method="spearman", exact=F)
   legend(pos, c(paste0("rho = ", round0(cor[[4]], digits=2)), paste0("p-value = ", scientific(cor[[3]], digits=2))), text.col=col, text.font=2, bty="n", cex=1.9)
 
   axis(side=1, at=seq(-0.5, 0.5, by=0.5), labels=c(-0.5, 0, 0.5), cex.axis=1.7)
   #axis(side=2, at=seq(6, 8, by=1), labels=c(6, 7, 8), cex.axis=1.7)   ## MARS
   #axis(side=2, at=seq(4, 6, by=1), labels=c(4, 5, 6), cex.axis=1.7)   ## GTPBP3
   mtext(ylab.text, side=2, line=2.74, cex=1.85)
   mtext(xlab.text, side=1, line=3.5, cex=1.9)
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