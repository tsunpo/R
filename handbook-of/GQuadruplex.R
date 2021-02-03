# =============================================================================
# Library      : G-quadruplex
# Name         : handbook-of/GQuadruplex.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 02/02/21
# =============================================================================

# -----------------------------------------------------------------------------
# G-quadruplex region skew (GRS)
# Last Modified: 02/02/21
# -----------------------------------------------------------------------------
getG4RS <- function(wd.meta, samples, bed.gc, nrds) {
   colnames <- c(paste0("chr", 1:22), "E", "L", "G4RS")   ##"ANOVA_P", "ANOVA_Q", 
   g4s <- toTable(0, length(colnames), nrow(samples), colnames)
   rownames(g4s) <- rownames(samples)
 
   for (s in 1:nrow(samples)) {
      sample <- samples$SAMPLE_ID[s]
      bed.g4 <- read.bed(file.path(wd.meta, paste0("_", sample, ".bed")))
      
      colnames <- c("BED", "SPLINE")
      tmp <- toTable(0, length(colnames), 0, colnames)
      for (c in 1:22) {
         chr <- chrs[c]
         bed.g4.chr <- subset(bed.g4, CHR == chr)
         bed.gc.chr <- subset(bed.gc, CHR == chr)
         
         bed.gc.chr$BED <- rownames(bed.gc.chr)
         overlaps <- intersect(nrds$BED, bed.gc.chr$BED)
         nrds.chr <- nrds[overlaps,]
         nrds.chr.RT <- setSpline(nrds.chr, bed.gc.chr, "RT")
         nrds.chr.RT <- cbind(nrds.chr.RT, bed.gc.chr[overlaps, c("START", "END")])
         
         ##
         tmp.chr <- toTable(0, length(colnames), 0, colnames)
         for (g in 1:nrow(bed.g4.chr)) {
            g4 <- bed.g4.chr[g,]
            
            nrds.chr.RT.end <- nrds.chr.RT[which(nrds.chr.RT$END >= g4$START),]   ## Important! Not >= seg$START if e.g. P2 chr1 10000 11000 0.646000
            nrds.chr.RT.end.start <- nrds.chr.RT.end[which(nrds.chr.RT.end$START <= g4$END),]
            tmp.chr <- rbind(tmp.chr, nrds.chr.RT.end.start[, c("BED", "SPLINE")])
         }
         
         e <- nrow(subset(tmp.chr, SPLINE > 0))
         l <- nrow(subset(tmp.chr, SPLINE < 0))
         g4s[s, c] <- (e - l)/(e + l)
         
         tmp <- rbind(tmp, tmp.chr)  
      }
      
      e <- nrow(subset(tmp, SPLINE > 0))
      l <- nrow(subset(tmp, SPLINE < 0))
      g4s[s, "E"] <- e
      g4s[s, "L"] <- l
      g4s[s, "G4RS"] <- (e - l)/(e + l)
   }
 
   return(g4s)
}

plotG4RS <- function(snr3, n3, file.name, main.text, xlab.text, ylab.text, col, col2, pos) {
   ylim <- c(-1.1, 1.1)
 
   pdf(paste0(file.name, ".pdf"), height=5, width=5)
   plot(snr3 ~ n3, ylim=ylim, ylab=ylab.text, xlab=xlab.text, main=main.text[1], col=col2, pch=15, cex=1.5, cex.axis=1.2, cex.lab=1.25, cex.main=1.35)
   abline(h=snr3[2], lty=5, lwd=1.3)
   
   lm.fit <- lm(snr3 ~ n3)
   abline(lm.fit, col=col, lwd=4)
 
   cor3 <- cor.test(snr3, n3, method="spearman", exact=F)
   legend("bottomright", c(paste0("rho = ", round0(cor3[[4]], digits=2)), paste0("p-value = ", scientific(cor3[[3]]))), text.col=cols, text.font=2, bty="n", cex=1.2)

   mtext(main.text[2], cex=1.25, line=0.3)
   dev.off()
}

plotBoxG4R <- function(wd.de.plots, file.name, tpm.1, tpm.2, main, names, cols, ylim, height=6, width=4.2) {
   trait <- rep(0, length(tpm.1))
   trait <- c(trait, rep(1, length(tpm.2)))
   trait <- as.factor(trait)
   expr <- as.numeric(c(tpm.1, tpm.2))
 
   pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=height, width=width)
   boxplot(expr ~ trait, outline=T, xaxt="n", ylab=paste0("Number of G4R [#]"), main=main, boxcol=cols, whiskcol=cols, outcol=cols, medcol=cols, staplecol=cols, ylim=ylim, cex.axis=1.2, cex.lab=1.25, cex.main=1.3)
 
   p <- testU(tpm.1, tpm.2)
   text(1.5, ylim[2], getPvalueSignificanceLevel(p), col="black", cex=2.5)
 
   ##
   axis(side=1, at=seq(1, 2, by=1), labels=names, font=2, cex.axis=1.2)
   #axis(side=1, at=2, labels="Total n=30,978", line=1.3, col=NA, cex.axis=1.2)
   axis(side=1, at=1, labels=paste0("n=", format(length(tpm.1), big.mark=",", scientific=F)), line=1.2, col=NA, cex.axis=1.25)
   axis(side=1, at=2, labels=paste0("n=", format(length(tpm.2), big.mark=",", scientific=F)), line=1.2, col=NA, cex.axis=1.25)
 
   mtext(paste0("p-value = ", scientific(wilcox.test(expr ~ trait, exact=F)$p.value)), cex=1.2, line=0.3)
   dev.off()
}

# =============================================================================
# Inner Class  : BED File Reader
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 02/02/21
# =============================================================================
read.bed <- function(bed.file) {
   bed <- readTable(bed.file, header=F, rownames=F, sep="")
   colnames(bed) <- c("CHR", "START", "END")
 
   return(bed)
}