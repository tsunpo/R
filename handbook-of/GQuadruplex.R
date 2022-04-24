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

plotG4RS <- function(rts, g4s, file.name, main.text, xlab.text, ylab.text, cols, pos) {
   ylim <- c(-1.1, 1.1)
   
   xlim <- c(0, 0)
   na <- which(is.na(g4s))
   if (length(na) != 0) {
      xlim <- c(min(g4s[-na]), max(g4s[-na]))
   } else {
      xlim <- c(min(g4s), max(g4s))
   }
   offset <- (xlim[2] - xlim[1])/20
   xlim <- c(xlim[1]-offset, xlim[2]+offset)
   
   pdf(paste0(file.name, ".pdf"), height=6, width=6)
   par(mar=c(5.1, 4.6, 4.1, 1.5))
   plot(rts ~ g4s, ylim=ylim, xlim=xlim, ylab=ylab.text, xlab=xlab.text, main=main.text[1], yaxt="n", xaxt="n", pch=19, col="white", cex=2, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
   abline(h=rts[2], lty=5, lwd=1.3)

   idx <- which(rts > rts[2])
   points(rts[idx] ~ g4s[idx], col=cols[1], pch=19, cex=2)
   idx <- which(rts < rts[2])
   points(rts[idx] ~ g4s[idx], col=cols[2], pch=19, cex=2)
   points(rts[2] ~ g4s[2], col=cols[3], pch=19, cex=2)
   if (length(na) != 0)
      points(rts[na] ~ g4s[na], col=cols[3], pch=1, cex=2)
   
   lm.fit <- lm(rts ~ g4s)
   abline(lm.fit, col=cols[4], lwd=5)
   
   cor3 <- cor.test(rts, g4s, method="spearman", exact=F)
   legend("bottomright", c(paste0("rho = ", round0(cor3[[4]], digits=2)), paste0("P = 1.00E-00")), text.col=c(cols[4], "white"), text.font=2, bty="n", cex=1.8)
   legend("bottomright", expression(bolditalic('P')~"                   "), text.col=cols[4], text.font=2, bty="n", cex=1.8)
   legend("bottomright", paste0("   = ", scientific(cor3[[3]])), text.col=cols[4], text.font=2, bty="n", cex=1.8)
   
   #text(max/4, 0.07, expression(italic('P')~"="~"                "), cex=1.7)
   #text(max/4, 0.07, paste0("      ", scientific(surv_pvalue(fit)$pval, digits=2)), cex=1.7)
   
   text(g4s[2], rts[2], "Chr2", col=cols[3], pos=3, cex=1.7)
   text(g4s[17], rts[17], "Chr17", col=cols[1], pos=3, cex=1.7)
   text(g4s[19], rts[19], "Chr19", col=cols[1], pos=3, cex=1.7)
   text(g4s[4], rts[4], "Chr4", col=cols[2], pos=1, cex=1.7)
   text(g4s[13], rts[13], "Chr13", col=cols[2], pos=1, cex=1.7)
   if (length(na) != 0)
      text(g4s[na], rts[na], paste0("Chr", na), col=cols[3], pos=1, cex=1.7) 
    
   legend("topleft", "Earlier than chr2", text.col=cols[1], pch=19, pt.cex=2.5, col=cols[1], cex=1.7)   
   legend("bottomleft", "Later than chr2", text.col=cols[2], pch=19, pt.cex=2.5, col=cols[2], cex=1.7)

   axis(side=2, at=seq(-1, 1, by=0.5), labels=c(-1, -0.5, 0, 0.5, 1), cex.axis=1.7)
   axis(side=1, at=seq(-0.4, 0.8, by=0.4), labels=c(-0.4, 0, 0.4, 0.8), cex.axis=1.7)
   axis(side=1, at=seq(-0.2, 1, by=0.4), labels=c(-0.2, 0.2, 0.6, 1), cex.axis=1.7)
   #mtext(main.text[2], cex=1.25, line=0.29)
   dev.off()
}

plotBoxG4R <- function(wd.de.plots, file.name, tpm.1, tpm.2, main.txt, names, cols, ylim, height=6, width=4.2) {
   trait <- rep(0, length(tpm.1))
   trait <- c(trait, rep(1, length(tpm.2)))
   trait <- as.factor(trait)
   expr <- as.numeric(c(tpm.1, tpm.2))
 
   pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=height, width=width)
   boxplot(expr ~ trait, outline=T, xaxt="n", ylab="", main=main.txt[1], col=cols, ylim=ylim, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
 
   p <- testU(tpm.1, tpm.2)
   #text(1.5, ylim[2]-1250, paste0("p-value = ", scientific(wilcox.test(expr ~ trait, exact=F)$p.value)), col="black", cex=1.5)
   #text(1.5, ylim[2], getPvalueSignificanceLevel(p), col="black", cex=3)
   text(1.5, ylim[2]-1250, paste0("p-value = ", scientific(wilcox.test(expr ~ trait, exact=F)$p.value)), col="black", cex=1.5)
   text(1.5, ylim[2], "*", col="black", cex=3)
   
   ##
   axis(side=1, at=seq(1, 2, by=1), labels=names, font=2, cex.axis=1.8)
   axis(side=1, at=1, labels=paste0("n=", format(sum(tpm.1), big.mark=",", scientific=F)), line=1.5, col=NA, cex.axis=1.6)
   axis(side=1, at=2, labels=paste0("n=", format(sum(tpm.2), big.mark=",", scientific=F)), line=1.5, col=NA, cex.axis=1.6)
 
   mtext("Number of G4R", side=2, line=2.9, cex=1.8)   
   mtext("", line=0.3, cex=1.25)
   dev.off()
}

plotG4RSvsIS <- function(n3, snr3, file.name, main.text, xlab.text, ylab.text, col, col2, pos) {
   #xlim <- c(0, xlim.max)
   ylim <- c(-0.35, 0.24)
 
   pdf(paste0(file.name, ".pdf"), height=6, width=6)
   plot(n3 ~ snr3, ylim=ylim, ylab="", xlab=xlab.text, main="", yaxt="n", xaxt="n", col=col2, pch=15, cex=2, lwd=0, cex.axis=1.5, cex.lab=1.6, cex.main=1.7)
   lm.fit <- lm(n3 ~ snr3)
   abline(lm.fit, col=col, lwd=4)
 
   cor3 <- cor.test(n3, snr3, method="spearman", exact=F)
   legend("bottomright", c(paste0("rho = ", round0(cor3[[4]], digits=2)), paste0("p-value = ", scientific(cor3[[3]]))), text.col=cols, text.font=2, bty="n", cex=1.5)
 
   axis(side=1, at=seq(-0.2, 1, by=0.2), labels=c(-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1), cex.axis=1.5)
   axis(side=2, at=seq(-0.2, 0.2, by=0.2), labels=c(-0.2, 0, 0.2), cex.axis=1.5)
   axis(side=2, at=seq(-0.3, 0.1, by=0.2), labels=c("", "", ""), cex.axis=1.5)
   mtext(ylab.text, side=2, line=2.75, cex=1.6)
 
   mtext(main.text[1], line=1.8, cex=1.7, font=2)
   mtext(main.text[2], line=0.33, cex=1.6)
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