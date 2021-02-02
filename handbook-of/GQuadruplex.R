# =============================================================================
# Library      : G-quadruplex
# Name         : handbook-of/GQuadruplex.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 02/02/21
# =============================================================================

# -----------------------------------------------------------------------------
# Methods: 
# Last Modified: //
# -----------------------------------------------------------------------------
getBEDFromSegment <- function(bed.gc.o, seg.gc) {
   #bed.gc.chr <- subset(bed.gc, CHR == seg.gc$CHR)                       ## Using subset() takes ~15 MINS on 1,040 segments over 2,861,558 entries (S00022)
   #bed.gc.chr.end <- subset(bed.gc.chr, END > seg.gc$START)              ## Important! Not >= seg$START if e.g. P2 chr1 10000 11000 0.646000
   #bed.gc.chr.end.start <- subset(bed.gc.chr.end, START <= seg.gc$END)
   bed.gc.chr <- bed.gc.o[which(bed.gc.o$CHR == seg.gc$CHR),]                 ## Using which() takes ~5 MINS on 1,040 segments over 2,861,558 entries (S00022)
   bed.gc.chr.end <- bed.gc.chr[which(bed.gc.chr$END >= seg.gc$START),]   ## Important! Not >= seg$START if e.g. P2 chr1 10000 11000 0.646000
   bed.gc.chr.end.start <- bed.gc.chr.end[which(bed.gc.chr.end$START <= seg.gc$END),]
 
   return(bed.gc.chr.end.start)
}

# -----------------------------------------------------------------------------
# G-quadruplex region skew (GRS)
# Last Modified: 02/02/21
# -----------------------------------------------------------------------------
getG4RS <- function(wd.meta, samples, nrds, bed.gc) {
   colnames <- c(paste0("chr", 1:22), "G4RS")   ##"ANOVA_P", "ANOVA_Q", 
   g4s <- toTable(0, length(colnames), nrow(samples), colnames)
   rownames(g4s) <- rownames(samples)
 
   for (s in 2:nrow(samples)) {
      sample <- samples$SAMPLE_ID[s]
      bed.g4 <- read.bed(file.path(wd.meta, paste0("_", sample, ".bed")))
      
      colnames <- c("BED", "SPLINE")
      tmp <- toTable(0, length(colnames), 0, colnames)
      for (c in 1:22) {
         chr <- chrs[c]
         bed.g4.chr <- subset(bed.g4, CHR == chr)
         bed.gc.chr <- subset(bed.gc, CHR == chr)
         
         bed.gc.chr$BED <- rownames(bed.gc.chr)
         overlaps <- intersect(nrds.lcl$BED, bed.gc.chr$BED)
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
      g4s[s, "G4RS"] <- (e - l)/(e + l)
   }
 
   return(g4s)
}

# -----------------------------------------------------------------------------
# Methods: 
# Last Modified: //
# -----------------------------------------------------------------------------
plotRTS <- function(sprs, file.name, main.text, cs=NULL, digits, unit, ylab.text, cex, chr2, offset) {
   xlab.text <- "Chromosome"
   cols <- c(red, blue, "black")
   #ylim <- getYlim(sprs$spr, unit)
   ylim <- c(-1.1, 1.1)
    
   pdf(paste0(file.name, ".pdf"), height=5, width=5)
   #plot(sprs$skew ~ sprs$chr, ylim=c(ymin, ymax), ylab=ylab.text, xlab=xlab.text, main=main.text, col=cols[3], xaxt="n", pch=19)   ## yaxt="n",
   plot(NULL, xlim=c(1, 22), ylim=ylim, xlab=xlab.text, ylab=ylab.text, main=main.text[1], col=cols[3], xaxt="n", yaxt="n", pch=19, cex.axis=1.2, cex.lab=1.25, cex.main=1.35)
   
   abline(h=sprs$spr[2], lty=5, lwd=1.3)
   lines(sprs$spr, y=NULL, type="l", lwd=2, col=cols[3])
   
   idx <- which(sprs$spr > sprs$spr[2])
   points(sprs$spr[idx] ~ sprs$chr[idx], col=cols[1], pch=19, cex=cex)
   idx <- which(sprs$spr < sprs$spr[2])
   points(sprs$spr[idx] ~ sprs$chr[idx], col=cols[2], pch=19, cex=cex)
   points(sprs$spr[2] ~ sprs$chr[2], col=cols[3], pch=19, cex=cex)
 
   text(sprs$chr[2], sprs$spr[2], paste0(offset, "Chr2 (", chr2, ")"), col=cols[3], pos=3, cex=1.2)
   #text(sprs$chr[2]+1.8, sprs$spr[2], paste0("Chr2 (", round0(sprs$spr[2], digits=digits), ")"), cex=1.1, col=cols[3], pos=3)
   if (!is.null(cs))
      for (c in 1:length(cs)) {
         c <- cs[c]
         if (sprs$spr[c] > sprs$spr[2])
            text(sprs$chr[c], sprs$spr[c], paste0("Chr", c), col="red", pos=3, cex=1.2)
            #text(sprs$chr[c]+1.8, sprs$spr[c], paste0("Chr", c, " (", round0(sprs$spr[c], digits=digits), ")"), cex=1.1, col=cols[1], pos=3)
         else
            text(sprs$chr[c], sprs$spr[c], paste0("Chr", c), col="blue", pos=1, cex=1.2)
            #text(sprs$chr[c]+1.8, sprs$spr[c], paste0("Chr", c, " (", round0(sprs$spr[c], digits=digits), ")"), cex=1.1, col=cols[2], pos=1)
      }
   legend("topleft", "Earlier than chr2", text.col=cols[1], pch=19, pt.cex=1.5, col=cols[1], cex=1.25)   
   legend("bottomleft", "Later than", text.col=cols[2], pch=19, pt.cex=1.5, col=cols[2], cex=1.25)

   axis(side=1, at=seq(2, 22, by=4), cex.axis=1.2)
   axis(side=1, at=seq(4, 20, by=4), cex.axis=1.2)
   axis(side=2, at=seq(-1, 1, by=0.5), labels=c(-1, -0.5, 0, 0.5, 1), cex.axis=1.2)
   mtext(main.text[2], line=0.3, cex=1.25)
   dev.off()
}

plotRTS2 <- function(sprs, means, file.name, main.text, cs, xlab.text, unit, ylab.text, cex) {
   cols <- c(red, blue, "black", green)
   ylim <- c(-1.1, 1.1)
   
   unit <- (max(means) - min(means))/15
   xlim <- c(min(means) - unit, max(means) + unit)
   
   pdf(paste0(file.name, ".pdf"), height=5, width=5)
   if (is.null(xlim))
      plot(sprs ~ means, ylim=ylim, ylab=ylab.text, xlab=xlab.text, main=main.text[1], col="white", yaxt="n", xaxt="n", pch=19, cex.axis=1.2, cex.lab=1.25, cex.main=1.35)
   else
      plot(sprs ~ means, ylim=ylim, xlim=xlim, ylab=ylab.text, xlab=xlab.text, main=main.text[1], col="white", yaxt="n", xaxt="n", pch=19, cex.axis=1.2, cex.lab=1.25, cex.main=1.35)
   
   abline(h=sprs[2], lty=5, lwd=1.3)
   lm.fit <- lm(sprs ~ means)
   abline(lm.fit, col=cols[4], lwd=3)
   
   idx <- which(sprs > sprs[2])
   points(sprs[idx] ~ means[idx], col=cols[1], pch=19, cex=cex)
   idx <- which(sprs < sprs[2])
   points(sprs[idx] ~ means[idx], col=cols[2], pch=19, cex=cex)
   points(sprs[2] ~ means[2], col=cols[3], pch=19, cex=cex)
 
   text(means[2], sprs[2], paste0("Chr", 2), col=cols[3], pos=3, cex=1.2)
   if (!is.null(cs))
      for (c in 1:length(cs)) {
         c <- cs[c]
         if (sprs[c] > sprs[2])
            #if (sprs[c] == sprs[22])
            #   text(means[c], sprs[c], paste0("Chr", c), col=cols[1], pos=1, cex=1.2)
            #else
               text(means[c], sprs[c], paste0("Chr", c), col="red", pos=3, cex=1.2)
         else
            text(means[c], sprs[c], paste0("Chr", c), col="blue", pos=1, cex=1.2)
      }
   
   cor <- cor.test(sprs, means, method="spearman", exact=F)
   legends <- c("topright", "bottomleft")
   if (cor[[4]] > 0) legends[1] <- "topleft"
   legend(legends[1], "Earlier than chr2", text.col=cols[1], pch=19, pt.cex=1.5, col=cols[1], cex=1.25)   ## bty="n"
   legend(legends[2], "Later than", text.col=cols[2], pch=19, pt.cex=1.5, col=cols[2], cex=1.25)
   
   #legend("bottomright", paste0("rho = ", round0(cor[[4]], digits=2)), text.col=cols[4], bty="n", cex=1.2)
   legend("bottomright", c(paste0("rho = ", round0(cor[[4]], digits=2)), paste0("p-value = ", scientific(cor[[3]]))), text.col=cols[4], text.font=2, bty="n", cex=1.25)
   #legend("bottomright", c(paste0("R^2 = ", round0(summary(lm.fit)$r.squared, digits=2)), paste0("p-value = ", scientific(summary(lm.fit)$coefficients[2, 4]))), text.col=cols[4], bty="n", cex=1.1)
   #axis(side=1, at=seq(2, 22, by=2))
   axis(side=1, at=seq(-0.4, 0.8, by=0.4), labels=c(-0.4, 0, 0.4, 0.8), cex.axis=1.2)
   axis(side=1, at=seq(-0.2, 0.6, by=0.4), labels=c(-0.2, 0.2, 0.6), cex.axis=1.2)
   axis(side=2, at=seq(-1, 1, by=0.5), labels=c(-1, -0.5, 0, 0.5, 1), cex.axis=1.2)
   mtext(main.text[2], cex=1.25, line=0.3)
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