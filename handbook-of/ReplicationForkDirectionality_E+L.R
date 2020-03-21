# =============================================================================
# Library      : Replication Fork Directionality
# Name         : handbook-of/ReplicationForkDirectionality.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 21/10/19; 13/09/19; 05/04/19; 21/02/19
# =============================================================================

# -----------------------------------------------------------------------------
# Read in bootstrap data (in 3b)
# Last Modified: 20/09/19; 31/10/18
# -----------------------------------------------------------------------------
getBSTRPS <- function(nrds.RT, colname, b) {
   nrds.RT <- nrds.RT[, c("BED", colname)]
   colnames(nrds.RT)[2] <- paste0(colname, "_", b)
 
   return(nrds.RT)
}

# -----------------------------------------------------------------------------
# Determine RFD from bootstrap data (in 3b)
# Last Modified: 20/09/19; 31/10/18
# -----------------------------------------------------------------------------
getPOS <- function(nrds.RT.BSTRPS) {   ## Left leadings are with a positive slope
   return(as.numeric(length(which(nrds.RT.BSTRPS > 0))))
}

getNEG <- function(nrds.RT.BSTRPS) {   ## Right leadings are with a negative slope
   return(as.numeric(length(which(nrds.RT.BSTRPS < 0))))
}

getRFD <- function(nrds.RT.BSTRPS) {   ## (R-L)/(R+L)
   return((nrds.RT.BSTRPS$NEG - nrds.RT.BSTRPS$POS)/(nrds.RT.BSTRPS$NEG + nrds.RT.BSTRPS$POS))
}

pipeBootstrap <- function(nrds.RT.BSTRPS, bstrps) {
   nrds.RT.BSTRPS$POS <- mapply(x = 1:nrow(nrds.RT.BSTRPS), function(x) as.numeric(getPOS(nrds.RT.BSTRPS[x, 1:bstrps])))   ## BUG FIX: 01/11/18
   nrds.RT.BSTRPS$NEG <- mapply(x = 1:nrow(nrds.RT.BSTRPS), function(x) as.numeric(getNEG(nrds.RT.BSTRPS[x, 1:bstrps])))   ## BUG FIX: 01/11/18

   return(nrds.RT.BSTRPS)
}

getBootstrap <- function(base, column) {
   load(file=file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.d.rt.RT.", column, "_", "chr1", ".RData")))
   nrds.RT.BSTRPS <- nrds.RT.BSTRPS.chr
   for (c in 2:22) {
      chr <- chrs[c]

      load(file=file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.d.rt.RT.", column, "_", chr, ".RData")))
      nrds.RT.BSTRPS <- rbind(nrds.RT.BSTRPS, nrds.RT.BSTRPS.chr)
   }
   
   return(nrds.RT.BSTRPS)
}

getRTNRFD <- function(nrds, nrds.RT.BSTRPS, bed.gc, kb) {
   nrds.RT <- NULL
   for (c in 1:22) {
      chr <- chrs[c]
      bed.gc.chr <- subset(bed.gc, CHR == chr)
  
      nrds.RT.chr <- setSpline(nrds, bed.gc.chr, "RT", returnAll=T)
      if (is.null(nrds.RT)) {
         nrds.RT <- nrds.RT.chr
      } else {
         nrds.RT <- rbind(nrds.RT, nrds.RT.chr)
      }
   }
   
   ##
   nrds.RT.BSTRPS$RFD <- NA
   nrds.RT.BSTRPS$RFD <- getRFD(nrds.RT.BSTRPS)
   colnames(nrds.RT.BSTRPS) <- c("L", "R", "RFD")

   overlaps <- intersect(rownames(nrds.RT), rownames(nrds.RT.BSTRPS))
   nrds.RT <- cbind(nrds.RT[overlaps,], nrds.RT.BSTRPS[overlaps,])

   nrds.RT.NRFD <- NULL
   for (c in 1:22) {
      chr <- chrs[c]
      bed.gc.chr <- subset(bed.gc, CHR == chr)
      
      overlaps <- intersect(rownames(bed.gc.chr), nrds.RT$BED)
      bed.gc.chr.o <- bed.gc.chr[overlaps,]
      nrds.chr.o <- nrds.RT[overlaps,]
   
      ## NRFD
      ## https://stackoverflow.com/questions/41061140/how-to-calculate-the-average-slope-within-a-moving-window-in-r
      ## https://www.rdocumentation.org/packages/zoo/versions/1.8-6/topics/rollapply
      rollingSlope.lm.fit <- function(vector) {
         return(coef(.lm.fit(cbind(1, seq(vector)), vector))[2])
      }
      Slope.lm.fit = rollapply(nrds.chr.o$RFD, width=kb, FUN=rollingSlope.lm.fit, fill=NA, partial=T)
      nrds.chr.o$NRFD <- Slope.lm.fit
      #nrow(subset(subset(subset(nrds.chr.o, RFD < 0.9), RFD > -0.9), NRFD == 0))
      
      if (is.null(nrds.RT.NRFD)) {
         nrds.RT.NRFD <- nrds.chr.o
      } else {
         nrds.RT.NRFD <- rbind(nrds.RT.NRFD, nrds.chr.o)
      }
   }
   
   return(nrds.RT.NRFD)
}

# -----------------------------------------------------------------------------
# Determine TTR and CTR
# Last Modified: 20/09/19; 31/10/18
# -----------------------------------------------------------------------------
getBootstrapCTR <- function(nrds.RFD, rfd) {
   return(subset(subset(nrds.RFD, RFD < rfd), RFD > -rfd))
}

getBootstrapTTR <- function(nrds.RFD, rfd) {
   nrds.RFD.ctr <- getBootstrapCTR(nrds.RFD, rfd)
 
   diff <- setdiff(rownames(nrds.RFD), rownames(nrds.RFD.ctr))
   return(nrds.RFD[diff,])
}

#getBootstrapCTR <- function(nrds.RFD, boundary.lower, boundary.upper) {
#   return(subset(subset(nrds.RFD,R > boundary.lower), R < boundary.upper))
#}
#
#getBootstrapTTR <- function(nrds.RFD, boundary.lower, boundary.upper) {
#   nrds.RFD.ctr <-  <- getBootstrapCTR(nrds.RFD, boundary.lower, boundary.upper) 
# 
#   diff <- setdiff(rownames(nrds.RFD), rownames(nrds.RFD.ctr <- ))
#   return(nrds.RFD[diff,])
#}

getBootstrapReport <- function(rfd, nrds.RT.RFD.1, nrds.RT.RFD.2, name.1, name.2) {
   colnames <- c("RFD", "SAMPLES", name.1, name.2, "Overlapping_N", "Overlapping_P", "RFD_N", "RFD_P", "Mappable")
   report <- toTable(NA, length(colnames), 6, colnames)
   
   ###
   ## |RFD| ≥ 0.9
   nrds.RT.RFD.1.t <- getBootstrapTTR(nrds.RT.RFD.1, rfd)
   nrds.RT.RFD.2.t <- getBootstrapTTR(nrds.RT.RFD.2, rfd)
   
   ## TTR
   report$RFD[1:2] <- "TTR"
   report$SAMPLES[1] <- name.1
   report$SAMPLES[2] <- name.2
   report$RFD_N[1] <- nrow(nrds.RT.RFD.1.t)
   report$RFD_N[2] <- nrow(nrds.RT.RFD.2.t)
   report$Mappable[1] <- nrow(nrds.RT.RFD.1)
   report$Mappable[2] <- nrow(nrds.RT.RFD.2)
   report$RFD_P[1] <- report$RFD_N[1]/report$Mappable[1]
   report$RFD_P[2] <- report$RFD_N[2]/report$Mappable[2]
   
   ##
   overlaps.t <- intersect(rownames(nrds.RT.RFD.1.t), rownames(nrds.RT.RFD.2.t))
   report$Overlapping_N[1] <- length(overlaps.t)
   report$Overlapping_N[2] <- length(overlaps.t)   
   report$Overlapping_P[1] <- length(overlaps.t)/nrow(nrds.RT.RFD.1)
   report$Overlapping_P[2] <- length(overlaps.t)/nrow(nrds.RT.RFD.2)
   
   nrds.1.RT.o <- nrds.RT.RFD.1.t[overlaps.t,]
   nrds.2.RT.o <- nrds.RT.RFD.2.t[overlaps.t,]

   nrds.1.RT.o$SIGN <- nrds.1.RT.o$SLOPE * nrds.2.RT.o$SLOPE
   report[2, 3] <- length(which(nrds.1.RT.o$SIGN > 0))/length(overlaps.t)
   
   ###
   ## |RFD| < 0.9
   nrds.RT.RFD.1.c <- getBootstrapCTR(nrds.RT.RFD.1, rfd)
   nrds.RT.RFD.2.c <- getBootstrapCTR(nrds.RT.RFD.2, rfd)
   
   nrds.RT.RFD.1.c.e <- subset(nrds.RT.RFD.1.c, NRFD > 0)
   nrds.RT.RFD.1.c.l <- subset(nrds.RT.RFD.1.c, NRFD < 0)
   nrds.RT.RFD.2.c.e <- subset(nrds.RT.RFD.2.c, NRFD > 0)
   nrds.RT.RFD.2.c.l <- subset(nrds.RT.RFD.2.c, NRFD < 0)
   
   ## CTR (IZ)
   report$RFD[3:4] <- "CTR-IZ"
   report$SAMPLES[3] <- name.1
   report$SAMPLES[4] <- name.2
   
   report$RFD_N[3] <- nrow(nrds.RT.RFD.1.c.e)
   report$RFD_N[4] <- nrow(nrds.RT.RFD.2.c.e)
   report$Mappable[3] <- nrow(nrds.RT.RFD.1)
   report$Mappable[4] <- nrow(nrds.RT.RFD.2)
   report$RFD_P[3] <- report$RFD_N[3]/report$Mappable[3]
   report$RFD_P[4] <- report$RFD_N[4]/report$Mappable[4]
   
   ##
   overlaps.c.e <- intersect(rownames(nrds.RT.RFD.1.c.e), rownames(nrds.RT.RFD.2.c.e))
   report$Overlapping_N[3] <- length(overlaps.c.e)
   report$Overlapping_N[4] <- length(overlaps.c.e)   
   report$Overlapping_P[3] <- length(overlaps.c.e)/nrow(nrds.RT.RFD.1)
   report$Overlapping_P[4] <- length(overlaps.c.e)/nrow(nrds.RT.RFD.2)
   
   nrds.1.RT.o <- nrds.RT.RFD.1.c.e[overlaps.c.e,]
   nrds.2.RT.o <- nrds.RT.RFD.2.c.e[overlaps.c.e,]

   nrds.1.RT.o$SIGN <- nrds.1.RT.o$SLOPE * nrds.2.RT.o$SLOPE
   report[4, 3] <- length(which(nrds.1.RT.o$SIGN > 0))/length(overlaps.c.e)

   ## CTR (TZ)
   report$RFD[5:6] <- "CTR-TZ"
   report$SAMPLES[5] <- name.1
   report$SAMPLES[6] <- name.2
   
   report$RFD_N[5] <- nrow(nrds.RT.RFD.1.c.l)
   report$RFD_N[6] <- nrow(nrds.RT.RFD.2.c.l)
   report$Mappable[5] <- nrow(nrds.RT.RFD.1)
   report$Mappable[6] <- nrow(nrds.RT.RFD.2)
   report$RFD_P[5] <- report$RFD_N[5]/report$Mappable[5]
   report$RFD_P[6] <- report$RFD_N[6]/report$Mappable[6]
   
   ##
   overlaps.c.l <- intersect(rownames(nrds.RT.RFD.1.c.l), rownames(nrds.RT.RFD.2.c.l))
   report$Overlapping_N[5] <- length(overlaps.c.l)
   report$Overlapping_N[6] <- length(overlaps.c.l)   
   report$Overlapping_P[5] <- length(overlaps.c.l)/nrow(nrds.RT.RFD.1)
   report$Overlapping_P[6] <- length(overlaps.c.l)/nrow(nrds.RT.RFD.2)
   
   nrds.1.RT.o <- nrds.RT.RFD.1.c.l[overlaps.c.l,]
   nrds.2.RT.o <- nrds.RT.RFD.2.c.l[overlaps.c.l,]
   
   nrds.1.RT.o$SIGN <- nrds.1.RT.o$SLOPE * nrds.2.RT.o$SLOPE
   report[6, 3] <- length(which(nrds.1.RT.o$SIGN > 0))/length(overlaps.c.l)
   
   return(report)
}

getReportRFD <- function(report, name) {
   return(subset(report, SAMPLES == name)$RFD_P)
}

getReportRFD12 <- function(report, name) {
   return(subset(report, SAMPLES == name)$Overlapping_P)
}

plotReportNRFD <- function(report.rfds, names, file.name, main.text) {
   titles <- c("TTR", "CTR_E", "CTR_L")
   cols <- c("black", "red", "blue")
   n <- length(names)
   
   colnames <- c("NAME", "X", "pch", "pos1", "pos2", "pos3", titles)
   rfds <- toTable(0, length(colnames), length(names), colnames)
   rfds$NAME <- names
   rfds$X    <- c(1, 3, 5, 7)
   rfds$pch  <- c(19, 17, 17, 17)
   rfds$pos1 <- c(3, 3, 3, 3)
   rfds$pos2 <- c(3, 3, 3, 1)
   rfds$pos3 <- c(1, 1, 1, 3)
   for (r in 1:length(names)) {
      rfds$TTR[r]   <- as.numeric(round0(report.rfds[[r]][1]*100, digit=1))
      rfds$CTR_E[r] <- as.numeric(round0(report.rfds[[r]][2]*100, digit=1))
      rfds$CTR_L[r] <- as.numeric(round0(report.rfds[[r]][3]*100, digit=1))
   }
   rfds$CTR_E <- c(7.4, 10.5, 13, 11.4)
   rfds$CTR_L <- c(4.7, 8.7, 10.1, 13.1)
   
   ##
   pdf(file.name, height=5, width=5.1)
   layout(matrix(c(1,2), 2, 1), widths=1, heights=c(1,1))           ## One figure each in row 1 and row 2; row 1 is 1/3 the height of row 2
   par(mar=c(1,4,3.6,1))
   plot(NULL, xlim=c(0.5, 9.1), ylim=c(55, 110), ylab="", main=main.text, col=cols[1], xaxt="n", yaxt="n", bty="n", pch=rfds$pch, cex.axis=1.1, cex.lab=1.2, cex.main=1.3)
   points(rfds$X, rfds$TTR, col=cols[1], pch=rfds$pch, cex=1.5)
   lines(x=rfds$X[1:4], y=rfds$TTR[1:4], type="l", lwd=3, col=cols[1])
   #lines(x=rfds$X[5:6], y=rfds$TTR[5:6], type="l", lwd=3, col=cols[1])
   
   text(8, rfds$TTR[n], " TTR ", cex=1.2, col="black")
   #abline(v=2, lty=5, lwd=1, col="black")
   #abline(v=8,  lty=5, lwd=1, col="black")
   for (n in 1:length(names))
      text(rfds$X[n], rfds$TTR[n], rfds$TTR[n], col=cols[1], pos=rfds$pos1[n], cex=1.2)
 
   ##
   axis(side=2, at=seq(60, 100, by=20), labels=c(60, 80, 100), cex.axis=1.1)
   legend("top", "Normal                              ", col="black", bty="n", pt.cex=1.4, pch=1, horiz=T, cex=1.2)
   legend("topright", "Tumour                    ", col="black", bty="n", pt.cex=1.2, pch=2, horiz=T, cex=1.2)
   #legend("top", "Primary bulks                            ", col="black", bty="n", pt.cex=1, pch=2, horiz=T, cex=1.2)
   #legend("topright", "Cell lines                ", col="black", bty="n", pt.cex=1, pch=0, horiz=T, cex=1.2)
   mtext("[%]                ", side=2, line=2.8, cex=1.2)
   
   ##
   par(mar=c(5,4,0,1))
   plot(NULL, xlim=c(0.5, 9.1), ylim=c(0, 15.8), ylab="", xlab="", col=cols[2], xaxt="n", bty="n", pch=rfds$pch, cex.axis=1.1, cex.lab=1.2, cex.main=1.3)
   points(rfds$X, rfds$CTR_L, col=cols[3], pch=rfds$pch, cex=1.5)
   lines(rfds$X[1:4], y=rfds$CTR_L[1:4], type="l", lwd=3, col=cols[3])
   lines(rfds$X[5:6], y=rfds$CTR_L[5:6], type="l", lwd=3, col=cols[3])
   
   points(rfds$X, rfds$CTR_E, col=cols[2], pch=rfds$pch, cex=1.5)
   lines(rfds$X[1:4], y=rfds$CTR_E[1:4], type="l", lwd=3, col=cols[2])
   lines(rfds$X[5:6], y=rfds$CTR_E[5:6], type="l", lwd=3, col=cols[2])
   
   text(8, rfds$CTR_E[n], "        CTR (E)", cex=1.2, col="red", pos=1) 
   text(8, rfds$CTR_L[n], "        CTR (L)", cex=1.2, col="blue", pos=3)
   #abline(v=2, lty=5, lwd=1, col="black")
   #abline(v=8,  lty=5, lwd=1, col="black")
   for (n in 1:length(names)) {
      text(rfds$X[n], rfds$CTR_E[n], rfds$CTR_E[n], col=cols[2], pos=rfds$pos2[n], cex=1.2)
      text(rfds$X[n], rfds$CTR_L[n], rfds$CTR_L[n], col=cols[3], pos=rfds$pos3[n], cex=1.2)
   }
   #axis(side=2, at=seq(0, 20, by=10), labels=c(0, 10, 20), cex.axis=1.1)

   ##
   #axis(side=1, at=seq(1, 1, by=2), labels=names[1], cex.axis=1.1)
   #axis(side=1, at=seq(1, 1, by=2), labels=c("n=92"), line=1.2, col=NA, cex.axis=1.1)
   axis(side=1, at=seq(1, 7, by=2), labels=names[1:4], cex.axis=1.2)
   axis(side=1, at=seq(1, 7, by=2), labels=c("n=92", "n=101", "n=56", "n=96"), line=1.2, col=NA, cex.axis=1.1)
   #axis(side=1, at=seq(9, 11, by=2), labels=names[5:6], cex.axis=1.1)
   #axis(side=1, at=seq(9, 11, by=2), labels=c("n=8", "n=14"), line=1.2, col=NA, cex.axis=1.1)
   
   mtext("                        Frequency", side=2, line=2.8, cex=1.2)
   dev.off()
}

plotReportNRFD12 <- function(report.rfds, names, file.name, main.text) {
   titles <- c("TTR", "CTR_E", "CTR_L")
   cols <- c("black", "red", "blue")
   n <- length(names)
 
   colnames <- c("NAME", "X", "pch", "pos1", "pos2", "pos3", titles)
   rfds <- toTable(0, length(colnames), length(names), colnames)
   rfds$NAME <- names
   rfds$X    <- c(1, 3, 5, 7)
   rfds$pch  <- c(19, 17, 17, 17)
   rfds$pos1 <- c(3, 3, 3, 3)
   rfds$pos2 <- c(3, 3, 3, 1)
   rfds$pos3 <- c(1, 1, 1, 3)
   #rfds$TTR0 <- c(0, 0, 0, 0)
   #for (r in 1:length(names)) {
   #   rfds$TTR0[r] <- as.numeric(round0((report.rfds[[r]][1] + report.rfds[[r]][2] + report.rfds[[r]][3])*100, digit=1))
   #}
   for (r in 1:length(names)) {
      rfds$TTR[r]   <- as.numeric(round0(report.rfds[[r]][1]*100, digit=1))
      rfds$CTR_E[r] <- as.numeric(round0(report.rfds[[r]][2]*100, digit=1))
      rfds$CTR_L[r] <- as.numeric(round0(report.rfds[[r]][3]*100, digit=1))
   }
   rfds$CTR_E <- c(7.2, 10.4, 13.6, 11.2)
   rfds$CTR_L <- c(5, 8.1, 10.7, 12.5)
   rfds$TTR0 <- rfds$TTR + rfds$CTR_E + rfds$CTR_L
    
   ##
   pdf(file.name, height=5, width=5.1)
   layout(matrix(c(1,2), 2, 1), widths=1, heights=c(1,1))           ## One figure each in row 1 and row 2; row 1 is 1/3 the height of row 2
   par(mar=c(1,4,3.6,1))
   plot(NULL, xlim=c(0.5, 9.1), ylim=c(55, 110), ylab="", main=main.text, col=cols[1], xaxt="n", yaxt="n", bty="n", pch=rfds$pch, cex.axis=1.1, cex.lab=1.2, cex.main=1.3)
   points(rfds$X, rfds$TTR, col=cols[1], pch=rfds$pch, cex=1.5)
   lines(x=rfds$X[1:4], y=rfds$TTR[1:4], type="l", lwd=3, col=cols[1])
   
   #rfds$TTR0 <- c(87.8, 80.7, 76.9, 75.5)
   points(rfds$X, rfds$TTR0, col="goldenrod2", pch=rfds$pch, cex=1.5)
   lines(x=rfds$X[1:4], y=rfds$TTR0[1:4], type="l", lty=5, lwd=2.5, col="goldenrod2")
   text(8, rfds$TTR0[n], "       TTR+CTR", cex=1.2, col="goldenrod2")
   for (n in 1:length(names))
      text(rfds$X[n], rfds$TTR0[n], rfds$TTR0[n], col="goldenrod2", pos=rfds$pos1[n], cex=1.2)
   
   text(8, rfds$TTR[n], " TTR ", cex=1.2, col="black")
   #abline(v=2, lty=5, lwd=1, col="black")
   for (n in 1:length(names))
      text(rfds$X[n], rfds$TTR[n], rfds$TTR[n], col=cols[1], pos=rfds$pos1[n], cex=1.2)

   ##
   axis(side=2, at=seq(60, 100, by=20), labels=c(60, 80, 100), cex.axis=1.1)
   legend("top", "Normal                              ", col="black", bty="n", pt.cex=1.4, pch=1, horiz=T, cex=1.2)
   legend("topright", "Tumour                    ", col="black", bty="n", pt.cex=1.2, pch=2, horiz=T, cex=1.2)
   mtext("[%]                ", side=2, line=2.8, cex=1.2)
   
   ##
   par(mar=c(5,4,0,1))
   plot(NULL, xlim=c(0.5, 9.1), ylim=c(0, 15.8), ylab="", xlab="", col=cols[2], xaxt="n", bty="n", pch=rfds$pch, cex.axis=1.1, cex.lab=1.2, cex.main=1.3)
   points(rfds$X, rfds$CTR_L, col=cols[3], pch=rfds$pch, cex=1.5)
   lines(rfds$X[1:4], y=rfds$CTR_L[1:4], type="l", lwd=3, col=cols[3])
 
   points(rfds$X, rfds$CTR_E, col=cols[2], pch=rfds$pch, cex=1.5)
   lines(rfds$X[1:4], y=rfds$CTR_E[1:4], type="l", lwd=3, col=cols[2])
   
   text(8, rfds$CTR_E[n], "        CTR (E)", cex=1.2, col="red", pos=1)
   text(8, rfds$CTR_L[n], "        CTR (L)", cex=1.2, col="blue", pos=3)
   #abline(v=2, lty=5, lwd=1, col="black")
   for (n in 1:length(names)) {
      text(rfds$X[n], rfds$CTR_E[n], rfds$CTR_E[n], col=cols[2], pos=rfds$pos2[n], cex=1.2)
      text(rfds$X[n], rfds$CTR_L[n], rfds$CTR_L[n], col=cols[3], pos=rfds$pos3[n], cex=1.2)
   }

   ##
   #axis(side=1, at=seq(1, 1, by=2), labels=names[1], cex.axis=1.1)
   #axis(side=1, at=seq(1, 1, by=2), labels=c("n=46+46"), line=1.2, col=NA, cex.axis=1.1)
   axis(side=1, at=seq(1, 7, by=2), labels=names[1:4], cex.axis=1.2)
   axis(side=1, at=seq(1, 7, by=2), labels=c("n=46,46", "n=50,51", "n=28,28", "n=48,48"), line=1.2, col=NA, cex.axis=1.2)

   mtext("                        Frequency", side=2, line=2.8, cex=1.2)
   dev.off()
}

# -----------------------------------------------------------------------------
# Visualisation of bootstrap re-sampling data (Histogram, RFD, and RT)
# Last Modified: 13/11/18
# -----------------------------------------------------------------------------
## https://www.r-graph-gallery.com/190-mirrored-histogram/
## http://www.r-graph-gallery.com/72-set-margin-size-with-par-mar-function/
## https://www.statmethods.net/advgraphs/layout.html
## https://stackoverflow.com/questions/6461209/how-to-round-up-to-the-nearest-10-or-100-or-x
## https://stackoverflow.com/questions/15717545/set-the-intervals-of-x-axis-using-r
plotBootstrapHist <- function(nrds.RT.BSTRPS, file.name, main.text, xlab.text, breaks, boundary.break) {
   adjustcolor.blue <- adjustcolor("steelblue1", alpha.f=0.15)
   adjustcolor.orange <- adjustcolor("sandybrown", alpha.f=0.15)
 
   cols <- rep("steelblue1", breaks)
   cols[(breaks/2 + boundary.break):breaks] <- "sandybrown"
   cols[(breaks/2 - boundary.break + 1):50] <- "white"
   cols[51:(breaks/2 + boundary.break)] <- "white"
   #cols[51:breaks] <- "sandybrown"
   
   h <- hist(nrds.RT.BSTRPS$NEG, breaks=breaks)
   ymax <- max(c(h$counts[2:4], h$counts[(breaks-3):(breaks-1)]))   ## Calculatte max frequency in row 2 before next line
   h$counts <- h$counts/1000                                        ## Change frequency scale to x1000 in row 1
   
   pdf(file.name, height=5, width=5)
   #par(mfrow=c(2,1))
   layout(matrix(c(1,2), 2, 1), widths=1, heights=c(1,2))           ## One figure each in row 1 and row 2; row 1 is 1/3 the height of row 2
   ylim <- sort(c(h$counts[1], h$counts[breaks]), decreasing=F)     ## Min and max frequencies in row 1 (in x1000 scale)
   if (ylim[1] > 500) {                                             ## Round up to the nearest power of 10 (in x1000 scale)
      ylim[1] <- floor(ylim[1]/100)*100
   } else if (ylim[1] > 10) {
      ylim[1] <- floor(ylim[1]/10)*10
   } else
      ylim[1] <- floor(ylim[1])
   par(mar=c(1,4,3.6,1))
   plot(h, main=main.text[1], ylab="Freq. (x1000)", xlab="", ylim=ylim, col=cols, xaxt="n", cex.axis=1.1, cex.lab=1.2, cex.main=1.25)
   #abline(v=500, lty=5, lwd=1, col="black")
   abline(v=950, lty=5, lwd=1, col="black")
   abline(v=50,  lty=5, lwd=1, col="black")

   ##
   par(mar=c(5,4,0,1))
   hist(nrds.RT.BSTRPS$NEG, main="", ylab="Frequency", xlab=xlab.text, ylim=c(0, ymax), breaks=breaks, col=cols, las=1, axes=F, cex.axis=1.1, cex.lab=1.2, cex.main=1.25)
   if (ymax < 1000) {
      axis(side=2, at=seq(0, ymax, by=250))
   } else if (ymax < 3000) {
      axis(side=2, at=seq(0, ymax, by=500))
   } else if (ymax > 20000) {
      axis(side=2, at=seq(0, ymax, by=10000))
   } else
      axis(side=2, at=seq(0, ymax, by=1000))
   axis(side=1, at=seq(0, 1000, by=250))
   #abline(v=500, lty=5, lwd=1, col="black")
   abline(v=950, lty=5, lwd=1, col="black")
   abline(v=50,  lty=5, lwd=1, col="black")
   text(500, ymax*4/5, "RFD = 0", cex=1.2, col="black") 
   
   mtext(main.text[2], line=4.8, cex=1.2)   ## separator(nrow(nrds.RT.BSTRPS)),
   dev.off()
}

#boundary.upper <- 950   ## 500-520 breaks
#boundary.lower <-  50   ## 480-500 breaks
#boundary.break <-  45   ## 1 breaks each centering 500
#file.name <- file.path(wd.rt.plots, paste0("hist_", base, "_rpkm_SLOPE_RFD>0.9_white.pdf"))
#main.text <- c(paste0(BASE, " bootstrap distribution"), "Chr1-22 (1-kbs)")   #paste0("Chr1-22 (1-kbs)"))
#xlab.text <- "Number of right-leading resamplings"
#plotBootstrapHist(nrds.RT.BSTRPS, file.name, main.text, xlab.text, 100, boundary.break)

# -----------------------------------------------------------------------------
# Visualisation of bootstrap re-sampling data (Histogram, RFD, and RT)
# Last Modified: 28/11/19
# -----------------------------------------------------------------------------
plotBootstrapRFD <- function(file.name, BASE, chr, xmin, xmax, nrds.RT.NRFD, bed.gc.chr, boundary.upper, boundary.lower, ext, width, kb, withUnclassified=F) {
   overlaps <- intersect(rownames(bed.gc.chr), nrds.RT.NRFD$BED)
   nrds.RT.NRFD.chr <- nrds.RT.NRFD[overlaps,]
   bed.gc.chr <- bed.gc.chr[overlaps,]
 
   adjustcolor.gray <- adjustcolor("darkgray", alpha.f=0.08)
     
   rights <- rownames(subset(nrds.RT.NRFD.chr, RFD >= 0.9))
   lefts  <- rownames(subset(nrds.RT.NRFD.chr, RFD <= -0.9))
   boundary.rights <- rownames(subset(subset(nrds.RT.NRFD.chr, RFD < 0.9), RFD >= 0))
   boundary.lefts <- rownames(subset(subset(nrds.RT.NRFD.chr,  RFD < 0),   RFD > -0.9))
   
   boundaries <- c(boundary.rights, boundary.lefts)
   initiations <- rownames(subset(nrds.RT.NRFD.chr[boundaries,], SPLINE > 0))
   terminations <- rownames(subset(nrds.RT.NRFD.chr[boundaries,], SPLINE < 0))
   #unclassified <- rownames(subset(nrds.RT.NRFD.chr[boundaries,], RT == 0))
   #unclassified <- c(unclassified, rownames(nrds.RT.NRFD.chr[boundaries,])[which(is.na(nrds.RT.NRFD.chr[boundaries,]$RT) == T)])
      
   if (width == 10) main.text <- paste0(BASE, " bootstrap replication fork directionality (RFD)")
   else main.text <- paste0(BASE, " bootstrap RFD")
   if (withUnclassified)
      main.text <- paste0(main.text, " (", kb, " kb)")
   
   if (!is.na(xmin) && !is.na(xmax)) file.name <- paste0(file.name, "_", xmin/1E6, "-", xmax/1E6, "Mb_", kb, "kb")
   if (withUnclassified) file.name <- paste0(file.name, "_with-un")
   if (is.na(xmin)) {
      start <- bed.gc.chr[rownames(nrds.RT.NRFD.chr)[1],]$START
      if (start < 5000000) xmin <- 0
      else xmin <- start
   }
   if (is.na(xmax)) xmax <- subset(chromInfo, chrom == chr)$size
   
   ## Initiation plot
   if (ext == "pdf") {
      pdf(paste0(file.name, ".pdf"), height=5, width=width)
   } else if (ext == "png")
      png(paste0(file.name, ".png"), height=5, width=width, units="in", res=300)   ## ADD 16/05/17: res=300
   
   ###
   ## Initiate RT plot
   layout(matrix(c(1,2), 2, 1), widths=1, heights=c(1,1))   ## One figure each in row 1 and row 2   ## See plotBootstrapsHist()
   par(mar=c(1,4,4,1))
   ylab.text <- "RT [log2]"

   plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=c(-2, 2), xlab="", ylab=ylab.text, main=main.text, xaxt="n", yaxt="n", cex.axis=1.1, cex.lab=1.2, cex.main=1.3)
   points(bed.gc.chr$START/1E6, nrds.RT.NRFD.chr$RT, col=adjustcolor.gray, pch=16, cex=0.35)
   
   axis(side=2, at=seq(-2, 2, by=4), labels=c("\u22122", 2), cex.axis=1.1)
   axis(side=2, at=seq(-1, 1, by=1), labels=c("\u22121", 0, 1), cex.axis=1.1)
   abline(h=0, lty=5, lwd=1, col="black")
   
   ## Plot cytobands (before smoothing spline)
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.4, col="lightgrey") 
   
   ## Plot smoothing spline with right-/left-leading positions
   points(bed.gc.chr[lefts,]$START/1E6,  nrds.RT.NRFD.chr[lefts,]$SPLINE, col="steelblue1", cex=0.4)
   points(bed.gc.chr[rights,]$START/1E6, nrds.RT.NRFD.chr[rights,]$SPLINE, col="sandybrown", cex=0.4)
   points(bed.gc.chr[terminations,]$START/1E6, nrds.RT.NRFD.chr[terminations,]$SPLINE, col="blue", pch=19, cex=0.4)
   points(bed.gc.chr[initiations,]$START/1E6,  nrds.RT.NRFD.chr[initiations,]$SPLINE,  col="red", pch=19, cex=0.4)
   if (withUnclassified && length(unclassified) != 0)
      points(bed.gc.chr[unclassified,]$START/1E6,  nrds.RT.NRFD.chr[unclassified,]$SPLINE, col="#01DF01", pch=19, cex=0.5)
   
   ## Plot legend
   legend("topright", "CTR (E)", col="red", bty="n", pt.cex=1, lty=1, lwd=3, pch=NA, horiz=T, cex=1.2)
   legend("bottomright", "CTR (L)", col="blue", bty="n", pt.cex=1, lty=1, lwd=3, pch=NA, horiz=T, cex=1.2)
   #if (withUnclassified) legend("bottomleft", "CTR (UN)", col="#01DF01", bty="n", pt.cex=1, lty=1, lwd=3, pch=NA, horiz=T, cex=1.2)
   mtext("", line=0.25, cex=1.2)   ## separator(nrow(nrds.RT.BSTRPS)),
      
   ###
   ## Initiate RFD plot
   par(mar=c(5.5,4,0,1))
   xlab.text <- paste0("Chromosome ", gsub("chr", "", chr), " position [Mb]")
   ylab.text <- "RFD"
   
   plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=c(-1.8, 1.8), xlab=xlab.text, ylab=ylab.text, main="", yaxt="n", cex.axis=1.1, cex.lab=1.2)
   axis(side=2, at=seq(-1, 1, by=1), labels=c("\u22121", 0, 1), cex.axis=1.1)
   #abline(h=0, lty=5, lwd=1, col="black")
   abline(h=0.9, lty=5, lwd=1, col="black")
   abline(h=-0.9, lty=5, lwd=1, col="black")
   
   ## Plot cytobands (before points)
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.4, col="lightgrey")

   points(bed.gc.chr[lefts,]$START/1E6,  nrds.RT.NRFD.chr[lefts,]$RFD,  col="steelblue1", cex=0.4)
   points(bed.gc.chr[rights,]$START/1E6, nrds.RT.NRFD.chr[rights,]$RFD, col="sandybrown", cex=0.4)
   points(bed.gc.chr[terminations,]$START/1E6, nrds.RT.NRFD.chr[terminations,]$RFD, col="blue", pch=19, cex=0.4)
   points(bed.gc.chr[initiations,]$START/1E6,  nrds.RT.NRFD.chr[initiations,]$RFD,  col="red",  pch=19, cex=0.4)
   #if (withUnclassified && length(unclassified) != 0)
   #   points(bed.gc.chr[unclassified,]$START/1E6, nrds.RT.NRFD.chr[unclassified,]$RFD, col="#01DF01", pch=19, cex=0.5)
   
   ## Plot legend
   legend("topright", "TTR (R)", col="sandybrown", bty="n", pt.cex=1, lty=1, lwd=3, pch=NA, horiz=T, cex=1.2)
   legend("bottomright", "TTR (L)", col="steelblue1", bty="n", pt.cex=1, lty=1, lwd=3, pch=NA, horiz=T, cex=1.2)
   dev.off()
}
