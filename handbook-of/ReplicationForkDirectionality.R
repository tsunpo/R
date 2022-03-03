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
   load(file=file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.m.rt.RT.", column, "_", "chr1", ".RData")))
   nrds.RT.BSTRPS <- nrds.RT.BSTRPS.chr
   for (c in 2:22) {
      chr <- chrs[c]

      load(file=file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.m.rt.RT.", column, "_", chr, ".RData")))
      nrds.RT.BSTRPS <- rbind(nrds.RT.BSTRPS, nrds.RT.BSTRPS.chr)
   }
   
   return(nrds.RT.BSTRPS)
}

getRTNRFD <- function(nrds, nrds.RT.BSTRPS, bed.gc, kb) {
   nrds.RT <- NULL
   for (c in 1:22) {
      chr <- chrs[c]
      bed.gc.chr <- subset(bed.gc, CHR == chr)
  
      nrds.RT.chr <- setSpline2(nrds, bed.gc.chr, "RT", returnAll=T)
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
   for (c in 2:22) {
      chr <- chrs[c]
      bed.gc.chr <- subset(bed.gc, CHR == chr)
      
      overlaps <- intersect(rownames(bed.gc.chr), nrds.RT$BED)
      bed.gc.chr.o <- bed.gc.chr[overlaps,]
      nrds.chr.o <- nrds.RT[overlaps,]
   
      ## NRFD
      ## https://stackoverflow.com/questions/41061140/how-to-calculate-the-average-slope-within-a-moving-window-in-r
      ## https://www.rdocumentation.org/packages/zoo/versions/1.8-6/topics/rollapply
      rollingSlope.lm <- function(vector) {
         return(coef(lm(vector ~ seq(vector)))[2])
      }
      rollingSlope.lm.fit <- function(vector) {
         return(coef(.lm.fit(cbind(1, seq(vector)), vector))[2])
      }
      
      Slope.lm.fit = rollapply(nrds.chr.o$RFD, width=kb, FUN=rollingSlope.lm, fill=NA, partial=T)
      nrds.chr.o$NRFD <- Slope.lm.fit
      #nrow(subset(subset(subset(nrds.chr.o, RFD < 0.9), RFD > -0.9), NRFD == 0))
      
      if (is.null(nrds.RT.NRFD)) {
         nrds.RT.NRFD <- nrds.chr.o
      } else {
         nrds.RT.NRFD <- rbind(nrds.RT.NRFD, nrds.chr.o)
      }
   }
   
   ## Test lm
   save(nrds.RT.NRFD, file=file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.d.rt.log2s.nrfd.", kb, "kb_", "m2-m1_fit", ".RData")))
   writeTable(nrds.RT.NRFD, gzfile(file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.d.rt.log2s.nrfd.", kb, "kb_", "m2-m1_fit", ".txt.gz"))), colnames=T, rownames=T, sep="\t")
   nrds.RT.NRFD.sclc <- nrds.RT.NRFD
   
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

getBootstrapReport3 <- function(rfd, nrds.RT.RFD.1, nrds.RT.RFD.2, nrds.RT.RFD.3, name.1) {
   colnames <- c("RFD", "SAMPLES", name.1, "ALL", "Overlapping_N", "Overlapping_P", "RFD_N", "RFD_P", "Mappable")
   report <- toTable(NA, length(colnames), 6, colnames)
 
   ###
   ## |RFD| ≥ 0.9
   nrds.RT.RFD.1.t <- getBootstrapTTR(nrds.RT.RFD.1, rfd)
   nrds.RT.RFD.2.t <- getBootstrapTTR(nrds.RT.RFD.2, rfd)
   nrds.RT.RFD.3.t <- getBootstrapTTR(nrds.RT.RFD.3, rfd)
   
   ## TTR
   report$RFD[1:2] <- "TTR"
   report$SAMPLES[1] <- name.1
   #report$SAMPLES[2] <- NA
   report$RFD_N[1] <- nrow(nrds.RT.RFD.1.t)
   #report$RFD_N[2] <- NA
   report$Mappable[1] <- nrow(nrds.RT.RFD.1)
   #report$Mappable[2] <- NA
   report$RFD_P[1] <- report$RFD_N[1]/report$Mappable[1]
   #report$RFD_P[2] <- NA
 
   ##
   overlaps.t <- intersect(intersect(rownames(nrds.RT.RFD.1.t), rownames(nrds.RT.RFD.2.t)), rownames(nrds.RT.RFD.3.t))
   report$Overlapping_N[1] <- length(overlaps.t)
   #report$Overlapping_N[2] <- NA
   report$Overlapping_P[1] <- length(overlaps.t)/nrow(nrds.RT.RFD.1)
   #report$Overlapping_P[2] <- NA
 
   nrds.1.RT.o <- nrds.RT.RFD.1.t[overlaps.t,]
   #nrds.2.RT.o <- nrds.RT.RFD.2.t[overlaps.t,]
 
   #nrds.1.RT.o$SIGN <- nrds.1.RT.o$SLOPE * nrds.2.RT.o$SLOPE
   #report[2, 3] <- length(which(nrds.1.RT.o$SIGN > 0))/length(overlaps.t)
 
   ###
   ## |RFD| < 0.9
   nrds.RT.RFD.1.c <- getBootstrapCTR(nrds.RT.RFD.1, rfd)
   nrds.RT.RFD.2.c <- getBootstrapCTR(nrds.RT.RFD.2, rfd)
   nrds.RT.RFD.3.c <- getBootstrapCTR(nrds.RT.RFD.3, rfd)
   
   nrds.RT.RFD.1.c.e <- subset(nrds.RT.RFD.1.c, NRFD > 0)
   nrds.RT.RFD.1.c.l <- subset(nrds.RT.RFD.1.c, NRFD < 0)
   nrds.RT.RFD.2.c.e <- subset(nrds.RT.RFD.2.c, NRFD > 0)
   nrds.RT.RFD.2.c.l <- subset(nrds.RT.RFD.2.c, NRFD < 0)
   nrds.RT.RFD.3.c.e <- subset(nrds.RT.RFD.3.c, NRFD > 0)
   nrds.RT.RFD.3.c.l <- subset(nrds.RT.RFD.3.c, NRFD < 0)
 
   ## CTR (IZ)
   report$RFD[3:4] <- "CTR-IZ"
   report$SAMPLES[3] <- name.1
   report$SAMPLES[4] <- NA
 
   report$RFD_N[3] <- nrow(nrds.RT.RFD.1.c.e)
   report$RFD_N[4] <- NA
   report$Mappable[3] <- nrow(nrds.RT.RFD.1)
   report$Mappable[4] <- NA
   report$RFD_P[3] <- report$RFD_N[3]/report$Mappable[3]
   report$RFD_P[4] <- NA
 
   ##
   overlaps.c.e <- intersect(intersect(rownames(nrds.RT.RFD.1.c.e), rownames(nrds.RT.RFD.2.c.e)), rownames(nrds.RT.RFD.3.c.e))
   report$Overlapping_N[3] <- length(overlaps.c.e)
   #report$Overlapping_N[4] <- length(overlaps.c.e)   
   report$Overlapping_P[3] <- length(overlaps.c.e)/nrow(nrds.RT.RFD.1)
   #report$Overlapping_P[4] <- length(overlaps.c.e)/nrow(nrds.RT.RFD.2)
 
   #nrds.1.RT.o <- nrds.RT.RFD.1.c.e[overlaps.c.e,]
   #nrds.2.RT.o <- nrds.RT.RFD.2.c.e[overlaps.c.e,]
 
   #nrds.1.RT.o$SIGN <- nrds.1.RT.o$SLOPE * nrds.2.RT.o$SLOPE
   #report[4, 3] <- length(which(nrds.1.RT.o$SIGN > 0))/length(overlaps.c.e)
 
   ## CTR (TZ)
   report$RFD[5:6] <- "CTR-TZ"
   report$SAMPLES[5] <- name.1
   report$SAMPLES[6] <- NA
 
   report$RFD_N[5] <- nrow(nrds.RT.RFD.1.c.l)
   report$RFD_N[6] <- NA
   report$Mappable[5] <- nrow(nrds.RT.RFD.1)
   report$Mappable[6] <- NA
   report$RFD_P[5] <- report$RFD_N[5]/report$Mappable[5]
   report$RFD_P[6] <- NA
 
   ##
   overlaps.c.l <- intersect(intersect(rownames(nrds.RT.RFD.1.c.l), rownames(nrds.RT.RFD.2.c.l)), rownames(nrds.RT.RFD.3.c.l))
   report$Overlapping_N[5] <- length(overlaps.c.l)
   #report$Overlapping_N[6] <- length(overlaps.c.l)   
   report$Overlapping_P[5] <- length(overlaps.c.l)/nrow(nrds.RT.RFD.1)
   #report$Overlapping_P[6] <- length(overlaps.c.l)/nrow(nrds.RT.RFD.2)
 
   #nrds.1.RT.o <- nrds.RT.RFD.1.c.l[overlaps.c.l,]
   #nrds.2.RT.o <- nrds.RT.RFD.2.c.l[overlaps.c.l,]
 
   #nrds.1.RT.o$SIGN <- nrds.1.RT.o$SLOPE * nrds.2.RT.o$SLOPE
   #report[6, 3] <- length(which(nrds.1.RT.o$SIGN > 0))/length(overlaps.c.l)
 
   return(report)
}

getReportRFD <- function(report, name) {
   return(subset(report, SAMPLES == name)$RFD_P)
}

getReportRFD12 <- function(report, name) {
   return(subset(report, SAMPLES == name)$Overlapping_P)
}

plotReportNRFD <- function(report.rfds, names, file.name, main.text) {
   titles <- c("TTR", "CTR_IZ", "CTR_TZ")
   cols <- c("black", red, blue)
   n <- length(names)
   
   colnames <- c("NAME", "X", "pch", "pos1", "pos2", "pos3", titles)
   rfds <- toTable(0, length(colnames), length(names), colnames)
   rfds$NAME <- names
   rfds$X    <- c(1, 3, 5, 7)
   rfds$pch  <- c(19, 17, 17, 17)
   rfds$cex  <- c(1.8, 1.6, 1.6, 1.6)
   rfds$pos1 <- c(3, 3, 3, 3)
   rfds$pos2 <- c(3, 1, 3, 1)
   rfds$pos3 <- c(1, 3, 1, 3)
   for (r in 1:length(names)) {
      rfds$TTR[r]   <- as.numeric(round0(report.rfds[[r]][1]*100, digit=1))
      rfds$CTR_E[r] <- as.numeric(round0(report.rfds[[r]][2]*100, digit=1))
      rfds$CTR_L[r] <- as.numeric(round0(report.rfds[[r]][3]*100, digit=1))
   }
   
   ##
   pdf(file.name, height=5, width=5.2)
   layout(matrix(c(1,2), 2, 1), widths=1, heights=c(1,1))           ## One figure each in row 1 and row 2; row 1 is 1/3 the height of row 2
   par(mar=c(1,3.6,3.6,0))
   plot(NULL, xlim=c(0.5, 9.1), ylim=c(55, 115), ylab="", main=main.text, col=cols[1], xaxt="n", yaxt="n", bty="n", pch=rfds$pch, cex.axis=1.25, cex.lab=1.3, cex.main=1.35)
   points(rfds$X, rfds$TTR, col=cols[1], pch=rfds$pch, cex=rfds$cex)
   lines(x=rfds$X[1:4], y=rfds$TTR[1:4], type="l", lwd=3, col=cols[1])
   #lines(x=rfds$X[5:6], y=rfds$TTR[5:6], type="l", lwd=3, col=cols[1])
   
   text(8, rfds$TTR[n], "   TTR ", cex=1.3, col="black", font=2)
   #abline(v=4, lty=5, lwd=1, col="black")
   #abline(v=8,  lty=5, lwd=1, col="black")
   for (n in 1:length(names))
      text(rfds$X[n], rfds$TTR[n], rfds$TTR[n], col=cols[1], pos=rfds$pos1[n], cex=1.3)
 
   ##
   axis(side=2, at=seq(60, 100, by=20), labels=c(60, 80, 100), cex.axis=1.25)
   legend("top", "Normal                              ", col="black", bty="n", pch=1, pt.cex=2, horiz=T, cex=1.3)
   legend("topright", "Tumor                    ", col="black", bty="n", pch=2, pt.cex=1.6, horiz=T, cex=1.3)
   #legend("top", "Primary bulks                            ", col="black", bty="n", pt.cex=1, pch=2, horiz=T, cex=1.2)
   #legend("topright", "Cell lines                ", col="black", bty="n", pt.cex=1, pch=0, horiz=T, cex=1.2)
   mtext("[%]          ", side=2, line=2.6, cex=1.3)
   
   ##
   par(mar=c(5,3.6,0,0))
   plot(NULL, xlim=c(0.5, 9.1), ylim=c(0, 15), ylab="", xlab="", col=cols[2], xaxt="n", bty="n", pch=rfds$pch, cex.axis=1.25, cex.lab=1.3, cex.main=1.35)
   points(rfds$X, rfds$CTR_L, col=cols[3], pch=rfds$pch, cex=rfds$cex)
   lines(rfds$X[1:4], y=rfds$CTR_L[1:4], type="l", lwd=3, col=cols[3])
   lines(rfds$X[5:6], y=rfds$CTR_L[5:6], type="l", lwd=3, col=cols[3])
   
   points(rfds$X, rfds$CTR_E, col=cols[2], pch=rfds$pch, cex=rfds$cex)
   lines(rfds$X[1:4], y=rfds$CTR_E[1:4], type="l", lwd=3, col=cols[2])
   lines(rfds$X[5:6], y=rfds$CTR_E[5:6], type="l", lwd=3, col=cols[2])
   
   text(8, rfds$CTR_E[n], "      IZ ", cex=1.3, col=red, pos=1, font=2) 
   text(8, rfds$CTR_L[n], "     TZ ", cex=1.3, col=blue, pos=3, font=2)
   #abline(v=4, lty=5, lwd=1, col="black")
   #abline(v=8,  lty=5, lwd=1, col="black")
   for (n in 1:length(names)) {
      text(rfds$X[n], rfds$CTR_E[n], rfds$CTR_E[n], col=red, pos=rfds$pos2[n], cex=1.3)
      text(rfds$X[n], rfds$CTR_L[n], rfds$CTR_L[n], col=blue, pos=rfds$pos3[n], cex=1.3)
   }
   #axis(side=2, at=seq(0, 20, by=10), labels=c(0, 10, 20), cex.axis=1.1)

   ##
   #axis(side=1, at=seq(1, 1, by=2), labels=names[1], cex.axis=1.1)
   #axis(side=1, at=seq(1, 1, by=2), labels=c("n=92"), line=1.2, col=NA, cex.axis=1.1)
   axis(side=1, at=seq(1, 7, by=4), labels=c(names[1], names[3]), cex.axis=1.3)
   axis(side=1, at=seq(3, 7, by=4), labels=c(names[2], names[4]), cex.axis=1.3)
   axis(side=1, at=seq(1, 7, by=2), labels=c("n=92", "n=101", "n=56", "n=96"), line=1.2, col=NA, cex.axis=1.3)
   #axis(side=1, at=seq(1, 7, by=2), labels=c("n=46", "n=51", "n=28", "n=48"), line=1.2, col=NA, cex.axis=1.25)
   #axis(side=1, at=seq(9, 11, by=2), labels=names[5:6], cex.axis=1.1)
   #axis(side=1, at=seq(9, 11, by=2), labels=c("n=8", "n=14"), line=1.2, col=NA, cex.axis=1.1)
   
   mtext("                      Frequency", side=2, line=2.6, cex=1.3)
   dev.off()
}

plotReportNRFD12 <- function(report.rfds, names, file.name, main.text) {
   titles <- c("TTR", "CTR_E", "CTR_L")
   cols <- c("black", red, blue)
   n <- length(names)
 
   colnames <- c("NAME", "X", "pch", "pos1", "pos2", "pos3", titles)
   rfds <- toTable(0, length(colnames), length(names), colnames)
   rfds$NAME <- names
   rfds$X    <- c(1, 3, 5, 7)
   rfds$pch  <- c(19, 17, 17, 17)
   rfds$cex  <- c(1.8, 1.6, 1.6, 1.6)
   rfds$pos1 <- c(3, 3, 3, 3)
   rfds$pos2 <- c(3, 1, 3, 1)
   rfds$pos3 <- c(1, 3, 1, 3)
   #rfds$TTR0 <- c(0, 0, 0, 0)
   #for (r in 1:length(names)) {
   #   rfds$TTR0[r] <- as.numeric(round0((report.rfds[[r]][1] + report.rfds[[r]][2] + report.rfds[[r]][3])*100, digit=1))
   #}
   for (r in 1:length(names)) {
      rfds$TTR[r]   <- as.numeric(round0(report.rfds[[r]][1]*100, digit=1))
      rfds$CTR_E[r] <- as.numeric(round0(report.rfds[[r]][2]*100, digit=1))
      rfds$CTR_L[r] <- as.numeric(round0(report.rfds[[r]][3]*100, digit=1))
   }
   rfds$TTR0 <- rfds$TTR + rfds$CTR_E + rfds$CTR_L
   
   ##
   pdf(file.name, height=5, width=5.2)
   layout(matrix(c(1,2), 2, 1), widths=1, heights=c(1,1))           ## One figure each in row 1 and row 2; row 1 is 1/3 the height of row 2
   par(mar=c(1,3.6,3.6,0))
   plot(NULL, xlim=c(0.5, 9.1), ylim=c(55, 115), ylab="", main=main.text, col=cols[1], xaxt="n", yaxt="n", bty="n", pch=rfds$pch, cex.axis=1.2, cex.lab=1.25, cex.main=1.3)
   points(rfds$X, rfds$TTR, col=cols[1], pch=rfds$pch, cex=rfds$cex)
   lines(x=rfds$X[1:4], y=rfds$TTR[1:4], type="l", lwd=3, col=cols[1])
   
   #rfds$TTR0 <- c(87.8, 80.7, 76.9, 75.5)
   #points(rfds$X, rfds$TTR0, col="#01DF01", pch=rfds$pch, cex=1.5)
   #lines(x=rfds$X[1:4], y=rfds$TTR0[1:4], type="l", lty=5, lwd=2.5, col="#01DF01")
   #text(8, rfds$TTR0[n], "       TTR+CTR", cex=1.2, col="#01DF01")
   #for (n in 1:length(names))
   #   text(rfds$X[n], rfds$TTR0[n], rfds$TTR0[n], col="#01DF01", pos=rfds$pos1[n], cex=1.2)

   text(8, rfds$TTR[n], " TTR ", cex=1.25, col="black")
   #abline(v=2, lty=5, lwd=1, col="black")
   for (n in 1:length(names))
      text(rfds$X[n], rfds$TTR[n], rfds$TTR[n], col=cols[1], pos=rfds$pos1[n], cex=1.25)

   ##
   axis(side=2, at=seq(60, 100, by=20), labels=c(60, 80, 100), cex.axis=1.2)
   legend("top", "Normal                              ", col="black", bty="n", pch=1, pt.cex=2, horiz=T, cex=1.3)
   legend("topright", "Tumour                    ", col="black", bty="n", pch=2, pt.cex=1.6, horiz=T, cex=1.3)
   mtext("[%]            ", side=2, line=2.6, cex=1.25)
   
   ##
   par(mar=c(5,3.6,0,0))
   plot(NULL, xlim=c(0.5, 9.1), ylim=c(0, 15), ylab="", xlab="", col=cols[2], xaxt="n", bty="n", pch=rfds$pch, cex.axis=1.2, cex.lab=1.25, cex.main=1.3)
   points(rfds$X, rfds$CTR_L, col=cols[3], pch=rfds$pch, cex=rfds$cex)
   lines(rfds$X[1:4], y=rfds$CTR_L[1:4], type="l", lwd=3, col=cols[3])
  
   points(rfds$X, rfds$CTR_E, col=cols[2], pch=rfds$pch, cex=rfds$cex)
   lines(rfds$X[1:4], y=rfds$CTR_E[1:4], type="l", lwd=3, col=cols[2])

   text(8, rfds$CTR_E[n], "        CTR (IZ)", cex=1.25, col=red, pos=1) 
   text(8, rfds$CTR_L[n], "         CTR (TZ)", cex=1.25, col=blue, pos=3)
   #abline(v=2, lty=5, lwd=1, col="black")
   for (n in 1:length(names)) {
      text(rfds$X[n], rfds$CTR_E[n], rfds$CTR_E[n], col=cols[2], pos=rfds$pos2[n], cex=1.25)
      text(rfds$X[n], rfds$CTR_L[n], rfds$CTR_L[n], col=cols[3], pos=rfds$pos3[n], cex=1.25)
   }

   ##
   #axis(side=1, at=seq(1, 1, by=2), labels=names[1], cex.axis=1.1)
   #axis(side=1, at=seq(1, 1, by=2), labels=c("n=46+46"), line=1.2, col=NA, cex.axis=1.1)
   axis(side=1, at=seq(1, 7, by=2), labels=names[1:4], cex.axis=1.25)
   axis(side=1, at=seq(1, 7, by=2), labels=c("n=46,46", "n=50,51", "n=28,28", "n=48,48"), line=1.2, col=NA, cex.axis=1.25)

   mtext("                        Frequency", side=2, line=2.6, cex=1.25)
   dev.off()
}

plotReportNRFDEL <- function(report.rfds, names, file.name, main.text) {
   titles <- c("IZ_E", "IZ_L", "TZ_E", "TZ_L")
   cols <- c(red, red, blue, blue)
   n <- length(names)
 
   colnames <- c("NAME", "X", "pch", "pos1", "pos2", "pos3", "pos4", titles)
   rfds <- toTable(0, length(colnames), length(names), colnames)
   rfds$NAME <- names
   rfds$X    <- c(1, 3, 5, 7)
   rfds$pch  <- c(19, 17, 17, 17)
   rfds$cex  <- c(1.8, 1.6, 1.6, 1.6)
   rfds$pos1 <- c(3, 3, 3, 3)
   rfds$pos2 <- c(3, 3, 3, 3)
   rfds$pos3 <- c(3, 1, 3, 3)
   rfds$pos4 <- c(1, 3, 1, 3)
   for (r in 1:length(names)) {
      rfds$IZ_E[r] <- as.numeric(round0(report.rfds[[r]][1]*100, digit=1))
      rfds$IZ_L[r] <- as.numeric(round0(report.rfds[[r]][2]*100, digit=1))
      rfds$TZ_E[r] <- as.numeric(round0(report.rfds[[r]][3]*100, digit=1))
      rfds$TZ_L[r] <- as.numeric(round0(report.rfds[[r]][4]*100, digit=1))
   }
 
   ##
   pdf(file.name, height=5, width=5.2)
   layout(matrix(c(1,2), ncol=1), widths=1, heights=c(1,1))           ## One figure each in row 1 and row 2; row 1 is 1/3 the height of row 2
   par(mar=c(1,3.6,3.6,0))
   
   plot(NULL, xlim=c(0.5, 9), ylim=c(1.5, 9.5), ylab="", xlab="", main=main.text, col=cols[1], xaxt="n", bty="n", pch=rfds$pch, cex.axis=1.2, cex.lab=1.25, cex.main=1.3)
   points(rfds$X, rfds$IZ_E, col=cols[1], pch=rfds$pch, cex=rfds$cex)
   lines(rfds$X[1:4], y=rfds$IZ_E[1:4], type="l", lwd=3, col=cols[1])
   points(rfds$X, rfds$IZ_L, col=cols[2], pch=rfds$pch, cex=rfds$cex)
   lines(rfds$X[1:4], y=rfds$IZ_L[1:4], type="l", lty=5, lwd=3, col=cols[2])
   
   text(8, rfds$IZ_E[n], "      IZ (Early)", cex=1.25, col=red) 
   text(8, rfds$IZ_L[n], "     IZ (Late)", cex=1.25, col=red)
   for (n in 1:length(names)) {
      text(rfds$X[n], rfds$IZ_E[n], rfds$IZ_E[n], col=cols[1], pos=rfds$pos1[n], cex=1.25)
      text(rfds$X[n], rfds$IZ_L[n], rfds$IZ_L[n], col=cols[2], pos=rfds$pos2[n], cex=1.25)
   }
   
   ##
   #axis(side=2, at=seq(60, 100, by=20), labels=c(60, 80, 100), cex.axis=1.1)
   #legend("top", "Normal                              ", col="black", bty="n", pt.cex=1.4, pch=1, horiz=T, cex=1.2)
   #legend("topright", "Tumour                    ", col="black", bty="n", pt.cex=1.2, pch=2, horiz=T, cex=1.2)
   
   ##
   mtext("[%]            ", side=2, line=2.6, cex=1.25)
   
   ##
   par(mar=c(5,3.6,0,0))
   plot(NULL, xlim=c(0.5, 9), ylim=c(1.5, 8.9), ylab="", xlab="", col=cols[3], xaxt="n", bty="n", pch=rfds$pch, cex.axis=1.2, cex.lab=1.25, cex.main=1.3)
   points(rfds$X, rfds$TZ_E, col=cols[3], pch=rfds$pch, cex=rfds$cex)
   lines(rfds$X[1:4], y=rfds$TZ_E[1:4], type="l", lwd=3, col=cols[3])
   points(rfds$X, rfds$TZ_L, col=cols[4], pch=rfds$pch, cex=rfds$cex)
   lines(rfds$X[1:4], y=rfds$TZ_L[1:4], type="l", lty=5, lwd=3, col=cols[4])
   
   text(8, rfds$TZ_E[n], "      TZ (Early)", cex=1.25, col=blue) 
   text(8, rfds$TZ_L[n], "     TZ (Late)", cex=1.25, col=blue)
   #abline(v=4, lty=5, lwd=1, col="black")
   #abline(v=8,  lty=5, lwd=1, col="black")
   for (n in 1:length(names)) {
      text(rfds$X[n], rfds$TZ_E[n], rfds$TZ_E[n], col=cols[3], pos=rfds$pos3[n], cex=1.25)
      text(rfds$X[n], rfds$TZ_L[n], rfds$TZ_L[n], col=cols[4], pos=rfds$pos4[n], cex=1.25)
   }
   #axis(side=2, at=seq(0, 20, by=10), labels=c(0, 10, 20), cex.axis=1.1)
 
   axis(side=1, at=seq(1, 7, by=2), labels=names[1:4], cex.axis=1.25)
   axis(side=1, at=seq(1, 7, by=2), labels=c("n=92", "n=101", "n=56", "n=96"), line=1.2, col=NA, cex.axis=1.25)
   mtext("                        Frequency", side=2, line=2.6, cex=1.25)
   dev.off()
}

###
## Distribution of expressed genes
plotReportEG <- function(report.rfds, names, file.name, main.text) {
   titles <- c("TTR", "CTR_IZ", "CTR_TZ")
   cols <- c("black", red, blue)
   n <- length(names)
 
   colnames <- c("NAME", "X", "pch", "pos1", "pos2", "pos3", titles)
   rfds <- toTable(0, length(colnames), length(names), colnames)
   rfds$NAME <- names
   rfds$X    <- c(1, 3, 5, 7)
   rfds$pch  <- c(19, 17, 17, 17)
   rfds$cex  <- c(1.8, 1.6, 1.6, 1.6)
   rfds$pos1 <- c(3, 3, 3, 3)
   rfds$pos2 <- c(3, 3, 3, 3)
   rfds$pos3 <- c(3, 3, 3, 3)
   for (r in 1:length(names)) {
      rfds$TTR[r]   <- as.numeric(round0(report.rfds[[r]][1]*100, digit=1))
      rfds$CTR_IZ[r] <- as.numeric(round0(report.rfds[[r]][2]*100, digit=1))
      rfds$CTR_TZ[r] <- as.numeric(round0(report.rfds[[r]][3]*100, digit=1))
   }
 
   ##
   pdf(file.name, height=5, width=5.2)
   layout(matrix(c(1,2), 2, 1), widths=1, heights=c(1,1))           ## One figure each in row 1 and row 2; row 1 is 1/3 the height of row 2
   par(mar=c(1,3.6,3.6,0))
   plot(NULL, xlim=c(0.5, 9.1), ylim=c(55, 115), ylab="", main=main.text, col=cols[1], xaxt="n", yaxt="n", bty="n", pch=rfds$pch, cex.axis=1.25, cex.lab=1.3, cex.main=1.35)
   points(rfds$X, rfds$TTR, col=cols[1], pch=rfds$pch, cex=rfds$cex)
   lines(x=rfds$X[2:4], y=rfds$TTR[2:4], type="l", lwd=3, col=cols[1])
   #lines(x=rfds$X[5:6], y=rfds$TTR[5:6], type="l", lwd=3, col=cols[1])
 
   text(8, rfds$TTR[n], "   TTR ", cex=1.3, col="black", font=2)
   #abline(v=4, lty=5, lwd=1, col="black")
   #abline(v=8,  lty=5, lwd=1, col="black")
   for (n in 2:length(names))
      text(rfds$X[n], rfds$TTR[n], rfds$TTR[n], col=cols[1], pos=rfds$pos1[n], cex=1.3)
 
   ##
   axis(side=2, at=seq(60, 100, by=20), labels=c(60, 80, 100), cex.axis=1.25)
   legend("top", "Normal                              ", col="white", bty="n", pch=1, pt.cex=2, horiz=T, cex=1.3, text.col="white")
   legend("topright", "Tumour                    ", col="black", bty="n", pch=2, pt.cex=1.6, horiz=T, cex=1.3)
   #legend("top", "Primary bulks                            ", col="black", bty="n", pt.cex=1, pch=2, horiz=T, cex=1.2)
   #legend("topright", "Cell lines                ", col="black", bty="n", pt.cex=1, pch=0, horiz=T, cex=1.2)
   mtext("[%]          ", side=2, line=2.6, cex=1.3)
 
   ##
   par(mar=c(5,3.6,0,0))
   plot(NULL, xlim=c(0.5, 9.1), ylim=c(8, 17.5), ylab="", xlab="", col=cols[2], xaxt="n", bty="n", pch=rfds$pch, cex.axis=1.25, cex.lab=1.3, cex.main=1.35)
   points(rfds$X, rfds$CTR_TZ, col=cols[3], pch=rfds$pch, cex=rfds$cex)
   lines(rfds$X[2:4], y=rfds$CTR_TZ[2:4], type="l", lwd=3, col=cols[3])
   #lines(rfds$X[5:6], y=rfds$CTR_L[5:6], type="l", lwd=3, col=cols[3])
 
   points(rfds$X, rfds$CTR_IZ, col=cols[2], pch=rfds$pch, cex=rfds$cex)
   lines(rfds$X[2:4], y=rfds$CTR_IZ[2:4], type="l", lwd=3, col=cols[2])
   #lines(rfds$X[5:6], y=rfds$CTR_E[5:6], type="l", lwd=3, col=cols[2])

   text(8, rfds$CTR_IZ[n], "      IZ ", cex=1.3, col=red, pos=1, font=2) 
   text(8, rfds$CTR_TZ[n], "     TZ ", cex=1.3, col=blue, pos=3, font=2)
   #abline(v=4, lty=5, lwd=1, col="black")
   #abline(v=8,  lty=5, lwd=1, col="black")
   for (n in 2:length(names)) {
      text(rfds$X[n], rfds$CTR_IZ[n], rfds$CTR_IZ[n], col=cols[2], pos=rfds$pos2[n], cex=1.3)
      text(rfds$X[n], rfds$CTR_TZ[n], rfds$CTR_TZ[n], col=cols[3], pos=rfds$pos3[n], cex=1.3)
   }
   #axis(side=2, at=seq(0, 20, by=10), labels=c(0, 10, 20), cex.axis=1.1)
 
   ##
   #axis(side=1, at=seq(1, 1, by=2), labels=names[1], cex.axis=1.1)
   #axis(side=1, at=seq(1, 1, by=2), labels=c("n=92"), line=1.2, col=NA, cex.axis=1.1)
   axis(side=1, at=seq(3, 7, by=2), labels=names[2:4], cex.axis=1.25)
   axis(side=1, at=seq(3, 7, by=2), labels=c("n=70", "n=53", "n=71"), line=1.2, col=NA, cex.axis=1.3)
   #axis(side=1, at=seq(9, 11, by=2), labels=names[5:6], cex.axis=1.1)
   #axis(side=1, at=seq(9, 11, by=2), labels=c("n=8", "n=14"), line=1.2, col=NA, cex.axis=1.1)
 
   mtext("                      Frequency", side=2, line=2.6, cex=1.3)
   dev.off()
}

plotReportEGEL <- function(report.rfds, names, file.name, main.text) {
   titles <- c("IZ_E", "IZ_L", "TZ_E", "TZ_L")
   cols <- c(red, red, blue, blue)
   n <- length(names)
 
   colnames <- c("NAME", "X", "pch", "pos1", "pos2", "pos3", "pos4", titles)
   rfds <- toTable(0, length(colnames), length(names), colnames)
   rfds$NAME <- names
   rfds$X    <- c(1, 3, 5, 7)
   rfds$pch  <- c(19, 17, 17, 17)
   rfds$cex  <- c(1.8, 1.6, 1.6, 1.6)
   rfds$pos1 <- c(3, 3, 3, 3)
   rfds$pos2 <- c(3, 3, 3, 3)
   rfds$pos3 <- c(3, 3, 3, 3)
   rfds$pos4 <- c(3, 3, 3, 3)
   for (r in 1:length(names)) {
      rfds$IZ_E[r] <- as.numeric(round0(report.rfds[[r]][1]*100, digit=1))
      rfds$IZ_L[r] <- as.numeric(round0(report.rfds[[r]][2]*100, digit=1))
      rfds$TZ_E[r] <- as.numeric(round0(report.rfds[[r]][3]*100, digit=1))
      rfds$TZ_L[r] <- as.numeric(round0(report.rfds[[r]][4]*100, digit=1))
   }
 
   ##
   pdf(file.name, height=5, width=5.2)
   layout(matrix(c(1,2), ncol=1), widths=1, heights=c(1,1))           ## One figure each in row 1 and row 2; row 1 is 1/3 the height of row 2
   par(mar=c(0,3.6,4.6,0))
 
   plot(NULL, xlim=c(0.5, 9), ylim=c(-1, 15), ylab="", xlab="", main=main.text, col=cols[1], xaxt="n", bty="n", pch=rfds$pch, cex.axis=1.2, cex.lab=1.25, cex.main=1.3)
   points(rfds$X, rfds$IZ_E, col=cols[1], pch=rfds$pch, cex=rfds$cex)
   lines(rfds$X[2:4], y=rfds$IZ_E[2:4], type="l", lwd=3, col=cols[1])
   points(rfds$X, rfds$IZ_L, col=cols[2], pch=rfds$pch, cex=rfds$cex)
   lines(rfds$X[2:4], y=rfds$IZ_L[2:4], type="l", lty=5, lwd=3, col=cols[2])
 
   text(8, rfds$IZ_E[n], "      IZ (Early)", cex=1.25, col=cols[1]) 
   text(8, rfds$IZ_L[n], "     IZ (Late)", cex=1.25, col=cols[2])
   for (n in 2:length(names)) {
      text(rfds$X[n], rfds$IZ_E[n], rfds$IZ_E[n], col=cols[1], pos=rfds$pos1[n], cex=1.25)
      text(rfds$X[n], rfds$IZ_L[n], rfds$IZ_L[n], col=cols[2], pos=rfds$pos2[n], cex=1.25)
   }
 
   ##
   #axis(side=2, at=seq(60, 100, by=20), labels=c(60, 80, 100), cex.axis=1.1)
   #legend("top", "Normal                              ", col="black", bty="n", pt.cex=1.4, pch=1, horiz=T, cex=1.2)
   #legend("topright", "Tumour                    ", col="black", bty="n", pt.cex=1.2, pch=2, horiz=T, cex=1.2)
 
   ##
   mtext("[%]     ", side=2, line=2.6, cex=1.25)
 
   ##
   par(mar=c(5,3.6,0,0))
   plot(NULL, xlim=c(0.5, 9.1), ylim=c(-1, 14), ylab="", xlab="", col=cols[3], xaxt="n", yaxt="n", bty="n", pch=rfds$pch, cex.axis=1.2, cex.lab=1.25, cex.main=1.3)
   points(rfds$X, rfds$TZ_E, col=cols[3], pch=rfds$pch, cex=rfds$cex)
   lines(rfds$X[2:4], y=rfds$TZ_E[2:4], type="l", lwd=3, col=cols[3])
   points(rfds$X, rfds$TZ_L, col=cols[4], pch=rfds$pch, cex=rfds$cex)
   lines(rfds$X[2:4], y=rfds$TZ_L[2:4], type="l", lty=5, lwd=3, col=cols[4])
 
   text(8, rfds$TZ_E[n], "       TZ (Early)", cex=1.25, col=blue) 
   text(8, rfds$TZ_L[n], "      TZ (Late)", cex=1.25, col=blue)
   #abline(v=4, lty=5, lwd=1, col="black")
   #abline(v=8,  lty=5, lwd=1, col="black")
   for (n in 2:length(names)) {
      text(rfds$X[n], rfds$TZ_E[n], rfds$TZ_E[n], col=cols[3], pos=rfds$pos3[n], cex=1.25)
      text(rfds$X[n], rfds$TZ_L[n], rfds$TZ_L[n], col=cols[4], pos=rfds$pos4[n], cex=1.25)
   }
   #axis(side=2, at=seq(0, 20, by=10), labels=c(0, 10, 20), cex.axis=1.1)
 
   axis(side=2, at=seq(0, 10, by=5), labels=c(0, 5, 10), cex.axis=1.2)
   axis(side=1, at=seq(3, 7, by=2), labels=names[2:4], cex.axis=1.2)
   axis(side=1, at=seq(3, 7, by=2), labels=c("n=70", "n=53", "n=71"), line=1.2, col=NA, cex.axis=1.25)
   mtext("                        Frequency", side=2, line=2.6, cex=1.25)
   dev.off()
}

plotSNR <- function(n, snr, file.name, main.text, xlab.text, ylab.text, col, pos) {
   unit <- (max(snr) - min(snr))/10
   xlim <- c(min(snr) - unit, max(snr) + unit)
   unit <- (max(n) - min(n))/10
   ylim <- c(min(n) - unit, max(n) + unit)
 
   pdf(paste0(file.name, ".pdf"), height=6, width=6)
   plot(n ~ snr, ylim=ylim, xlim=xlim, ylab="", xlab=xlab.text, main=main.text[1], col=col, pch=c(19, 17, 17, 17), cex=c(2.1, 2, 2, 2), cex.axis=1.7, cex.lab=1.7, cex.main=1.8)

   samples <- c("SCLC-NL", "SCLC", "NBL", "CLL")
   text(snr, n, samples, col=col, pos=c(3,3,3,3), cex=1.75)
   
   lm.fit <- lm(n ~ snr)
   abline(lm.fit, col=col, lwd=3)
 
   cor <- cor.test(n, snr, method="spearman", exact=F)
   legend(pos, c(paste0("rho = ", round0(cor[[4]], digits=1)), paste0("p-value = ", scientific(cor[[3]], digits=1))), text.col=col, bty="n", cex=1.75)

   mtext(ylab.text, side=2, line=2.85, cex=1.7)
   #mtext(main.text[2], cex=1.2, line=0.3)
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
   adjustcolor.blue <- lightblue   ## adjustcolor(lightblue, alpha.f=0.15)
   adjustcolor.orange <- orange    ## adjustcolor(orange, alpha.f=0.15)
 
   cols <- rep(lightblue, breaks)
   cols[(breaks/2 + boundary.break):breaks] <- orange
   cols[(breaks/2 - boundary.break + 1):50] <- "white"
   cols[51:(breaks/2 + boundary.break)] <- "white"
   #cols[51:breaks] <- "sandybrown"
   #cols[50:50]     <- "steelblue1"
   
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
   plot(h, main=main.text[1], ylab="Fr. [x1000]", xlab="", ylim=ylim, border="darkgray", col=cols, xaxt="n", cex.axis=1.2, cex.lab=1.3, cex.main=1.35)
   #abline(v=500, lty=5, lwd=1, col="black")
   abline(v=950, lty=5, lwd=1.5, col="black")
   abline(v=50,  lty=5, lwd=1.5, col="black")

   ##
   par(mar=c(5,4,0,1))
   hist(nrds.RT.BSTRPS$NEG, main="", ylab="Frequency", xlab=xlab.text, ylim=c(0, ymax), breaks=breaks, border="darkgray", col=cols, las=1, axes=F, cex.axis=1.2, cex.lab=1.3, cex.main=1.4)
   if (ymax < 1000) {
      axis(side=2, at=seq(0, ymax, by=250), cex.axis=1.2)
   } else if (ymax < 3000) {
      axis(side=2, at=seq(0, ymax, by=500), cex.axis=1.2)
   } else if (ymax > 20000) {
      axis(side=2, at=seq(0, ymax, by=10000), cex.axis=1.2)
   } else
      axis(side=2, at=seq(0, ymax, by=1000), cex.axis=1.2)
   axis(side=1, at=seq(0, 1000, by=250), cex.axis=1.2)
   #abline(v=500, lty=5, lwd=1, col="black")
   abline(v=950, lty=5, lwd=1.5, col="black")
   abline(v=50,  lty=5, lwd=1.5, col="black")
   ctr <- (nrow(subset(nrds.RT.BSTRPS, RFD >= 0.9)) + nrow(subset(nrds.RT.BSTRPS, RFD <= -0.9))) / nrow(nrds.RT.BSTRPS)
   text(500, ymax*4/5, paste0("|RFD| > 0.9 (", round0(ctr*100, digits=1), "%)"), cex=1.3, col="black", font=2) 
   
   mtext(main.text[2], line=4.7, cex=1.3)   ## separator(nrow(nrds.RT.BSTRPS)),
   dev.off()
}

boundary.upper <- 950   ## 500-520 breaks
boundary.lower <-  50   ## 480-500 breaks
boundary.break <-  45   ## 1 breaks each centering 500
file.name <- file.path(wd.rt.plots, paste0("hist_", base, "_rpkm_SLOPE_RFD>0.9.pdf"))
main.text <- c(paste0(BASE, " bootstrap distribution"), "")
xlab.text <- "Number of rightward forks per kb window"
plotBootstrapHist(nrds.RT.BSTRPS, file.name, main.text, xlab.text, 100, boundary.break)

# -----------------------------------------------------------------------------
# Visualisation of bootstrap re-sampling data (Histogram, RFD, and RT)
# Last Modified: 28/11/19
# -----------------------------------------------------------------------------
plotBootstrapRFD <- function(file.name, BASE, chr, xmin, xmax, nrds.RT.NRFD, bed.gc.chr, boundary.upper, boundary.lower, ext, width, kb, withUnclassified=F, gene="") {
   overlaps <- intersect(rownames(bed.gc.chr), nrds.RT.NRFD$BED)
   nrds.RT.NRFD.chr <- nrds.RT.NRFD[overlaps,]
   bed.gc.chr <- bed.gc.chr[overlaps,]

   adjustcolor.gray <- adjustcolor("darkgray", alpha.f=0.08)
     
   rights <- rownames(subset(nrds.RT.NRFD.chr, RFD >= 0.9))
   lefts  <- rownames(subset(nrds.RT.NRFD.chr, RFD <= -0.9))
   boundary.rights <- rownames(subset(subset(nrds.RT.NRFD.chr, RFD < 0.9), RFD >= 0))
   boundary.lefts <- rownames(subset(subset(nrds.RT.NRFD.chr,  RFD < 0),   RFD > -0.9))
   
   boundaries <- c(boundary.rights, boundary.lefts)
   initiations <- rownames(subset(nrds.RT.NRFD.chr[boundaries,], NRFD > 0))
   terminations <- rownames(subset(nrds.RT.NRFD.chr[boundaries,], NRFD < 0))
   unclassified <- rownames(subset(nrds.RT.NRFD.chr[boundaries,], NRFD == 0))
   unclassified <- c(unclassified, rownames(nrds.RT.NRFD.chr[boundaries,])[which(is.na(nrds.RT.NRFD.chr[boundaries,]$NRFD) == T)])
      
   if (width == 10) main.text <- paste0(BASE, " bootstrap-based replication fork directionality (RFD)")
   else main.text <- paste0(BASE, " bootstrap-based RFD")
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

   plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=c(-2, 2), xlab="", ylab=ylab.text, main=main.text, xaxt="n", yaxt="n", cex.axis=1.2, cex.lab=1.3, cex.main=1.35)
   points(bed.gc.chr$START/1E6, nrds.RT.NRFD.chr$RT, col=adjustcolor.gray, pch=16, cex=0.3)
   
   axis(side=2, at=seq(-2, 2, by=4), labels=c("\u22122", 2), cex.axis=1.2)
   axis(side=2, at=seq(-1, 1, by=1), labels=c("\u22121", 0, 1), cex.axis=1.2)
   abline(h=0, lty=5, lwd=1.5, col="black")

   ##
   if (gene != "") {
      g <- getGene(gene)
      abline(v=g$start_position/1E6, lty=5, lwd=1.5, col="#01DF01")
      abline(v=g$end_position/1E6, lty=5, lwd=1.5, col="#01DF01")
   }
   #abline(v=27548716/1E6, lty=2, lwd=1, col="#01DF01")
   #abline(v=27579868/1E6, lty=2, lwd=1, col="#01DF01")

   ## Plot cytobands (before smoothing spline)
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.5, col="lightgrey") 
   
   ## Plot smoothing spline with right-/left-leading positions
   points(bed.gc.chr[lefts,]$START/1E6,  nrds.RT.NRFD.chr[lefts,]$SPLINE, col=lightblue, pch=19, cex=0.5)
   points(bed.gc.chr[rights,]$START/1E6, nrds.RT.NRFD.chr[rights,]$SPLINE, col=orange, pch=19, cex=0.5)
   points(bed.gc.chr[terminations,]$START/1E6, nrds.RT.NRFD.chr[terminations,]$SPLINE, col=blue, pch=19, cex=0.5)
   points(bed.gc.chr[initiations,]$START/1E6,  nrds.RT.NRFD.chr[initiations,]$SPLINE,  col=red, pch=19, cex=0.5)
   if (withUnclassified && length(unclassified) != 0)
      points(bed.gc.chr[unclassified,]$START/1E6,  nrds.RT.NRFD.chr[unclassified,]$SPLINE, col="#01DF01", pch=19, cex=0.5)
   
   ## Plot legend
   legend("topright", "IZ", col=red, bty="n", pt.cex=1, lty=1, lwd=3, pch=NA, horiz=T, cex=1.3)
   legend("bottomright", "TZ", col=blue, bty="n", pt.cex=1, lty=1, lwd=3, pch=NA, horiz=T, cex=1.3)
   if (withUnclassified) legend("bottomleft", "UN", col="#01DF01", bty="n", pt.cex=1, lty=1, lwd=3, pch=NA, horiz=T, cex=1.3)
   mtext("", line=0.25, cex=1.2)   ## separator(nrow(nrds.RT.BSTRPS)),
   
   #if (width == 10) {
      legend("topleft", "Early", bty="n", text.col="black", cex=1.3, x.intersp=0)   
      legend("bottomleft", "Late", bty="n", text.col="black", cex=1.3, x.intersp=0)
   #}
   
   ###
   ## Initiate RFD plot
   par(mar=c(5.5,4,0,1))
   xlab.text <- paste0("Chromosome ", gsub("chr", "", chr), " position [Mb]")
   ylab.text <- "RFD"
   
   plot(NULL, xlim=c(xmin/1E6, xmax/1E6), ylim=c(-1.8, 1.8), xlab=xlab.text, ylab=ylab.text, main="", yaxt="n", cex.axis=1.2, cex.lab=1.3)
   axis(side=2, at=seq(-1, 1, by=1), labels=c("\u22121", 0, 1), cex.axis=1.2)
   #abline(h=0, lty=5, lwd=1, col="black")
   abline(h=0.9, lty=5, lwd=1.5, col="black")
   abline(h=-0.9, lty=5, lwd=1.5, col="black")
   
   ##
   if (gene != "") {
      g <- getGene(gene)
      abline(v=g$start_position/1E6, lty=5, lwd=1.5, col="#01DF01")
      abline(v=g$end_position/1E6, lty=5, lwd=1.5, col="#01DF01")
   }
   #abline(v=27548716/1E6, lty=2, lwd=1, col="#01DF01")
   #abline(v=27579868/1E6, lty=2, lwd=1, col="#01DF01")
   #arrows(65107831/1E6, -2, 65117867/1E6, -2, length=0.15, angle=90, lty=1, lwd=2, col="#01DF01")
   
   ## Plot cytobands (before points)
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.5, col="lightgrey")

   points(bed.gc.chr[lefts,]$START/1E6,  nrds.RT.NRFD.chr[lefts,]$RFD,  col=lightblue, pch=19, cex=0.5)
   points(bed.gc.chr[rights,]$START/1E6, nrds.RT.NRFD.chr[rights,]$RFD, col=orange, pch=19, cex=0.5)
   points(bed.gc.chr[terminations,]$START/1E6, nrds.RT.NRFD.chr[terminations,]$RFD, col=blue, pch=19, cex=0.5)
   points(bed.gc.chr[initiations,]$START/1E6,  nrds.RT.NRFD.chr[initiations,]$RFD,  col=red,  pch=19, cex=0.5)
   if (withUnclassified && length(unclassified) != 0)
      points(bed.gc.chr[unclassified,]$START/1E6, nrds.RT.NRFD.chr[unclassified,]$RFD, col="#01DF01", pch=19, cex=0.5)
  
   ## Plot legend
   legend("topright", "Rightward", col=orange, bty="n", pt.cex=1, lty=1, lwd=3, pch=NA, horiz=T, cex=1.3)
   legend("bottomright", "Leftward", col=lightblue, bty="n", pt.cex=1, lty=1, lwd=3, pch=NA, horiz=T, cex=1.3)
   dev.off()
}

# -----------------------------------------------------------------------------
# Report (between NL and Ts)
# Last Modified: 23/03/20
# -----------------------------------------------------------------------------
getBootstrapSummary <- function(report.sclc.nl.vs.sclc) {
   colnames <- c("Conserved", "Preserved", "Fold")
   summary <- toTable(0, length(colnames), 3, colnames)
   rownames(summary) <- c("TTR", "IZ", "TZ")
 
   summary$Conserved[1] <- report.sclc.nl.vs.sclc$Overlapping_P[1]
   summary$Conserved[2] <- report.sclc.nl.vs.sclc$Overlapping_P[3]
   summary$Conserved[3] <- report.sclc.nl.vs.sclc$Overlapping_P[5]
 
   summary$Preserved[1] <- (report.sclc.nl.vs.sclc$RFD_N[1] - report.sclc.nl.vs.sclc$Overlapping_N[1])/report.sclc.nl.vs.sclc$Mappable[1]
   summary$Preserved[2] <- (report.sclc.nl.vs.sclc$RFD_N[3] - report.sclc.nl.vs.sclc$Overlapping_N[3])/report.sclc.nl.vs.sclc$Mappable[3]
   summary$Preserved[3] <- (report.sclc.nl.vs.sclc$RFD_N[5] - report.sclc.nl.vs.sclc$Overlapping_N[5])/report.sclc.nl.vs.sclc$Mappable[5]
 
   for (c in 1:3) {
      if (summary$Conserved[c] < summary$Preserved[c])
         summary$Fold[c] <- summary$Preserved[c]/summary$Conserved[c]
      else
         summary$Fold[c] <- summary$Conserved[c]/summary$Preserved[c]
   }
   return(summary)
}

## https://stackoverflow.com/questions/40919759/stacked-barplot-using-r-base-how-to-add-values-inside-each-stacked-bar
plotBootstrapSummary <- function(summary, file.name, main.text) {
   xlab.text <- "Percentage [%]"
   #cols <- c(adjustcolor(green, alpha.f=0.7), adjustcolor(yellow, alpha.f=0.7))
   cols <- c(green, yellow)
   summary$Conserved <- summary$Conserved*100
   summary$Preserved <- summary$Preserved*100
   summary <- t(as.matrix(summary[,-3]))
   summary <- summary[,3:1]
 
   pdf(paste0(file.name, ".pdf"), height=3.5, width=10)
   par(mar=c(5.1, 1.1, 4.1, 1.1), xpd=TRUE)
   #bp <- barplot(summary, col=cols, xlim=c(0, 90), xlab=xlab.text, names.arg=c("","",""), main=main.text[1], horiz=T, cex.axis=1.5, cex.lab=1.6, cex.main=1.7)
   bp <- barplot(summary, col=cols, xlim=c(0, 75), xlab=xlab.text, names.arg=c("","",""), main=main.text[1], horiz=T, cex.axis=1.5, cex.lab=1.6, cex.main=1.7)
 
   h <- summary
   h[1,] <- as.numeric(round0(summary[1,], digits=1))
   h[2,] <- as.numeric(round0(summary[2,], digits=1))
   H <- apply(summary, 2L, cumsum) - summary/2
   text(H, rep(bp, each = nrow(H)), labels=h, cex=1.5, col="black")
 
   legend("bottomright", c("Shared", "Specific"), pt.cex=2.5, cex=1.6, fill=cols, horiz=T, bty="n")
   #mtext(ylab.text, side=2, line=2.75, cex=1.6)
   dev.off()
}

plotBootstrapSummaryTotal <- function(summary, file.name, main.text) {
   xlab.text <- "Percentage [%]"
   #cols <- c(adjustcolor(green, alpha.f=0.7), adjustcolor(yellow, alpha.f=0.7))
   cols <- c(green, yellow)
   summary$Conserved <- summary$Conserved*100
   summary$Preserved <- summary$Preserved*100
   summary <- as.data.frame(t(summary[,-3]))
   summary <- summary[,3:1]
   summary$Total <- 0
   summary$Total[1] <- sum(summary[1,])
   summary$Total[2] <- sum(summary[2,])
   summary <- as.matrix(summary[,4])
 
   pdf(paste0(file.name, "_total.pdf"), height=3.5, width=4)
   par(mar=c(5.1, 1.1, 4.1, 1.1), xpd=TRUE)
   bp <- barplot(summary, col=cols, xlim=c(0, 100), xlab=xlab.text, names.arg="", main=main.text[1], horiz=T, cex.axis=1.5, cex.lab=1.6, cex.main=1.7)
 
   h <- summary
   h[1,] <- as.numeric(round0(summary[1,], digits=1))
   h[2,] <- as.numeric(round0(summary[2,], digits=1))
   H <- apply(summary, 2L, cumsum) - summary/2
   text(H, rep(bp, each = nrow(H)), labels=h, cex=1.5, col="black")
 
   #legend("bottomright", rownames(summary), cex=1.6, fill=cols, horiz=T, bty="n")
   #mtext(ylab.text, side=2, line=2.75, cex=1.6)
   dev.off()
}
