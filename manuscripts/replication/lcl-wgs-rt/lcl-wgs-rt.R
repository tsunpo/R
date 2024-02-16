# =============================================================================
# Manuscript   : 
# Chapter      : Chromosome replication timing of the human genome
# Name         : manuscripts/replication/lcl-wgs-rt.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 16/11/23; 16/06/19; 25/02/19; 16/05/18
# =============================================================================
#wd.src <- "/projects/cangen/tyang2/dev/R"            ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"                 ## tyang2@gauss
wd.src <- "/Users/ty2/Work/dev/R"                  ## tpyang@localhost

wd.src.handbook <- file.path(wd.src, "handbook-of")   ## Required handbooks/libraries for the manuscript
handbooks <- c("Commons.R", "ReplicationTiming.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.handbook, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")     ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.bed.gc.1kb.RData"))
load(file.path(wd.src.ref, "hg19.lcl.koren.woodfine.RData"))

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 30/01/18
# -----------------------------------------------------------------------------
#wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/ty2/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE  <- "LCL"
BASE1 <- "T"
BASE0 <- "N"
PAIR1 <- "S"
PAIR0 <- "G1"
base  <- tolower(BASE)
method <- "rpkm"

wd.ngs   <- file.path(wd, BASE, "ngs/WGS")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.rt       <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

samples1 <- readTable(file.path(wd.ngs, "lcl_wgs_n7.list"), header=F, rownames=F, sep="")
samples0 <- readTable(file.path(wd.ngs, "lcl_wgs_n7.list"), header=F, rownames=F, sep="")
n1 <- length(samples1)
n0 <- length(samples0)

# -----------------------------------------------------------------------------
# Plot RD and RT (see ReplicationTiming.R)
# Last Modified: 07/08/19; 28/05/19; 14/02/19; 10/01/19; 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
nrds <- getLog2ScaledRT(wd.rt.data, base, method, BASE1, BASE0, n1, n0, chrs, bed.gc)
save(nrds, file=file.path(wd.rt.data, paste0(base, "_", method, ".cn.m.rt.log2s_", "s-g1", ".RData")))
#load(file.path(wd.rt.data, paste0(base, "_", method, ".cn.m.rt.log2s_", "s-g1", ".RData")))
# [1] 2684771 (All)
# [1] 2654359 (M; RT != NA; cn.m.rt.log2s)
# [1] 2582940 (D)
# [1] 2582940 - 22
nrds.lcl <- nrds

ymax <- 0.6
ymin <- 0.15
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]
   lcl.rt.chr <- subset(lcl.rt, CHR == chr)   ## Koren 2012
   
   ## Plot RT
   main.text <- c(paste0(BASE, " S to G1 replication timing (RT)"), "")
   file.name <- file.path(wd.rt.plots, paste0("RT_", BASE, "_", method, ".m.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, "_3_spline"))   
   plotRT(file.name, main.text, chr, NA, NA, nrds.chr, bed.gc.chr, c(red, blue), c("S phase", "G1 phase"), c(red, blue), c("S", "G1"), "png", width=13, peaks=c(), ylim=c(ymin, ymax), lcl.rt.chr=NULL, nrds.lcl.chr=NULL, legend="bottomright")
   
   file.name <- file.path(wd.rt.plots, "with-koren", paste0("RT_", BASE, "_", method, ".m.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, "_3_spline"))   
   plotRT(file.name, main.text, chr, NA, NA, nrds.chr, bed.gc.chr, c(red, blue, green), c("S phase", "G1 phase"), c(red, blue), c("S", "G1"), "png", width=13, peaks=c(), ylim=c(ymin, ymax), lcl.rt.chr=lcl.rt.chr, nrds.lcl.chr=NULL, legend="bottomright")
   #plotRT(file.name, main.text, chr, NA, NA, nrds.chr, bed.gc.chr, c(red, blue, green), c("S phase", "G1 phase"), c(red, blue), c("S", "G1"), "png", width=10, peaks=c(), ylim=c(ymin, ymax), lcl.rt.chr=lcl.rt.chr, nrds.lcl.chr=NULL, legend="bottomright")

   ## chr2
   #main.text <- paste0(BASE, " S/G1 replication timing")  
   #file.name <- file.path(wd.rt.plots, paste0("RT_", BASE, "_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ""))
   #plotRTAbstract(file.name, main.text, chr,  13000000,  17000000, nrds.chr, bed.gc.chr, c(red, blue, green), c("S", "G1"), c(red, blue), c("S", "G1"), "png", width=5, peaks=c(), ylim=c(ymin, ymax), NULL, NULL, 0.1)

   #plotRTAbstract(file.name, main.text, chr,  69500000,  81000000, nrds.chr, bed.gc.chr, c(red, blue, green), c("S", "G1"), c(red, blue), c("S", "G1"), "png", width=5, peaks=c(), ylim=c(ymin, ymax), NULL, NULL, 0.05)
   #plotRTAbstract(file.name, main.text, chr,  70000000,  80000000, nrds.chr, bed.gc.chr, c(red, blue, green), c("S", "G1"), c(red, blue), c("S", "G1"), "png", width=5, peaks=c(), ylim=c(ymin, ymax), NULL, NULL, 0.05)
   #plotRTAbstract(file.name, main.text, chr, 160000000, 170000000, nrds.chr, bed.gc.chr, c(red, blue, green), c("S", "G1"), c(red, blue), c("S", "G1"), "png", width=5, peaks=c(), ylim=c(ymin, ymax), NULL, NULL, 0.05)
}

# -----------------------------------------------------------------------------
# RD vs RT (RDS and SPR)
# Last Modified: 11/07/19; 27/05/19
# -----------------------------------------------------------------------------
sprs <- getSPR(nrds, bed.gc)
save(sprs, file=file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-s-g1_spline_spearman.RData")))
writeTable(sprs, file=file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-s-g1_spline_spearman.txt")), colnames=T, rownames=F, sep="\t")
#load(file.path("/Users/ty2/Work/uni-koeln/tyang2/LCL/analysis/replication/lcl-wgs-rt/data", paste0("rd-vs-rt_", base, "-s-g1_spline_spearman.RData")))
sprs.lcl <- sprs
#rts.lcl <- sprs
#colnames(rts.lcl)[7] <- "rts"

for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   
   nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]
   nrds.chr.T  <- setSpline(nrds.chr, bed.gc.chr, "T")
   nrds.chr.N  <- setSpline(nrds.chr, bed.gc.chr, "N")
   nrds.chr.RT <- setSpline(nrds.chr, bed.gc.chr, "RT")

   ## Figure 1
   xlab.text <- expression("RT [log" * ""[2] * "]")
   ylab.text <- "Read depth [RPKM]"
   main.text <- c(paste0("Correlation to RT (Chr", c, ")"), "")   #, paste0("rho = ", round0(sprs$cor[c], digits=2), " (S vs. G1)"))
   file.name <- file.path(wd.rt.plots, "", paste0("RD-vs-RT_LCL-S-G1_chr", c, "_Fig. 1_P_S+G1_0.01_NEW_spline0"))
   #plotSvsRT(nrds.chr.T$SPLINE, nrds.chr.N$SPLINE, nrds.chr.RT$SPLINE, file.name, main.text, ylab.text, xlab.text, c(red, adjustcolor(red.lighter, alpha.f=0.01)), c("S", "G1"), method="spearman", ylim=c(0.15, 0.6))
   #plotRD2vsRT(nrds.chr.T$SPLINE, nrds.chr.N$SPLINE, nrds.chr.RT$SPLINE, file.name, main.text, ylab.text, xlab.text, c(red, blue), c("S", "G1"), method="spearman", ylim=c(0.15, 0.6))
   plotRD2vsRT(nrds.chr.T$SPLINE, nrds.chr.N$SPLINE, nrds.chr.RT$SPLINE, file.name, main.text, ylab.text, xlab.text, c(red, blue), c(adjustcolor(red.lighter, alpha.f=0.01), adjustcolor(blue.lighter, alpha.f=0.01)), c("S", "G1"), method="spearman", ylim=c(0.15, 0.6))
   #plotRD2vsRT(nrds.chr.T$SPLINE, nrds.chr.N$SPLINE, nrds.chr$RT, file.name, main.text, ylab.text, xlab.text, c(red, blue), c(adjustcolor(red, alpha.f=0.01), adjustcolor(blue, alpha.f=0.01)), c("S", "G1"), method="spearman", ylim=c(0.15, 0.6))

   file.name <- file.path(wd.rt.plots, "", paste0("RD-vs-RT_LCL-S-G1_chr", c, "_Fig. 1_P_S+G1_0.01_NEW_T+N-RT"))
   plotRD2vsRT(nrds.chr.T$T, nrds.chr.N$N, nrds.chr.RT$RT, file.name, main.text, ylab.text, xlab.text, c(red, blue), c(adjustcolor(red.lighter, alpha.f=0.01), adjustcolor(blue.lighter, alpha.f=0.01)), c("S", "G1"), method="spearman", ylim=c(0.15, 0.6))
   
   file.name <- file.path(wd.rt.plots, "", paste0("RD-vs-RT_GC-RT_chr", c, "_spline_spearman"))
   overlaps <- intersect(rownames(bed.gc.chr), rownames(nrds.chr.RT))
   plotRDvsRT(bed.gc.chr[overlaps,]$GC, nrds.chr.RT[overlaps,]$SPLINE, file.name, main.text, ylab.text, xlab.text, c(red, adjustcolor(red.lighter, alpha.f=0.01)), c("GC", ""), method="spearman")
   
   
   ## SFigure 1
   xlab.text <- expression("RT [Log" * ""[2] * "]")
   ylab.text <- "Read depth [RPKM]"
   main.text <- c(paste0("Chr", c), "")   #, paste0("rho = ", round0(sprs$cor[c], digits=2), " (S vs. G1)"))
   file.name <- file.path(wd.rt.plots, "chrs", paste0("RD-vs-RT_LCL-S-G1_chr", c, "_spline_spearman"))
   #plotRD2vsRT(nrds.chr.T$SPLINE, nrds.chr.N$SPLINE, nrds.chr.RT$SPLINE, file.name, main.text, ylab.text, xlab.text, c(red, blue), c("S", "G1"), method="spearman")
   
   main.text <- c(paste0("Chr", c), "")
   file.name <- file.path(wd.rt.plots, "chrs", paste0("RD-vs-RT_LCL-G1_chr", c, "_spline_spearman"))
   #plotRDvsRT(nrds.chr.N$SPLINE, nrds.chr.RT$SPLINE, file.name, main.text, ylab.text, xlab.text, c(blue, adjustcolor(blue.lighter, alpha.f=0.01)), c("S", "G1"), method="spearman")
   
   main.text <- c(paste0("Chr", c), "")
   file.name <- file.path(wd.rt.plots, "chrs", paste0("RD-vs-RT_LCL-S_chr", c, "_spline_spearman"))
   #plotRDvsRT(nrds.chr.T$SPLINE, nrds.chr.RT$SPLINE, file.name, main.text, ylab.text, xlab.text, c(red, adjustcolor(red.lighter, alpha.f=0.01)), c("S", "G1"), method="spearman")
}

###
##
ylim <- c(-0.05, 0.95)
file.name <- file.path(wd.rt.plots, paste0("DNS_LCL_COR_cex=2_-0.05"))
main.text <- c("LCL read depth correlation", "")
ylab.text <- "Correlation to RT"
xlab.text <- "Chromosome"
legends <- c("S vs. RT", "|G1 vs. RT|")
cols <- c(red, blue)

pdf(paste0(file.name, ".pdf"), height=4.5, width=9.8)
par(mar=c(5.1, 4.6, 4.1, 1.5))
plot(NULL, xlim=c(1, 22), ylim=ylim, xlab=xlab.text, ylab=ylab.text, main=main.text, yaxt="n", xaxt="n", pch=19, cex.axis=1.8, cex.lab=1.9, cex.main=2)
#abline(v=7, lty=3, lwd=3)

points(abs(sprs.order$cor2) ~ rownames(sprs.order), col=cols[2], pch=19, cex=2)
lines(rownames(sprs.order), y=abs(sprs.order$cor2), lty=5, lwd=2.5, col=cols[2])

points(abs(sprs.order$cor1) ~ rownames(sprs.order), col=cols[1], pch=19, cex=2)
lines(rownames(sprs.order), y=abs(sprs.order$cor1), lty=5, lwd=2.5, col=cols[1])

axis(side=2, at=seq(-0.2, 1, by=0.2), labels=c("", 0, "", 0.4, "", 0.8, ""), cex.axis=1.8)	  
axis(side=2, at=0.4, cex.axis=1.8)	
axis(side=2, at=0.8, cex.axis=1.8)	
axis(side=1, at=seq(1, 22, by=2), labels=insilico$CHR[seq(1, 22, by=2)], cex.axis=1.8)
axis(side=1, at=seq(2, 22, by=2), labels=insilico$CHR[seq(2, 22, by=2)], cex.axis=1.8)
legend("bottomleft", legend=legends, col=cols, lty=5, lwd=3, pt.cex=2, cex=1.9)
dev.off()

###
## GC-corrected
load("/Users/ty2/Work/uni-koeln/tyang2/LCL/analysis/replication/lcl-wgs-rt/data/rd-vs-rt_lcl-s-g1_spline_spearman.RData")
sprs.order <- sprs[order(abs(sprs$cor1), decreasing=T),]
rownames(sprs.order) <- 1:22

ylim <- c(-0.23, 0.95)
file.name <- file.path(wd.rt.plots, paste0("DNS_LCL_COR_-1_3_simple"))
main.text <- c("LCL read depth correlation", "")
ylab.text <- "Correlation to RT"
xlab.text <- "Chromosome"
legends <- c("S phase", "G1 phase (Inverted)")
cols <- c(red, blue)

pdf(paste0(file.name, ".pdf"), height=4.5, width=9.8)
par(mar=c(5.1, 4.6, 4.1, 1.5))
plot(NULL, xlim=c(1, 22), ylim=ylim, xlab=xlab.text, ylab=ylab.text, main=main.text, yaxt="n", xaxt="n", pch=19, cex.axis=1.8, cex.lab=1.9, cex.main=2)
#abline(v=7, lty=3, lwd=3)

points(sprs.order$cor1 ~ rownames(sprs.order), col=cols[1], pch=19, cex=2.5)
lines(rownames(sprs.order), y=sprs.order$cor1, lty=5, lwd=2.5, col=cols[1])

points(sprs.order$cor2*-1 ~ rownames(sprs.order), col=cols[2], pch=19, cex=2.5)
lines(rownames(sprs.order), y=sprs.order$cor2 * -1, lty=5, lwd=2.5, col=cols[2])

axis(side=2, at=seq(-0.2, 1, by=0.2), labels=c("", 0, "", 0.4, "", 0.8, ""), cex.axis=1.8)	  
axis(side=2, at=0.4, cex.axis=1.8)	
axis(side=2, at=0.8, cex.axis=1.8)	
axis(side=1, at=seq(1, 22, by=2), labels=sprs.order$chr[seq(1, 22, by=2)], cex.axis=1.8)
axis(side=1, at=seq(2, 22, by=2), labels=sprs.order$chr[seq(2, 22, by=2)], cex.axis=1.8)
legend("bottomright", legend=legends, col=cols, pch=19, lty=5, lwd=3, pt.cex=3, cex=1.9)
dev.off()



# -----------------------------------------------------------------------------
# RD vs RT (RDS and SPR) -- OLD GC-corrected S and G1 and RT
# Last Modified: 15/02/24; 11/07/19; 27/05/19
# -----------------------------------------------------------------------------
###
## GC-corrected
load("/Users/ty2/Work/uni-koeln/tyang2/LCL/analysis/replication/lcl-wgs-rt/data_gc/rd-vs-rt_lcl-s-g1_spline_spearman.RData")
sprs.order <- sprs[order(abs(sprs$cor2), decreasing=T),]
rownames(sprs.order) <- 1:22

ylim <- c(-0.23, 0.95)
file.name <- file.path(wd.rt.plots, paste0("DNS_LCL-GC_COR_-1_3_simple_NEW"))
main.text <- c("LCL read depth correlation (GC-corrected)", "")
ylab.text <- "Correlation to RT"
xlab.text <- "Chromosome"
legends <- c("S phase", "G1 phase (Inverted)")
cols <- c(red, blue)

pdf(paste0(file.name, ".pdf"), height=4.5, width=9.8)
par(mar=c(5.1, 4.6, 4.1, 1.5))
plot(NULL, xlim=c(1, 22), ylim=ylim, xlab=xlab.text, ylab=ylab.text, main=main.text, yaxt="n", xaxt="n", pch=19, cex.axis=1.8, cex.lab=1.9, cex.main=2)
#abline(v=7, lty=3, lwd=3)

points(sprs.order$cor1 ~ rownames(sprs.order), col=cols[1], pch=19, cex=2.5)
lines(rownames(sprs.order), y=sprs.order$cor1, lty=5, lwd=2.5, col=cols[1])

points(sprs.order$cor2*-1 ~ rownames(sprs.order), col=cols[2], pch=19, cex=2.5)
lines(rownames(sprs.order), y=sprs.order$cor2 * -1, lty=5, lwd=2.5, col=cols[2])

axis(side=2, at=seq(-0.2, 1, by=0.2), labels=c("", 0, "", 0.4, "", 0.8, ""), cex.axis=1.8)	  
axis(side=2, at=0.4, cex.axis=1.8)	
axis(side=2, at=0.8, cex.axis=1.8)	
axis(side=1, at=seq(1, 22, by=2), labels=sprs.order$chr[seq(1, 22, by=2)], cex.axis=1.8)
axis(side=1, at=seq(2, 22, by=2), labels=sprs.order$chr[seq(2, 22, by=2)], cex.axis=1.8)
legend("bottomright", legend=legends[2:1], col=cols[2:1], pch=19, lty=5, lwd=3, pt.cex=3, cex=1.9)
dev.off()






# -----------------------------------------------------------------------------
# RD vs RT (RDS and SPR)
# Last Modified: 11/07/19; 27/05/19
# -----------------------------------------------------------------------------
sprs.gc <- getGCSPR(nrds, bed.gc)
sprs.gc.order <- sprs.gc[order(sprs.gc$cor2),]
rownames(sprs.gc.order) <- 1:22
save(sprs.gc.order, file=file.path(wd.rt.data, paste0("rd-vs-gc_", base, "-s-g1_spline_spearman.RData")))

##
ylim <- c(-0.23, 0.95)
file.name <- file.path(wd.rt.plots, paste0("DNS_LCL_GC"))
main.text <- c("LCL read depth correlation", "")
ylab.text <- "Correlation to GC"
xlab.text <- "Chromosome"
legends <- c("S phase", "G1 phase (Inverted)")
cols <- c(red, blue)

pdf(paste0(file.name, ".pdf"), height=4.5, width=9.8)
par(mar=c(5.1, 4.6, 4.1, 1.5))
plot(NULL, xlim=c(1, 22), ylim=ylim, xlab=xlab.text, ylab=ylab.text, main=main.text, yaxt="n", xaxt="n", pch=19, cex.axis=1.8, cex.lab=1.9, cex.main=2)
#abline(v=7, lty=3, lwd=3)

points(sprs.gc.order$cor1 ~ rownames(sprs.gc.order), col=cols[1], pch=19, cex=2.5)
lines(rownames(sprs.gc.order), y=sprs.gc.order$cor1, lty=5, lwd=2.5, col=cols[1])

points(sprs.gc.order$cor2*-1 ~ rownames(sprs.gc.order), col=cols[2], pch=19, cex=2.5)
lines(rownames(sprs.gc.order), y=sprs.gc.order$cor2 * -1, lty=5, lwd=2.5, col=cols[2])

axis(side=2, at=seq(-0.2, 1, by=0.2), labels=c("", 0, "", 0.4, "", 0.8, ""), cex.axis=1.8)	  
axis(side=2, at=0.4, cex.axis=1.8)	
axis(side=2, at=0.8, cex.axis=1.8)	
axis(side=1, at=seq(1, 22, by=2), labels=sprs.gc.order$chr[seq(1, 22, by=2)], cex.axis=1.8)
axis(side=1, at=seq(2, 22, by=2), labels=sprs.gc.order$chr[seq(2, 22, by=2)], cex.axis=1.8)
legend("bottomright", legend=legends[2:1], col=cols[2:1], pch=19, lty=5, lwd=3, pt.cex=3, cex=1.9)
dev.off()







sprs.gc <- toTable(0, 4, 22, c("chr", "cor0", "cor1", "cor2"))
sprs.gc$chr <- 1:22
for (c in 1:22) {
	  chr <- chrs[c]
	  bed.gc.chr <- subset(bed.gc, CHR == chr)
	  bed.gc.chr$BED <- rownames(bed.gc.chr)
	  nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]
	  nrds.chr.T  <- setSpline(nrds.chr, bed.gc.chr, "T")
	  nrds.chr.N  <- setSpline(nrds.chr, bed.gc.chr, "N")
	  bed.gc.chr <- bed.gc.chr[intersect(nrds.chr.T$BED, bed.gc.chr$BED),]
	  nrds.chr.GC <- setSpline(bed.gc.chr, bed.gc.chr, "GC")
	
	  sprs.gc$cor1[c] <- getCor(nrds.chr.T$SPLINE, nrds.chr.GC$SPLINE, method="spearman")[[4]]
	  sprs.gc$cor2[c] <- getCor(nrds.chr.N$SPLINE, nrds.chr.GC$SPLINE, method="spearman")[[4]]
	  
	  ##
	  xlab.text <- expression("GC content")
	  ylab.text <- "Read depth [RPKM]"
	  main.text <- c(paste0("Chr", c), "")
	  file.name <- file.path(wd.rt.plots, "GC", paste0("RD-vs-GC_G1_chr", c, "_spline_spearman"))
	  plotRDvsRT(nrds.chr.N$SPLINE, nrds.chr.GC$SPLINE, file.name, main.text, ylab.text, xlab.text, c(blue, adjustcolor(blue.lighter, alpha.f=0.01)), c("S", "G1"), method="spearman")
	  
	  main.text <- c(paste0("Chr", c), "")
	  file.name <- file.path(wd.rt.plots, "GC", paste0("RD-vs-GC_S_chr", c, "_spline_spearman"))
	  plotRDvsRT(nrds.chr.T$SPLINE, nrds.chr.GC$SPLINE, file.name, main.text, ylab.text, xlab.text, c(red, adjustcolor(red.lighter, alpha.f=0.01)), c("S", "G1"), method="spearman")
}





for (c in 1:22) {
	  chr <- chrs[c]
	  bed.gc.chr <- subset(bed.gc, CHR == chr)
	
	  nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]
	  nrds.chr.T  <- setSpline(nrds.chr, bed.gc.chr, "T")
	  nrds.chr.N  <- setSpline(nrds.chr, bed.gc.chr, "N")
	  nrds.chr.RT <- setSpline(nrds.chr, bed.gc.chr, "RT")
	
	  ## Figure 1
	  xlab.text <- expression("RT [log" * ""[2] * "]")
	  ylab.text <- "Read depth [RPKM]"
	  main.text <- c(paste0("Chr", c), "")   #, paste0("rho = ", round0(sprs$cor[c], digits=2), " (S vs. G1)"))
	  file.name <- file.path(wd.rt.plots, "chrs_test", paste0("RD-vs-RT_LCL-S-G1_chr", c, "_Fig. 1_P_S+G1_0.01"))
	  plotRD2vsRT(nrds.chr.T$SPLINE, nrds.chr.N$SPLINE, nrds.chr$RT, file.name, main.text, ylab.text, xlab.text, c(red, blue), c(adjustcolor(red, alpha.f=0.01), adjustcolor(blue, alpha.f=0.01)), c("S", "G1"), method="spearman", ylim=c(0.15, 0.6))

	}














# -----------------------------------------------------------------------------
# GC correlation (GCC) 
# Last Modified: 25/11/23; 11/07/19; 27/05/19
# -----------------------------------------------------------------------------
rdc <- getGCC(nrds, bed.gc)
save(rdc, file=file.path(wd.rt.data, paste0("gc-vs-rd_", base, "-s-g1_spearman.RData")))
#writeTable(sprs, file=file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-s-g1_spearman.txt")), colnames=T, rownames=F, sep="\t")
#load(file.path("/Users/ty2/Work/uni-koeln/tyang2/LCL/analysis/replication/lcl-wgs-rt/data", paste0("rd-vs-rt_", base, "-s-g1_spearman.RData")))

for (c in 1:22) {
	chr <- chrs[c]
	bed.gc.chr <- subset(bed.gc, CHR == chr)
	
	nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]
	#nrds.chr.T  <- setSpline(nrds.chr, bed.gc.chr, "T")
	#nrds.chr.N  <- setSpline(nrds.chr, bed.gc.chr, "N")
	#nrds.chr.RT <- setSpline(nrds.chr, bed.gc.chr, "RT")
	
	## Figure 1
	xlab.text <- expression("RT [log" * ""[2] * "]")
	ylab.text <- "Read depth [RPKM]"
	main.text <- c(paste0("Correlation to RT"), "")   #, paste0("rho = ", round0(sprs$cor[c], digits=2), " (S vs. G1)"))
	file.name <- file.path(wd.rt.plots,  paste0("RD-vs-RATIO_LCL-S-G1_chr", c, "_Fig. 1"))
	#plotRD2vsRT(nrds.chr.T$SPLINE, nrds.chr.N$SPLINE, nrds.chr.RT$SPLINE, file.name, main.text, ylab.text, xlab.text, c(red, blue), c("S", "G1"), method="spearman", ylim=c(0.15, 0.6))
	
	## SFigure 1
	xlab.text <- expression("RT [Log" * ""[2] * "]")
	ylab.text <- "Read depth [RPKM]"
	main.text <- c(paste0("Chr", c), "")   #, paste0("rho = ", round0(sprs$cor[c], digits=2), " (S vs. G1)"))
	file.name <- file.path(wd.rt.plots, "RDC", paste0("RD-vs-RT_LCL-S-G1_chr", c, "_SPR"))
	plotRD2vsRT(nrds.chr$T, nrds.chr$N, nrds.chr$RT, file.name, main.text, ylab.text, xlab.text, c(red, blue), c("S", "G1"), method="spearman")
	
	main.text <- c(paste0("Chr", c), "")
	file.name <- file.path(wd.rt.plots, "RDC", paste0("RD-vs-RT_LCL-G1_chr", c, "_SPR"))
	plotRDvsRT(nrds.chr$N, nrds.chr$RT, file.name, main.text, ylab.text, xlab.text, c(blue, adjustcolor(blue.lighter, alpha.f=0.01)), c("S", "G1"), method="spearman")
	
	main.text <- c(paste0("Chr", c), "")
	file.name <- file.path(wd.rt.plots, "RDC", paste0("RD-vs-RT_LCL-S_chr", c, "_SPR"))
	plotRDvsRT(nrds.chr$T, nrds.chr$RT, file.name, main.text, ylab.text, xlab.text, c(red, adjustcolor(red.lighter, alpha.f=0.01)), c("S", "G1"), method="spearman")
}

###
##
ylim <- c(-0.05, 0.95)
file.name <- file.path(wd.rt.plots, paste0("DNS_LCL_COR_cex=2_-0.05"))
main.text <- c("LCL read depth correlation", "")
ylab.text <- "Correlation to RT"
xlab.text <- "Chromosome"
legends <- c("S vs. RT", "|G1 vs. RT|")
cols <- c(red, blue)

pdf(paste0(file.name, ".pdf"), height=4.5, width=9.8)
par(mar=c(5.1, 4.6, 4.1, 1.5))
plot(NULL, xlim=c(1, 22), ylim=ylim, xlab=xlab.text, ylab=ylab.text, main=main.text, yaxt="n", xaxt="n", pch=19, cex.axis=1.8, cex.lab=1.9, cex.main=2)
#abline(v=7, lty=3, lwd=3)

points(abs(sprs.order$cor2) ~ rownames(sprs.order), col=cols[2], pch=19, cex=2)
lines(rownames(sprs.order), y=abs(sprs.order$cor2), lty=5, lwd=2.5, col=cols[2])

points(abs(sprs.order$cor1) ~ rownames(sprs.order), col=cols[1], pch=19, cex=2)
lines(rownames(sprs.order), y=abs(sprs.order$cor1), lty=5, lwd=2.5, col=cols[1])

axis(side=2, at=seq(-0.2, 1, by=0.2), labels=c("", 0, "", 0.4, "", 0.8, ""), cex.axis=1.8)	  
axis(side=2, at=0.4, cex.axis=1.8)	
axis(side=2, at=0.8, cex.axis=1.8)	
axis(side=1, at=seq(1, 22, by=2), labels=insilico$CHR[seq(1, 22, by=2)], cex.axis=1.8)
axis(side=1, at=seq(2, 22, by=2), labels=insilico$CHR[seq(2, 22, by=2)], cex.axis=1.8)
legend("bottomleft", legend=legends, col=cols, lty=5, lwd=3, pt.cex=2, cex=1.9)
dev.off()

# -----------------------------------------------------------------------------
# Read depth correlation (RDC) 
# Last Modified: 25/11/23; 11/07/19; 27/05/19
# -----------------------------------------------------------------------------
rdc <- getRDC(nrds, bed.gc)
save(rdc, file=file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-s-g1_spearman.RData")))
#writeTable(sprs, file=file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-s-g1_spearman.txt")), colnames=T, rownames=F, sep="\t")
#load(file.path("/Users/ty2/Work/uni-koeln/tyang2/LCL/analysis/replication/lcl-wgs-rt/data", paste0("rd-vs-rt_", base, "-s-g1_spearman.RData")))

for (c in 1:22) {
	chr <- chrs[c]
	bed.gc.chr <- subset(bed.gc, CHR == chr)
	
	nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]
	#nrds.chr.T  <- setSpline(nrds.chr, bed.gc.chr, "T")
	#nrds.chr.N  <- setSpline(nrds.chr, bed.gc.chr, "N")
	#nrds.chr.RT <- setSpline(nrds.chr, bed.gc.chr, "RT")
	
	## Figure 1
	xlab.text <- expression("RT [log" * ""[2] * "]")
	ylab.text <- "Read depth [RPKM]"
	main.text <- c(paste0("Correlation to RT"), "")   #, paste0("rho = ", round0(sprs$cor[c], digits=2), " (S vs. G1)"))
	file.name <- file.path(wd.rt.plots,  paste0("RD-vs-RATIO_LCL-S-G1_chr", c, "_Fig. 1"))
	#plotRD2vsRT(nrds.chr.T$SPLINE, nrds.chr.N$SPLINE, nrds.chr.RT$SPLINE, file.name, main.text, ylab.text, xlab.text, c(red, blue), c("S", "G1"), method="spearman", ylim=c(0.15, 0.6))
	
	## SFigure 1
	xlab.text <- expression("RT [Log" * ""[2] * "]")
	ylab.text <- "Read depth [RPKM]"
	main.text <- c(paste0("Chr", c), "")   #, paste0("rho = ", round0(sprs$cor[c], digits=2), " (S vs. G1)"))
	file.name <- file.path(wd.rt.plots, "RDC", paste0("RD-vs-RT_LCL-S-G1_chr", c, "_SPR"))
	plotRD2vsRT(nrds.chr$T, nrds.chr$N, nrds.chr$RT, file.name, main.text, ylab.text, xlab.text, c(red, blue), c("S", "G1"), method="spearman")
	
	main.text <- c(paste0("Chr", c), "")
	file.name <- file.path(wd.rt.plots, "RDC", paste0("RD-vs-RT_LCL-G1_chr", c, "_SPR"))
	plotRDvsRT(nrds.chr$N, nrds.chr$RT, file.name, main.text, ylab.text, xlab.text, c(blue, adjustcolor(blue.lighter, alpha.f=0.01)), c("S", "G1"), method="spearman")
	
	main.text <- c(paste0("Chr", c), "")
	file.name <- file.path(wd.rt.plots, "RDC", paste0("RD-vs-RT_LCL-S_chr", c, "_SPR"))
	plotRDvsRT(nrds.chr$T, nrds.chr$RT, file.name, main.text, ylab.text, xlab.text, c(red, adjustcolor(red.lighter, alpha.f=0.01)), c("S", "G1"), method="spearman")
}

###
##
ylim <- c(-0.05, 0.95)
file.name <- file.path(wd.rt.plots, paste0("DNS_LCL_COR_cex=2_-0.05"))
main.text <- c("LCL read depth correlation", "")
ylab.text <- "Correlation to RT"
xlab.text <- "Chromosome"
legends <- c("S vs. RT", "|G1 vs. RT|")
cols <- c(red, blue)

pdf(paste0(file.name, ".pdf"), height=4.5, width=9.8)
par(mar=c(5.1, 4.6, 4.1, 1.5))
plot(NULL, xlim=c(1, 22), ylim=ylim, xlab=xlab.text, ylab=ylab.text, main=main.text, yaxt="n", xaxt="n", pch=19, cex.axis=1.8, cex.lab=1.9, cex.main=2)
#abline(v=7, lty=3, lwd=3)

points(abs(sprs.order$cor2) ~ rownames(sprs.order), col=cols[2], pch=19, cex=2)
lines(rownames(sprs.order), y=abs(sprs.order$cor2), lty=5, lwd=2.5, col=cols[2])

points(abs(sprs.order$cor1) ~ rownames(sprs.order), col=cols[1], pch=19, cex=2)
lines(rownames(sprs.order), y=abs(sprs.order$cor1), lty=5, lwd=2.5, col=cols[1])

axis(side=2, at=seq(-0.2, 1, by=0.2), labels=c("", 0, "", 0.4, "", 0.8, ""), cex.axis=1.8)	  
axis(side=2, at=0.4, cex.axis=1.8)	
axis(side=2, at=0.8, cex.axis=1.8)	
axis(side=1, at=seq(1, 22, by=2), labels=insilico$CHR[seq(1, 22, by=2)], cex.axis=1.8)
axis(side=1, at=seq(2, 22, by=2), labels=insilico$CHR[seq(2, 22, by=2)], cex.axis=1.8)
legend("bottomleft", legend=legends, col=cols, lty=5, lwd=3, pt.cex=2, cex=1.9)
dev.off()











###
##
file.name <- file.path(wd.rt.plots, "RD2_LCL-S-G1-vs-LCL-RT")
main.text <- c("Correlation with RT (Chr1-22)", "Spearman's rho")
plotRD2(sprs, file.name, main.text, 0.1, 0.75)




## Replication timing skew (RTS)
file.name <- file.path(wd.rt.plots, "RTS_LCL-S-G1_spline_spearman_chr2_NEW_lwd=3")
main.text <- c("Replication timing skew (RTS)", "")
ylab.text <- "RTS = (E-L)/(E+L)"
plotRTS(sprs.lcl, file.name, main.text, c(4, 13, 17, 19), digits=3, unit=5, ylab.text, cex=2)

file.name <- file.path(wd.rt.plots, "RTS_Koren_16_Fig")
main.text <- c("Mitotic RT (Koren 2014)", "")
ylab.text <- "Replication timing skew"
plotRTS(sprs.lcl, file.name, main.text, c(4, 13, 16, 17), digits=3, unit=5, ylab.text, cex=2, size=6)




## S-phase progression rate (SPR)
ylab.text <- "SPR"
file.name <- file.path(wd.rt.plots, "SPR_LCL-S-G1_spline_spearman")
main.text <- c(paste0(BASE, " S-phase progression rate"), "SPR = (E-L)/(E+L)")
plotSPR(sprs, file.name, main.text, c(13, 17), digits=3, unit=5, ylab.text)

## SPR vs Read depth correlation
file.name <- file.path(wd.rt.plots, "SPR-RDC_LCL-S-G1_spline_spearman")
main.text <- c(paste0(BASE, " SPR vs. Read depths correlation"), "")
xlab.text <- "S vs. G1 [rho]"
plotSPRRDC(sprs$spr, sprs$cor, file.name, main.text, c(4, 13, 17, 19, 22), xlab.text, unit=5, ylab.text)

## SPR vs Woodfine 2004
file.name <- file.path(wd.rt.plots, "SPR-Woodfine_LCL-S-G1_spline_spearman")
main.text <- c(paste0(BASE, " SPR vs. Woodfine et al. 2004"), "")
xlab.text <- "Mean replication timing ratio"
plotSPRRDC(sprs$spr, lcl.mean$Mean, file.name, main.text, c(4, 13, 17, 19, 22), xlab.text, unit=5, ylab.text)

## Fighte 2C
file.name <- file.path(wd.rt.plots, "RD-vs-RT_LCL_spline_spearman_nasa.blue_google.red_lwd=2_cex=1.5_text.font=2_lwd=3")
main.text <- paste0("LCL read depth vs. RT")
ymin <- -0.8773492
ymax <- 0.8392611
plotRD2vsRTALL(sprs, file.name, main.text, ymin, ymax, cols=c(red, blue), cols2=c("red", "blue"), c("S", "G1"), c=2)


# -----------------------------------------------------------------------------
# Mapping replication origins
# Last Modified: 22/10/19
# -----------------------------------------------------------------------------










# -----------------------------------------------------------------------------
# Plot RD and RT (see ReplicationTiming.R)
# Last Modified: 04/08/19; 28/05/19; 14/02/19; 10/01/19; 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
rts <- toTable(NA, 4, 0, c("BED", "T", "N", "RT"))
for (c in 1:22) {
   chr <- chrs[c]
 
   rts.chr <- readTable(file.path(wd.rt.data, paste0(base, "_rpkb.gc.cn.d.rt_", chr, "_", BASE1, "-", BASE0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rts.chr <- setScaledRT(rts.chr, pseudocount=0.01, recaliRT=T, scaledRT=F)
   rts <- rbind(rts, rts.chr)
}
rts$RT <- scale(rts$RT)

for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
 
   nrds.chr.rt <- readTable(file.path(wd.rt.data, paste0(base, "_rpkb.gc.cn.d.rt_", chr, "_", BASE1, "-", BASE0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   nrds.chr.rt <- setScaledRT(nrds.chr.rt, pseudocount=0, recaliRT=T, scaledRT=T)

   ## Plot RT
   main.text <- paste0(BASE, " S/G1 read depth ratio between S phase (n=", n1, ") and G1 phase (n=", n0, ") cells")  
   file.name <- file.path(wd.rt.plots, paste0("RT_", base, "_rpkb.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, "_scale"))   
   plotRT(file.name, main.text, chr, NA, NA, nrds.chr.rt, bed.gc.chr, c("red", "blue"), c("S phase", "G1 phase"), c(adjustcolor.red, adjustcolor.blue), c("S", "G1"), "png", width=10, peaks=c(), NA, NA, 3, 3, isKoren=T)
}
# > 9.5 - 7.25
# [1] 2.25

#spline <- smooth.spline(x=bed.gc.chr$START, y=rpkms.chr.rt$RT)
#idx.1 <- which(spline$x > 73353001)
#idx.2 <- which(spline$x < 75353001)
#idx <- intersect(idx.1, idx.2)

#which(spline$y == max(spline$y[idx]))
#spline$x[72636]

# -----------------------------------------------------------------------------
# Plot RT for individual genes (see ReplicationTiming.R)
# Last Modified: 26/04/19; 14/02/19; 10/01/19; 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
genes <- c("MYC", "MYCN", "MYCL1")
genes <- c("NRG1", "DDX55", "KNTC1")
genes <- c("CCAR1", "SMARCD1")

plotWholeChr <- T
ranges <- c(50000, 500000, 5000000)
for (g in 1:length(genes)) {
   gene  <- getGene(genes[g])
   chr   <- gene$chromosome_name
   start <- gene$start_position
   end   <- gene$end_position
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   
   rpkms.chr.rt <- readTable(file.path(wd.rt.data, paste0(base, "_rpkb.gc.cn.d.rt_", chr, "_", BASE, "-", BASE, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, scaledRT=T) 

   ## RD & RT 
   main.text <- paste0(BASE, " S/G1 read depth ratio between S phase (n=", n1, ") and G1 phase (n=", n0, ") cells")   
   file.name <- file.path(wd.rt.plots, "genes", paste0("RT_", base, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, "_", genes[g]))
   if (plotWholeChr)
      plotRT(file.name, main.text, chr, NA, NA, rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("S phase", "G1 phase"), c(adjustcolor.gray, adjustcolor.gray), c("S", "G1"), "png", width=10, peaks=c(start, end), 7.75, 9, 3, 3)

   for (r in 1:length(ranges))
      plotRT(file.name, paste0(BASE, " S/G1 read depth ratio"), chr, start-ranges[r],	end+ranges[r], rpkms.chr.rt, bed.gc.chr, c("red", "blue"), c("S", "G1"), c(adjustcolor.red, adjustcolor.blue), c("S", "G1"), "png", width=5, peaks=c(start, end), 7.75, 9, 3, 3)
}

# -----------------------------------------------------------------------------
# RD vs RT (RDS and SPR)
# Last Modified: 11/07/19; 27/05/19
# -----------------------------------------------------------------------------
cors <- toTable(0, 7, 22, c("chr", "length", "cor", "cor1", "cor2", "intercept1", "intercept2"))
cors$chr <- 1:22

for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
 
   nrds.chr.rt <- readTable(file.path(wd.rt.data, paste0(base, "_rpkb.gc.cn.d.rt_", chr, "_", BASE1, "-", BASE0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   nrds.chr.rt <- setScaledRT(nrds.chr.rt, pseudocount=0, recaliRT=T, scaledRT=T)
   nrds.chr.rt.T  <- setSpline(nrds.chr.rt, bed.gc.chr, "T")
   nrds.chr.rt.N  <- setSpline(nrds.chr.rt, bed.gc.chr, "N")
   nrds.chr.rt.RT <- setSpline(nrds.chr.rt, bed.gc.chr, "RT")
   cors$length[c] <- nrow(nrds.chr.rt.RT)
   
   cor <- getCor(nrds.chr.rt.T$SPLINE, nrds.chr.rt.N$SPLINE, method="spearman")
   cors$cor[c] <- cor
   
   main.text <- c(paste0("LCL read depth correlation (", "Chr", c, ")"), paste0("rho = ", round0(cor, digits=2), " (S vs. G1)"))
   xlab.text <- "LCL S/G1"
   ylab.text <- "LCL read depth [log2]"
   file.name <- file.path(wd.rt.plots, "chrs", paste0("RD-vs-RT_LCL-S-G1_chr", c, "_spline_spearman"))
   plotRD2vsRT(nrds.chr.rt.T$SPLINE, nrds.chr.rt.N$SPLINE, nrds.chr.rt.RT$SPLINE, file.name, main.text, ylab.text, xlab.text, c("red", "blue"), c("S", "G1"), method="spearman")
   
   cors$cor1[c] <- getCor(nrds.chr.rt.T$SPLINE, nrds.chr.rt.RT$SPLINE, method="spearman")
   cors$cor2[c] <- getCor(nrds.chr.rt.N$SPLINE, nrds.chr.rt.RT$SPLINE, method="spearman")
   cors$intercept1[c] <- lm(nrds.chr.rt.T$SPLINE ~ nrds.chr.rt.RT$SPLINE)[[1]][1]
   cors$intercept2[c] <- lm(nrds.chr.rt.N$SPLINE ~ nrds.chr.rt.RT$SPLINE)[[1]][1]
   
   ## Read depth skew (RDS)
   cors$skew <- (cors$intercept1 - cors$intercept2) / (cors$intercept1 + cors$intercept2)
}
save(cors, file=file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-s-g1_spline_spearman.RData")))
writeTable(cors, file=file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-s-g1_spline_spearman.txt")), colnames=T, rownames=F, sep="\t")

#load(file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-s-g1_spline_spearman.RData")))
file.name <- file.path(wd.rt.plots, "RD-vs-RT_LCL_spline_spearman")
main.text <- paste0("LCL read depth vs. LCL S/G1")
ymin <- -0.8
ymax <- 0.8
plotRD2vsRTALL(cors, file.name, main.text, ymin, ymax, cols=c("red", "blue"), c("S", "G1"), c=2)

## Read depth skew (RDS)
file.name <- file.path(wd.rt.plots, "RDS_LCL-S-G1_spline_spearman")
main.text <- c(paste0(BASE, " read depth imbalance"), "Y-axis intercept")
plotRDS(cors, file.name, main.text, ymin=8, ymax=9, cols=c("red", "blue"), c("S-phase intercept (S)", "G1-phase intercept (G1)"), c(2, 13, 17), digits=3)

## S-phase progression rate (SPR)
file.name <- file.path(wd.rt.plots, "RDS-SPR_LCL-S-G1_spline_spearman")
main.text <- c(paste0(BASE, " S-phase progression rate"), "SPR = (S-G1)/(S+G1)")
plotSPR(cors, file.name, main.text, c(13, 17), digits=3, unit=5)

## SPR vs Read depth correlation (RDC)
file.name <- file.path(wd.rt.plots, "RDS-SPR-RDC_LCL-S-G1_spline_spearman")
main.text <- c(paste0(BASE, " SPR vs. Read depth correlation"), "SPR = (S-G1)/(S+G1)")
xlab.text <- "Read depth correlation [rho]"
plotSPRRDC(cors, file.name, main.text, c(4, 13, 17, 19, 21, 22), xlab.text, unit=5)

# -----------------------------------------------------------------------------
# Test
# Last Modified: 01/08/19
# -----------------------------------------------------------------------------
timings <- toTable(0, 2, 22, c("chr", "RT"))
timings$chr <- 1:22
timings.log2 <- toTable(0, 2, 22, c("chr", "RT"))
timings.log2$chr <- 1:22
for (c in 22:1) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
 
   rpkms.chr.rt <- readTable(file.path(wd.rt.data, paste0(base, "_nrd.cn.d.rt_", chr, "_", BASE1, "-", BASE0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   timings$RT[c] <- sum(rpkms.chr.rt$T)/sum(rpkms.chr.rt$N)
   
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, scaledRT=F)
   timings.log2$RT[c] <- sum(rpkms.chr.rt$T)/sum(rpkms.chr.rt$N)
}
timings <- timings[order(timings$RT, decreasing=T),]
timings.log2 <- timings.log2[order(timings.log2$RT, decreasing=T),]
