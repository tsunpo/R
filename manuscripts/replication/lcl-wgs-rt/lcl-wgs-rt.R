# =============================================================================
# Manuscript   : 
# Chapter      : Chromosome replication timing of the human genome
# Name         : manuscripts/replication/lcl-wgs-rt.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 16/06/19; 25/02/19; 16/05/18
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"            ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"                 ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"                  ## tpyang@localhost

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
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
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
save(nrds, file=file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d.rt.log2s_", "s-g1", ".RData")))
#load(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d.rt.log2s_", "s-g1", ".RData")))
# [1] 2684771 (All)
# [1] 2654359 (M; RT != NA)
# [1] 2582940 (D)
# [1] 2582940 - 22
nrds.lcl <- nrds

load(file.path(wd.rt.data, paste0("lcl", "_", "rpkm", ".gc.cn.d.rt.log2s_", "s-g1", ".RData")))
nrds.lcl <- nrds

ymax <- 0.6
ymin <- 0.15
for (c in 2:2) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]
   #lcl.rt.chr <- subset(lcl.rt, CHR == chr)   ## Koren 2012
   
   ## Plot RT
   main.text <- paste0(BASE, " S/G1 read depth ratio between S phase (n=", n1, ") and G1 phase (n=", n0, ") cells")  
   file.name <- file.path(wd.rt.plots, paste0("RT_", BASE, "_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, "_1kb"))   
   plotRT(file.name, main.text, chr, NA, NA, nrds.chr, bed.gc.chr, c(red, blue, green), c("S phase", "G1 phase"), c(red, blue), c("S", "G1"), "png", width=10, peaks=c(), ylim=c(ymin, ymax), lcl.rt.chr=NULL, nrds.lcl.chr=NULL, legend="bottomright")
   #plotRT(file.name, main.text, chr, NA, NA, nrds.chr, bed.gc.chr, c(red, blue, green), c("S phase", "G1 phase"), c(red, blue), c("S", "G1"), "png", width=10, peaks=c(), ylim=c(ymin, ymax), lcl.rt.chr=NULL, nrds.lcl.chr=NULL, legend="bottomright")
   #plotRT(file.name, main.text, chr, NA, NA, nrds.chr, bed.gc.chr, c(red, blue, green), c("S phase", "G1 phase"), c(red, blue), c("S", "G1"), "png", width=10, peaks=c(), ylim=c(ymin, ymax), lcl.rt.chr=lcl.rt.chr, nrds.lcl.chr=NULL, legend="bottomright")

   ## chr2
   main.text <- paste0(BASE, " S/G1 replication timing")  
   file.name <- file.path(wd.rt.plots, paste0("RT_", BASE, "_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ""))
   plotRTAbstract(file.name, main.text, chr,  13000000,  17000000, nrds.chr, bed.gc.chr, c(red, blue, green), c("S", "G1"), c(red, blue), c("S", "G1"), "png", width=5, peaks=c(), ylim=c(ymin, ymax), NULL, NULL, 0.1)

   
   plotRTAbstract(file.name, main.text, chr,  69500000,  81000000, nrds.chr, bed.gc.chr, c(red, blue, green), c("S", "G1"), c(red, blue), c("S", "G1"), "png", width=5, peaks=c(), ylim=c(ymin, ymax), NULL, NULL, 0.05)
   plotRTAbstract(file.name, main.text, chr,  70000000,  80000000, nrds.chr, bed.gc.chr, c(red, blue, green), c("S", "G1"), c(red, blue), c("S", "G1"), "png", width=5, peaks=c(), ylim=c(ymin, ymax), NULL, NULL, 0.05)
   plotRTAbstract(file.name, main.text, chr, 160000000, 170000000, nrds.chr, bed.gc.chr, c(red, blue, green), c("S", "G1"), c(red, blue), c("S", "G1"), "png", width=5, peaks=c(), ylim=c(ymin, ymax), NULL, NULL, 0.05)
}

# -----------------------------------------------------------------------------
# RD vs RT (RDS and SPR)
# Last Modified: 11/07/19; 27/05/19
# -----------------------------------------------------------------------------
sprs <- getSPR(nrds, bed.gc)
save(sprs, file=file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-s-g1_spline_spearman.RData")))
writeTable(sprs, file=file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-s-g1_spline_spearman.txt")), colnames=T, rownames=F, sep="\t")
#load(file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-s-g1_spline_spearman.RData")))
sprs.lcl <- sprs

for (c in 2:2) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   
   nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]
   nrds.chr.T  <- setSpline(nrds.chr, bed.gc.chr, "T")
   nrds.chr.N  <- setSpline(nrds.chr, bed.gc.chr, "N")
   nrds.chr.RT <- setSpline(nrds.chr, bed.gc.chr, "RT")

   ## Figure 1
   xlab.text <- "RT [log2]"
   ylab.text <- "Read depth [RPKM]"
   main.text <- c(paste0("Correlation (", "Chr", c, ")"))   #, paste0("rho = ", round0(sprs$cor[c], digits=2), " (S vs. G1)"))
   file.name <- file.path(wd.rt.plots, "chrs", paste0("RD-vs-RT_LCL-S-G1_chr", c, ""))
   plotRD2vsRT(nrds.chr.T$SPLINE, nrds.chr.N$SPLINE, nrds.chr.RT$SPLINE, file.name, main.text, ylab.text, xlab.text, c(red, blue), c("S", "G1"), method="spearman")
   
   ## SFigure 1
   xlab.text <- "RT [log2]"
   ylab.text <- "Read depth [RPKM]"
   #main.text <- c(paste0("Chr", c))   #, paste0("rho = ", round0(sprs$cor[c], digits=2), " (S vs. G1)"))
   #file.name <- file.path(wd.rt.plots, "chrs", paste0("RD-vs-RT_LCL-S-G1_chr", c, "_spline_spearman"))
   #plotRD2vsRT(nrds.chr.T$SPLINE, nrds.chr.N$SPLINE, nrds.chr.RT$SPLINE, file.name, main.text, ylab.text, xlab.text, c(red, blue), c("S", "G1"), method="spearman")
   
   main.text <- c(paste0("Chr", c), "")
   file.name <- file.path(wd.rt.plots, "chrs", paste0("RD-vs-RT_LCL-G1_chr", c, "_spline_spearman"))
   plotRDvsRT(nrds.chr.N$SPLINE, nrds.chr.RT$SPLINE, file.name, main.text, ylab.text, xlab.text, c(blue, adjustcolor(lighterblue, alpha.f=0.01)), c("S", "G1"), method="spearman")
   
   main.text <- c(paste0("Chr", c), "")
   file.name <- file.path(wd.rt.plots, "chrs", paste0("RD-vs-RT_LCL-S_chr", c, "_spline_spearman"))
   plotRDvsRT(nrds.chr.T$SPLINE, nrds.chr.RT$SPLINE, file.name, main.text, ylab.text, xlab.text, c(red, adjustcolor(lighterred, alpha.f=0.01)), c("S", "G1"), method="spearman")
}

###
##
file.name <- file.path(wd.rt.plots, "RD2_LCL-S-G1-vs-LCL-S-G1")
main.text <- c("Correlation (Chr1-22)", "Spearman's rho")
plotRD2(sprs, file.name, main.text, 0.1, 0.75)




## Replication timing skew (RTS)
file.name <- file.path(wd.rt.plots, "RTS_LCL-S-G1_spline_spearman_chr2")
main.text <- c("Replication timing skew", "RTS = (E-L)/(E+L)")
ylab.text <- "LCL RTS"
plotRTS(sprs.lcl, file.name, main.text, c(4, 13, 17, 19), digits=3, unit=5, ylab.text, cex=1.2, chr2="0.08", offset="           ")




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
