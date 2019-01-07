# =============================================================================
# Manuscript   : Using mutational strand asymmetry to study replication-transcription conflicts through space and time
# Chapter II   : Reconstruction of tissue-specific replicaiton timing profile in tumours
# Name         : manuscripts/replicaiton/sclc-wgs-rt-bstrps.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 12/11/18
# =============================================================================
#wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/re/home/tyang2/dev/R"                ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Common.R", "DifferentialExpression.R", "ReplicationTiming.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))

# -----------------------------------------------------------------------------
# Variations in bootstrapped data (Ensembl genes)
# Last Modified: 31/10/18
# -----------------------------------------------------------------------------
BASE <- "LCL"
base <- tolower(BASE)
bstrps       <- 1000
origin.upper <- 510   ## 500-505, 505-510 breaks
origin.lower <- 490   ## 490-495, 495-500 breaks
origin.break <- 2     ## 2 breaks each centering 500

#wd <- "/projects/cangen/tyang2/"             ## tyang2@cheops
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@local
wd.anlys <- file.path(wd, BASE, "analysis")
wd.rt    <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data/bstrps")
wd.rt.plots <- file.path(wd.rt, "plots/bstrps")
wd.de    <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

# -----------------------------------------------------------------------------
# Distributions in bootstrapped data
# Last Modified: 02/11/18
# -----------------------------------------------------------------------------
bed.gc.rt <- NULL
for (c in 8:8) {
   chr <- chrs[c]
 
   load(file=file.path(wd.rt.data, paste0("bed.gc.rt_", base, "_bstrps", bstrps, "_", chr, ".RData")))
   if (is.null(bed.gc.rt))
      bed.gc.rt <- bed.gc.rt.chr
   bed.gc.rt <- rbind(bed.gc.rt, bed.gc.rt.chr)

   file.name <- file.path(wd.rt.plots, paste0("hist_", base, "_bed.gc.rt_bstrps", bstrps, "_", chr, ".pdf"))
   plotHistBootstraps(bed.gc.rt.chr, file.name, paste0("Chr", c), BASE, 200, origin.break)   ## See ReplicationTiming.R
}
#save(bed.gc.rt, file=file.path(wd.rt.data, paste0("bed.gc.rt_", base,"_bstrps", bstrps, ".RData")))

file.name <- file.path(wd.rt.plots, paste0("hist_", base, "_bed.gc.rt_bstrps", bstrps, ".pdf"))
plotHistBootstraps(bed.gc.rt, file.name, "Chr1-22", BASE, 200, origin.break)                 ## See ReplicationTiming.R

# -----------------------------------------------------------------------------
# Plot RO and RT from bootstrapped data
# Last Modified: 04/11/18
# -----------------------------------------------------------------------------
for (c in 1:22) {
   c <- 17
   chr <- chrs[c]

   ## Replication origins   
   load(file=file.path(wd.rt.data, paste0("bed.gc.rt_", base, "_bstrps", bstrps, "_", chr, ".RData")))
   right.idx  <- which(bed.gc.rt.chr$RIGHT_LEADING > origin.upper)
   left.idx   <- which(bed.gc.rt.chr$RIGHT_LEADING < origin.lower)
   origin.idx <- setdiff(c(1:nrow(bed.gc.rt.chr)), c(right.idx, left.idx))
   bed.gc.chr <- bed.gc[rownames(bed.gc.rt.chr),]
   
   file.name <- file.path(wd.rt.plots, paste0(base, "_RO_bstrps1000_", chr))
   plotRO(file.name, BASE, chr, 37000000, 40000000, bed.gc.chr, bed.gc.rt.chr, right.idx, left.idx, origin.idx, "png")      ## See ReplicationTiming.R
   
   ## Replication timing
   rt.chr <- readTable(file.path(wd.rt, "data", paste0(base, "_rpkm.corr.gc.d.rt_", chr, "_", BASE, "-", BASE, "_n101-92.txt.gz")), header=T, rownames=T, sep="\t") 
   rt.chr <- rt.chr[rownames(bed.gc.rt.chr),]
   
   file.name  <- file.path(wd.rt.plots, paste0(base, "_RT_bstrps1000_", chr))
   plotRT(file.name, BASE, chr, 37000000, 40000000, rt.chr, bed.gc.chr, right.idx, left.idx, origin.idx, ymax=1.5, "png")   ## See ReplicationTiming.R
}

# -----------------------------------------------------------------------------
# Plot RO and RT for these 30 candidate genes (sclc.tx.5050.5050)
# Last Modified: 05/12/18
# -----------------------------------------------------------------------------
wd.rt.plots.5050 <- file.path(wd.rt.plots, "5050")
for (g in 1:length(sclc.tx.5050.5050)) {
   #g <- 19   #27
   #gene <- ensGene[sclc.tx.5050.5050[g],]
   gene <- ensGene["ENSG00000136997",]   ## MYC, ATM
   chr <- gene$chromosome_name
   length <- gene$end_position - gene$start_position
   #start <- gene$start_position - 500000
   #end <- gene$end_position + 500000
 
   ## Replication fork directionality (RFD)  
   load(file=file.path(wd.rt.data, paste0("bed.gc.rt_", base, "_bstrps", bstrps, "_", chr, ".RData")))
   bed.gc.rt.chr$RFD <- mapply(x = 1:nrow(bed.gc.rt.chr), function(x) as.numeric(getRFD(bed.gc.rt.chr[x, ])))   ## ADD 05/12/18
   right.idx  <- which(bed.gc.rt.chr$RIGHT_LEADING > origin.upper)
   left.idx   <- which(bed.gc.rt.chr$RIGHT_LEADING < origin.lower)
   origin.idx <- setdiff(c(1:nrow(bed.gc.rt.chr)), c(right.idx, left.idx))
   bed.gc.chr <- bed.gc[rownames(bed.gc.rt.chr),]
 
   start <- 125000000
   end   <- 135000000
   file.name <- file.path(wd.rt.plots.5050, paste0(base, "_RFD_bstrps1000_", chr))
   plotBootstrapsRFD(file.name, BASE, chr, start, end, bed.gc.chr, bed.gc.rt.chr, right.idx, left.idx, origin.idx, "png")       ## see ReplicationTiming.R
 
   ## Replication timing
   rt.chr <- readTable(file.path(wd.rt, "data", paste0(base, "_rpkm.corr.gc.d.rt_", chr, "_", "T", "-", "N", "_n7.txt.gz")), header=T, rownames=T, sep="\t") 
   rt.chr <- rt.chr[rownames(bed.gc.rt.chr),]
   rt.chr$RT <- rt.chr$MEDIAN
   
   file.name  <- file.path(wd.rt.plots.5050, paste0(base, "_RT_bstrps1000_", chr))
   plotBootstrapsRT(file.name, gene=gene$external_gene_name, BASE, chr, start, end, rt.chr, bed.gc.chr, right.idx, left.idx, origin.idx, ymax=1, "png")   ## see ReplicationTiming.R
}

# -----------------------------------------------------------------------------
# Transcription vs. replication time
# Last Modified: 08/11/18
# -----------------------------------------------------------------------------
load("/Users/tpyang/Work/uni-koeln/tyang2/SCLC/analysis/expression/kallisto/sclc-tpm-de/data/sclc_kallisto_0.43.1_tpm.gene_r5p47.RData")
#load(file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
tpm.gene      <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=T, proteinCodingOnly=F)
tpm.gene.log2 <- getLog2andMedian(tpm.gene, pseudocount=1)
# > nrow(tpm.gene.log2)
# [1] 18440

load(file=file.path(wd.rt.data, paste0("ensGene.rt_", base, "_bstrps", bstrps, ".RData")))   ## Load objects ensGene.rt.start and ensGene.rt.end
#load("/Users/tpyang/Work/uni-koeln/tyang2/LCL/analysis/replication/lcl-wgs-rt/data/bstrps/ensGene.rt_lcl_bstrps1000.RData")   ## ADD 04/12/18
ensGene.rt.tx.lcl <- getEnsGeneTxRFD(ensGene, ensGene.rt.start, ensGene.rt.end, tpm.gene.log2)
# > nrow(ensGene.rt.tx.lcl)
# [1] 15948
# > nrow(ensGene.rt.tx.sclc)
# [1] 18078

##
file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_TSS+TES_sclc"))
plotEnsGeneTxRFD(file.name, BASE, ensGene.rt.tx.lcl, origin.upper, "pdf")

leading.ratio <- log2(origin.upper/500)
lcl.tx.right.0 <- rownames(subset(ensGene.rt.tx.lcl, TSS >  leading.ratio))
lcl.tx.left.0  <- rownames(subset(ensGene.rt.tx.lcl, TSS < -leading.ratio))
lcl.tx.5050    <- setdiff(rownames(ensGene.rt.tx.lcl), c(lcl.tx.right.0, lcl.tx.left.0))
lcl.tx.inconsist.0 <- rownames(subset(ensGene.rt.tx.lcl, CONSIST < 0))   ## This may include "50/50" if cutoff for leading-count ratio change
# length(intersect(lcl.tx.5050, lcl.tx.inconsist.0))
# [1] 61

## ADD 20/11/18: Remove inconsistent genes
lcl.tx.right   <- setdiff(lcl.tx.right.0, lcl.tx.inconsist.0)
lcl.tx.left    <- setdiff(lcl.tx.left.0,  lcl.tx.inconsist.0)
lcl.tx.consist <- c(lcl.tx.right, lcl.tx.left)

length(intersect(sclc.tx.consist, lcl.tx.consist))
# [1] 14082
length(intersect(sclc.tx.right, lcl.tx.right))
# [1] 5110   ## 60.16% (5110/8494)
length(intersect(sclc.tx.right, lcl.tx.left))
# [1] 1944   ## 22.89% (1944/8494)
length(intersect(sclc.tx.left, lcl.tx.left))
# [1] 4985   ## 59.12% (4985/8432)
length(intersect(sclc.tx.left, lcl.tx.right))
# [1] 2043   ## 24.23% (2043/8432)

sclc.tx.consist.lcl <- intersect(sclc.tx.consist, lcl.tx.consist)
sclc.tx.right.lcl <- intersect(sclc.tx.right, lcl.tx.right)
sclc.tx.left.lcl  <- intersect(sclc.tx.left,  lcl.tx.left)

lcl.tx.inconsist <- setdiff(lcl.tx.inconsist.0, lcl.tx.5050)
lcl.tx.inconsist.right <- rownames(subset(ensGene.rt.tx.lcl[lcl.tx.inconsist,], TSS > leading.ratio))
lcl.tx.inconsist.left  <- rownames(subset(ensGene.rt.tx.lcl[lcl.tx.inconsist,], TSS < -leading.ratio))
testWbyEnsGeneRTTx(ensGene.rt.tx.lcl, lcl.tx.inconsist, lcl.tx.consist)
# [1] 0.05760913
## ADD 20/11/18: How many origins are inconsistent (which will be visulised in the boxplot later)
# > 7558+7305+934+151   ## Right + Left + Inconsistent + 50/50
# [1] 15948

length(intersect(sclc.tx.inconsist.0, lcl.tx.inconsist.0))
# [1] 228
length(intersect(sclc.tx.inconsist, lcl.tx.inconsist))
# [1] 223
length(intersect(sclc.tx.inconsist.right, lcl.tx.inconsist.right))
# [1] 100   ??
length(intersect(sclc.tx.inconsist.right, lcl.tx.inconsist.left))
# [1] 15
length(intersect(sclc.tx.inconsist.left, lcl.tx.inconsist.left))
# [1] 94    ??
length(intersect(sclc.tx.inconsist.left, lcl.tx.inconsist.right))
# [1] 14

sclc.tx.inconsist.lcl <- intersect(sclc.tx.inconsist, lcl.tx.inconsist)
sclc.tx.inconsist.right.lcl <- intersect(sclc.tx.inconsist.right, lcl.tx.inconsist.right)
sclc.tx.inconsist.left.lcl  <- intersect(sclc.tx.inconsist.left,  lcl.tx.inconsist.left)
sclc.tx.5050.lcl      <- intersect(sclc.tx.5050,      lcl.tx.5050)

# > length(sclc.tx.consist)
# [1] 16926
# > length(lcl.tx.consist)
# [1] 14863

sclc.tx.lcl <- c(sclc.tx.consist.lcl, sclc.tx.inconsist.lcl, sclc.tx.5050.lcl)
length(sclc.tx.consist.lcl)
# [1] 14082
length(sclc.tx.inconsist.lcl)
# [1] 223
length(sclc.tx.5050.lcl)
# [1] 1       ## 14306=14082+223+1

###
##
file.name  <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_RFD_TSS+TTS_lcl_consist.pdf"))
main.text <- "Replication fork directionality (SCLC vs. LCL)"
xlab.text <- c("RFD = (R-L)/(R+L)", "SCLC")
ylab.text <- "LCL (Koren, 2012)"
cols <- rep("gray", length(sclc.tx.consist.lcl))

pdf(file.name, height=6, width=6)
plot(ensGene.rt.tx.lcl[sclc.tx.consist.lcl,]$RFD_TSS ~ ensGene.rt.tx[sclc.tx.consist.lcl,]$RFD_TSS, xlab=xlab.text, ylab=ylab.text, main=main.text, col=cols)

points(ensGene.rt.tx.lcl[sclc.tx.right.lcl,]$RFD_TSS ~ ensGene.rt.tx[sclc.tx.right.lcl,]$RFD_TSS, col="sandybrown")
points(ensGene.rt.tx.lcl[sclc.tx.left.lcl,]$RFD_TSS  ~ ensGene.rt.tx[sclc.tx.left.lcl,]$RFD_TSS,  col="steelblue1")
abline(v=0, lty=5, lwd=0.85, col="black")
text( 0.5,  0.5, "72.44% (5,110)", cex=0.85, col="black")
text(-0.5, -0.5, "70.93% (4,985)", cex=0.85, col="black") 

#legend("bottom", c("Left-leading in SCLC", paste0("L/R (n=", length(sclc.tx.origin.o.c), ") in SCLC"), "Right-leading in SCLC"), col=c("steelblue1", "red", "sandybrown"), pch=1, cex=0.8, horiz=T)
mtext(paste0("Overlapping consist genes (n=", separator(length(sclc.tx.consist.lcl)), ")"), cex=1, line=0.4)
dev.off()

###
##
file.name  <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_TSS+TES_lcl_inconsist.pdf"))
main.text <- "Replication time (SCLC vs. LCL)"
xlab.text <- c("Leading ratio (log2)", "SCLC")
ylab.text <- c("LCL (Koren, 2012)", "Leading ratio (log2)")
cols <- rep("gray", length(sclc.tx.inconsist.lcl))

pdf(file.name, height=6, width=6)
plot(ensGene.rt.tx.lcl[sclc.tx.inconsist.lcl,]$TSS ~ ensGene.rt.tx.sclc[sclc.tx.inconsist.lcl,]$TSS, xlab=xlab.text, ylab=ylab.text, main=main.text, col=cols)

points(ensGene.rt.tx.lcl[sclc.tx.inconsist.right.lcl,]$TSS ~ ensGene.rt.tx.sclc[sclc.tx.inconsist.right.lcl,]$TSS, col="purple1")
points(ensGene.rt.tx.lcl[sclc.tx.inconsist.left.lcl,]$TSS  ~ ensGene.rt.tx.sclc[sclc.tx.inconsist.left.lcl,]$TSS,  col="purple1")

#legend("bottom", c("Left-leading in SCLC", paste0("L/R (n=", length(sclc.tx.origin.o.c), ") in SCLC"), "Right-leading in SCLC"), col=c("steelblue1", "red", "sandybrown"), pch=1, cex=0.8, horiz=T)
mtext(paste0("Overlapping inconsist genes (n=", separator(length(sclc.tx.inconsist.lcl)), ")"), cex=1, line=0.4)
dev.off()


# -----------------------------------------------------------------------------
# Inconsistent genes have lower expression levels
# Last Modified: 12/11/18
# -----------------------------------------------------------------------------
file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_length_inconsist.pdf"))
main.text <- c("Inconsistent genes in SCLC", paste0("n=", separator(length(sclc.tx.inconsist.lcl))))
plotSRCTxLength(file.name, main.text, ensGene.rt.tx.lcl, sclc.tx.inconsist.lcl, "purple1")

file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_length_right.pdf"))
main.text <- c("Right-leading genes in SCLC", paste0("n=", separator(length(sclc.tx.right.lcl))))
plotSRCTxLength(file.name, main.text, ensGene.rt.tx.lcl, sclc.tx.right.lcl, "sandybrown")

file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_length_left.pdf"))
main.text <- c("Left-leading genes in SCLC", paste0("n=", separator(length(sclc.tx.left.lcl))))
plotSRCTxLength(file.name, main.text, ensGene.rt.tx.lcl, sclc.tx.left.lcl, "steelblue1")

file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_length_5050.pdf"))
main.text <- c("Left-/Right-leading genes in SCLC", paste0("n=", separator(length(sclc.tx.50.50.lcl))))
plotSRCTxLength(file.name, main.text, ensGene.rt.tx.lcl, sclc.tx.50.50.lcl, "red")







xlab.text <- "Gene length (log10)"
ylab.text <- "Expression (log2)"

pdf(file.name, height=6, width=6)
plot(  MEDIAN ~ log10(LENGTH), data=subset(ensGene.rt.start.sclc.inconsist, RATIO > 0), ylab=ylab.text, xlab=xlab.text, main=main.text, col="sandybrown")
points(MEDIAN ~ log10(LENGTH), data=subset(ensGene.rt.start.sclc.inconsist, RATIO < 0), col="steelblue1", pch=1)
abline(lm(MEDIAN ~ log10(LENGTH), data=subset(ensGene.rt.start.sclc.inconsist, RATIO > 0)), col="sandybrown")
abline(lm(MEDIAN ~ log10(LENGTH), data=subset(ensGene.rt.start.sclc.inconsist, RATIO < 0)), col="steelblue1")
legend("topright", c("Left-leading", "Right-leading"), col=c("steelblue1", "sandybrown"), pch=1, cex=0.8, horiz=F)
dev.off()

##
xlab.text <- "Gene length (log10)"
ylab.text <- "Expression (log2)"
ymin <- min(ensGene.rt.tx$MEDIAN)
ymax <- max(ensGene.rt.tx$MEDIAN)
xmin <- min(ensGene.rt.tx$MEDIAN)
xmax <- max(ensGene.rt.tx$MEDIAN)

file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_length_inconsistent.pdf"))
main.text <- c("Right-leading genes in SCLC", paste0("n=", separator(length(sclc.tx.right))))
pdf(file.name, height=6, width=6)
plot(MEDIAN ~ log10(LENGTH), data=ensGene.rt.tx[sclc.tx.inconsist,], ylab=ylab.text, xlab=xlab.text, main=main.text, col="purple1")
abline(lm(MEDIAN ~ log10(LENGTH), data=ensGene.rt.tx[sclc.tx.inconsist,]), ylim=c(ymin, ymax), col="black")
dev.off()
cor.test(ensGene.rt.tx[sclc.tx.inconsist,]$MEDIAN, log10(ensGene.rt.tx[sclc.tx.inconsist,]$LENGTH), method="spearman", exact=F)[[4]]
cor.test(ensGene.rt.tx[sclc.tx.inconsist,]$MEDIAN, log10(ensGene.rt.tx[sclc.tx.inconsist,]$LENGTH), method="spearman", exact=F)[[3]]


file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_length_right.pdf"))
main.text <- c("Right-leading genes in SCLC", paste0("n=", separator(length(sclc.tx.right))))
pdf(file.name, height=6, width=6)
plot(MEDIAN ~ log10(LENGTH), data=ensGene.rt.tx[sclc.tx.right,], ylab=ylab.text, xlab=xlab.text, main=main.text, col="sandybrown")
abline(lm(MEDIAN ~ log10(LENGTH), data=ensGene.rt.tx[sclc.tx.right,]), col="black")
dev.off()
cor.test(ensGene.rt.tx[sclc.tx.right,]$MEDIAN, log10(ensGene.rt.tx[sclc.tx.right,]$LENGTH), method="spearman", exact=F)[[4]]
cor.test(ensGene.rt.tx[sclc.tx.right,]$MEDIAN, log10(ensGene.rt.tx[sclc.tx.right,]$LENGTH), method="spearman", exact=F)[[3]]


##
file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_length_left.pdf"))
main.text <- c("Left-leading genes in SCLC", paste0("n=", separator(length(sclc.tx.left))))
pdf(file.name, height=6, width=6)
plot(MEDIAN ~ log10(LENGTH), data=ensGene.rt.tx[sclc.tx.left, ], ylab=ylab.text, xlab=xlab.text, main=main.text, col="steelblue1")
abline(lm(MEDIAN ~ log10(LENGTH), data=ensGene.rt.tx[sclc.tx.left,]), col="black")
dev.off()
cor.test(ensGene.rt.tx[sclc.tx.left,]$MEDIAN, log10(ensGene.rt.tx[sclc.tx.left,]$LENGTH), method="spearman", exact=F)[[4]]
cor.test(ensGene.rt.tx[sclc.tx.left,]$MEDIAN, log10(ensGene.rt.tx[sclc.tx.left,]$LENGTH), method="spearman", exact=F)[[3]]


##
file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_length_5050.pdf"))
main.text <- c("Left/Right genes in SCLC", paste0("n=", separator(length(sclc.tx.left))))
pdf(file.name, height=6, width=6)
plot(MEDIAN ~ log10(LENGTH), data=ensGene.rt.tx[sclc.tx.50.50,], ylab=ylab.text, xlab=xlab.text, main=main.text, col="red")
abline(lm(MEDIAN ~ log10(LENGTH), data=ensGene.rt.tx[sclc.tx.50.50,]), col="black")
dev.off()
cor.test(ensGene.rt.tx[sclc.tx.50.50,]$MEDIAN, log10(ensGene.rt.tx[sclc.tx.50.50,]$LENGTH), method="spearman", exact=F)[[4]]
cor.test(ensGene.rt.tx[sclc.tx.50.50,]$MEDIAN, log10(ensGene.rt.tx[sclc.tx.50.50,]$LENGTH), method="spearman", exact=F)[[3]]



##
file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.tx_bstrps1000_Start+End.png"))
plotEnsGeneRTTxRO(file.name, ensGene.rt.tx, origin.upper, "png")

sclc.tx.right  <- rownames(subset(ensGene.rt.start.sclc, RIGHT_LEADING > cutoff.upper))
sclc.tx.left   <- rownames(subset(ensGene.rt.start.sclc, RIGHT_LEADING < cutoff.lower))
sclc.tx.origin <- setdiff(rownames(ensGene.rt.start.sclc), c(sclc.tx.right, sclc.tx.left))
sclc.tx.inconsist <- rownames(subset(ensGene.rt.start.sclc, SIGN < 0))
## ADD 20/11/18: Remove inconsistent genes
sclc.tx.inconsist <- setdiff(sclc.tx.inconsist, sclc.tx.origin)
sclc.tx.right     <- setdiff(sclc.tx.right, sclc.tx.inconsist)
sclc.tx.left      <- setdiff(sclc.tx.left,  sclc.tx.inconsist)
## ADD 20/11/18: How many origins are inconsistent (which will be visulised in the boxplot later)
length(intersect(sclc.tx.origin, sclc.tx.inconsist))
# [1] 45
# > 8486+8442+1098+97-45   ## Right+Left+Inconsistent+Origin-Overlaps
# [1] 18078

file.name  <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_Start+End.png"))
main.text  <- "Transcription vs. replication time in SCLC"
mtext.text <- paste0("Expressed genes (n=", separator(nrow(ensGene.rt.start.sclc)), ")")
xlab.text <- "Right-leading ratio (log2)"
ylab.text <- "log2(TPM+1)"

pdf(file.name, height=6, width=6)   #, units="in", res=300)
cols <- rep("steelblue1", nrow(ensGene.rt.start.sclc))
plot(  MEDIAN ~ RATIO, data=ensGene.rt.start.sclc, xlab=xlab.text, ylab=ylab.text, main=main.text, col=cols)
points(MEDIAN ~ RATIO, data=ensGene.rt.start.sclc[sclc.tx.right,],     col="sandybrown")
points(MEDIAN ~ RATIO, data=ensGene.rt.start.sclc[sclc.tx.inconsist,], col="purple1")
points(MEDIAN ~ RATIO, data=ensGene.rt.start.sclc[sclc.tx.origin,],    col="red")
legend("top", c("Left-leading", paste0("L/R (n=", separator(length(sclc.tx.origin)), ")"), "Inconsistent", "Right-leading"), col=c("steelblue1", "red", "purple1", "sandybrown"), pch=1, cex=0.75, horiz=T)
mtext(mtext.text, cex=1, line=0.5)
dev.off()

# -----------------------------------------------------------------------------
# Inconsistent genes have lower expression levels
# Last Modified: 12/11/18
# -----------------------------------------------------------------------------
sclc.tx.inconsist <- subset(ensGene.rt.start.sclc, SIGN == -1)$MEDIAN
sclc.tx.consist   <- subset(ensGene.rt.start.sclc, SIGN == 1)$MEDIAN
testW(sclc.tx.inconsist, sclc.tx.consist)
#[1] 0.007379227

file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_inconsistent_beeswarm.pdf"))
main.text <- c("Inconsistent genes in SCLC", paste0("n=", separator(length(sclc.tx.inconsist))))
xlab.text <- "Replication time"
ylab.text <- "log2(TPM+1)"
ymin <- min(ensGene.rt.start.sclc$MEDIAN)
ymax <- max(ensGene.rt.start.sclc$MEDIAN)
overlaps <- intersect(sclc.tx.inconsist, sclc.tx.origin)
names <- c(paste0("Consistent (n=", separator(length(sclc.tx.consist)),")"), paste0("Inconsistent (n=", separator(length(sclc.tx.inconsist)),")"))

pdf(file.name, height=6, width=3.5)   #, units="in", res=300)
boxplot( MEDIAN ~ RT, data=ensGene.rt.start.sclc.inconsist, ylim=c(ymin, ymax), ylab=ylab.text, xlab=xlab.text, names=names, main=main.text, border="purple1", outline=T)
beeswarm(MEDIAN ~ RT, data=ensGene.rt.start.sclc.inconsist[overlaps,], col="red", pch=1, add=T)
legend("topleft", paste0("L/R (n=", separator(length(overlaps)), ")"), col="red", pch=1, cex=0.8)
dev.off()



# -----------------------------------------------------------------------------
# Inconsistent genes have lower expression levels (TO-DO)
# Last Modified: 12/11/18
# -----------------------------------------------------------------------------
ensGene.rt.start.sclc.inconsist <- ensGene.rt.start.sclc[sclc.tx.inconsist,]
sclc.tx.inconsist.right <- subset(ensGene.rt.start.sclc.inconsist, RATIO > 0)$MEDIAN
sclc.tx.inconsist.left  <- subset(ensGene.rt.start.sclc.inconsist, RATIO < 0)$MEDIAN
testW(sclc.tx.inconsist.right, sclc.tx.inconsist.left)
# [1] 1.235423e-15

file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_inconsistent_beeswarm.pdf"))
main.text <- c("Inconsistent genes in SCLC", paste0("n=", separator(length(sclc.tx.inconsist))))
xlab.text <- "Replication time"
ylab.text <- "log2(TPM+1)"
ymin <- min(ensGene.rt.start.sclc$MEDIAN)
ymax <- max(ensGene.rt.start.sclc$MEDIAN)
overlaps <- intersect(sclc.tx.inconsist, sclc.tx.origin)
names <- c(paste0("L (n=", separator(length(sclc.tx.inconsist.left)),")"), paste0("R (n=", separator(length(sclc.tx.inconsist.right)),")"))

pdf(file.name, height=6, width=3.5)   #, units="in", res=300)
boxplot( MEDIAN ~ RT, data=ensGene.rt.start.sclc.inconsist, ylim=c(ymin, ymax), ylab=ylab.text, xlab=xlab.text, names=names, main=main.text, border="purple1", outline=T)
beeswarm(MEDIAN ~ RT, data=ensGene.rt.start.sclc.inconsist[overlaps,], col="red", pch=1, add=T)
legend("topleft", paste0("L/R (n=", separator(length(overlaps)), ")"), col="red", pch=1, cex=0.8)
dev.off()

##
sclc.tx.all <- ensGene.rt.start.sclc$MEDIAN
testW(sclc.tx.all, c(sclc.tx.inconsist.right, sclc.tx.inconsist.left))
# [1] 0.01134909
testW(sclc.tx.all, sclc.tx.inconsist.right)
# [1] 8.092788e-05
testW(sclc.tx.all, sclc.tx.inconsist.left)
# [1] 3.372877e-11

##
sclc.tx.inconsist.right.length <- subset(ensGene.rt.start.sclc.inconsist, RATIO > 0)$LENGTH
sclc.tx.inconsist.left.length  <- subset(ensGene.rt.start.sclc.inconsist, RATIO < 0)$LENGTH
testW(sclc.tx.inconsist.right.length, sclc.tx.inconsist.left.length)
# [1] 2.826203e-13

##
file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_length_inconsistent.pdf"))
main.text <- c("Inconsistent genes in SCLC", paste0("n=", separator(length(sclc.tx.inconsist))))
xlab.text <- "Gene length (log10)"
ylab.text <- "log2(TPM+1)"

pdf(file.name, height=6, width=6)
plot(  MEDIAN ~ log10(LENGTH), data=subset(ensGene.rt.start.sclc.inconsist, RATIO > 0), ylab=ylab.text, xlab=xlab.text, main=main.text, col="sandybrown")
points(MEDIAN ~ log10(LENGTH), data=subset(ensGene.rt.start.sclc.inconsist, RATIO < 0), col="steelblue1", pch=1)
abline(lm(MEDIAN ~ log10(LENGTH), data=subset(ensGene.rt.start.sclc.inconsist, RATIO > 0)), col="sandybrown")
abline(lm(MEDIAN ~ log10(LENGTH), data=subset(ensGene.rt.start.sclc.inconsist, RATIO < 0)), col="steelblue1")
legend("topright", c("Left-leading", "Right-leading"), col=c("steelblue1", "sandybrown"), pch=1, cex=0.8, horiz=F)
dev.off()

writeTable(rownames(ensGene.rt.start.sclc.inconsist), file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_inconsistent.list")), colnames=F, rownames=F, sep="")
writeTable(rownames(subset(ensGene.rt.start.sclc.inconsist, RATIO > 0)), file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_inconsistent_right.list")), colnames=F, rownames=F, sep="")
writeTable(rownames(subset(ensGene.rt.start.sclc.inconsist, RATIO < 0)), file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_inconsistent_left.list")), colnames=F, rownames=F, sep="")

# -----------------------------------------------------------------------------
# 
# Last Modified: 12/11/18
# -----------------------------------------------------------------------------
ensGene.rt.start.sclc <- cbind(ensGene[rownames(ensGene.rt.start.sclc), 1:7], ensGene.rt.start.sclc)
ensGene.rt.start.sclc$LENGTH <- abs(ensGene.rt.start.sclc$start_position - ensGene.rt.start.sclc$end_position)

right <- subset(ensGene.rt.start.sclc, RIGHT_LEADING > cutoff.upper)$MEDIAN
left  <- subset(ensGene.rt.start.sclc, RIGHT_LEADING < cutoff.lower)$MEDIAN
testW(right, left)
# [1] 0.1335412

##
xlab.text <- "Gene length (log10)"
ylab.text <- "log2(TPM+1)"

file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_length_right.pdf"))
main.text <- c("Right-leading genes in SCLC", paste0("n=", separator(length(sclc.tx.right))))
pdf(file.name, height=6, width=6)
plot(MEDIAN ~ log10(LENGTH), data=subset(ensGene.rt.start.sclc, RIGHT_LEADING > cutoff.upper), ylab=ylab.text, xlab=xlab.text, main=main.text, col="sandybrown")
abline(lm(MEDIAN ~ log10(LENGTH), data=subset(ensGene.rt.start.sclc, RIGHT_LEADING > cutoff.upper)), col="orange")
dev.off()

##
file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_length_left.pdf"))
main.text <- c("Left-leading genes in SCLC", paste0("n=", separator(length(sclc.tx.left))))
pdf(file.name, height=6, width=6)
plot(MEDIAN ~ log10(LENGTH), data=subset(ensGene.rt.start.sclc, RIGHT_LEADING < cutoff.lower), ylab=ylab.text, xlab=xlab.text, main=main.text, col="steelblue1")
abline(lm(MEDIAN ~ log10(LENGTH), data=subset(ensGene.rt.start.sclc, RIGHT_LEADING < cutoff.lower)), col="steelblue2")
dev.off()

# -----------------------------------------------------------------------------
# Transcription vs. replication time (Cell cycle genes; Dominguez et al 2016)
# Last Modified: 08/11/18
# -----------------------------------------------------------------------------
file.name  <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_periodic.png"))
main.text <- "Transcription vs. replication time in SCLC"
xlab.text <- "Right-leading ratio (log2)"
ylab.text <- "log2(TPM+1)"
cols <- rep("grey", nrow(ensGene.rt.start.sclc))

png(file.name, height=6, width=6, units="in", res=300)
plot(  MEDIAN ~ RATIO, data=ensGene.rt.start.sclc, xlab=xlab.text, ylab=ylab.text, main=main.text, col=cols)
points(MEDIAN ~ RATIO, data=ensGene.rt.start.sclc[periodic.G2M,], col="forestgreen")
points(MEDIAN ~ RATIO, data=ensGene.rt.start.sclc[periodic.G1S,], col="orange")
legend("bottom", c("G1/S", "G2/M"), col=c("orange", "forestgreen"), pch=1, cex=1, horiz=T)
mtext("Periodic gene lists (Dominguez 2016)", cex=1, line=0.5)
dev.off()

