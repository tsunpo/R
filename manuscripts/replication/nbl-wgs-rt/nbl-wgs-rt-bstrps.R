# =============================================================================
# Manuscript   : Using mutational strand asymmetry to study replication-transcription conflicts through space and time
# Chapter II   : Reconstruction of tissue-specific replicaiton timing profile in tumours
# Name         : manuscripts/replicaiton/nbl-wgs-rt-bstrps.R
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
BASE <- "NBL"
base <- tolower(BASE)
bstrps       <- 1000
origin.upper <- 510   ## 500-505, 505-510 breaks
origin.lower <- 490   ## 490-495, 495-500 breaks
origin.break <- 2     ## 2 breaks each centering 500

#wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@local
wd.anlys <- file.path(wd, BASE, "analysis")
wd.rt    <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data/bstrps")
wd.rt.plots <- file.path(wd.rt, "plots/bstrps")
wd.de    <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")

# -----------------------------------------------------------------------------
# Distributions in bootstrapped data
# Last Modified: 02/11/18
# -----------------------------------------------------------------------------
bed.gc.rt <- NULL
for (c in 1:22) {
   chr <- chrs[c]
 
   load(file=file.path(wd.rt.data, paste0("bed.gc.rt_", base, "_bstrps", bstrps, "_", chr, ".RData")))
   if (is.null(bed.gc.rt))
      bed.gc.rt <- bed.gc.rt.chr
   bed.gc.rt <- rbind(bed.gc.rt, bed.gc.rt.chr)

   file.name <- file.path(wd.rt.plots, paste0("hist_", base, "_bed.gc.rt_bstrps", bstrps, "_", chr, ".png"))
   plotHistBootstraps(bed.gc.rt.chr, file.name, paste0("Chr", c), BASE, 200, origin.break)   ## see ReplicationTiming.R
}
#save(bed.gc.rt, file=file.path(wd.rt.data, paste0("bed.gc.rt_", base,"_bstrps", bstrps, ".RData")))

file.name <- file.path(wd.rt.plots, paste0("hist_", base, "_bed.gc.rt_bstrps", bstrps, ".png"))
plotHistBootstraps(bed.gc.rt, file.name, "Chr1-22", BASE, 200, origin.break)                 ## see ReplicationTiming.R

# -----------------------------------------------------------------------------
# Plot RO and RT from bootstrapped data
# Last Modified: 04/11/18
# -----------------------------------------------------------------------------
for (c in 1:22) {
   c <- 5
   chr <- chrs[c]
 
   ## Replication origins   
   load(file=file.path(wd.rt.data, paste0("bed.gc.rt_", base, "_bstrps", bstrps, "_", chr, ".RData")))
   right.idx  <- which(bed.gc.rt.chr$RIGHT_LEADING > origin.upper)
   left.idx   <- which(bed.gc.rt.chr$RIGHT_LEADING < origin.lower)
   origin.idx <- setdiff(c(1:nrow(bed.gc.rt.chr)), c(right.idx, left.idx))
   bed.gc.chr <- bed.gc[rownames(bed.gc.rt.chr),]
   
   file.name <- file.path(wd.rt.plots, paste0(base, "_RO_bstrps1000_", chr))
   plotRO(file.name, BASE, chr, 1, 50000000, bed.gc.chr, bed.gc.rt.chr, right.idx, left.idx, origin.idx, "png")      ## see ReplicationTiming.R
   
   ## Replication timing
   rt.chr <- readTable(file.path(wd.rt, "data", paste0(base, "_rpkm.corr.gc.d.rt_", chr, "_", BASE, "-", BASE, "_n56-56.txt.gz")), header=T, rownames=T, ymax=1.5, sep="\t") 
   rt.chr <- rt.chr[rownames(bed.gc.rt.chr),]
   rt.chr$RT <- rt.chr$MEDIAN   ## Only for older version
   
   file.name  <- file.path(wd.rt.plots, paste0(base, "_RT_bstrps1000_", chr))
   plotRT(file.name, BASE, chr, 1, 50000000, rt.chr, bed.gc.chr, right.idx, left.idx, origin.idx, ymax=1.5, "png")   ## see ReplicationTiming.R
}

# -----------------------------------------------------------------------------
# Transcription vs. replication time
# Last Modified: 08/11/18
# -----------------------------------------------------------------------------
load(file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
tpm.gene      <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=T, proteinCodingOnly=F)
tpm.gene.log2 <- getLog2andMedian(tpm.gene, pseudocount=1)

load(file=file.path(wd.rt.data, paste0("ensGene.rt_", base, "_bstrps", bstrps, ".RData")))         ## Load objects ensGene.rt.start and ensGene.rt.end
ensGene.rt.start.nbl <- getEnsGeneRTTx(ensGene.rt.start, ensGene.rt.end, ensGene, tpm.gene.log2)   ## Gene start information; see ReplicationTiming.R
ensGene.rt.end.nbl   <- ensGene.rt.end[rownames(ensGene.rt.start.nbl),]                            ## Gene end information

##
nbl.tx.right  <- rownames(subset(ensGene.rt.start.nbl, RIGHT_LEADING > cutoff.upper))
nbl.tx.left   <- rownames(subset(ensGene.rt.start.nbl, RIGHT_LEADING < cutoff.lower))
nbl.tx.origin <- setdiff(rownames(ensGene.rt.start.nbl), c(nbl.tx.right, nbl.tx.left))
nbl.tx.inconsist <- rownames(subset(ensGene.rt.start.nbl, SIGN < 0))

file.name  <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_TSS+TES.pdf"))
main.text  <- "Transcription vs. replication time in NBL"
mtext.text <- paste0("Expressed genes (n=", separator(nrow(ensGene.rt.start.nbl)), ")")
xlab.text <- "Right-leading ratio (log2)"
ylab.text <- "log2(TPM+1)"

pdf(file.name, height=6, width=6)
cols <- rep("steelblue1", nrow(ensGene.rt.start.nbl))
plot(  MEDIAN ~ RATIO, data=ensGene.rt.start.nbl, xlab=xlab.text, ylab=ylab.text, main=main.text, col=cols)
points(MEDIAN ~ RATIO, data=ensGene.rt.start.nbl[nbl.tx.right,],     col="sandybrown")
points(MEDIAN ~ RATIO, data=ensGene.rt.start.nbl[nbl.tx.inconsist,], col="purple1")
points(MEDIAN ~ RATIO, data=ensGene.rt.start.nbl[nbl.tx.origin,],    col="red")
legend("top", c("Left-leading", paste0("L/R (n=", separator(length(nbl.tx.origin)), ")"), "Inconsistent", "Right-leading"), col=c("steelblue1", "red", "purple1", "sandybrown"), pch=1, cex=0.75, horiz=T)
mtext(mtext.text, cex=1, line=0.5)
dev.off()

# -----------------------------------------------------------------------------
# Inconsistent genes have lower expression levels
# Last Modified: 08/11/18
# -----------------------------------------------------------------------------
ensGene.rt.start.nbl.inconsist <- ensGene.rt.start.nbl[nbl.tx.inconsist,]
nbl.tx.inconsist.right <- subset(ensGene.rt.start.nbl.inconsist, RATIO > 0)$MEDIAN
nbl.tx.inconsist.left  <- subset(ensGene.rt.start.nbl.inconsist, RATIO < 0)$MEDIAN
testW(nbl.tx.inconsist.right, nbl.tx.inconsist.left)
# [1] 6.486166e-11

file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_inconsistent_beeswarm.png"))
main.text <- c("Inconsistent genes in NBL", paste0("n=", length(nbl.tx.inconsist)))
xlab.text <- "Replication time"
ylab.text <- "log2(TPM+0.01)"
ymin <- min(ensGene.rt.start.sclc$MEDIAN)
ymax <- max(ensGene.rt.start.sclc$MEDIAN)
overlaps <- intersect(nbl.tx.inconsist, nbl.tx.origin)
names <- c(paste0("L (n=", length(nbl.tx.inconsist.left),")"), paste0("R (n=", length(nbl.tx.inconsist.right),")"))

png(file.name, height=6, width=3.5, units="in", res=300)
 boxplot(MEDIAN ~ RT, data=ensGene.rt.start.nbl.inconsist, ylim=c(ymin, ymax), ylab=ylab.text, xlab=xlab.text, names=names, main=main.text, border="purple1", outline=T)
beeswarm(MEDIAN ~ RT, data=ensGene.rt.start.nbl.inconsist[overlaps,], col="red", pch=1, add=T)
legend("bottomleft", paste0("L/R (n=", length(overlaps), ")"), col="red", pch=1, cex=0.8)
dev.off()

##
nbl.tx.all <- ensGene.rt.start.nbl$MEDIAN
testW(nbl.tx.all, c(nbl.tx.inconsist.right, nbl.tx.inconsist.left))
# [1] 0.1960944
testW(nbl.tx.all, nbl.tx.inconsist.right)
# [1] 9.194107e-05
testW(nbl.tx.all, nbl.tx.inconsist.left)
# [1] 1.548922e-06

##
nbl.tx.inconsist.right.length <- subset(ensGene.rt.start.nbl.inconsist, RATIO > 0)$LENGTH
nbl.tx.inconsist.left.length  <- subset(ensGene.rt.start.nbl.inconsist, RATIO < 0)$LENGTH
testW(nbl.tx.inconsist.right.length, nbl.tx.inconsist.left.length)
# [1] 2.128381e-09

##
file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_inconsistent_length.png"))
main.text <- c("Inconsistent genes in NBL", paste0("n=", length(nbl.tx.inconsist)))
xlab.text <- "Gene length (log10)"
ylab.text <- "log2(TPM+0.01)"

png(file.name, height=6, width=6, units="in", res=300)
  plot(MEDIAN ~ log10(LENGTH), data=subset(ensGene.rt.start.nbl.inconsist, RATIO > 0), ylab=ylab.text, xlab=xlab.text, main=main.text, col="sandybrown")
points(MEDIAN ~ log10(LENGTH), data=subset(ensGene.rt.start.nbl.inconsist, RATIO < 0), col="steelblue1", pch=1)
abline(lm(MEDIAN ~ log10(LENGTH), data=subset(ensGene.rt.start.nbl.inconsist, RATIO > 0)), col="sandybrown")
abline(lm(MEDIAN ~ log10(LENGTH), data=subset(ensGene.rt.start.nbl.inconsist, RATIO < 0)), col="steelblue1")
legend("bottomleft", c("Left-leading", "Right-leading"), col=c("steelblue1", "sandybrown"), pch=1, cex=0.8, horiz=F)
dev.off()

writeTable(rownames(ensGene.rt.start.nbl.inconsist), file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_inconsistent.list")), colnames=F, rownames=F, sep="")
writeTable(rownames(subset(ensGene.rt.start.nbl.inconsist, RATIO > 0)), file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_inconsistent_right.list")), colnames=F, rownames=F, sep="")
writeTable(rownames(subset(ensGene.rt.start.nbl.inconsist, RATIO < 0)), file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_inconsistent_left.list")), colnames=F, rownames=F, sep="")

# -----------------------------------------------------------------------------
# Transcription vs. replication time (Cell cycle genes; Dominguez et al 2016)
# Last Modified: 08/11/18
# -----------------------------------------------------------------------------
file.name  <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_periodic.png"))
main.text <- "Transcription vs. replication time in SCLC"
xlab.text <- "Right-leading ratio (log2)"
ylab.text <- "Expression log2(TPM+0.01)"
cols <- rep("grey", nrow(ensGene.rt.start.nbl))

png(file.name, height=6, width=6, units="in", res=300)
  plot(MEDIAN ~ RATIO, data=ensGene.rt.start.nbl, xlab=xlab.text, ylab=ylab.text, main=main.text, col=cols)
points(MEDIAN ~ RATIO, data=ensGene.rt.start.nbl[periodic.G2M,], col="forestgreen")
points(MEDIAN ~ RATIO, data=ensGene.rt.start.nbl[periodic.G1S,], col="orange")
legend("bottom", c("G1/S", "G2/M"), col=c("orange", "forestgreen"), pch=1, cex=1, horiz=T)
mtext("Periodic gene lists (Dominguez 2016)", cex=1, line=0.5)
dev.off()








## RB1
genes.rb1.up   <- rownames(subset(subset(de.tpm.gene, FDR <= 0.05), LOG2_FC > 0))
genes.rb1.down <- rownames(subset(subset(de.tpm.gene, FDR <= 0.05), LOG2_FC < 0))

file.name  <- file.path(wd.rt.plots, paste0(tolower(BASE), "_ensGene.rt.start.tx_bstrps1000_rb1_q0.05.png"))
main.text <- "Transcription vs. replication time in SCLC"
xlab.text <- "Right-leading ratio (log2)"
ylab.text <- "log2(TPM+0.01) expression"
cols <- rep("grey", nrow(ensGene.rt.start.tx))

png(file.name, height=6, width=6, units="in", res=300)
plot(ensGene.rt.start.tx$MEDIAN ~ ensGene.rt.start.tx$RATIO, xlab=xlab.text, ylab=ylab.text, main=main.text, col=cols)
points(ensGene.rt.start.tx[genes.rb1.down,]$MEDIAN ~ ensGene.rt.start.tx[genes.rb1.down,]$RATIO, col="dodgerblue")
points(ensGene.rt.start.tx[genes.rb1.up,]$MEDIAN ~ ensGene.rt.start.tx[genes.rb1.up,]$RATIO, col="red")
legend("bottom", c("Downregulated", "Upregulated in RB1"), col=c("dodgerblue", "red"), pch=1, cex=1, horiz=T)
mtext("RB1-loss D.E. in LCNEC (FDR=0.05)", cex=1, line=0.5)
dev.off()

###
##
inconsist.right <- subset(ensGene.rt.start.tx[sclc.tx.inconsist,], RATIO > 0)$MEDIAN
inconsist.left  <- subset(ensGene.rt.start.tx[sclc.tx.inconsist,], RATIO < 0)$MEDIAN
testW(inconsist.right, inconsist.left)
# [1] 1.235423e-15

ensGene.rt.start.tx <- cbind(ensGene[rownames(ensGene.rt.start.tx), 1:7], ensGene.rt.start.tx)
ensGene.rt.start.tx$SIGN_RTC <- ensGene.rt.start.tx$RT * ensGene.rt.start.tx$strand
inconsist.right.cd <- subset(subset(ensGene.rt.start.tx[sclc.tx.inconsist,], RATIO > 0), SIGN_RTC > 0)$MEDIAN
inconsist.right.ho <- subset(subset(ensGene.rt.start.tx[sclc.tx.inconsist,], RATIO > 0), SIGN_RTC < 0)$MEDIAN
inconsist.left.cd <- subset(subset(ensGene.rt.start.tx[sclc.tx.inconsist,], RATIO < 0), SIGN_RTC > 0)$MEDIAN
inconsist.left.ho <- subset(subset(ensGene.rt.start.tx[sclc.tx.inconsist,], RATIO < 0), SIGN_RTC < 0)$MEDIAN
testW(inconsist.right.cd, inconsist.right.ho)
# [1] 0.6957948
testW(inconsist.left.cd, inconsist.left.ho)
# [1] 0.5029064
testW(inconsist.right.cd, inconsist.left.cd)
# [1] 1.451833e-08
testW(inconsist.right.ho, inconsist.left.ho)
# [1] 1.83927e-08

expressions <- c(inconsist.left, inconsist.right)
groups <- rep("left", length(inconsist.left))
groups <- c(groups, rep("right", length(inconsist.right)))
x <- as.data.frame(cbind(expressions, groups))
x$expressions <- as.vector(x$expressions)
x$expressions <- 0
x$expressions <- expressions
x$groups <- as.factor(x$groups)

file.name <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_inconsistent_beeswarm.pdf"))
main.text <- c("Inconsistent genes in SCLC", paste0("n=", length(expressions)))
xlab.text <- "Replication time"
ylab.text <- "log2(TPM+0.01) expression"
x <- subset(ensGene.rt.start.tx, SIGN < 0)
overlaps <- intersect(rownames(x), sclc.tx.origin)
names <- c(paste0("L (n=", length(inconsist.left),")"), paste0("R (n=", length(inconsist.right),")"))

pdf(file.name, height=6, width=3.5)
boxplot(MEDIAN ~ RT, data=x, ylab=ylab.text, xlab=xlab.text, names=names, main=main.text, border="purple1", outline=T)
beeswarm(MEDIAN ~ RT, data=x[overlaps,], col="red", pch=1, add=T)
#beeswarm(expressions ~ groups, data=subset(x, groups == "right"), col="purple1", pch=1, add=T)
legend("bottomleft", paste0("L/R (n=", length(overlaps), ")"), col="red", pch=1, cex=0.8)
dev.off()

##
origin.right <- subset(ensGene.rt.start.tx[sclc.tx.origin,], RATIO > 0)$MEDIAN
origin.left  <- subset(ensGene.rt.start.tx[sclc.tx.origin,], RATIO < 0)$MEDIAN
testW(origin.right, origin.left)
# [1] 0.4519312

origin.right.inconsist <- subset(subset(ensGene.rt.start.tx[sclc.tx.origin,], SIGN < 0), RATIO > 0)$MEDIAN
origin.left.inconsist  <- subset(subset(ensGene.rt.start.tx[sclc.tx.origin,], SIGN < 0), RATIO < 0)$MEDIAN
testW(origin.right.inconsist, origin.left.inconsist)
# [1] 0.3232271

origin.right.consist <- subset(subset(ensGene.rt.start.tx[sclc.tx.origin,], SIGN > 0), RATIO > 0)$MEDIAN
origin.left.consist  <- subset(subset(ensGene.rt.start.tx[sclc.tx.origin,], SIGN > 0), RATIO < 0)$MEDIAN
testW(origin.right.consist, origin.left.consist)

## FDR 0.05
genes.rb1.up   <- rownames(subset(subset(de.tpm.gene, FDR <= 0.1), LOG2_FC > 0))
genes.rb1.down <- rownames(subset(subset(de.tpm.gene, FDR <= 0.1), LOG2_FC < 0))
#genes.rb1.up   <- rownames(subset(de.tpm.gene, LOG2_FC > 0))
#genes.rb1.down <- rownames(subset(de.tpm.gene, LOG2_FC < 0))

up.right <- subset(ensGene.rt.start.tx[genes.rb1.up,], RATIO > 0)$MEDIAN
up.left  <- subset(ensGene.rt.start.tx[genes.rb1.up,], RATIO < 0)$MEDIAN
testW(up.right, up.left)
# [1] 0.7866975
down.right <- subset(ensGene.rt.start.tx[genes.rb1.down,], RATIO > 0)$MEDIAN
down.left  <- subset(ensGene.rt.start.tx[genes.rb1.down,], RATIO < 0)$MEDIAN
testW(down.right, down.left)
# [1] 0.9586084

## FDR 0.1
up.right <- subset(ensGene.rt.start.tx[genes.rb1.up,], RATIO > 0)$MEDIAN
up.left  <- subset(ensGene.rt.start.tx[genes.rb1.up,], RATIO < 0)$MEDIAN
testW(up.right, up.left)
# [1] 0.7700956
down.right <- subset(ensGene.rt.start.tx[genes.rb1.down,], RATIO > 0)$MEDIAN
down.left  <- subset(ensGene.rt.start.tx[genes.rb1.down,], RATIO < 0)$MEDIAN
testW(down.right, down.left)
# [1] 0.1839742

## FDR 0.5
up.right <- subset(ensGene.rt.start.tx[genes.rb1.up,], RATIO > 0)$MEDIAN
up.left  <- subset(ensGene.rt.start.tx[genes.rb1.up,], RATIO < 0)$MEDIAN
testW(up.right, up.left)
# [1] 0.4962725
down.right <- subset(ensGene.rt.start.tx[genes.rb1.down,], RATIO > 0)$MEDIAN
down.left  <- subset(ensGene.rt.start.tx[genes.rb1.down,], RATIO < 0)$MEDIAN
testW(down.right, down.left)
# [1] 0.09845262

##
g1s.right <- subset(ensGene.rt.start.tx[periodic.G1S,], RATIO > 0)$MEDIAN
g1s.left  <- subset(ensGene.rt.start.tx[periodic.G1S,], RATIO < 0)$MEDIAN
testW(g1s.right, g1s.left)
# [1] 0.2413578

g2m.right <- subset(ensGene.rt.start.tx[periodic.G2M,], RATIO > 0)$MEDIAN
g2m.left  <- subset(ensGene.rt.start.tx[periodic.G2M,], RATIO < 0)$MEDIAN
testW(g2m.right, g2m.left)
# [1] 0.3858542

##
test <- ensGene.rt.start.tx[genes.rb1.down,]
subset(subset(subset(test, MEDIAN < 5), RATIO > -0.5), RATIO < 0.5)
#ENSG00000029559 ENSG00000029559            chr4      1       88720733     88733074 protein_coding               IBSP
#ENSG00000105696 ENSG00000105696           chr19      1       18718240     18731849 protein_coding            TMEM59L
#ENSG00000154328 ENSG00000154328            chr8      1       11627148     11644855 protein_coding              NEIL2

###
## Compare SCLC RT with LCL RT
load("/Users/tpyang/Work/uni-koeln/tyang2/LCL/analysis/replication/lcl-wgs-rt/data/bstrps/ensGene.rt_lcl_bstrps1000.RData")
ensGene.rt.start.lcl <- ensGene.rt.start
ensGene.rt.start.lcl$RATIO <- mapply(x = 1:nrow(ensGene.rt.start.lcl), function(x) as.numeric(getLeadingRatio(ensGene.rt.start.lcl[x,])))
ensGene.rt.end.lcl <- ensGene.rt.end

sclc.o.lcl <- intersect(rownames(ensGene.rt.start.tx), rownames(ensGene.rt.start.lcl))
ensGene.rt.start.sclc.o.lcl  <- ensGene.rt.start.tx[sclc.o.lcl,]
ensGene.rt.start.lcl.o.sclc <- ensGene.rt.start.lcl[sclc.o.lcl,]
periodic.G1S.sclc.o.lcl <- intersect(periodic.G1S, sclc.o.lcl)
periodic.G2M.sclc.o.lcl <- intersect(periodic.G2M, sclc.o.lcl)
sclc.tx.right.o.lcl  <- intersect(sclc.tx.right, sclc.o.lcl)
sclc.tx.left.o.lcl   <- intersect(sclc.tx.left, sclc.o.lcl)
sclc.tx.origin.o.lcl <- intersect(sclc.tx.origin, sclc.o.lcl)
sclc.tx.inconsist.o.lcl <- intersect(sclc.tx.inconsist, sclc.o.lcl)

ensGene.rt.start.lcl.o.sclc$SIGN_SCLC <- ensGene.rt.start.lcl.o.sclc$RT * ensGene.rt.start.sclc.o.lcl$RT
sclc.o.c.lcl <- rownames(subset(ensGene.rt.start.lcl.o.sclc, SIGN_SCLC > 0))
periodic.G1S.sclc.o.c.lcl <- intersect(periodic.G1S, sclc.o.c.lcl)
periodic.G2M.sclc.o.c.lcl <- intersect(periodic.G2M, sclc.o.c.lcl)
sclc.tx.right.o.c.lcl  <- intersect(sclc.tx.right, sclc.o.c.lcl)
sclc.tx.left.o.c.lcl   <- intersect(sclc.tx.left, sclc.o.c.lcl)
sclc.tx.origin.o.c.lcl <- intersect(sclc.tx.origin, sclc.o.c.lcl)
sclc.tx.inconsist.o.c.lcl <- intersect(sclc.tx.inconsist, sclc.o.c.lcl)

length(sclc.o.c.lcl)/length(sclc.o.lcl)
# [1] 0.6967831
length(periodic.G1S.sclc.o.c.lcl)/length(periodic.G1S.sclc.o.lcl)
# [1] 0.720524
length(periodic.G2M.sclc.o.c.lcl)/length(periodic.G2M.sclc.o.lcl)
# [1] 0.6657303

length(sclc.tx.right.o.c.lcl)/length(sclc.tx.right.o.lcl)
# [1] 0.7014773
length(sclc.tx.left.o.c.lcl)/length(sclc.tx.left.o.lcl)
# [1] 0.6933915
length(sclc.tx.origin.o.c.lcl)/length(sclc.tx.origin.o.lcl)
# [1] 0.4772727
length(sclc.tx.inconsist.lcl.o.c)/length(sclc.tx.inconsist.lcl.o)
# [1] 0.592397

##
file.name  <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_sclc.o.c.lcl_tx.png"))
main.text <- "Transcription vs. replication time"
xlab.text <- c("Right-leading ratio (log2)", "LCL")
ylab.text <- c("SCLC", "log2(TPM+0.01)")
cols <- rep("gray", nrow(ensGene.rt.start.sclc.o.lcl))

png(file.name, height=6, width=6, units="in", res=300)
plot(ensGene.rt.start.sclc.o.lcl$MEDIAN ~ ensGene.rt.start.lcl.o.sclc$RATIO, xlab=xlab.text, ylab=ylab.text, main=main.text, col=cols)
points(ensGene.rt.start.sclc.o.lcl[sclc.tx.left.o.c,]$MEDIAN ~ ensGene.rt.start.lcl.o.sclc[sclc.tx.left.o.c,]$RATIO, col="steelblue1")
points(ensGene.rt.start.sclc.o.lcl[sclc.tx.right.o.c,]$MEDIAN ~ ensGene.rt.start.lcl.o.sclc[sclc.tx.right.o.c,]$RATIO, col="sandybrown")
points(ensGene.rt.start.sclc.o.lcl[sclc.tx.origin.o.c,]$MEDIAN ~ ensGene.rt.start.lcl.o.sclc[sclc.tx.origin.o.c,]$RATIO, col="red")
legend("bottom", c("Left-leading in SCLC", paste0("L/R (n=", length(sclc.tx.origin.o.c), ") in SCLC"), "Right-leading in SCLC"), col=c("steelblue1", "red", "sandybrown"), pch=1, cex=0.8, horiz=T)
mtext(paste0("Overlapping expressed genes (n=", nrow(ensGene.rt.start.tx.o), ")"), cex=1, line=0.5)
dev.off()

##
file.name  <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_sclc.o.c.lcl.png"))
main.text <- "Replication time (SCLC vs. LCL)"
xlab.text <- c("Right-leading ratio (log2)", "SCLC")
ylab.text <- c("LCL (Koren, 2012)", "Right-leading ratio (log2)")
cols <- rep("gray", nrow(ensGene.rt.start.sclc.o.lcl))

png(file.name, height=6, width=6, units="in", res=300)
plot(ensGene.rt.start.lcl.o.sclc$RATIO ~ ensGene.rt.start.sclc.o.lcl$RATIO, xlab=xlab.text, ylab=ylab.text, main=main.text, col=cols)
points(ensGene.rt.start.lcl.o.sclc[sclc.tx.left.o.c,]$RATIO ~ ensGene.rt.start.sclc.o.lcl[sclc.tx.left.o.c,]$RATIO, col="steelblue1")
points(ensGene.rt.start.lcl.o.sclc[sclc.tx.right.o.c,]$RATIO ~ ensGene.rt.start.sclc.o.lcl[sclc.tx.right.o.c,]$RATIO, col="sandybrown")
points(ensGene.rt.start.lcl.o.sclc[sclc.tx.origin.o.c,]$RATIO ~ ensGene.rt.start.sclc.o.lcl[sclc.tx.origin.o.c,]$RATIO, col="red")
legend("bottom", c("Left-leading in SCLC", paste0("L/R (n=", length(sclc.tx.origin.o.c), ") in SCLC"), "Right-leading in SCLC"), col=c("steelblue1", "red", "sandybrown"), pch=1, cex=0.8, horiz=T)
mtext(paste0("Overlapping expressed genes (n=", nrow(ensGene.rt.start.sclc.o.lcl), ")"), cex=1, line=0.5)
dev.off()

###
## Compare SCLC RT with NBL RT
load("/Users/tpyang/Work/uni-koeln/tyang2/NBL/analysis/replication/nbl-wgs-rt/data/bstrps/ensGene.rt_nbl_bstrps1000.RData")
ensGene.rt.start.nbl <- ensGene.rt.start
ensGene.rt.start.nbl$RATIO <- mapply(x = 1:nrow(ensGene.rt.start.nbl), function(x) as.numeric(getLeadingRatio(ensGene.rt.start.nbl[x,])))
ensGene.rt.end.nbl <- ensGene.rt.end

sclc.o.nbl <- intersect(rownames(ensGene.rt.start.tx), rownames(ensGene.rt.start.nbl))
ensGene.rt.start.sclc.o.nbl  <- ensGene.rt.start.tx[sclc.o.nbl,]
ensGene.rt.start.nbl.o.sclc <- ensGene.rt.start.nbl[sclc.o.nbl,]
periodic.G1S.sclc.o.nbl <- intersect(periodic.G1S, sclc.o.nbl)
periodic.G2M.sclc.o.nbl <- intersect(periodic.G2M, sclc.o.nbl)
sclc.tx.right.o.nbl  <- intersect(sclc.tx.right, sclc.o.nbl)
sclc.tx.left.o.nbl   <- intersect(sclc.tx.left, sclc.o.nbl)
sclc.tx.origin.o.nbl <- intersect(sclc.tx.origin, sclc.o.nbl)
sclc.tx.inconsist.o.nbl <- intersect(sclc.tx.inconsist, sclc.o.nbl)

ensGene.rt.start.nbl.o.sclc$SIGN <- ensGene.rt.start.nbl.o.sclc$RT * ensGene.rt.start.sclc.o.nbl$RT
sclc.o.c.nbl <- rownames(subset(ensGene.rt.start.nbl.o.sclc, SIGN > 0))
periodic.G1S.sclc.o.c.nbl <- intersect(periodic.G1S, sclc.o.c.nbl)
periodic.G2M.sclc.o.c.nbl <- intersect(periodic.G2M, sclc.o.c.nbl)
sclc.tx.right.o.c.nbl  <- intersect(sclc.tx.right, sclc.o.c.nbl)
sclc.tx.left.o.c.nbl   <- intersect(sclc.tx.left, sclc.o.c.nbl)
sclc.tx.origin.o.c.nbl <- intersect(sclc.tx.origin, sclc.o.c.nbl)
sclc.tx.inconsist.o.c.nbl <- intersect(sclc.tx.inconsist, sclc.o.c.nbl)

length(sclc.o.c.nbl)/length(sclc.o.nbl)
# [1] 0.7343118
length(periodic.G1S.sclc.o.c.nbl)/length(periodic.G1S.sclc.o.nbl)
# [1] 0.6920152
length(periodic.G2M.sclc.o.c.nbl)/length(periodic.G2M.sclc.o.nbl)
# [1] 0.7319202

length(sclc.tx.right.o.c.nbl)/length(sclc.tx.right.o.nbl)
# [1] 0.7214198
length(sclc.tx.left.o.c.nbl)/length(sclc.tx.left.o.nbl)
# [1] 0.7480975
length(sclc.tx.origin.o.c.nbl)/length(sclc.tx.origin.o.nbl)
# [1] 0.5614035
length(sclc.tx.inconsist.o.c.nbl)/length(sclc.tx.inconsist.o.nbl)
# [1] 0.5992714

##
file.name  <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_sclc.o.c.nbl_tx.png"))
main.text <- "Transcription vs. replication time"
xlab.text <- c("Right-leading ratio (log2)", "NBL")
ylab.text <- c("SCLC", "log2(TPM+0.01)")
cols <- rep("gray", nrow(ensGene.rt.start.sclc.o.lcl))

png(file.name, height=6, width=6, units="in", res=300)
plot(ensGene.rt.start.sclc.o.nbl$MEDIAN ~ ensGene.rt.start.nbl.o.sclc$RATIO, xlab=xlab.text, ylab=ylab.text, main=main.text, col=cols)
points(ensGene.rt.start.sclc.o.nbl[sclc.tx.left.o.c.nbl,]$MEDIAN ~ ensGene.rt.start.nbl.o.sclc[sclc.tx.left.o.c.nbl,]$RATIO, col="steelblue1")
points(ensGene.rt.start.sclc.o.nbl[sclc.tx.right.o.c.nbl,]$MEDIAN ~ ensGene.rt.start.nbl.o.sclc[sclc.tx.right.o.c.nbl,]$RATIO, col="sandybrown")
points(ensGene.rt.start.sclc.o.nbl[sclc.tx.origin.o.c.nbl,]$MEDIAN ~ ensGene.rt.start.nbl.o.sclc[sclc.tx.origin.o.c.nbl,]$RATIO, col="red")
legend("bottom", c("Left-leading in SCLC", paste0("L/R (n=", length(sclc.tx.origin.o.c.nbl), ") in SCLC"), "Right-leading in SCLC"), col=c("steelblue1", "red", "sandybrown"), pch=1, cex=0.8, horiz=T)
mtext(paste0("Overlapping expressed genes (n=", nrow(ensGene.rt.start.sclc.o.nbl), ")"), cex=1, line=0.5)
dev.off()


###
## Compare NBL RT with LCL RT (TO-DO)
lcl.o.nbl <- intersect(rownames(ensGene.rt.start.lcl), rownames(ensGene.rt.start.nbl))
ensGene.rt.start.lcl.o.nbl <- ensGene.rt.start.lcl[lcl.o.nbl,]
ensGene.rt.start.nbl.o.lcl <- ensGene.rt.start.nbl[lcl.o.nbl,]
periodic.G1S.lcl.o.nbl <- intersect(periodic.G1S, lcl.o.nbl)
periodic.G2M.lcl.o.nbl <- intersect(periodic.G2M, lcl.o.nbl)

ensGene.rt.start.lcl.o.nbl$SIGN <- ensGene.rt.start.lcl.o.nbl$RT * ensGene.rt.start.nbl.o.lcl$RT
lcl.o.c.nbl <- rownames(subset(ensGene.rt.start.lcl.o.nbl, SIGN > 0))
periodic.G1S.lcl.o.c.nbl <- intersect(periodic.G1S, lcl.o.c.nbl)
periodic.G2M.lcl.o.c.nbl <- intersect(periodic.G2M, lcl.o.c.nbl)

length(lcl.o.c.nbl)/length(lcl.o.nbl)
# [1] 0.6583011
length(periodic.G1S.lcl.o.c.nbl)/length(periodic.G1S.lcl.o.nbl)
# [1] 0.6130435
length(periodic.G2M.lcl.o.c.nbl)/length(periodic.G2M.lcl.o.nbl)
# [1] 0.6596639



##
file.name  <- file.path(wd.rt.plots, paste0(base, "_ensGene.rt.start.tx_bstrps1000_nbl.png"))
main.text <- "Transcription vs. replication time"
xlab.text <- c("Right-leading ratio (log2)", "NBL")
ylab.text <- c("SCLC", "Expression (log2 TPM+0.01)")
cols <- rep("steelblue1", nrow(ensGene.rt.start.tx.o))

png(file.name, height=6, width=6, units="in", res=300)
plot(ensGene.rt.start.tx.o$MEDIAN ~ ensGene.rt.start.nbl.o$RATIO, xlab=xlab.text, ylab=ylab.text, main=main.text, col=cols)
points(ensGene.rt.start.tx.o[sclc.tx.right.o,]$MEDIAN ~ ensGene.rt.start.nbl.o[sclc.tx.right.o,]$RATIO, col="sandybrown")
points(ensGene.rt.start.tx.o[sclc.tx.origin.o,]$MEDIAN ~ ensGene.rt.start.nbl.o[sclc.tx.origin.o,]$RATIO, col="red")
legend("bottom", c("Left-leading in SCLC", paste0("L/R (n=", length(sclc.tx.origin.o), ") in SCLC"), "Right-leading in SCLC"), col=c("steelblue1", "red", "sandybrown"), pch=1, cex=0.8, horiz=T)
mtext(paste0("Overlapping expressed genes (n=", nrow(ensGene.rt.start.tx.o), ")"), cex=1, line=0.5)
dev.off()




subset(subset(subset(subset(test, MEDIAN > 0), MEDIAN < 5), RATIO > -0.5), RATIO < 0.5)

###
## Dominguez et al 2016
file.name  <- file.path(wd.rt.plots, paste0(tolower(BASE), "_ensGene.rt.start.tx_bstrps1000_periodic.png"))
main.text <- "Transcription vs. replication time in NBL"
xlab.text <- "Right-leading ratio (log2)"
ylab.text <- "Expression log2(TPM+0.01)"
cols <- rep("grey", nrow(ensGene.rt.start.tx))

png(file.name, height=6, width=6, units="in", res=300)
plot(ensGene.rt.start.tx$MEDIAN ~ ensGene.rt.start.tx$RATIO, xlab=xlab.text, ylab=ylab.text, main=main.text, col=cols)
points(ensGene.rt.start.tx[periodic.G2M,]$MEDIAN ~ ensGene.rt.start.tx[periodic.G2M,]$RATIO, col="forestgreen")
points(ensGene.rt.start.tx[periodic.G1S,]$MEDIAN ~ ensGene.rt.start.tx[periodic.G1S,]$RATIO, col="orange")
legend("bottom", c("G1/S", "G2/M"), col=c("orange", "forestgreen"), pch=1, cex=1, horiz=T)
mtext("Periodic gene lists (Dominguez 2016)", cex=1, line=0.5)
dev.off()

## c("red","deepskyblue","forestgreen","purple3","blue","gold","lightsalmon","turquoise1","limegreen")
## Tirosh et al 2016 
core.G1S.idx  <- which(rownames(ensGene.rt.start.tx) %in% core.G1S)
core.G2M.idx  <- which(rownames(ensGene.rt.start.tx) %in% core.G2M)
core.Stemness.idx  <- which(rownames(ensGene.rt.start.tx) %in% core.Stemness)

file.name  <- file.path(wd.rt.plots, paste0(tolower(BASE), "_ensGene.rt.start.tx_bstrps1000_core.png"))
main.text <- "Transcription vs. replication time"
xlab.text <- c("Right-leading ratio (log2)", BASE)
ylab.text <- c("SCLC", "Expression (log2(TPM+0.01))")
cols <- rep("grey", nrow(ensGene.rt.start.tx))

png(file.name, height=6, width=6, units="in", res=300)
plot(ensGene.rt.start.tx$MEDIAN ~ ensGene.rt.start.tx$RATIO, xlab=xlab.text, ylab=ylab.text, main=main.text, col=cols)
points(ensGene.rt.start.tx$MEDIAN[core.Stemness.idx] ~ ensGene.rt.start.tx$RATIO[core.Stemness.idx], col="purple3")
points(ensGene.rt.start.tx$MEDIAN[core.G2M.idx] ~ ensGene.rt.start.tx$RATIO[core.G2M.idx], col="gold")
points(ensGene.rt.start.tx$MEDIAN[core.G1S.idx] ~ ensGene.rt.start.tx$RATIO[core.G1S.idx], col="lightsalmon")
legend("bottom", c("G1/S", "G2/M", "Stemness"), col=c("lightsalmon", "gold", "purple3"), pch=1, cex=1, horiz=T)
dev.off()






bed.gc.rt.chr$RATIO <- mapply(x = 1:nrow(bed.gc.rt.chr), function(x) as.numeric(getLeadingRatio(bed.gc.rt.chr[x,])))
right.idx  <- which(bed.gc.rt.chr$RIGHT_LEADING > 505)
left.idx   <- which(bed.gc.rt.chr$RIGHT_LEADING < 495)
origin.idx <- setdiff(c(1:nrow(bed.gc.rt.chr)), c(right.idx, left.idx))





# -----------------------------------------------------------------------------
# Replication origins from bootstrapped data
# Last Modified: 04/11/18
# -----------------------------------------------------------------------------
getLeadingRatio <- function(bed.gc.rt) {
   if (bed.gc.rt$RT == 1) return(log2(bed.gc.rt$RIGHT_LEADING/500))
   else if (bed.gc.rt$RT == -1) return(log2(bed.gc.rt$LEFT_LEADING/500) * -1)
   else return(0)
}

plotRO <- function(wd.rt.plots, BASE, chr, xmin, xmax, bed.gc.rt.chr, bed.gc.chr, ext) {
   file.name  <- file.path(wd.rt.plots, paste0(tolower(BASE), "_RO_bstrps1000_", chr))
   main.text <- paste0("Bootstrapped right-leading ratio in ", BASE)
   xlab.text <- paste0("Chromosome ", gsub("chr", "", chr), " coordinate (Mb)")
   ylab.text <- "Right-leading ratio (log2)"
   if (!is.na(xmin) && !is.na(xmax)) file.name <- paste0(file.name, "_", xmin/1E6, "-", xmax/1E6, "Mb")
   if (is.na(xmin)) xmin <- 0
   if (is.na(xmax)) xmax <- subset(chromInfo, chrom == chr)$size
   ymin <- -2
   ymax <- 1.5
   
   ## Initiation plot
   if (ext == "pdf") {
      pdf(paste0(file.name, ".pdf"), height=3.5, width=10)
   } else if (ext == "png")
      png(paste0(file.name, ".png"), height=3.5, width=10, units="in", res=300)   ## ADD 16/05/17: res=300
   plot(NULL, ylim=c(ymin, ymax), xlim=c(xmin/1E6, xmax/1E6), xlab=xlab.text, ylab=ylab.text, main=main.text)
   abline(h=0, lwd=0.5, col="grey")
   
   ## Plot right-/left-leading positions
   bed.gc.rt.chr$RATIO <- mapply(x = 1:nrow(bed.gc.rt.chr), function(x) as.numeric(getLeadingRatio(bed.gc.rt.chr[x,])))
   right.idx  <- which(bed.gc.rt.chr$RIGHT_LEADING > 505)
   left.idx   <- which(bed.gc.rt.chr$RIGHT_LEADING < 495)
   origin.idx <- setdiff(c(1:nrow(bed.gc.rt.chr)), c(right.idx, left.idx))
   
   points(bed.gc.chr$START[left.idx]/1E6, bed.gc.rt.chr$RATIO[left.idx], col="steelblue1", cex=0.2)
   points(bed.gc.chr$START[right.idx]/1E6, bed.gc.rt.chr$RATIO[right.idx], col="sandybrown", cex=0.2)
   points(bed.gc.chr$START[origin.idx]/1E6, bed.gc.rt.chr$RATIO[origin.idx], col="red", cex=0.5)
   
   spline <- smooth.spline(x=bed.gc.chr$START, y=bed.gc.rt.chr$RATIO)
   lines(spline$x/1E6, spline$y)

   ## Plot cytobands and legend
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.4, col="lightgrey")
   
   legend("bottomright", c("Left-leading", "Left / Right", "Right-leading"), col=c("steelblue1", "red", "sandybrown"), pch=16, cex=0.8, horiz=T)
   dev.off()
}

plotRT <- function(wd.rt.plots, BASE, chr, xmin, xmax, rt.chr, bed.gc.chr, bed.gc.rt.chr, ext) {
   file.name  <- file.path(wd.rt.plots, paste0(tolower(BASE), "_RT_bstrps1000_", chr))
   main.text <- paste0("Read depth (CN-, GC-corrected) ratio (T/N) in ", BASE)
   xlab.text <- paste0("Chromosome ", gsub("chr", "", chr), " coordinate (Mb)")
   ylab.text <- "Replication time (log2)"
   if (!is.na(xmin) && !is.na(xmax)) file.name <- paste0(file.name, "_", xmin/1E6, "-", xmax/1E6, "Mb")
   if (is.na(xmin)) xmin <- 0
   if (is.na(xmax)) xmax <- subset(chromInfo, chrom == chr)$size
   ymin <- -0.15
   ymax <- 0.15
   
   ## Initiation plot
   if (ext == "pdf") {
      pdf(paste0(file.name, ".pdf"), height=4, width=10)
   } else if (ext == "png")
      png(paste0(file.name, ".png"), height=4, width=10, units="in", res=300)   ## ADD 16/05/17: res=300
   plot(NULL, ylim=c(ymin, ymax), xlim=c(xmin/1E6, xmax/1E6), xlab=xlab.text, ylab=ylab.text, main=main.text)
   points(bed.gc.chr$START/1E6, rt.chr$RT, col="grey", cex=0.3)
   abline(h=0, lwd=0.5, col="grey")
   
   ## Plot right-/left-leading positions
   spline <- smooth.spline(x=bed.gc.chr$START, y=rt.chr$RT)
   right.idx  <- which(bed.gc.rt.chr$RIGHT_LEADING > 505)
   left.idx   <- which(bed.gc.rt.chr$RIGHT_LEADING < 495)
   origin.idx <- setdiff(c(1:nrow(bed.gc.rt.chr)), c(right.idx, left.idx))
   
   points(spline$x[left.idx]/1E6, spline$y[left.idx], col="steelblue1", pch=16, cex=0.3)
   points(spline$x[right.idx]/1E6, spline$y[right.idx], col="sandybrown", pch=16, cex=0.3)
   points(spline$x[origin.idx]/1E6, spline$y[origin.idx], col="red", pch=16, cex=0.3)

   ## Plot cytobands and legend
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.4, col="lightgrey")
 
   legend("bottomright", c("Left-leading", "Left / Right", "Right-leading"), col=c("steelblue1", "red", "sandybrown"), pch=16, cex=0.8, horiz=T)
   dev.off()
}

for (c in 2:2) {
   chr <- chrs[c]

   ## Replication origins   
   load(file=file.path(wd.rt.data, paste0("bed.gc.rt_", base, "_bstrps", bstrps, "_", chr, ".RData")))
   bed.gc.rt.chr <- bed.gc.rt.chr[,1001:1003]
   bed.gc.chr <- bed.gc[rownames(bed.gc.rt.chr),]
   plotRO(wd.rt.plots, BASE, chr, NA, NA, bed.gc.rt.chr, bed.gc.chr, "png")
   
   ## Replication timing
   rt.chr <- readTable(file.path(wd.rt, "data", paste0(base, "_rpkm.corr.gc.d.rt_", chr, "_", BASE, "-", BASE, "_n101-92.txt.gz")), header=T, rownames=T, sep="\t") 
   rt.chr <- rt.chr[rownames(bed.gc.rt.chr),]
   plotRT(wd.rt.plots, BASE, chr, NA, NA, rt.chr, bed.gc.chr, bed.gc.rt.chr, "png")
}










load(file.path(wd.rt.data, "ensGene.rt_sclc_bstrp1000.RData"))
ensGene.rt.start.right <- subset(ensGene.rt.start, RT == 1)
ensGene.rt.start.left  <- subset(ensGene.rt.start, RT == -1)
ensGene.rt.start.zero  <- subset(ensGene.rt.start, RT == 0)
ensGene.rt.end.right   <- subset(ensGene.rt.end, RT == 1)
ensGene.rt.end.left    <- subset(ensGene.rt.end, RT == -1)
ensGene.rt.end.zero    <- subset(ensGene.rt.end, RT == 0)

overlaps.right <- intersect(rownames(ensGene.rt.start.right), rownames(ensGene.rt.end.right))
overlaps.left  <- intersect(rownames(ensGene.rt.start.left), rownames(ensGene.rt.end.left))

## TSS & TES
file.name <- file.path(wd.rt.plots, "hist_ensGene.rt.bstrp1000_Ensembl_TSS.png")
png(file.name, height=6, width=6, units="in", res=300)
hist(ensGene.rt.start$RIGHT_LEADING, xlab="Right-leading", main=c("TSS (Ensembl)", paste0("n=", nrow(ensGene.rt.start))), breaks=200)
dev.off()

file.name <- file.path(wd.rt.plots, "hist_ensGene.rt.bstrp1000_Ensembl_TES.png")
png(file.name, height=6, width=6, units="in", res=300)
hist(ensGene.rt.start$RIGHT_LEADING, xlab="Right-leading", main=c("TES (Ensembl)", paste0("n=", nrow(ensGene.rt.start))), breaks=200)
dev.off()

## Right-leading
file.name <- file.path(wd.rt.plots, "hist_ensGene.rt.bstrp1000_Ensembl_Right.png")
png(file.name, height=6, width=6, units="in", res=300)
hist(ensGene.rt.start[overlaps.right,]$RIGHT_LEADING, xlab="Right-leading", main=c("Right-leading (Ensembl)", paste0("n=", length(overlaps.right))), breaks=200)
dev.off()

file.name <- file.path(wd.rt.plots, "hist_ensGene.rt.bstrp1000_Ensembl_Left.png")
png(file.name, height=6, width=6, units="in", res=300)
hist(ensGene.rt.start[overlaps.left,]$LEFT_LEADING, xlab="Left-leading", main=c("Left-leading (Ensembl)", paste0("n=", length(overlaps.left))), breaks=200)
dev.off()

file.name <- file.path(wd.rt.plots, "hist_ensGene.rt.bstrp1000_Ensembl_TSS+TES.png")
png(file.name, height=6, width=6, units="in", res=300)
hist(ensGene.rt.start[c(overlaps.right, overlaps.left),]$RIGHT_LEADING, xlab="Right-leading", main=c("TSS+TES (Ensembl)", paste0("n=", length(c(overlaps.right, overlaps.left)))), breaks=200)
dev.off()














##
BASE <- "SCLC"
base <- "sclc"
wd.anlys <- file.path(wd, BASE, "analysis")
load(file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))

tpm.gene.input <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=F, proteinCodingOnly=F, proteinCodingNonRedundantOnly=F)
file.name <- file.path(wd.rt.plots, "hist_ensGene.rt.bstrp1000_Ensembl_SCLC_TSS_breaks200.png")
png(file.name, height=6, width=6, units="in", res=300)
hist(ensGene.rt.start[rownames(tpm.gene.input),]$RIGHT_LEADING, xlab="Right-leading", main=c("TSS", paste0("SCLC (n=", nrow(tpm.gene.input), ")")), breaks=200)
dev.off()

tpm.gene.input <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=T, proteinCodingOnly=F, proteinCodingNonRedundantOnly=F)
file.name <- file.path(wd.rt.plots, "hist_ensGene.rt.bstrp1000_Ensembl_SCLC_AUTO_TSS_breaks200.png")
png(file.name, height=6, width=6, units="in", res=300)
hist(ensGene.rt.start[rownames(tpm.gene.input),]$RIGHT_LEADING, xlab="Right-leading", main=c("TSS", paste0("SCLC autosome (n=", nrow(tpm.gene.input), ")")), breaks=200)
dev.off()

##
BASE <- "CLL"
base <- "cll"
wd.anlys <- file.path(wd, BASE, "analysis")
load(file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))

tpm.gene.input <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=F, proteinCodingOnly=F, proteinCodingNonRedundantOnly=F)
file.name <- file.path(wd.rt.plots, "hist_ensGene.rt.bstrp1000_Ensembl_CLL_TSS_breaks200.png")
png(file.name, height=6, width=6, units="in", res=300)
hist(ensGene.rt.start[rownames(tpm.gene.input),]$RIGHT_LEADING, xlab="Right-leading", main=c("TSS", paste0("CLL (n=", nrow(tpm.gene.input), ")")), breaks=200)
dev.off()

tpm.gene.input <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=T, proteinCodingOnly=F, proteinCodingNonRedundantOnly=F)
file.name <- file.path(wd.rt.plots, "hist_ensGene.rt.bstrp1000_Ensembl_CLL_AUTO_TSS_breaks200.png")
png(file.name, height=6, width=6, units="in", res=300)
hist(ensGene.rt.start[rownames(tpm.gene.input),]$RIGHT_LEADING, xlab="Right-leading", main=c("TSS", paste0("CLL autosome (n=", nrow(tpm.gene.input), ")")), breaks=200)
dev.off()



## TSS & TES
file.name <- file.path(wd.rt.plots, "ensGene.rt.bstrp1000_TSS.png")
#d <- getDensityCount(log10(ensGene.rt.start$RIGHT_LEADING/1000))
d <- density(ensGene.rt.start$RIGHT_LEADING, bw=0.05)
d$y <- d$y/sum(d$y)*d$n   ## Convert to counts
ymax <- max(d$y)
png(file.name, height=6, width=6, units="in", res=300)
plot(d, ylab="Frequency", xlab="Right-leading (%)", main="TSS (Ensembl)", ylim=c(0, ymax))
rug(jitter(d$y))
mtext(paste0("n=", nrow(ensGene.rt.start)), cex=0.8, font=3, line=0.5)
dev.off()



file.name <- file.path(wd.rt.plots, "ensGene.rt.bstrp1000_TSS_1000.png")
d <- getDensityCount(ensGene.rt.end$RIGHT_LEADING/10)
ymax <- max(d$y)
png(file.name, height=6, width=6, units="in", res=300)
plot(d, ylab="Frequency", xlab="Right-leading (%)", main="TES (Ensembl)", ylim=c(0, ymax))
dev.off()


## Right-leading (Negative slope)
file.name <- file.path(wd.rt.plots, "ensGene.rt.bstrp1000_TSS.png")
#d <- getDensityCount(log10(ensGene.tx.rt.start.right$RATIO))
d <- density(-log10(ensGene.rt.start.right$RATIO), adjust=1)
d$y <- d$y/sum(d$y)*d$n   ## Convert to counts
ymax <- max(d$y)
png(file.name, height=6, width=6, units="in", res=300)
plot(d, ylab="Frequency", xlab="Ratio (-log10)", main="TSS (Right-leading)", ylim=c(0, ymax))
dev.off()

file.name <- "/Users/tpyang/Work/uni-koeln/tyang2/SCLC/analysis/replication/sclc-wgs-rt/data/bstrp/ensGene.tx.rt.bstrp1000_right_TES.png"
d <- getDensityCount(-log10(ensGene.tx.rt.end.right$RATIO))
#ymax <- max(d$y)
png(file.name, height=6, width=6, units="in", res=300)
plot(d, ylab="Frequency", xlab="Ratio (-log10)", main="TES (Right-leading)", ylim=c(0, ymax))
dev.off()

file.name <- "/Users/tpyang/Work/uni-koeln/tyang2/SCLC/analysis/replication/sclc-wgs-rt/data/bstrp/ensGene.tx.rt.bstrp1000_right_TSS-TES.png"
png(file.name, height=6, width=6, units="in", res=300)
overlaps.right <- intersect(rownames(ensGene.tx.rt.start.right), rownames(ensGene.tx.rt.end.right))
plot(-log10(ensGene.tx.rt.start.right[overlaps.right,]$RATIO), -log10(ensGene.tx.rt.end.right[overlaps.right,]$RATIO), ylab="TES (-log10)", xlab="TSS (-log10)", main="TSS vs TES (Right-leading)")
dev.off()

##
file.name <- "/Users/tpyang/Work/uni-koeln/tyang2/SCLC/analysis/replication/sclc-wgs-rt/data/bstrp/ensGene.tx.rt.bstrp1000_left_TSS.png"
d <- getDensityCount(-log10(ensGene.tx.rt.start.left$RATIO))
#ymax <- max(d$y)
png(file.name, height=6, width=6, units="in", res=300)
plot(d, ylab="Frequency", xlab="Ratio (-log10)", main="TSS (Left-leading)", ylim=c(0, ymax))
dev.off()

file.name <- "/Users/tpyang/Work/uni-koeln/tyang2/SCLC/analysis/replication/sclc-wgs-rt/data/bstrp/ensGene.tx.rt.bstrp1000_left_TES.png"
d <- getDensityCount(-log10(ensGene.tx.rt.end.left$RATIO))
#ymax <- max(d$y)
png(file.name, height=6, width=6, units="in", res=300)
plot(d, ylab="Frequency", xlab="Ratio (-log10)", main="TES (Left-leading)", ylim=c(0, ymax))
dev.off()

file.name <- "/Users/tpyang/Work/uni-koeln/tyang2/SCLC/analysis/replication/sclc-wgs-rt/data/bstrp/ensGene.tx.rt.bstrp1000_left_TSS-TES.png"
png(file.name, height=6, width=6, units="in", res=300)
overlaps.left <- intersect(rownames(ensGene.tx.rt.start.left), rownames(ensGene.tx.rt.end.left))
plot(-log10(ensGene.tx.rt.start.left[overlaps.left,]$RATIO), -log10(ensGene.tx.rt.end.left[overlaps.left,]$RATIO), ylab="TES (-log10)", xlab="TSS (-log10)", main="TSS vs TES (Left-leading)")
dev.off()
