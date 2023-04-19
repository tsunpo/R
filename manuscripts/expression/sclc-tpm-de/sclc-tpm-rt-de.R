# =============================================================================
# Manuscript   : 
# Chapter      : 
# Name         : manuscripts/expression/sclc-tpm-rt-de.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 27/05/20
# =============================================================================
wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Commons.R", "Transcription.R", "TranscriptionReplicationInteraction.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.ensGene.bed.1kb.RData"))

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/ty2/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE <- "SCLC"
base <- tolower(BASE)

wd.wgs   <- file.path(wd, BASE, "ngs/WGS")
wd.rna   <- file.path(wd, BASE, "ngs/RNA")
wd.anlys <- file.path(wd, BASE, "analysis")
wd.meta  <- file.path(wd, BASE, "metadata", "George 2015")

wd.de       <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.gsea  <- file.path(wd.de, "gsea")
wd.de.plots <- file.path(wd.de, "plots")

samples.rna <- readTable(file.path(wd.rna, "sclc_rna_n81.list"), header=F, rownames=T, sep="")
samples.wgs <- readTable(file.path(wd.wgs, "sclc_wgs_n101.txt"), header=T, rownames=T, sep="")
overlaps <- intersect(rownames(samples.rna), rownames(samples.wgs))

samples.tpm.sclc  <- samples.wgs[overlaps,]
samples.tpm.sclc$M2 <- samples.tpm.sclc$M2 + 1
#samples.sclc.tpm$M2 <- as.factor(samples.sclc.tpm$M2)
# > length(which(samples.sclc.tpm$M2 == 1))
# [1] 35
# > length(which(samples.sclc.tpm$M2 == 0))
# [1] 35

#load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.RData")))
#load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median1.RData")))
#expressed <- rownames(tpm.gene)
#load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.r5p47.RData")))
tpm.gene <- tpm.gene[, rownames(samples.sclc.tpm)]   ## VERY VERY VERY IMPORTANT!!!
tpm.gene.log2   <- log2(tpm.gene + 1)
tpm.gene.log2.m <- getLog2andMedian(tpm.gene, 1)
nrow(tpm.gene.log2)
# [1] 34908
# [1] 23572
# [1] 19131
# [1] 15803

q4 <- quantile(as.numeric(tpm.gene.log2["ENSG00000140451",]))

ylim <- c(0, quantile(as.numeric(tpm.gene.log2["ENSG00000140451",]))[5])
ylim <- c(0, 6.5)
file.name <- paste0("tpm.gene.log2_TCGA_SCLC_PIF1_0.5_1.3.4.5_3.5")
plotBox2(wd.de.plots, file.name, as.numeric(tpm.gene.lung.log2["ENSG00000140451",]), as.numeric(tpm.gene.log2["ENSG00000140451",]), main="PIF1 expression", names=c("TCGA normal lung", "Lung-SCLC"), cols=c("black", red), ylim, ylab=text.Log2.TPM.1)
# [1] 1.871352e-18

# -----------------------------------------------------------------------------
# GTEx
# Last Modified: 22/03/23
# -----------------------------------------------------------------------------
gtex <- readTable(file.path("/Users/ty2/Work/uni-koeln/tyang2/ICGC/metadata/GTEx", "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"), header=T, rownames=T, sep="\t")
gtex <- gtex["ENSG00000140451.12", -c(1,2)]
gtex <- gtex[, order(as.numeric(gtex))]
gtex <- gtex[, colnames(gtex)[38:54]]

tpm <- readTable(file.path("/Users/ty2/Work/uni-koeln/tyang2/ICGC/metadata/GTEx", "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct"), header=T, rownames=T, sep="\t")
tpm <- tpm["ENSG00000140451.12", -c(1,2)]

samples <- readTable(file.path("/Users/ty2/Work/uni-koeln/tyang2/ICGC/metadata/GTEx", "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"), header=T, rownames=T, sep="\t")
samples <- samples[, c("SMTS", "SMTSD")]

names <- c()
labels <- c()
trait <- c()
expr  <- c()
for (h in 1:ncol(gtex)) {
	  hist <- unlist(strsplit(colnames(gtex)[h], "\\.\\.\\."))
	  if (length(hist) != 1) {
	  	  hist[1] <- gsub("\\.", " ", hist[1])
	  	  hist[2] <- gsub("\\.", " ", hist[2])
	  	  hist <- paste0(hist, collapse=" - ")
	  }
	  if (hist == "Cells - EBV transformed lymphocytes")
	  	  hist <- "Cells - EBV-transformed lymphocytes"
	  if (hist == "Fallopian.Tube")
	  	  hist <- "Fallopian Tube"
	  samples.1 <- gsub("-", ".", rownames(subset(samples, SMTSD == hist)))
	  samples.1 <- intersect(samples.1, colnames(tpm))

	  if (hist == "Esophagus - Gastroesophageal Junction")
	    	hist <- "Esophagus - Gastroesophageal junction"
	  if (hist == "Small Intestine - Terminal Ileum")
	  	  hist <- "Small intestine - Terminal ileum"
	  if (hist == "Fallopian Tube")
	  	  hist <- "Fallopian tube"
	  names <- c(names, hist)
	  labels <- c(labels, "")
	  
	  tpm.1 <- tpm[, samples.1]
	  if (h == 1)
	     trait <- rep(h-1, ncol(tpm.1))
	  else
	     trait <- c(trait, rep(h-1, ncol(tpm.1)))
	  
	  expr <- c(expr, log2(as.numeric(tpm.1) + 1))
}
trait <- as.factor(trait)

cols <- "black"
ylab <-	text.Log2.TPM.1
pdf(file.path(wd.de.plots, paste0("GTEx_top17_1_1.3.4.5_abline.pdf")), height=5, width=10)
par(mar=c(16, 4.5, 2.1, 1.4))
boxplot(expr ~ trait, outline=F, xlab="", xaxt="n", col="white", boxcol=cols, whiskcol=cols, outcol=cols, medcol=cols, staplecol=cols, ylab=ylab, main="PIF1 expression in top 17 GTEx tissues", ylim=ylim, cex.axis=1.3, cex.lab=1.4, cex.main=1.5)

abline(h=q4[3], lty=5, lwd=2, col=red)

axis(side=1, at=seq(1, 17, by=1), labels=labels, cex.axis=1.4)
text(labels=names, x=1:17, y=par("usr")[3] - 1, srt=45, adj=0.965, xpd=NA, cex=1.4)
dev.off()

# -----------------------------------------------------------------------------
# RFD vs. TPM
# Last Modified: 16/02/23; 31/10/18
# -----------------------------------------------------------------------------
nrds.RT.NRFD <- nrds.RT.NRFD.sclc   ## 17/05/20 WHY nrds.RT.NRFD.sclc.nl? ## MUY MUY IMPORTANTE!!
load(file=file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.d.rt.log2s.nrfd.", kb, "kb_", "m2-m1", ".RData")))
nrow(nrds.RT.NRFD)
# [1] 2650083

tpm.gene.log2.m <- tpm.gene.log2.m
nrow(tpm.gene.log2.m)
# [1] 34908
# [1] 23572
# [1] 15803

tpm.gene.log2.m.rfd <- getTRC(tpm.gene.log2.m, nrds.RT.NRFD)
tpm.gene.log2.m.rfd$length <- abs(tpm.gene.log2.m.rfd$end_position - tpm.gene.log2.m.rfd$start_position)
nrow(tpm.gene.log2.m.rfd)
# [1] 31616
# [1] 22023
# [1] 18007
# [1] 14986
nrow(subset(tpm.gene.log2.m.rfd, TRC == 0))
# [1] 3
# [1] 1
# [1] 1
# [1] 1

###
## TTR and CTR (IZ +TZ)
#save(tpm.gene.log2.m.rfd, file=file.path(wd.de.data, "tpm_gene_log2_m_rfd0+1.RData"))
#save(tpm.gene.log2.m.rfd, file=file.path(wd.de.data, "tpm_gene_median0_log2+1_m_rfd0.RData"))
#save(tpm.gene.log2.m.rfd, file=file.path(wd.de.data, "tpm_gene_r5p47_log2_m_rfd0.RData"))
save(tpm.gene.log2.m.rfd, file=file.path(wd.de.data, "tpm_gene_median1_log2+1_m_rfd0.RData"))

#load(file=file.path(wd.de.data, "tpm_gene_log2+1_m_rfd0.RData"))
#load(file=file.path(wd.de.data, "tpm_gene_median0_log2+1_m_rfd0.RData"))
#load(file=file.path(wd.de.data, "tpm_gene_r5p47_log2_m_rfd0.RData"))
load(file=file.path(wd.de.data, "tpm_gene_median1_log2+1_m_rfd0.RData"))

#file.name <- file.path(wd.de.data, "tpm_gene_log2+1_m_rfd.RData")
#file.name <- file.path(wd.de.data, "tpm_gene_median0_log2+1_m_rfd.RData")
#file.name <- file.path(wd.de.data, "tpm_gene_r5p47_log2_m_rfd.RData")
file.name <- file.path(wd.de.data, "tpm_gene_median1_log2+1_m_rfd.RData")
setTRC(tpm.gene.log2.m.rfd, rfd=0.9, file.name)
load(file.name)

## CTR
length(which(tpm.gene.log2.m.rfd.ctr.iz$GENE_NRFD < 0))
# [1] 84
# [1] 71
# [1] 60
# [1] 50
length(which(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD > 0))
# [1] 95
# [1] 85
# [1] 80
# [1] 70

# -----------------------------------------------------------------------------
# RFD vs. TPM (median1)
# Last Modified: 16/02/23; 11/02/21; 27/11/20; 01/09/20; 29/05/20
# -----------------------------------------------------------------------------
file.name <- file.path(wd.de.data, "tpm_gene_median1_log2+1_m_rfd.RData")
load(file.name)

ylim <- c(min(tpm.gene.log2.m.rfd$MEDIAN), 11)
file.name <- paste0("boxplot3_sclc_tpm.gene_log2+1_RFD_test")
plotBox3(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr, tpm.gene.log2.m.rfd.ctr.tz, tpm.gene.log2.m.rfd.ctr.iz, main="SCLC expression", names=c("TTR", "TZ", "IZ"), cols=c("black", blue, red), ylim)

# -----------------------------------------------------------------------------
# DeepMap SCLC-CL
# Last Modified: 23/02/23; 17/02/23
# -----------------------------------------------------------------------------
subtypes <- read.csv("/Users/ty2/Work/sanger/ty2/Pan-Pody/metadata/DepMap/22Q2/sample_info.csv", header=T)
rownames(subtypes) <- subtypes$DepMap_ID
subtypes <- subset(subtypes, Subtype == "Small Cell Lung Cancer (SCLC)")
subtypes <- subset(subtypes, primary_or_metastasis != "Metastasis")
sclcs <- sort(rownames(subtypes))
writeTable(sclcs, "/Users/ty2/Work/sanger/ty2/Pan-Pody/metadata/DepMap/22Q2/sample_info_SCLC.txt", colnames=F, rownames=F, sep="\t")
#commons <- readTable("/Users/ty2/Work/sanger/ty2/Pan-Pody/metadata/DepMap/22Q2/common_essentials.txt", header=F, rownames=F, sep="")[, 1]
#nons    <- readTable("/Users/ty2/Work/sanger/ty2/Pan-Pody/metadata/DepMap/22Q2/nonessentials.txt", header=F, rownames=F, sep="")[, 1]

effects <- read.csv("/Users/ty2/Work/sanger/ty2/Pan-Pody/metadata/DepMap/22Q4/CRISPRGeneEffect.csv", header=T)
rownames(effects) <- effects$X
#effects <- read.csv("/Users/ty2/Work/sanger/ty2/Pan-Pody/metadata/DepMap/22Q2/CRISPR_gene_effect.csv", header=T)
#rownames(effects) <- effects$DepMap_ID
effects <- effects[,-1]
# > dim(effects)
# [1]  1078 17453
# > dim(effects)
# [1]  1086 17386
colnames(effects) <- mapply(x = 1:length(colnames(effects)), function(x) strsplit0(toString(colnames(effects)[x]), "\\.\\.")[1])
effects.sclc <- effects[intersect(sclcs, rownames(effects)), ]
# > dim(effects.sclc)
# [1]    20 17453
# > dim(effects.sclc)
# [1]    20 17386
effects.sclc.t <- as.data.frame(t(effects.sclc))
effects.sclc.t.m <- effects.sclc.t
effects.sclc.t.m$MEDIAN <- mapply(x = 1:nrow(effects.sclc.t), function(x) median(as.numeric(effects.sclc.t[x,])))

depends <- read.csv("/Users/ty2/Work/sanger/ty2/Pan-Pody/metadata/DepMap/22Q4/CRISPRGeneDependency.csv", header=T)
rownames(depends) <- depends$ModelID
depends <- depends[,-1]
colnames(depends) <- mapply(x = 1:length(colnames(depends)), function(x) strsplit0(toString(colnames(depends)[x]), "\\.\\.")[1])
# > dim(depends)
# [1]  1078 17453
depends.sclc <- depends[intersect(sclcs, rownames(depends)), ]
# > dim(effects.sclc)
# [1]    20 17453
depends.sclc.t <- as.data.frame(t(depends.sclc))
depends.sclc.t.m <- depends.sclc.t
depends.sclc.t.m$MEDIAN <- mapply(x = 1:nrow(depends.sclc.t), function(x) median(as.numeric(depends.sclc.t[x,])))

# -----------------------------------------------------------------------------
# 
# Last Modified: 23/02/23
# -----------------------------------------------------------------------------
#overlaps <- intersect(rownames(effects.sclc.t.m), rownames(depends.sclc.t.m))
#sclc.nl.rfd <- c(sclc.nl.ttr, sclc.nl.ctr.iz, sclc.nl.ctr.tz)
#sclc.rfd <- c(sclc.ttr, sclc.ctr.iz, sclc.ctr.tz)
#nbl.rfd  <- c(nbl.ttr, nbl.ctr.iz, nbl.ctr.tz)
#cll.rfd  <- c(cll.ttr, cll.ctr.iz, cll.ctr.tz)
#overlaps <- intersect(overlaps, intersect(intersect(intersect(sclc.nl.rfd, sclc.rfd), nbl.rfd), cll.rfd)) 
# > length(overlaps)
# [1] 14669
load("/Users/ty2/Work/uni-koeln/tyang2/SCLC/analysis/expression/kallisto/sclc-tpm-de/data/tpm_gene_median1_log2+1_m_rfd.RData")

sclc.ttr <- intersect(rownames(effects.sclc.t.m), tpm.gene.log2.m.rfd.ttr$external_gene_name)
sclc.ctr.tz <- intersect(rownames(effects.sclc.t.m), tpm.gene.log2.m.rfd.ctr.tz$external_gene_name)
sclc.ctr.iz <- intersect(rownames(effects.sclc.t.m), tpm.gene.log2.m.rfd.ctr.iz$external_gene_name)
length(sclc.ctr.iz)
# [1] 2573   ## median0
# [1] 1790   ## median1

#sclc.ttr <- intersect(overlaps, sclc.ttr)
#sclc.ctr.iz <- intersect(overlaps, sclc.ctr.iz)
#sclc.ctr.tz <- intersect(overlaps, sclc.ctr.tz)
# > length(sclc.ctr.iz)
# [1] 1969

# quantile(effects.sclc.t.m[sclc.ctr.iz,]$MEDIAN)
#          0%         25%         50%         75%        100% 
# -3.82070041 -0.16319426 -0.05735207  0.01382089  0.23845697
ylim <- c(-0.67, 0.4)
file.name <- paste0("DepMap_effects_median1_SCLC-CL_RFD-vs-SCLC")
plotBoxDepMap(wd.de.plots, file.name, effects.sclc.t.m[sclc.ttr,], effects.sclc.t.m[sclc.ctr.tz,], effects.sclc.t.m[sclc.ctr.iz,], main="DepMap SCLC-CL", names=c("TTR", "TZ", "IZ"), cols=c("white", blue.lighter, red.lighter), ylim, ylab="Gene effect")

# quantile(depends.sclc.t.m[sclc.ctr.iz,]$MEDIAN)
#          0%         25%         50%         75%        100% 
# 0.003388999 0.021971653 0.042265738 0.104133349 1.000000000
#ylim <- c(0, 0.23)
#file.name <- paste0("DepMap_depends_boxplot3_SCLC-CL_RFD")
#plotBoxDepMap(wd.de.plots, file.name, depends.sclc.t.m[sclc.ttr,], depends.sclc.t.m[sclc.ctr.tz,], depends.sclc.t.m[sclc.ctr.iz,], main="DepMap SCLC-CL", names=c("TTR", "TZ", "IZ"), cols=c("black", blue, red), ylim, ylab="Gene dependency")

testU(effects.sclc.t.m[sclc.ttr,]$MEDIAN, effects.sclc.t.m[sclc.ctr.iz,]$MEDIAN)
testU(effects.sclc.t.m[sclc.nl.ttr,]$MEDIAN, effects.sclc.t.m[sclc.nl.ctr.iz,]$MEDIAN)
testU(effects.sclc.t.m[cll.ttr,]$MEDIAN, effects.sclc.t.m[cll.ctr.iz,]$MEDIAN)
testU(effects.sclc.t.m[nbl.ttr,]$MEDIAN, effects.sclc.t.m[nbl.ctr.iz,]$MEDIAN)

testU(effects.sclc.t.m[sclc.ctr.tz,]$MEDIAN, effects.sclc.t.m[sclc.ctr.iz,]$MEDIAN)
testU(effects.sclc.t.m[sclc.nl.ctr.tz,]$MEDIAN, effects.sclc.t.m[sclc.nl.ctr.iz,]$MEDIAN)
testU(effects.sclc.t.m[cll.ctr.tz,]$MEDIAN, effects.sclc.t.m[cll.ctr.iz,]$MEDIAN)
testU(effects.sclc.t.m[nbl.ctr.tz,]$MEDIAN, effects.sclc.t.m[nbl.ctr.iz,]$MEDIAN)

save(effects.sclc.t.m, effects.nbl.t.m, effects.cll.t.m, sclc.ttr, sclc.ctr.iz, sclc.ctr.tz, sclc.nl.ttr, sclc.nl.ctr.iz, sclc.nl.ctr.tz, nbl.ttr, nbl.ctr.iz, nbl.ctr.tz, cll.ttr, cll.ctr.iz, cll.ctr.tz, file=file.path(wd.de.data, "DepMap_tpm_gene_median1_log2+1_m_rfd.RData"))

## Gene effect and gene expression (IZ)
ylim <- c(-3.830851, 0.3436195)

duplicates <- as.data.frame(table(subset(tpm.gene.log2.m.rfd.ctr.iz, external_gene_name %in% sclc.ctr.iz)$external_gene_name))
duplicates <- as.vector(subset(duplicates, Freq > 1)$Var1)
genes <- setdiff(sclc.ctr.iz, duplicates)

ens <- subset(tpm.gene.log2.m.rfd.ctr.iz, external_gene_name %in% genes)
rownames(ens) <- ens$external_gene_name
ens <- ens[genes,]$ensembl_gene_id

file.name <- file.path(wd.de.plots, "TPM-vs-Effect_SCLC-IZ.pdf")
xlab.text <- expression("log" * ""[2] * "(TPM + 1)")
ylab.text <- "Gene effect"
x <- as.numeric(tpm.gene.log2.m[ens,]$MEDIAN)
y <- effects.sclc.t.m[genes,]$MEDIAN
plotCorrelation(file.name, "SCLC IZ", xlab.text, ylab.text, x, y, pos="bottomright", cols=c(red.lighter, red), 6, ylim)

## Gene effect and gene expression (TZ)
duplicates <- as.data.frame(table(subset(tpm.gene.log2.m.rfd.ctr.tz, external_gene_name %in% sclc.ctr.tz)$external_gene_name))
duplicates <- as.vector(subset(duplicates, Freq > 1)$Var1)
genes <- setdiff(sclc.ctr.tz, duplicates)

ens <- subset(tpm.gene.log2.m.rfd.ctr.tz, external_gene_name %in% genes)
rownames(ens) <- ens$external_gene_name
ens <- ens[genes,]$ensembl_gene_id

file.name <- file.path(wd.de.plots, "TPM-vs-Effect_SCLC-TZ.pdf")
xlab.text <- expression("log" * ""[2] * "(TPM + 1)")
ylab.text <- "Gene effect"
x <- as.numeric(tpm.gene.log2.m[ens,]$MEDIAN)
y <- effects.sclc.t.m[genes,]$MEDIAN
plotCorrelation(file.name, "SCLC TZ", xlab.text, ylab.text, x, y, pos="bottomright", cols=c(blue.lighter, blue), 6, ylim)

## Gene effect and gene expression (TTR)
duplicates <- as.data.frame(table(subset(tpm.gene.log2.m.rfd.ttr, external_gene_name %in% sclc.ttr)$external_gene_name))
duplicates <- as.vector(subset(duplicates, Freq > 1)$Var1)
genes <- setdiff(sclc.ttr, duplicates)

ens <- subset(tpm.gene.log2.m.rfd.ttr, external_gene_name %in% genes)
rownames(ens) <- ens$external_gene_name
ens <- ens[genes,]$ensembl_gene_id

file.name <- file.path(wd.de.plots, "TPM-vs-Effect_SCLC-TTR.pdf")
xlab.text <- expression("log" * ""[2] * "(TPM + 1)")
ylab.text <- "Gene effect"
x <- as.numeric(tpm.gene.log2.m[ens,]$MEDIAN)
y <- effects.sclc.t.m[genes,]$MEDIAN
plotCorrelation(file.name, "SCLC TTR", xlab.text, ylab.text, x, y, pos="bottomright", cols=c("lightgray", "black"), 6, ylim)


# -----------------------------------------------------------------------------
# Remove shared genes
# Last Modified: 01/03/23
# -----------------------------------------------------------------------------
length(sclc.ctr.iz)
# [1] 2163
length(sclc.nl.ctr.iz)
# [1] 1547
length(cll.ctr.iz)
# [1] 2537
length(nbl.ctr.iz)
# [1] 2679

length(sclc.ctr.iz)   ## median1
# [1] 1790
length(sclc.nl.ctr.iz)
# [1] 1262
length(cll.ctr.iz)
# [1] 1801
length(nbl.ctr.iz)
# [1] 2177

sclc.shared <- intersect(sclc.ctr.iz, sclc.nl.ctr.iz)
length(sclc.shared)
# [1] 1067
# [1] 878
sclc.specific <- setdiff(sclc.ctr.iz, intersect(sclc.ctr.iz, sclc.nl.ctr.iz))
length(sclc.specific)
# [1] 1096
# [1] 912
testU(effects.sclc.t.m[sclc.shared,]$MEDIAN, effects.sclc.t.m[sclc.specific,]$MEDIAN)
# [1] 0.6670427
# [1] 0.6670149

sclc.specific <- setdiff(sclc.specific, intersect(sclc.specific, cll.ctr.iz))
length(sclc.specific)
# [1] 721
# [1] 650

sclc.specific <- setdiff(sclc.specific, intersect(sclc.specific, nbl.ctr.iz))
length(sclc.specific)
length(sclc.specific)
# [1] 382
# [1] 362

sclc.shared <- setdiff(sclc.ctr.iz, sclc.specific)
length(sclc.shared)
# [1] 1781
# [1] 1428
testU(effects.sclc.t.m[sclc.shared,]$MEDIAN, effects.sclc.t.m[sclc.specific,]$MEDIAN)
# [1] 0.3692903
# [1] 0.0997264

# -----------------------------------------------------------------------------
# Shared genes (IZ)
# Last Modified: 205/03/23; 4/02/23
# -----------------------------------------------------------------------------
shared.3 <- intersect(intersect(sclc.ctr.iz, nbl.ctr.iz), cll.ctr.iz)
length(shared.3)
# [1] 299

shared.2 <- intersect(sclc.ctr.iz, sclc.nl.ctr.iz)
length(shared.2)
# [1] 878

sclc.specific <- setdiff(sclc.ctr.iz, unique(c(shared.2, shared.3)))
length(sclc.specific)
# [1] 785

ylim <- c(-0.8, 0.4)
file.name <- paste0("DepMap_effects_median1_SCLC-CL_IZ_SCLC-specific_0.8")
plotBoxDepMapSolid(wd.de.plots, file.name, effects.sclc.t.m[shared.3,], effects.sclc.t.m[shared.2,], effects.sclc.t.m[sclc.specific,], main="DepMap SCLC-CL", names=c("SCLC & NB & CLL", "SCLC & SCLC-NL", "SCLC-specific IZ"), cols=c(red.lighter, red.light, red), ylim, ylab="Gene effect")
testU(effects.sclc.t.m[shared.3,]$MEDIAN, effects.sclc.t.m[sclc.specific,]$MEDIAN)
# [1] 0.2012056
testU(effects.sclc.t.m[shared.2,]$MEDIAN, effects.sclc.t.m[sclc.specific,]$MEDIAN)
# [1] 0.3338531

#ylim <- c(-0.8, 0.4)
#file.name <- paste0("DepMap_effects_median1_specific_NBL-CL_RFD_SCLC-specific")
#plotBoxDepMapSolid(wd.de.plots, file.name, effects.nbl.t.m[shared.3,], effects.nbl.t.m[shared.2,], effects.nbl.t.m[sclc.specific,], main="DepMap NB-CL", names=c("SCLC & NB & CLL", "SCLC & SCLC-NL", "SCLC-specific IZ"), cols=c(red.lighter, red.light, red), ylim, ylab="Gene effect")
#testU(effects.nbl.t.m[shared.3,]$MEDIAN, effects.nbl.t.m[sclc.specific,]$MEDIAN)
# [1] 0.01341057
#testU(effects.nbl.t.m[shared.2,]$MEDIAN, effects.nbl.t.m[sclc.specific,]$MEDIAN)
# [1] 0.1439723

#ylim <- c(-0.8, 0.4)
#file.name <- paste0("DepMap_effects_median1_specific_CLL-CL_RFD_SCLC-specific")
#plotBoxDepMapSolid(wd.de.plots, file.name, effects.cll.t.m[shared.3,], effects.cll.t.m[shared.2,], effects.cll.t.m[sclc.specific,], main="DepMap CLL-CL", names=c("SCLC & NB & CLL", "SCLC & SCLC-NL", "SCLC-specific IZ"), cols=c(red.lighter, red.light, red), ylim, ylab="Gene effect")
#testU(effects.cll.t.m[shared.3,]$MEDIAN, effects.cll.t.m[sclc.specific,]$MEDIAN)
# [1] 0.002959214
#testU(effects.cll.t.m[shared.2,]$MEDIAN, effects.cll.t.m[sclc.specific,]$MEDIAN)
# [1] 0.07217747

# -----------------------------------------------------------------------------
# Shared genes (TZ)
# Last Modified: 14/03/23; 4/02/23
# -----------------------------------------------------------------------------
shared.3 <- intersect(intersect(sclc.ctr.tz, nbl.ctr.tz), cll.ctr.tz)
length(shared.3)
# [1] 129

shared.2 <- intersect(sclc.ctr.tz, sclc.nl.ctr.tz)
length(shared.2)
# [1] 457

sclc.specific.tz <- setdiff(sclc.ctr.tz, unique(c(shared.2, shared.3)))
length(sclc.specific.tz)
# [1] 541

ylim <- c(-0.6, 0.4)
file.name <- paste0("DepMap_effects_median1_SCLC-CL_SCLC-specific-TZ")
plotBoxDepMapSolid(wd.de.plots, file.name, effects.sclc.t.m[shared.3,], effects.sclc.t.m[shared.2,], effects.sclc.t.m[sclc.specific.tz,], main="", names=c("SCLC & NB & CLL", "SCLC & SCLC-NL", "SCLC-specific TZ"), cols=c(blue.lighter, blue.light, blue), ylim, ylab="Gene effect")
testU(effects.sclc.t.m[shared.3,]$MEDIAN, effects.sclc.t.m[sclc.specific.tz,]$MEDIAN)
# [1] 0.1550376
testU(effects.sclc.t.m[shared.2,]$MEDIAN, effects.sclc.t.m[sclc.specific.tz,]$MEDIAN)
# [1] 0.9143408

median(effects.sclc.t.m[shared.3,]$MEDIAN)
# [1] -0.05958306
median(effects.sclc.t.m[shared.2,]$MEDIAN)
# [1] -0.04270429
median(effects.sclc.t.m[sclc.specific.tz,]$MEDIAN)
# [1] -0.04658705

# -----------------------------------------------------------------------------
# Shared genes (TTR)
# Last Modified: 14/03/23; 4/02/23
# -----------------------------------------------------------------------------
shared.3 <- intersect(intersect(sclc.ttr, nbl.ttr), cll.ttr)
length(shared.3)
# [1] 4717

shared.2 <- intersect(sclc.ttr, sclc.nl.ttr)
length(shared.2)
# [1] 9047

sclc.specific.ttr <- setdiff(sclc.ttr, unique(c(shared.2, shared.3)))
length(sclc.specific.ttr)
# [1] 424

ylim <- c(-0.6, 0.4)
file.name <- paste0("DepMap_effects_median1_SCLC-CL_SCLC-specific-TTR")
plotBoxDepMapSolid(wd.de.plots, file.name, effects.sclc.t.m[shared.3,], effects.sclc.t.m[shared.2,], effects.sclc.t.m[sclc.specific.ttr,], main="", names=c("SCLC & NB & CLL", "SCLC & SCLC-NL", "SCLC-specific TTR"), cols=c("white", "lightgray", "dimgray"), ylim, ylab="Gene effect")
testU(effects.sclc.t.m[shared.3,]$MEDIAN, effects.sclc.t.m[sclc.specific.ttr,]$MEDIAN)
# [1] 0.9881853
testU(effects.sclc.t.m[shared.2,]$MEDIAN, effects.sclc.t.m[sclc.specific.ttr,]$MEDIAN)
# [1] 0.04832639 
testU(effects.sclc.t.m[shared.3,]$MEDIAN, effects.sclc.t.m[shared.2,]$MEDIAN)
# [1] 2.518918e-07

median(effects.sclc.t.m[shared.3,]$MEDIAN)
# [1] -0.05555021
median(effects.sclc.t.m[shared.2,]$MEDIAN)
# [1] -0.04275472
median(effects.sclc.t.m[sclc.specific.ttr,]$MEDIAN)
# [1] -0.05887152

# -----------------------------------------------------------------------------
# Specific TTR genes
# Last Modified: 14/03/23; 4/02/23
# -----------------------------------------------------------------------------
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median1.RData")))

ids <- subset(ensGene, external_gene_name %in% sclc.specific.ttr)$ensembl_gene_id
length(ids)
# [1] 429
ids <- intersect(rownames(tpm.gene), ids)
length(ids)
# [1] 425

tpm.gene <- tpm.gene[ids, rownames(samples.sclc.tpm)]   ## VERY VERY VERY IMPORTANT!!!
tpm.gene.log2   <- log2(tpm.gene + 1)
tpm.gene.log2.m <- getLog2andMedian(tpm.gene, 1)
nrow(tpm.gene.log2.m)
# [1] 425

# -----------------------------------------------------------------------------
# Correlation bwteen TPM and in-silico sorting
# Last Modified: 19/02/23; 24/04/22; 12/12/19; 08/01/19; 17/08/17
# -----------------------------------------------------------------------------
colnames <- c("RHO", "P", "Q", "G1", "S", "LOG2_FC")   ##"ANOVA_P", "ANOVA_Q", 
de <- toTable(0, length(colnames), nrow(tpm.gene.log2), colnames)
rownames(de) <- rownames(tpm.gene.log2)

## SRC
de$RHO <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples.sclc.tpm$COR, method="spearman", exact=F)[[4]])
de$P   <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples.sclc.tpm$COR, method="spearman", exact=F)[[3]])
de <- de[!is.na(de$P),]

## Log2 fold change
de$G1 <- median00(tpm.gene.log2, rownames(subset(samples.sclc.tpm, M2 == 0)))
de$S  <- median00(tpm.gene.log2, rownames(subset(samples.sclc.tpm, M2 == 1)))
de$LOG2_FC <- de$S - de$G1

## FDR
library(qvalue)
de$Q <- qvalue(de$P)$qvalue
de <- de[order(de$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)   ## BE EXTRA CAREFUL!!

writeTable(de.tpm.gene, file.path(wd.de.data, "src_sclc_tpm-gene-median1_SORTING_q_n70_specific-425TTR.txt"), colnames=T, rownames=F, sep="\t")
save(de.tpm.gene, samples, file=file.path(wd.de.data, "src_sclc_tpm-gene-median1_SORTING_q_n70_specific-425TTR.RData"))
nrow(de.tpm.gene)
# [1] 425

# -----------------------------------------------------------------------------
# Shared TTR genes
# Last Modified: 14/03/23; 4/02/23
# -----------------------------------------------------------------------------
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median1.RData")))

ids <- subset(ensGene, external_gene_name %in% shared.2)$ensembl_gene_id
length(ids)
# [1] 9087
ids <- intersect(rownames(tpm.gene), ids)
length(ids)
# [1] 9057

tpm.gene <- tpm.gene[ids, rownames(samples.sclc.tpm)]   ## VERY VERY VERY IMPORTANT!!!
tpm.gene.log2   <- log2(tpm.gene + 1)
tpm.gene.log2.m <- getLog2andMedian(tpm.gene, 1)
nrow(tpm.gene.log2.m)
# [1] 15803   ## TPM > 1
# [1] 786     ## SCLC-specific IZ
# [1] 9057    ## SCLC+SCLC-NL TTR

# -----------------------------------------------------------------------------
# Correlation bwteen TPM and in-silico sorting
# Last Modified: 19/02/23; 24/04/22; 12/12/19; 08/01/19; 17/08/17
# -----------------------------------------------------------------------------
colnames <- c("RHO", "P", "Q", "G1", "S", "LOG2_FC")   ##"ANOVA_P", "ANOVA_Q", 
de <- toTable(0, length(colnames), nrow(tpm.gene.log2), colnames)
rownames(de) <- rownames(tpm.gene.log2)

## SRC
de$RHO <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples.sclc.tpm$COR, method="spearman", exact=F)[[4]])
de$P   <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples.sclc.tpm$COR, method="spearman", exact=F)[[3]])
de <- de[!is.na(de$P),]

## Log2 fold change
de$G1 <- median00(tpm.gene.log2, rownames(subset(samples.sclc.tpm, M2 == 0)))
de$S  <- median00(tpm.gene.log2, rownames(subset(samples.sclc.tpm, M2 == 1)))
de$LOG2_FC <- de$S - de$G1

## FDR
library(qvalue)
de$Q <- qvalue(de$P)$qvalue
de <- de[order(de$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)   ## BE EXTRA CAREFUL!!

writeTable(de.tpm.gene, file.path(wd.de.data, "src_sclc_tpm-gene-median1_SORTING_q_n70_shared-9057TTR.txt"), colnames=T, rownames=F, sep="\t")
save(de.tpm.gene, samples, file=file.path(wd.de.data, "src_sclc_tpm-gene-median1_SORTING_q_n70_shared-9057TTR.RData"))
nrow(de.tpm.gene)
# [1] 9057









# -----------------------------------------------------------------------------
# Shared genes
# Last Modified: 24/02/23
# -----------------------------------------------------------------------------
shared.3 <- intersect(intersect(sclc.ctr.iz, nbl.ctr.iz), cll.ctr.iz)
length(shared.3)
# [1] 454
# [1] 420
shared.2.1 <- setdiff(intersect(sclc.ctr.iz, nbl.ctr.iz), shared.3)
shared.2.2 <- setdiff(intersect(sclc.ctr.iz, cll.ctr.iz), shared.3)
shared.2 <- c(shared.2.1, shared.2.2)
length(shared.2)
# [1] 1060
# [1] 860

sclc.specific <- setdiff(setdiff(sclc.ctr.iz, shared.3), shared.2)
length(sclc.specific)
# [1] 649
# [1] 631

ylim <- c(-0.67, 0.4)
file.name <- paste0("DepMap_effects_median1_specific_SCLC-CL_RFD_SCLC-specific")
plotBoxDepMapSolid(wd.de.plots, file.name, effects.sclc.t.m[shared.3,], effects.sclc.t.m[shared.2,], effects.sclc.t.m[sclc.specific,], main="DepMap SCLC-CL", names=c("SCLC + NB + CLL", "SCLC + (NB/CLL)", "SCLC-specific IZ"), cols=c(red.lighter, red.light, red), ylim, ylab="Gene effect")
testU(effects.sclc.t.m[shared.3,]$MEDIAN, effects.sclc.t.m[sclc.specific,]$MEDIAN)
# [1] 0.7139454
# [1] 0.07118778
testU(effects.sclc.t.m[shared.2,]$MEDIAN, effects.sclc.t.m[sclc.specific,]$MEDIAN)
# [1] 0.2757434
# [1] 0.02013834

ylim <- c(-0.67, 0.4)
file.name <- paste0("DepMap_effects_median1_specific_NBL-CL_RFD_SCLC-specific")
plotBoxDepMapSolid(wd.de.plots, file.name, effects.nbl.t.m[shared.3,], effects.nbl.t.m[shared.2,], effects.nbl.t.m[sclc.specific,], main="DepMap NB-CL", names=c("SCLC + NB + CLL", "SCLC + (NB/CLL)", "SCLC-specific IZ"), cols=c(red.lighter, red.light, red), ylim, ylab="Gene effect")
testU(effects.nbl.t.m[shared.3,]$MEDIAN, effects.nbl.t.m[sclc.specific,]$MEDIAN)
# [1] 0.0008233387
testU(effects.nbl.t.m[shared.2,]$MEDIAN, effects.nbl.t.m[sclc.specific,]$MEDIAN)
# [1] 0.001626442

ylim <- c(-0.67, 0.4)
file.name <- paste0("DepMap_effects_median1_specific_CLL-CL_RFD_SCLC-specific")
plotBoxDepMapSolid(wd.de.plots, file.name, effects.cll.t.m[shared.3,], effects.cll.t.m[shared.2,], effects.cll.t.m[sclc.specific,], main="DepMap CLL-CL", names=c("SCLC + NB + CLL", "SCLC + (NB/CLL)", "SCLC-specific IZ"), cols=c(red.lighter, red.light, red), ylim, ylab="Gene effect")
testU(effects.cll.t.m[shared.3,]$MEDIAN, effects.cll.t.m[sclc.specific,]$MEDIAN)
# [1] 0.0004397824
testU(effects.cll.t.m[shared.2,]$MEDIAN, effects.cll.t.m[sclc.specific,]$MEDIAN)
# [1] 0.01201821













# -----------------------------------------------------------------------------
# 
# Last Modified: 03/03/23
# -----------------------------------------------------------------------------
sclc.specific <- setdiff(sclc.specific, intersect(sclc.specific, sclc.nl.ctr.iz))
length(sclc.specific)
# [1] 362

sclc.shared <- setdiff(sclc.ctr.iz, sclc.specific)
length(sclc.shared)
# [1] 1428
testU(effects.sclc.t.m[sclc.shared,]$MEDIAN, effects.sclc.t.m[sclc.specific,]$MEDIAN)
# [1] 0.0997264



sclc.iz.shared   <- intersect(rownames(depends.sclc.t.m), sclc.iz.shared.all$external_gene_name)
sclc.iz.specific <- intersect(rownames(depends.sclc.t.m), sclc.iz.specific.all$external_gene_name)
# > testU(effects.sclc.t.m[sclc.iz.shared, ]$MEDIAN, effects.sclc.t.m[sclc.iz.specific, ]$MEDIAN)
# [1] 0.1621926
# testU(depends.sclc.t.m[sclc.iz.shared, ]$MEDIAN, depends.sclc.t.m[sclc.iz.specific, ]$MEDIAN)
# [1] 0.2081016

# -----------------------------------------------------------------------------
# HO vs. CD
# Last Modified: 16/02/23; 11/02/21; 27/11/20; 01/09/20; 29/05/20
# -----------------------------------------------------------------------------
file.name <- paste0("boxplot2_sclc_tpm.gene.median1_rfd_IZ_E+L")
plotBox2(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.l, tpm.gene.log2.m.rfd.ctr.iz.e, main="SCLC IZ genes", names=c("Late", "Early"), cols=c(red, red), ylim)

file.name <- paste0("boxplot2_sclc_tpm.gene.median1_rfd_IZ_E_HO+CD")
plotBox2(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.e.ho, tpm.gene.log2.m.rfd.ctr.iz.e.cd, main="Early IZ genes", names=c("HO", "CD"), cols=c(red, red), ylim)

##
file.name <- paste0("boxplot2_sclc_tpm.gene.median1_rfd_TZ_E+L")
plotBox2(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.l, tpm.gene.log2.m.rfd.ctr.tz.e, main="SCLC TZ genes", names=c("Late", "Early"), cols=c(blue, blue), ylim)

file.name <- paste0("boxplot2_sclc_tpm.gene.median1_rfd_TZ_E_HO+CD")
plotBox2(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.e.ho, tpm.gene.log2.m.rfd.ctr.tz.e.cd, main="Early TZ genes", names=c("HO", "CD"), cols=c(blue, blue), ylim)

# -----------------------------------------------------------------------------
# SCLC-specific genes
# Last Modified: 21/02/23; 19/02/23
# -----------------------------------------------------------------------------
nrds.RT.NRFD <- nrds.RT.NRFD.sclc   ## 17/05/20 WHY nrds.RT.NRFD.sclc.nl? ## MUY MUY IMPORTANTE!!
load(file=file.path(wd.rt.data, paste0(base, "_rpkm.gc.cn.d.rt.log2s.nrfd.", kb, "kb_", "m2-m1", ".RData")))
nrow(nrds.RT.NRFD)
# [1] 2650083

nrds.RT.RFD.1 <- nrds.RT.NRFD.sclc
nrds.RT.RFD.2 <- nrds.RT.NRFD.sclc.nl

nrds.RT.RFD.1.c <- getBootstrapCTR(nrds.RT.RFD.1, rfd)
nrds.RT.RFD.2.c <- getBootstrapCTR(nrds.RT.RFD.2, rfd)
nrds.RT.RFD.1.c.iz <- subset(nrds.RT.RFD.1.c, NRFD > 0)
nrds.RT.RFD.1.c.tz <- subset(nrds.RT.RFD.1.c, NRFD < 0)
nrds.RT.RFD.2.c.iz <- subset(nrds.RT.RFD.2.c, NRFD > 0)
nrds.RT.RFD.2.c.tz <- subset(nrds.RT.RFD.2.c, NRFD < 0)

sclc.iz.shared <- intersect(rownames(nrds.RT.RFD.1.c.iz), rownames(nrds.RT.RFD.2.c.iz))
sclc.iz.specific <- setdiff(rownames(nrds.RT.RFD.1.c.iz), sclc.iz.shared)
sclc.tz.shared <- intersect(rownames(nrds.RT.RFD.1.c.tz), rownames(nrds.RT.RFD.2.c.tz))
sclc.tz.specific <- setdiff(rownames(nrds.RT.RFD.1.c.tz), sclc.tz.shared)

###
## 19/02/23
# > nrow(tpm.gene.log2.m)
# [1] 15803
sclc.iz.shared.median1 <- getTRC(tpm.gene.log2.m, nrds.RT.NRFD[sclc.iz.shared,])
# > nrow(sclc.iz.shared.all)
# [1] 1625
# > nrow(sclc.iz.shared.median1)
# [1] 732
sclc.iz.specific.median1 <- getTRC(tpm.gene.log2.m, nrds.RT.NRFD[sclc.iz.specific,])
# > nrow(sclc.iz.specific.all)
# [1] 1574
# > nrow(sclc.iz.specific.median1)
# [1] 701

###
## 21/02/23
# > nrow(tpm.gene.log2.m)
# [1] 15803
sclc.tz.shared.median1 <- getTRC(tpm.gene.log2.m, nrds.RT.NRFD[sclc.tz.shared,])
# > nrow(sclc.tz.shared.all)
# [1] 
# > nrow(sclc.tz.shared.median1)
# [1] 358
sclc.tz.specific.median1 <- getTRC(tpm.gene.log2.m, nrds.RT.NRFD[sclc.tz.specific,])
# > nrow(sclc.tz.specific.all)
# [1] 
# > nrow(sclc.tz.specific.median1)
# [1] 466

# -----------------------------------------------------------------------------
# 
# Last Modified: 19/02/23
# -----------------------------------------------------------------------------
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median1.RData")))
tpm.gene <- tpm.gene[rownames(sclc.iz.specific.median1), rownames(samples.sclc.tpm)]   ## VERY VERY VERY IMPORTANT!!!
tpm.gene.log2   <- log2(tpm.gene + 1)
tpm.gene.log2.m <- getLog2andMedian(tpm.gene, 1)
dim(tpm.gene.log2.m)
# > dim(tpm.gene.log2.m)
# [1] 701  71

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median1.RData")))
tpm.gene <- tpm.gene[rownames(sclc.tz.specific.median1), rownames(samples.sclc.tpm)]   ## VERY VERY VERY IMPORTANT!!!
tpm.gene.log2   <- log2(tpm.gene + 1)
tpm.gene.log2.m <- getLog2andMedian(tpm.gene, 1)
dim(tpm.gene.log2.m)
# > dim(tpm.gene.log2.m)
# [1] 466  71























#deepmaps <- read.csv("/Users/ty2/Work/sanger/ty2/Pan-Pody/metadata/DepMap/22Q4/CRISPRGeneEffect.csv", header=T)
#deepmaps <- mapply(x = 1:length(deepmaps), function(x) strsplit0(toString(deepmaps[x]), " ")[1])
#essentials <- readTable("/Users/ty2/Work/sanger/ty2/Pan-Pody/metadata/DeepMap/22Q4/CRISPRInferredCommonEssentials.txt", header=F, rownames=F, sep="")[, 1]

effects <- read.csv("/Users/ty2/Work/sanger/ty2/Pan-Pody/metadata/DepMap/22Q2/Achilles_gene_effect.csv", header=T)
rownames(effects) <- effects$DepMap_ID
effects <- effects[,-1]
colnames(effects) <- mapply(x = 1:length(colnames(effects)), function(x) strsplit0(toString(colnames(effects)[x]), "\\.\\.")[1])
# > dim(effects)
# [1]   957 18017
effects.sclc <- effects[intersect(sclcs, rownames(effects)), ]
#effects.sclc <- effects[, intersect(colnames(effects), tpm.gene.log2.m.rfd.ctr.iz$external_gene_name)]

depends <- read.csv("/Users/ty2/Work/sanger/ty2/Pan-Pody/metadata/DepMap/22Q2/Achilles_gene_dependency.csv", header=T)
rownames(depends) <- depends$DepMap_ID
depends <- depends[,-1]
colnames(depends) <- mapply(x = 1:length(colnames(depends)), function(x) strsplit0(toString(colnames(depends)[x]), "\\.\\.")[1])
# > dim(depends)
# [1]   957 18017
depends.sclc <- depends[intersect(sclcs, rownames(depends)), ]

##
median0col <- function(expr) {
	  return(mapply(x = 1:ncol(expr), function(x) median(as.numeric(expr[,x])))) 
}

##
effects.sclc.iz.shared <- effects[, intersect(colnames(effects), sclc.iz.shared$external_gene_name)]
effects.sclc.iz.specific <- effects[, intersect(colnames(effects), sclc.iz.specific$external_gene_name)]

ylab="DepMap gene effect"
ylim <- c(min(c(median0col(effects.sclc.iz.shared), median0col(effects.sclc.iz.specific))), max(c(median0col(effects.sclc.iz.shared), median0col(effects.sclc.iz.specific))))
file.name <- paste0("boxplot2_sclc_tpm.gene.median1_rfd_IZ_Shared-Specific_effects")
plotBox2length(wd.de.plots, file.name, median0col(effects.sclc.iz.shared), median0col(effects.sclc.iz.specific), main="SCLC IZ genes", names=c("Shared", "Specific"), cols=c(red, red), ylim, ylab)

##
depends.sclc.iz.shared <- depends[, intersect(colnames(depends), sclc.iz.shared$external_gene_name)]
depends.sclc.iz.specific <- depends[, intersect(colnames(depends), sclc.iz.specific$external_gene_name)]

ylab="DepMap gene dependency"
ylim <- c(min(c(median0col(depends.sclc.iz.shared), median0col(depends.sclc.iz.specific))), max(c(median0col(depends.sclc.iz.shared), median0col(depends.sclc.iz.specific))))
file.name <- paste0("boxplot2_sclc_tpm.gene.median1_rfd_IZ_Shared-Specific_depends")
plotBox2length(wd.de.plots, file.name, median0col(depends.sclc.iz.shared), median0col(depends.sclc.iz.specific), main="SCLC IZ genes", names=c("Shared", "Specific"), cols=c(red, red), ylim, ylab)










##
test <- toTable(0, 2, 2, c("HO", "CD"))
rownames(test) <- c("ESSENTIAL", "NON")
test[1, 1] <- length(intersect(essentials, tpm.gene.log2.m.rfd.ctr.iz.e.ho$external_gene_name))
test[1, 2] <- length(intersect(essentials, tpm.gene.log2.m.rfd.ctr.iz.e.cd$external_gene_name))
test[2, 1] <- length(intersect(nons,       tpm.gene.log2.m.rfd.ctr.iz.e.ho$external_gene_name))
test[2, 2] <- length(intersect(nons,       tpm.gene.log2.m.rfd.ctr.iz.e.cd$external_gene_name))

#fisher.test(test)[[1]]
chisq.test(test)[[3]]

##
test <- toTable(0, 2, 2, c("HO", "CD"))
rownames(test) <- c("ESSENTIAL", "NON")
test[1, 1] <- length(intersect(essentials, tpm.gene.log2.m.rfd.ctr.tz.e.ho$external_gene_name))
test[1, 2] <- length(intersect(essentials, tpm.gene.log2.m.rfd.ctr.tz.e.cd$external_gene_name))
test[2, 1] <- length(intersect(nons,       tpm.gene.log2.m.rfd.ctr.tz.e.ho$external_gene_name))
test[2, 2] <- length(intersect(nons,       tpm.gene.log2.m.rfd.ctr.tz.e.cd$external_gene_name))

#fisher.test(test)[[1]]
chisq.test(test)[[3]]

##
test <- toTable(0, 2, 2, c("HO", "CD"))
rownames(test) <- c("ESSENTIAL", "NON")
test[1, 1] <- length(intersect(essentials, tpm.gene.log2.m.rfd.ttr.e.ho$external_gene_name))
test[1, 2] <- length(intersect(essentials, tpm.gene.log2.m.rfd.ttr.e.cd$external_gene_name))
test[2, 1] <- length(intersect(nons,       tpm.gene.log2.m.rfd.ttr.e.ho$external_gene_name))
test[2, 2] <- length(intersect(nons,       tpm.gene.log2.m.rfd.ttr.e.cd$external_gene_name))

#fisher.test(test)[[1]]
chisq.test(test)[[3]]



# -----------------------------------------------------------------------------
# RFD vs. TPM (All)
# Last Modified: 16/02/23; 27/11/20; 01/09/20; 29/05/20
# -----------------------------------------------------------------------------
file.name <- file.path(wd.de.data, "tpm_gene_log2_m_rfd+1.RData")
load(file.name)

ylim <- c(min(tpm.gene.log2.m.rfd$MEDIAN), 11)

file.name <- paste0("boxplot3_sclc_tpm.gene+1_RFD_")
plotBox3(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr, tpm.gene.log2.m.rfd.ctr.tz, tpm.gene.log2.m.rfd.ctr.iz, main="SCLC expression", names=c("TTR", "TZ", "IZ"), cols=c("black", blue, red), ylim, text="")

# -----------------------------------------------------------------------------
# 
# Last Modified: 16/02/23; 29/05/20
# -----------------------------------------------------------------------------
genes <- c("DHX36", "DNAJC2", "TERT")
genes <- c("RAD9A", "PIF1", "E2F3", "BLM", "KIF18B", "GTPBP3", "BRCA2")
for (g in 1:length(genes)) {
	id <- subset(ensGene, external_gene_name == genes[g])$ensembl_gene_id
	plotSRC(genes[g], as.numeric(tpm.gene.log2[id,]), samples.sclc.tpm$COR, 1, red, "bottomright", xlab.text="")
}





# -----------------------------------------------------------------------------
# Purities
# Last Modified: 10/07/22
# -----------------------------------------------------------------------------
samples.wgs$purity <- NA
samples.wgs$ploidy <- NA
samples.wgs$purity2 <- NA
samples.wgs$ploidy2 <- NA

for (s in 1:nrow(samples.wgs)) {
   sample <- samples.wgs$SAMPLE_ID[s]
   purities <- readTable(file.path(wd.wgs, "peiflyne", sample, paste0(sample, "_ANALYSIS"), paste0(sample, "_pupl.txt")), header=F, rownames=F, sep="\t")
   purities2 <- readTable(file.path(wd.wgs, "peiflyne/2015", sample, paste0(sample, "_ANALYSIS"), paste0(sample, "_pupl.txt")), header=F, rownames=F, sep="\t")
   
   samples.wgs$purity[s]  <- purities$V2
   samples.wgs$ploidy[s]  <- purities$V3
   samples.wgs$purity2[s] <- purities2$V2
   samples.wgs$ploidy2[s] <- purities2$V3
}

file.name <- file.path(wd.rt.plots, paste0("Correlation_SCLC_purity"))
x <- samples.wgs$purity
y <- samples.wgs$purity2
plotCorrelation(file.name, "SCLC purity", "BWA-MEM", "BWA", x, y, pos="bottomright", cols=c("dimgray", "black"), size=5)

file.name <- file.path(wd.rt.plots, paste0("Correlation_SCLC_ploidy"))
x <- samples.wgs$ploidy
y <- samples.wgs$ploidy2
plotCorrelation(file.name, "SCLC ploidy", "BWA-MEM", "BWA", x, y, pos="bottomright", cols=c("dimgray", "black"), size=5)

save(samples.wgs, samples.sclc.tpm, file=file.path(wd.rt.data, "SCLC_purities_n101.RData"))

# -----------------------------------------------------------------------------
# Correlation bwteen TPM and in-silico sorting
# Last Modified: 12/12/19; 08/01/19; 17/08/17
# -----------------------------------------------------------------------------
samples.tpm.sclc$Purity <- samples.tpm.sclc$purity2

xlab.text <- "Expression vs. SCF index [rho]"
ylab.text <- "Expression vs. Purity [rho]"
pvalue <- 0.01

colnames <- c("GENE", "ADJ_1", "ADJ_2")
genes0 <- c("PIF1", "E2F3", "RAD9A", "BLM", "RFC5", "RAD54L")
genes <- toTable(NA, length(colnames), length(genes0), colnames)
genes$GENE <- genes0
#genes[2, 2] <- -0.15
#genes[2, 3] <- 1.2

for (q in 1:4) {
   samples.sclc.tpm.Q4 <- subset(samples.sclc.tpm, Q4 == q)
   n <- nrow(samples.sclc.tpm.Q4)
 
   load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
   tpm.gene <- tpm.gene[, rownames(samples.sclc.tpm.Q4)]   ## VERY VERY VERY IMPORTANT!!!
   tpm.gene.log2   <- log2(tpm.gene + 1)
 
   #getSRC(wd.de.data, BASE, tpm.gene.log2, samples.sclc.tpm.Q4, "COR", n, "Q", q)
   #getSRC(wd.de.data, BASE, tpm.gene.log2, samples.sclc.tpm.Q4, "Purity", n, "Q", q)
 
   ##
   expressed <- rownames(tpm.gene)
   de <- getCannoli(wd.de.data, BASE, n, expressed, TEST="COR", TEST2="Purity", M2=paste0("_Q", q))
   plot.de <- file.path(wd.de.plots, paste0("cannoliplot_SRC_", BASE, "_TPM-COR-Purity_P1E02_MEDIAN0_Q", q))
   #genes <- getCannoliGenes(de, pvalue, genes0)
   #genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
   file.de <- paste0(plot.de, ".pdf")
   file.main <- c(paste0("Lung-SCLC Q", q, " (n=", n, ")"), "")
   plotCannoli(de, pvalue, genes, file.de, file.main, xlab.text, ylab.text, "bottomleft", c("", ""), c(red.lighter, blue.lighter), c(red, blue), fold=0, pos="bottomright")
}

for (m in 1:2) {
   samples.tpm.sclc.M2 <- subset(samples.tpm.sclc, M2 == m)
   n <- nrow(samples.tpm.sclc.M2)
 
   load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
   tpm.gene <- tpm.gene[, rownames(samples.tpm.sclc.M2)]   ## VERY VERY VERY IMPORTANT!!!
   tpm.gene.log2   <- log2(tpm.gene + 1)
 
   #getSRC(wd.de.data, BASE, tpm.gene.log2, samples.tpm.sclc.M2, "COR", n, "M", m)
   #getSRC(wd.de.data, BASE, tpm.gene.log2, samples.tpm.sclc.M2, "Purity", n, "M", m)
 
   ##
   expressed <- rownames(tpm.gene)
   de <- getCannoli(wd.de.data, BASE, n, expressed, TEST="COR", TEST2="Purity", M2=paste0("_M", m))
   plot.de <- file.path(wd.de.plots, paste0("cannoliplot_SRC_", BASE, "_TPM-COR-Purity_P1E02_MEDIAN0_M", m))
   #genes <- getCannoliGenes(de, pvalue, genes0)
   #genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
   file.de <- paste0(plot.de, ".pdf")
   file.main <- c(paste0("Lung-SCLC M", m, " (n=", n, ")"), "")
   plotCannoli(de, pvalue, genes, file.de, file.main, xlab.text, ylab.text, "bottomleft", c("", ""), c(red.lighter, blue.lighter), c(red, blue), fold=0, pos="bottomright")
}

n <- nrow(samples.tpm.sclc)

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
tpm.gene <- tpm.gene[, rownames(samples.tpm.sclc)]   ## VERY VERY VERY IMPORTANT!!!
tpm.gene.log2   <- log2(tpm.gene + 1)

#getSRC0(wd.de.data, BASE, tpm.gene.log2, samples.tpm.sclc, "COR", n)
#getSRC0(wd.de.data, BASE, tpm.gene.log2, samples.tpm.sclc, "Purity", n)

##
expressed <- rownames(tpm.gene)
de <- getCannoli(wd.de.data, BASE, n, expressed, TEST="COR", TEST2="Purity", M2="")
plot.de <- file.path(wd.de.plots, paste0("cannoliplot_SRC_", BASE, "_TPM-COR-Purity_P1E02_MEDIAN0"))
#genes <- getCannoliGenes(de, pvalue, genes0)
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c(paste0("Lung-SCLC (n=", n, ")"), "")
plotCannoli(de, pvalue, genes, file.de, file.main, xlab.text, ylab.text, "bottomleft", c("", ""), c(red.lighter, blue.lighter), c(red, blue), fold=0, pos="bottomright")

# -----------------------------------------------------------------------------
# 
# Last Modified: 19/10/22
# -----------------------------------------------------------------------------
genes.G1S <- readTable(file.path(wd.meta, "Dominguez_G1-S.list"), header=F, rownames=F, sep="")
genes.G2M <- readTable(file.path(wd.meta, "Dominguez_G2-M.list"), header=F, rownames=F, sep="")

for (q in 1:4) {
   n <- nrow(subset(samples.tpm.sclc, Q4 == q))
   expressed <- rownames(tpm.gene)
   de <- getCannoli(wd.de.data, BASE, n, expressed, TEST="COR", TEST2="Purity", M2=paste0("_Q", q))
 
   plot.de <- file.path(wd.de.plots, paste0("cannoliplot_SRC_", BASE, "_TPM-COR-Purity_MEDIAN0_G1-S_Q", q))
   file.de <- paste0(plot.de, ".pdf")
   file.main <- c(paste0("Lung-SCLC Q", q, " (n=", n, ")"), "")
   plotCannoliGenes(de, genes.G1S, file.de, file.main, xlab.text, ylab.text, col=red, pos="bottomright")
 
   plot.de <- file.path(wd.de.plots, paste0("cannoliplot_SRC_", BASE, "_TPM-COR-Purity_MEDIAN0_G2-M_Q", q))
   file.de <- paste0(plot.de, ".pdf")
   file.main <- c(paste0("Lung-SCLC Q", q, " (n=", n, ")"), "")
   plotCannoliGenes(de, genes.G2M, file.de, file.main, xlab.text, ylab.text, col=green, pos="bottomright")
}

##
file.name <- file.path(wd.de.plots, paste0("Correlation_SCLC_purity_SCF"))
x <- samples.tpm.sclc$COR
y <- samples.tpm.sclc$purity2
plotCorrelation(file.name, "SCLC", "SCF index", "Purity", x, y, pos="bottomright", cols=c("dimgray", "black"), size=5)

for (q in 1:4) {
   samples.tpm.sclc.Q4 <- subset(samples.tpm.sclc, Q4 == q)
 
   file.name <- file.path(wd.de.plots, paste0("Correlation_SCLC_purity_SCF_Q", q))
   x <- samples.tpm.sclc.Q4$COR
   y <- samples.tpm.sclc.Q4$purity2
   plotCorrelation(file.name, paste0("SCLC (Q", q, ")"), "SCF index", "Purity", x, y, pos="bottomright", cols=c("dimgray", "black"), size=5)
}

##
file.name <- file.path(wd.de.plots, paste0("Correlation_SCLC_ploidy_SCF"))
x <- samples.tpm.sclc$COR
y <- samples.tpm.sclc$ploidy2
plotCorrelation(file.name, "SCLC", "SCF index", "Ploidy", x, y, pos="bottomright", cols=c("dimgray", "black"), size=5)

for (q in 1:4) {
   samples.tpm.sclc.Q4 <- subset(samples.tpm.sclc, Q4 == q)
 
   file.name <- file.path(wd.de.plots, paste0("Correlation_SCLC_ploidy_SCF_Q", q))
   x <- samples.tpm.sclc.Q4$COR
   y <- samples.tpm.sclc.Q4$ploidy2
   plotCorrelation(file.name, paste0("SCLC (Q", q, ")"), "SCF index", "Ploidy", x, y, pos="bottomright", cols=c("dimgray", "black"), size=5)
}

# -----------------------------------------------------------------------------
# PCA
# Last Modified: 26/03/21; 23/08/20
# -----------------------------------------------------------------------------
dim(tpm.gene)
# [1] 23572    70

overlaps <- intersect(rownames(samples.sclc.tpm), rownames(samples.wgs))
test <- tpm.gene[, rownames(samples.sclc.tpm)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))
save(pca.de, file=file.path(wd.de.data, paste0("PCA_SCLC_n70.RData")))

samples.sclc.tpm$GROUP_ID <- "Proliferative"
samples.sclc.tpm[rownames(subset(samples.sclc.tpm, M2 == 0)),]$GROUP_ID <- "Resting"

trait <- samples.sclc.tpm$GROUP_ID
file.main <- c("SCLC expression", "")
plotPCA(1, 2, pca.de, trait, wd.de.plots, "PCA_NBL_S-G1_n70", size=5, file.main, "topleft", c("Proliferative", "Resting"), c(red, blue), flip.x=1, flip.y=1)

scores.iCN <- pcaScores(pca.de)
x <- scores.iCN[overlaps,]$PC1
y <- wgds[overlaps,]$decision_value
file.name <- file.path(wd.de.plots, paste0("correlation_PCA-SCLC-TPM-MEDIAN0_PC1-vs-WGD"))
plotCorrelation(file.name, "SCLC expression", paste0("PC", 1, " (", pcaProportionofVariance(pca.iCN, 1), "%)"), "WGD", x, y, size=5)

x <- scores.iCN[overlaps,]$PC1
y <- samples.sclc[overlaps,]$COR
file.name <- file.path(wd.de.plots, paste0("correlation_PCA-SCLC-TPM-MEDIAN0_PC1-vs-SORTING"))
plotCorrelation(file.name, "SCLC expression", paste0("PC", 1, " (", pcaProportionofVariance(pca.iCN, 1), "%)"), "Proliferation rate", x, y, size=5)

x <- scores.iCN[overlaps,]$PC1
y <- samples.wgs[overlaps,]$purity
file.name <- file.path(wd.de.plots, paste0("correlation_PCA-SCLC-TPM-MEDIAN0_PC1-vs-purity"))
plotCorrelation(file.name, "SCLC expression", paste0("PC", 1, " (", pcaProportionofVariance(pca.iCN, 1), "%)"), "Purity", x, y, size=5)

x <- scores.iCN[overlaps,]$PC1
y <- samples.wgs[overlaps,]$ploidy
file.name <- file.path(wd.de.plots, paste0("correlation_PCA-SCLC-TPM-MEDIAN0_PC1-vs-ploidy"))
plotCorrelation(file.name, "SCLC expression", paste0("PC", 1, " (", pcaProportionofVariance(pca.iCN, 1), "%)"), "Ploidy", x, y, size=5)

# -----------------------------------------------------------------------------
# Correlation bwteen TPM and in-silico sorting
# Last Modified: 12/12/19; 08/01/19; 17/08/17
# -----------------------------------------------------------------------------
colnames <- c("RHO", "P", "Q", "M1", "M2", "LOG2_FC")   ##"ANOVA_P", "ANOVA_Q", 
src <- toTable(0, length(colnames), nrow(tpm.gene.log2), colnames)
rownames(src) <- rownames(tpm.gene.log2)

## SRC
src$RHO <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples.sclc.tpm$purity2, method="spearman", exact=F)[[4]])
src$P   <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples.sclc.tpm$purity2, method="spearman", exact=F)[[3]])
src <- src[!is.na(src$P),]

## Log2 fold change
#de$M1 <- median00(tpm.gene.log2, rownames(subset(samples, M2 == 0)))
#de$M2 <- median00(tpm.gene.log2, rownames(subset(samples, M2 == 1)))
#de$LOG2_FC <- de$M2 - de$M1

## FDR
library(qvalue)
src$Q <- qvalue(src$P)$qvalue
src <- src[order(src$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
src.tpm.gene <- cbind(annot[rownames(src),], src)   ## BE EXTRA CAREFUL!!

writeTable(src.tpm.gene, file.path(wd.de.data, "SRC_SCLC_tpm-gene_Purity-vs-TPM_q_n70.txt"), colnames=T, rownames=F, sep="\t")
save(src.tpm.gene, samples.sclc.tpm, file=file.path(wd.de.data, "SRC_SCLC_tpm-gene_Purity-vs-TPM_q_n70.RData"))
nrow(src.tpm.gene)
# [1] 23572
# [1] 19131

###
##
genes <- c("NSUN2", "BRD9", "HIC2", "PIF1", "BLM", "E2F3", "TONSL", "SULT1A3")
genes <- c("RELA", "REL", "TONSL")
for (g in 1:length(genes)) {
   id <- subset(ensGene, external_gene_name == genes[g])$ensembl_gene_id
   plotSRC(genes[g], as.numeric(tpm.gene.log2[id,]), samples.sclc.tpm$COR, 1, pos="bottomright", col=red, size=5)
}

genes <- c("CMPK2", "RSAD2", "CCDC85A", "GPR37")
genes <- c("NLRP3", "NFKB1", "NFKB2", "RELB")
genes <- c("SLC25A46", "CDC37L1")
genes <- c("NDUFA2", "PDHB", "ISCA1", "NDUFS4", "NDUFA8")
genes <- c("CDC37L1")
for (g in 1:length(genes)) {
   id <- subset(ensGene, external_gene_name == genes[g])$ensembl_gene_id
   plotSRC(file.path(wd.de.plots, ""), genes[g], as.numeric(tpm.gene.log2[id,]), samples.sclc.tpm$COR, 1, pos="bottomright", col=blue, size=5)
}

# -----------------------------------------------------------------------------
# Correlation bwteen TPM and CNAs
# Last Modified: 26/06/22; 12/12/19; 08/01/19; 17/08/17
# -----------------------------------------------------------------------------
dim(cna.gene.nona.tpm)
# [1] 34895   101
cna.gene.nona.tpm.src <- cna.gene.nona.tpm[, samples.sclc.tpm$SAMPLE_ID]
tpm.gene.log2.cna <- tpm.gene.log2[rownames(cna.gene.nona.tpm.src), rownames(samples.sclc.tpm)]
#colnames(tpm.gene.log2.cna) <- samples.sclc.tpm$SAMPLE_ID
dim(cna.gene.nona.tpm.src)
# [1] 34895    70
dim(tpm.gene.log2.cna)
# [1] 34895    70

colnames <- c("RHO", "P", "Q", "DEL", "AMP", "LOG2_FC")   ##"ANOVA_P", "ANOVA_Q", 
src <- toTable(0, length(colnames), nrow(tpm.gene.log2.cna), colnames)
rownames(src) <- rownames(tpm.gene.log2.cna)

## SRC
src$RHO <- mapply(x = 1:nrow(tpm.gene.log2.cna), function(x) cor.test(as.numeric(tpm.gene.log2.cna[x,]), as.numeric(cna.gene.nona.tpm.src[x,]), method="spearman", exact=F)[[4]])
src$P   <- mapply(x = 1:nrow(tpm.gene.log2.cna), function(x) cor.test(as.numeric(tpm.gene.log2.cna[x,]), as.numeric(cna.gene.nona.tpm.src[x,]), method="spearman", exact=F)[[3]])

## Log2 fold change
#src$DEL <- mapply(x = 1:nrow(tpm.gene.log2.cna), function(x) median(as.numeric(cna.gene.nona.tpm.src[x, colnames(cna.gene.nona.tpm.src)[which(cna.gene.nona.tpm.src[x,] < 2)]])))
#src$AMP <- mapply(x = 1:nrow(tpm.gene.log2.cna), function(x) median(as.numeric(cna.gene.nona.tpm.src[x, colnames(cna.gene.nona.tpm.src)[which(cna.gene.nona.tpm.src[x,] > 2)]])))
#src$LOG2_FC <- src$AMP - src$DEL
src <- src[!is.na(src$P),]

## FDR
library(qvalue)
src$Q <- qvalue(src$P)$qvalue
src <- src[order(src$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
src.tpm.gene <- cbind(annot[rownames(src),], src)   ## BE EXTRA CAREFUL!!

writeTable(src.tpm.gene, file.path(wd.de.data, "SRC_SCLC_tpm-gene_CNA-vs-TPM_q_n70.txt"), colnames=T, rownames=F, sep="\t")
save(src.tpm.gene, samples.sclc.tpm, file=file.path(wd.de.data, "SRC_SCLC_tpm-gene_CNA-vs-TPM_q_n70.RData"))
#writeRNKformat(src.tpm.gene, wd.de.gsea, "SRC_SCLC_tpm-gene-median0_CNA-vs-TPM_q_n70")   ## GSEA

###
## Volcano plots
#xlab.text <- "CNA difference"
#ylab.text <- expression("Significance " * "[-log" * ""[10] * "(" * italic("P") * ")]")

#plot.de <- file.path(wd.de.plots, "volcanoplot_SRC_NBL_r5p47_SORTING_TERT")
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
#file.de <- paste0(plot.de, ".pdf")
#file.main <- c("S (n=25) vs. G1 (n=29)", "")
#plotVolcano(src.tpm.gene, 1, genes, file.de, file.main, xlab.text, ylab.text, "topright", c("S > G1", "S < G1"), c(adjustcolor.red, adjustcolor.blue), c(red, blue), fold=0)

###
##
genes <- c("NSUN2", "BRD9", "HIC2", "PIF1", "BLM", "E2F3", "TONSL", "SULT1A3")
genes <- c("NLRP3", "TONSL", "NFKB1", "RELB", "RELA", "REL")
for (g in 1:length(genes)) {
   id <- subset(ensGene, external_gene_name == genes[g])$ensembl_gene_id
   plotCNA(genes[g], as.numeric(tpm.gene.log2.cna[id,]), as.numeric(cna.gene.nona.tpm.src[id,]), 1, pos="bottomright", col="black", size=5)
}

genes <- c("CMPK2", "RSAD2", "CCDC85A", "GPR37")
genes <- c("NFKB2")
genes <- c("SLC25A46", "CDC37L1")
genes <- c("CDC37L1")
genes <- c("YWHAE")
for (g in 1:length(genes)) {
   id <- subset(ensGene, external_gene_name == genes[g])$ensembl_gene_id
   plotCNA(genes[g], as.numeric(tpm.gene.log2.cna[id,]), as.numeric(cna.gene.nona.tpm.src[id,]), 1, pos="bottomright", col="black", size=5)
}

# -----------------------------------------------------------------------------
# Cannoli plot (P < 0.001)
# Last Modified: 26/06/22; 11/10/17
# -----------------------------------------------------------------------------
xlab.text <- "Expression vs. SCF index [rho]"
#ylab.text <- expression("Expression vs."~italic('in silico')~"sorting")
ylab.text <- "Expression vs. Purity [rho]"
pvalue <- 0.001

colnames <- c("GENE", "ADJ_1", "ADJ_2")
genes0 <- c("CMPK2", "BLM", "E2F3", "BRD9")
genes <- toTable(NA, length(colnames), length(genes0), colnames)
genes$GENE <- genes0
genes[2, 2] <- -0.15
#genes[2, 3] <- 1.2

## Total
de <- getCannoli(wd.de.data, BASE, 70, NULL, TEST="SORTING", TEST2="Purity")
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
expressed <- rownames(tpm.gene)

de <- getCannoli(wd.de.data, BASE, 70, expressed, TEST="SORTING", TEST2="Purity")
plot.de <- file.path(wd.de.plots, "cannoliplot_SRC_SCLC_TPM-Purity-SORTING_P1E03_MEDIAN0")
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c(paste0("Lung-SCLC (n=", 70, ")"), "")
plotCannoli(de, pvalue, genes, file.de, file.main, xlab.text, ylab.text, "bottomleft", c("", ""), c(red.lighter, blue.lighter), c(red, blue), fold=0, pos="bottomright")





10290/23559
# [1] 0.4367758
8267/23559
# [1] 0.3509062
3600/23559
# [1] 0.1528078
1402/23559
# [1] 0.05951017

#plot.de <- file.path(wd.de.plots, "cannoliplot_SRC_SCLC_TPM-CNA-SORTING_P1E03_E2F3")
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
#file.de <- paste0(plot.de, ".pdf")
#file.main <- c(paste0(BASE, " total genes"), "")
#plotCannoli(de, pvalue, genes, file.de, file.main, xlab.text, ylab.text, "topleft", c("", ""), c(red, blue), c(red, blue), fold=0)

## Expressed
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
expressed <- rownames(tpm.gene)
#load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.RData")))

de <- getCannoli(wd.de.data, BASE, 70, expressed, TEST="SORTING", TEST2="Purity")
plot.de <- file.path(wd.de.plots, "cannoliplot_SRC_SCLC_TPM-Purity-SORTING_P1E03_MEDIAN0_GSEA")
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c(paste0("Lung-SCLC (n=", 70, ")"), "")
plotCannoli(de, pvalue, genes, file.de, file.main, xlab.text, ylab.text, "bottomright", c("", ""), c(red.lighter, blue.lighter), c(red, blue), fold=0, pos="bottomright")




plot.de <- file.path(wd.de.plots, "cannoliplot_SRC_SCLC_TPM-CNA-SORTING_P1E03_MEDIAN0_POS")
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c(paste0("Lung-SCLC (n=", 70, ")"), "")
plotCannoli(de, 1, genes, file.de, file.main, xlab.text, ylab.text, "topleft", c("", ""), c(red.lighter, "lightgray"), c(red.lighter, "lightgray"), fold=0)

plot.de <- file.path(wd.de.plots, "cannoliplot_SRC_SCLC_TPM-CNA-SORTING_P1E03_MEDIAN0_NEG")
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c(paste0("Lung-SCLC (n=", 70, ")"), "")
plotCannoli(de, 1, genes, file.de, file.main, xlab.text, ylab.text, "topleft", c("", ""), c("lightgray", blue.lighter), c("lightgray", blue.lighter), fold=0)

## GSEA
de.pos.pos <- subset(subset(de, Effect1 > 0), Effect2 > 0)
de.pos.pos.sig <- subset(subset(de.pos.pos, P1 <= 0.001), P2 <= 0.001)

de.pos.neg <- subset(subset(de, Effect1 > 0), Effect2 < 0)
de.pos.neg.sig <- subset(subset(de.pos.neg, P1 <= 0.001), P2 <= 0.001)

de.neg.neg <- subset(subset(de, Effect1 < 0), Effect2 < 0)
de.neg.neg.sig <- subset(subset(de.neg.neg, P1 <= 0.001), P2 <= 0.001)

de.neg.pos <- subset(subset(de, Effect1 < 0), Effect2 > 0)
de.neg.pos.sig <- subset(subset(de.neg.pos, P1 <= 0.001), P2 <= 0.001)

writeTable(de.positive.sig, file.path(wd.de.data, "SRC_SCLC_tpm-gene_SORTING-Purity-TPM_q_n70_pos16.txt"), colnames=T, rownames=T, sep="\t")

writeRNKformatCNA(rbind(de.pos.pos, de.neg.neg), wd.de.gsea, "SRC_SCLC_tpm-gene-median0_SORTING-Purity-TPM_q_n70_POS-SELECTION")   ## GSEA
writeRNKformatCNA(rbind(de.neg.pos, de.pos.neg), wd.de.gsea, "SRC_SCLC_tpm-gene-median0_SORTING-Purity-TPM_q_n70_NEG-SELECTION")   ## GSEA

writeRNKformatCNA(rbind(de.pos.pos, de.pos.neg), wd.de.gsea, "SRC_SCLC_tpm-gene-median0_SORTING-Purity-TPM_q_n70_GAIN")   ## GSEA
writeRNKformatCNA(rbind(de.neg.pos, de.neg.neg), wd.de.gsea, "SRC_SCLC_tpm-gene-median0_SORTING-Purity-TPM_q_n70_LOSS")   ## GSEA

# -----------------------------------------------------------------------------
# GSEA (11/07/22)
# -----------------------------------------------------------------------------
## POS
gsea.pos <- readTable(file.path(wd.de.gsea, "HALLMARK", "gsea_report_for_na_pos_1657521569538.tsv"), header=T, rownames=F, sep="\t")
gsea.pos$NAME <- mapply(x = 1:nrow(gsea.pos), function(x) unlist(strsplit(gsea.pos$NAME[x], "HALLMARK_"))[2])
gsea.pos$NAME <- mapply(x = 1:nrow(gsea.pos), function(x) gsub("_", " ", gsea.pos$NAME[x]))
gsea.pos <- gsea.pos[5:1, c("NAME", "NES")]
gsea.pos[5, 1] <- "E2F targets"
gsea.pos[4, 1] <- "G2M checkpoint"
gsea.pos[3, 1] <- "MYC targets v1"
gsea.pos[2, 1] <- firstup(gsea.pos[2, 1])
gsea.pos[1, 1] <- "MYC targets v2"
 
main.text <- c("Hallmark enrichment score                         ", "")
file.de <- file.path(wd.de.gsea, "HALLMARK_POS.pdf")
pdf(file.de, height=2, width=4)
par(mar=c(2,8.5,2,4))   # increase y-axis margin.
barplot(gsea.pos$NES, main=main.text[1], las=1, horiz=T, xlim=c(0, 6), xaxt="n", names.arg=gsea.pos$NAME, col=red.lighter, xlab="")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
#abline(v=2, lty=5)
axis(side=1, at=seq(0, 6, by=2), labels=F)
axis(side=1, at=seq(0, 6, by=2), pos=0.5, tick=F)
#mtext(main.text[2], line=0.3)
dev.off()

## NEG
gsea.neg <- readTable(file.path(wd.de.gsea, "HALLMARK", "gsea_report_for_na_neg_1657521569538.tsv"), header=T, rownames=F, sep="\t")
gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) unlist(strsplit(gsea.neg$NAME[x], "HALLMARK_"))[2])
gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) gsub("_", " ", gsea.neg$NAME[x]))
gsea.neg <- gsea.neg[5:1, c("NAME", "NES")]
gsea.neg[5, 1] <- "Epithelial mesenchymal trans"
gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) firstup(gsea.neg$NAME[x]))

main.text <- c("", "")
file.de <- file.path(wd.de.gsea, "HALLMARK_NEG.pdf")
pdf(file.de, height=1.6, width=4)
par(mar=c(2,1,0,11.5))   # increase y-axis margin.
posbar <- barplot(gsea.neg$NES, main=main.text[1], las=1, horiz=T, xlim=c(-6, 0), xaxt="n", names.arg="", col=blue.lighter, xlab="")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=gsea.neg$NAME, xpd=NA)
#abline(v=-2, lty=5)
axis(side=1, at=seq(-6, 0, by=2), labels=F)
axis(side=1, at=seq(-6, 0, by=2), pos=0.5, tick=F)
#mtext(main.text[2], line=0.3)
dev.off()

# -----------------------------------------------------------------------------
# CNA GAIN or LOSS
# -----------------------------------------------------------------------------
overlaps <- intersect(rownames(cna.gene.nona), expressed)
cna.gene.nona.median0 <- cna.gene.nona[overlaps, ]
dim(cna.gene.nona.median0)
# [1] 23559   101

cna.median0.sclc <- toTable(0, 4, nrow(cna.gene.nona.median0), c("N_LOSS", "N_2N", "N_GAIN", "MEDIAN"))
rownames(cna.median0.sclc) <- overlaps
cna.median0.sclc$N_LOSS <- mapply(x = 1:nrow(cna.gene.nona.median0), function(x) length(which(cna.gene.nona.median0[x,] < 2)))
cna.median0.sclc$N_2N   <- mapply(x = 1:nrow(cna.gene.nona.median0), function(x) length(which(cna.gene.nona.median0[x,] == 2)))
cna.median0.sclc$N_GAIN <- mapply(x = 1:nrow(cna.gene.nona.median0), function(x) length(which(cna.gene.nona.median0[x,] > 2)))
cna.median0.sclc$MEDIAN <- mapply(x = 1:nrow(cna.gene.nona.median0), function(x) median(as.numeric(cna.gene.nona.median0[x,])))

nrow(subset(cna.median0.sclc, MEDIAN < 2))
# [1] 3360
nrow(subset(cna.median0.sclc, MEDIAN == 2))
# [1] 13007
nrow(subset(cna.median0.sclc, MEDIAN > 2))
# [1] 7192

length(which(cna.median0.sclc$N_LOSS > cna.median0.sclc$N_GAIN))
length(which(cna.median0.sclc$N_LOSS == cna.median0.sclc$N_GAIN))
length(which(cna.median0.sclc$N_LOSS < cna.median0.sclc$N_GAIN))

loss <- unique(c(rownames(subset(cna.median0.sclc, MEDIAN < 2)), rownames(cna.median0.sclc[which(cna.median0.sclc$N_LOSS > cna.median0.sclc$N_GAIN),])))
gain <- intersect(rownames(subset(cna.median0.sclc, MEDIAN > 2)), rownames(cna.median0.sclc[which(cna.median0.sclc$N_LOSS < cna.median0.sclc$N_GAIN),]))
normal <- setdiff(overlaps, c(loss, gain))             
save(expressed, overlaps, loss, normal, gain, file=file.path(wd.de.data, "SCLC-MEDIAN0-CNA_n101.RData"))
               
# -----------------------------------------------------------------------------
# Cannoli plot: CNA NORMAL
# Last Modified: 26/06/22; 11/10/17
# -----------------------------------------------------------------------------
xlab.text <- "Expression vs. proliferation [rho]"
#ylab.text <- expression("Expression vs."~italic('in silico')~"sorting")
ylab.text <- "Expression vs. CNA [rho]"
pvalue <- 0.001

## CN = 2
colnames <- c("GENE", "ADJ_1", "ADJ_2")
genes0 <- c("PIF1", "BLM", "NFKB2", "RELA")
genes <- toTable(NA, length(colnames), length(genes0), colnames)
genes$GENE <- genes0
genes[2, 2] <- -0.15
#genes[3, 2] <- -0.1
#genes[3, 3] <- 1

de <- getCannoli(wd.de.data, BASE, 70, normal)
plot.de <- file.path(wd.de.plots, "cannoliplot_SRC_SCLC_TPM-CNA-SORTING_P1E03_MEDIAN0_CNA_NORMAL_white_NFKB2")
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c(paste0("Lung-SCLC (CN = 2)"), "")
plotCannoli(de, pvalue, genes, file.de, file.main, xlab.text, ylab.text, "bottomleft", c("", ""), c(red.lighter, blue.lighter), c(red, blue), fold=0, pos="bottomright")

## GSEA
de.pos.pos <- subset(subset(de, Effect2 > 0), Effect1 > 0)
de.pos.neg <- subset(subset(de, Effect2 > 0), Effect1 < 0)
writeRNKformatCNA(rbind(de.pos.pos, de.pos.neg), wd.de.gsea, "SRC_SCLC_tpm-gene-median0_SORTING-CNA-TPM_q_n70_CNA-NORMAL-POS")   ## GSEA

## POS
gsea.pos <- readTable(file.path(wd.de.gsea, "CNA_NORMAL_POS", "HALLMARK", "gsea_report_for_na_pos_1657718053844.tsv"), header=T, rownames=F, sep="\t")
gsea.pos$NAME <- mapply(x = 1:nrow(gsea.pos), function(x) unlist(strsplit(gsea.pos$NAME[x], "HALLMARK_"))[2])
gsea.pos$NAME <- mapply(x = 1:nrow(gsea.pos), function(x) gsub("_", " ", gsea.pos$NAME[x]))
gsea.pos <- gsea.pos[5:1, c("NAME", "NES")]
gsea.pos[5, 1] <- "E2F targets"
gsea.pos[4, 1] <- "G2M checkpoint"
gsea.pos[3, 1] <- "MYC targets v1"
gsea.pos[2, 1] <- firstup(gsea.pos[2, 1])
gsea.pos[1, 1] <- "MYC targets v2"

main.text <- c("Hallmark enrichment score                         ", "")
file.de <- file.path(wd.de.gsea, "CNA_NORMAL_POS", "HALLMARK_POS.pdf")
pdf(file.de, height=2, width=4)
par(mar=c(2,8.5,2,4))   # increase y-axis margin.
barplot(gsea.pos$NES, main=main.text[1], las=1, horiz=T, xlim=c(0, 4), xaxt="n", names.arg=gsea.pos$NAME, col=red.lighter, xlab="")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
#abline(v=2, lty=5)
axis(side=1, at=seq(0, 4, by=2), labels=F)
axis(side=1, at=seq(0, 4, by=2), pos=0.5, tick=F)
#mtext(main.text[2], line=0.3)
dev.off()

## NEG
gsea.neg <- readTable(file.path(wd.de.gsea, "CNA_NORMAL_POS", "HALLMARK", "gsea_report_for_na_neg_1657718053844.tsv"), header=T, rownames=F, sep="\t")
gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) unlist(strsplit(gsea.neg$NAME[x], "HALLMARK_"))[2])
gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) gsub("_", " ", gsea.neg$NAME[x]))
gsea.neg <- gsea.neg[5:1, c("NAME", "NES")]
gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) firstup(gsea.neg$NAME[x]))

main.text <- c("", "")
file.de <- file.path(wd.de.gsea, "CNA_NORMAL_POS", "HALLMARK_NEG.pdf")
pdf(file.de, height=1.6, width=4)
par(mar=c(2,1,0,11.5))   # increase y-axis margin.
posbar <- barplot(gsea.neg$NES, main=main.text[1], las=1, horiz=T, xlim=c(-4, 0), xaxt="n", names.arg="", col=blue.lighter, xlab="")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=gsea.neg$NAME, xpd=NA)
#abline(v=-2, lty=5)
axis(side=1, at=seq(-4, 0, by=2), labels=F)
axis(side=1, at=seq(-4, 0, by=2), pos=0.5, tick=F)
#mtext(main.text[2], line=0.3)
dev.off()

# -----------------------------------------------------------------------------
# Cannoli plot: CNA GAIN
# Last Modified: 26/06/22; 11/10/17
# -----------------------------------------------------------------------------
colnames <- c("GENE", "ADJ_1", "ADJ_2")
genes0 <- c("TONSL", "E2F3", "BRD9", "RELB")
genes <- toTable(NA, length(colnames), length(genes0), colnames)
genes$GENE <- genes0
genes[1, 2] <- 0.13
genes[1, 3] <- 1.5

de <- getCannoli(wd.de.data, BASE, 70, gain)
plot.de <- file.path(wd.de.plots, "cannoliplot_SRC_SCLC_TPM-CNA-SORTING_P1E03_MEDIAN0_CNA_GAIN_white_RELB")
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c(paste0("Lung-SCLC (CN > 2)"), "")
plotCannoli(de, pvalue, genes, file.de, file.main, xlab.text, ylab.text, "bottomleft", c("", ""), c(red.lighter, blue.lighter), c(red, blue), fold=0, pos="bottomright")

## GSEA
de.pos.pos <- subset(subset(de, Effect1 > 0), Effect2 > 0)
de.pos.neg <- subset(subset(de, Effect1 > 0), Effect2 < 0)
writeRNKformatCNA(rbind(de.pos.pos, de.pos.neg), wd.de.gsea, "SRC_SCLC_tpm-gene-median0_SORTING-CNA-TPM_q_n70_CNA-GAIN-POS")   ## GSEA

## POS
gsea.pos <- readTable(file.path(wd.de.gsea, "CNA_GAIN_POS", "HALLMARK", "gsea_report_for_na_pos_1657718679317.tsv"), header=T, rownames=F, sep="\t")
gsea.pos$NAME <- mapply(x = 1:nrow(gsea.pos), function(x) unlist(strsplit(gsea.pos$NAME[x], "HALLMARK_"))[2])
gsea.pos$NAME <- mapply(x = 1:nrow(gsea.pos), function(x) gsub("_", " ", gsea.pos$NAME[x]))
gsea.pos <- gsea.pos[5:1, c("NAME", "NES")]
gsea.pos[5, 1] <- "G2M checkpoint"
gsea.pos[4, 1] <- "E2F targets"
gsea.pos[3, 1] <- firstup(gsea.pos[3, 1])
gsea.pos[2, 1] <- "MYC targets v1"
gsea.pos[1, 1] <- "Unfolded protein"

main.text <- c("Hallmark enrichment score                         ", "")
file.de <- file.path(wd.de.gsea, "CNA_GAIN_POS", "HALLMARK_POS.pdf")
pdf(file.de, height=2, width=4)
par(mar=c(2,8.5,2,4))   # increase y-axis margin.
barplot(gsea.pos$NES, main=main.text[1], las=1, horiz=T, xlim=c(0, 4), xaxt="n", names.arg=gsea.pos$NAME, col=red.lighter, xlab="")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
#abline(v=2, lty=5)
axis(side=1, at=seq(0, 4, by=2), labels=F)
axis(side=1, at=seq(0, 4, by=2), pos=0.5, tick=F)
#mtext(main.text[2], line=0.3)
dev.off()

## NEG
gsea.neg <- readTable(file.path(wd.de.gsea, "CNA_GAIN_POS", "HALLMARK", "gsea_report_for_na_neg_1657718679317.tsv"), header=T, rownames=F, sep="\t")
gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) unlist(strsplit(gsea.neg$NAME[x], "HALLMARK_"))[2])
gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) gsub("_", " ", gsea.neg$NAME[x]))
gsea.neg <- gsea.neg[5:1, c("NAME", "NES")]
gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) firstup(gsea.neg$NAME[x]))

main.text <- c("", "")
file.de <- file.path(wd.de.gsea, "CNA_GAIN_POS", "HALLMARK_NEG.pdf")
pdf(file.de, height=1.6, width=4)
par(mar=c(2,1,0,11.5))   # increase y-axis margin.
posbar <- barplot(gsea.neg$NES, main=main.text[1], las=1, horiz=T, xlim=c(-4, 0), xaxt="n", names.arg="", col=blue.lighter, xlab="")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=gsea.neg$NAME, xpd=NA)
#abline(v=-2, lty=5)
axis(side=1, at=seq(-4, 0, by=2), labels=F)
axis(side=1, at=seq(-4, 0, by=2), pos=0.5, tick=F)
#mtext(main.text[2], line=0.3)
dev.off()

# -----------------------------------------------------------------------------
# Cannoli plot: CNA LOSS
# Last Modified: 26/06/22; 11/10/17
# -----------------------------------------------------------------------------
colnames <- c("GENE", "ADJ_1", "ADJ_2")
genes0 <- c("NFKB1")
genes <- toTable(NA, length(colnames), length(genes0), colnames)
genes$GENE <- genes0
#genes[1, 2] <- -0.1

de <- getCannoli(wd.de.data, BASE, 70, loss)
plot.de <- file.path(wd.de.plots, "cannoliplot_SRC_SCLC_TPM-CNA-SORTING_P1E03_MEDIAN0_CNA_LOSS_white_NFKB1")
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c(paste0("Lung-SCLC (CN < 2)"), "")
plotCannoli(de, pvalue, genes, file.de, file.main, xlab.text, ylab.text, "bottomleft", c("", ""), c(red.lighter, blue.lighter), c(red, blue), fold=0, pos="bottomright")

## GSEA
de.pos.pos <- subset(subset(de, Effect2 > 0), Effect1 > 0)
de.pos.neg <- subset(subset(de, Effect2 > 0), Effect1 < 0)
writeRNKformatCNA(rbind(de.pos.pos, de.pos.neg), wd.de.gsea, "SRC_SCLC_tpm-gene-median0_SORTING-CNA-TPM_q_n70_CNA-LOSS-POS")   ## GSEA

## POS
gsea.pos <- readTable(file.path(wd.de.gsea, "CNA_LOSS_POS", "HALLMARK", "gsea_report_for_na_pos_1657719534051.tsv"), header=T, rownames=F, sep="\t")
gsea.pos$NAME <- mapply(x = 1:nrow(gsea.pos), function(x) unlist(strsplit(gsea.pos$NAME[x], "HALLMARK_"))[2])
gsea.pos$NAME <- mapply(x = 1:nrow(gsea.pos), function(x) gsub("_", " ", gsea.pos$NAME[x]))
gsea.pos <- gsea.pos[5:1, c("NAME", "NES")]
gsea.pos[5, 1] <- "G2M checkpoint"
gsea.pos[4, 1] <- "E2F targets"
gsea.pos[3, 1] <- firstup(gsea.pos[3, 1])
gsea.pos[2, 1] <- firstup(gsea.pos[2, 1])
gsea.pos[1, 1] <- firstup(gsea.pos[1, 1])

#main.text <- c("", "")
main.text <- c("Hallmark enrichment score                         ", "")
file.de <- file.path(wd.de.gsea, "CNA_LOSS_POS", "HALLMARK_POS.pdf")
pdf(file.de, height=2, width=4)
par(mar=c(2,9.5,2,3))   # increase y-axis margin.
barplot(gsea.pos$NES, main=main.text[1], las=1, horiz=T, xlim=c(0, 4), xaxt="n", names.arg=gsea.pos$NAME, col=red.lighter, xlab="")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
#abline(v=2, lty=5)
axis(side=1, at=seq(0, 4, by=2), labels=F)
axis(side=1, at=seq(0, 4, by=2), pos=0.5, tick=F)
#mtext(main.text[2], line=0.3)
dev.off()

## NEG
gsea.neg <- readTable(file.path(wd.de.gsea, "CNA_LOSS_POS", "HALLMARK", "gsea_report_for_na_neg_1657719534051.tsv"), header=T, rownames=F, sep="\t")
gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) unlist(strsplit(gsea.neg$NAME[x], "HALLMARK_"))[2])
gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) gsub("_", " ", gsea.neg$NAME[x]))
gsea.neg <- gsea.neg[5:1, c("NAME", "NES")]
gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) firstup(gsea.neg$NAME[x]))
gsea.neg[1, 1] <- "UV response up"
gsea.neg[4, 1] <- "MYC targets v1"

#main.text <- c("                              Hallmark enrichment score", "")
main.text <- c("", "")
file.de <- file.path(wd.de.gsea, "CNA_LOSS_POS", "HALLMARK_NEG.pdf")
pdf(file.de, height=1.6, width=4)
par(mar=c(2,2,0,10.5))   # increase y-axis margin.
posbar <- barplot(gsea.neg$NES, main=main.text[1], las=1, horiz=T, xlim=c(-4, 0), xaxt="n", names.arg="", col=blue.lighter, xlab="")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=gsea.neg$NAME, xpd=NA)
#abline(v=-2, lty=5)
axis(side=1, at=seq(-4, 0, by=2), labels=F)
axis(side=1, at=seq(-4, 0, by=2), pos=0.5, tick=F)
#mtext(main.text[2], line=0.3)
dev.off()




















# -----------------------------------------------------------------------------
# Cannoli plot (RT)
# Last Modified: 02/07/22
# -----------------------------------------------------------------------------

## Not expressed
de <- getCannoli(wd.de.data, BASE, 70, setdiff(rownames(src.tpm.gene), expressed))
plot.de <- file.path(wd.de.plots, "cannoliplot_SRC_SCLC_TPM-CNA-SORTING_P1E03_MEDIAN0-0")
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c(paste0(BASE, " not expressed genes"), "")
plotCannoli(de, pvalue, toTable(NA, length(colnames), 0, colnames), file.de, file.main, xlab.text, ylab.text, "topleft", c("", ""), c(red, blue), c(red, blue), fold=0)

## Early-replicating
de <- getCannoli(wd.de.data, BASE, 70, gene.sclc.e)
plot.de <- file.path(wd.de.plots, "cannoliplot_SRC_SCLC_TPM-CNA-SORTING_P1E03_EARLY")
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c(paste0(BASE, " early-repli genes (n=", separator(nrow(de)), ")"), "")
plotCannoli(de, pvalue, genes, file.de, file.main, xlab.text, ylab.text, "topleft", c("", ""), c(red, blue), c(red, blue), fold=0)

## Late-replicating
de <- getCannoli(wd.de.data, BASE, 70, gene.sclc.l)
plot.de <- file.path(wd.de.plots, "cannoliplot_SRC_SCLC_TPM-CNA-SORTING_P1E03_LATE")
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c(paste0(BASE, " late-repli genes (n=", separator(nrow(de)), ")"), "")
plotCannoli(de, pvalue, genes, file.de, file.main, xlab.text, ylab.text, "topleft", c("", ""), c(red, blue), c(red, blue), fold=0)


# -----------------------------------------------------------------------------
# Cannoli plot (Late)
# Last Modified: 26/06/22; 11/10/17
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, "SRC_SCLC_tpm-gene-median0_CNA-vs-TPM_q_n70.RData"))
overlaps <- intersect(rownames(src.tpm.gene), gene.sclc.l)
de1 <- src.tpm.gene[overlaps,]
de1$Effect <- de1$RHO
load(file=file.path(wd.de.data, "SRC_SCLC_tpm-gene-median0_SORTING_q_n70.RData"))
de2 <- src.tpm.gene[overlaps,]
de2$Effect <- de2$RHO

overlaps <- intersect(rownames(de1), rownames(de2))
de <- as.data.frame(cbind(de1[overlaps, c("P", "Q", "Effect")], de2[overlaps, c("P", "Q", "Effect")]))
rownames(de) <- overlaps
colnames(de) <- c("P1", "Q1", "Effect1", "P2", "Q2", "Effect2")
de <- cbind(de1[overlaps,]$external_gene_name, de)
colnames(de) <- c("external_gene_name", "P1", "Q1", "Effect1", "P2", "Q2", "Effect2")

de.positive <- subset(subset(de, Effect1 > 0), Effect2 > 0)
de.positive.sig <- subset(subset(de.positive, P1 <= 0.001), P2 <= 0.001)


###
## Cannoli plots
xlab.text <- "Expression vs. CNA"
ylab.text <- expression("Expression vs."~italic('in silico')~"sorting")

plot.de <- file.path(wd.de.plots, "cannoliplot_SRC_SCLC_median0_TPM-CNA-SORTING_P1E03_LATE")
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("SCLC late-repli genes (n=4,009)", "")
plotCannoli(de1, de2, 0.001, genes, file.de, file.main, xlab.text, ylab.text, "topleft", c("", ""), c(red, blue), c(red, blue), fold=0)

# -----------------------------------------------------------------------------
# Cannoli plot (Early)
# Last Modified: 26/06/22; 11/10/17
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, "SRC_SCLC_tpm-gene-median0_CNA-vs-TPM_q_n70.RData"))
overlaps <- intersect(rownames(src.tpm.gene), )
de1 <- src.tpm.gene[overlaps,]
de1$Effect <- de1$RHO
load(file=file.path(wd.de.data, "SRC_SCLC_tpm-gene-median0_SORTING_q_n70.RData"))
de2 <- src.tpm.gene[overlaps,]
de2$Effect <- de2$RHO

overlaps <- intersect(rownames(de1), rownames(de2))
de <- as.data.frame(cbind(de1[overlaps, c("P", "Q", "Effect")], de2[overlaps, c("P", "Q", "Effect")]))
rownames(de) <- overlaps
colnames(de) <- c("P1", "Q1", "Effect1", "P2", "Q2", "Effect2")
de <- cbind(de1[overlaps,]$external_gene_name, de)
colnames(de) <- c("external_gene_name", "P1", "Q1", "Effect1", "P2", "Q2", "Effect2")

de.positive <- subset(subset(de, Effect1 > 0), Effect2 > 0)
de.positive.sig <- subset(subset(de.positive, P1 <= 0.001), P2 <= 0.001)

###
## Cannoli plots
xlab.text <- "Expression vs. CNA"
ylab.text <- expression("Expression vs."~italic('in silico')~"sorting")


# -----------------------------------------------------------------------------
# Cannoli plot (IZ)
# Last Modified: 26/06/22; 11/10/17
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, "SRC_SCLC_tpm-gene-median0_CNA-vs-TPM_q_n70.RData"))
overlaps <- intersect(rownames(src.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.iz))
de1 <- src.tpm.gene[overlaps,]
de1$Effect <- de1$RHO
load(file=file.path(wd.de.data, "SRC_SCLC_tpm-gene-median0_SORTING_q_n70.RData"))
de2 <- src.tpm.gene[overlaps,]
de2$Effect <- de2$RHO

overlaps <- intersect(rownames(de1), rownames(de2))
de <- as.data.frame(cbind(de1[overlaps, c("P", "Q", "Effect")], de2[overlaps, c("P", "Q", "Effect")]))
rownames(de) <- overlaps
colnames(de) <- c("P1", "Q1", "Effect1", "P2", "Q2", "Effect2")
de <- cbind(de1[overlaps,]$external_gene_name, de)
colnames(de) <- c("external_gene_name", "P1", "Q1", "Effect1", "P2", "Q2", "Effect2")

de.positive <- subset(subset(de, Effect1 > 0), Effect2 > 0)
de.positive.sig <- subset(subset(de.positive, P1 <= 0.001), P2 <= 0.001)

###
## Cannoli plots
xlab.text <- "Expression vs. CNA"
ylab.text <- expression("Expression vs."~italic('in silico')~"sorting")

plot.de <- file.path(wd.de.plots, "cannoliplot_SRC_SCLC_median0_TPM-CNA-SORTING_P1E03_IZ")
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("SCLC IZ genes (n=3,075)", "")
plotCannoli(de1, de2, 0.001, genes, file.de, file.main, xlab.text, ylab.text, "topleft", c("", ""), c(red, blue), c(red, blue), fold=0)

# -----------------------------------------------------------------------------
# Cannoli plot (TZ)
# Last Modified: 26/06/22; 11/10/17
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, "SRC_SCLC_tpm-gene-median0_CNA-vs-TPM_q_n70.RData"))
overlaps <- intersect(rownames(src.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.tz))
de1 <- src.tpm.gene[overlaps,]
de1$Effect <- de1$RHO
load(file=file.path(wd.de.data, "SRC_SCLC_tpm-gene-median0_SORTING_q_n70.RData"))
de2 <- src.tpm.gene[overlaps,]
de2$Effect <- de2$RHO

overlaps <- intersect(rownames(de1), rownames(de2))
de <- as.data.frame(cbind(de1[overlaps, c("P", "Q", "Effect")], de2[overlaps, c("P", "Q", "Effect")]))
rownames(de) <- overlaps
colnames(de) <- c("P1", "Q1", "Effect1", "P2", "Q2", "Effect2")
de <- cbind(de1[overlaps,]$external_gene_name, de)
colnames(de) <- c("external_gene_name", "P1", "Q1", "Effect1", "P2", "Q2", "Effect2")

de.positive <- subset(subset(de, Effect1 > 0), Effect2 > 0)
de.positive.sig <- subset(subset(de.positive, P1 <= 0.001), P2 <= 0.001)

de.negative <- subset(subset(de, Effect1 < 0), Effect2 < 0)
de.negative.sig <- subset(subset(de.negative, P1 <= 0.001), P2 <= 0.001)

###
## Cannoli plots
xlab.text <- "Expression vs. CNA"
ylab.text <- expression("Expression vs."~italic('in silico')~"sorting")

plot.de <- file.path(wd.de.plots, "cannoliplot_SRC_SCLC_median0_TPM-CNA-SORTING_P1E03_TZ")
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("SCLC TZ genes (n=1,944)", "")
plotCannoli(de1, de2, 0.001, genes, file.de, file.main, xlab.text, ylab.text, "topleft", c("", ""), c(red, blue), c(red, blue), fold=0)

# -----------------------------------------------------------------------------
# Cannoli plot (TTR)
# Last Modified: 26/06/22; 11/10/17
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, "SRC_SCLC_tpm-gene-median0_CNA-vs-TPM_q_n70.RData"))
overlaps <- intersect(rownames(src.tpm.gene), rownames(tpm.gene.log2.m.rfd.ttr))
de1 <- src.tpm.gene[overlaps,]
de1$Effect <- de1$RHO
load(file=file.path(wd.de.data, "SRC_SCLC_tpm-gene-median0_SORTING_q_n70.RData"))
de2 <- src.tpm.gene[overlaps,]
de2$Effect <- de2$RHO

overlaps <- intersect(rownames(de1), rownames(de2))
de <- as.data.frame(cbind(de1[overlaps, c("P", "Q", "Effect")], de2[overlaps, c("P", "Q", "Effect")]))
rownames(de) <- overlaps
colnames(de) <- c("P1", "Q1", "Effect1", "P2", "Q2", "Effect2")
de <- cbind(de1[overlaps,]$external_gene_name, de)
colnames(de) <- c("external_gene_name", "P1", "Q1", "Effect1", "P2", "Q2", "Effect2")

de.positive <- subset(subset(de, Effect1 > 0), Effect2 > 0)
de.positive.sig <- subset(subset(de.positive, P1 <= 0.001), P2 <= 0.001)

de.negative <- subset(subset(de, Effect1 < 0), Effect2 < 0)
de.negative.sig <- subset(subset(de.negative, P1 <= 0.001), P2 <= 0.001)

###
## Cannoli plots
xlab.text <- "Expression vs. CNA"
ylab.text <- expression("Expression vs."~italic('in silico')~"sorting")

plot.de <- file.path(wd.de.plots, "cannoliplot_SRC_SCLC_median0_TPM-CNA-SORTING_P1E03_TTR")
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("SCLC TTR genes (n=17,001)", "")
plotCannoli(de1, de2, 0.001, genes, file.de, file.main, xlab.text, ylab.text, "topleft", c("", ""), c(red, blue), c(red, blue), fold=0)

# -----------------------------------------------------------------------------
# Cannoli plot (P < 0.01)
# Last Modified: 26/06/22; 11/10/17
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, "SRC_SCLC_tpm-gene-median0_CNA-vs-TPM_q_n70.RData"))
de1 <- src.tpm.gene
de1$Effect <- de1$RHO
load(file=file.path(wd.de.data, "SRC_SCLC_tpm-gene-median0_SORTING_q_n70.RData"))
de2 <- src.tpm.gene
de2$Effect <- de2$RHO

###
## Cannoli plots
xlab.text <- "Expression vs. CNA"
ylab.text <- expression("Expression vs."~italic('in silico')~"sorting")

plot.de <- file.path(wd.de.plots, "cannoliplot_SRC_SCLC_median0_TPM-CNA-SORTING_P1E02")
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("SCLC expressed genes (n=23,572)", "")
plotCannoli(de1, de2, 0.01, genes, file.de, file.main, xlab.text, ylab.text, "topleft", c("", ""), c(red, blue), c(red, blue), fold=0)

# -----------------------------------------------------------------------------
# Cannoli plot (Late)
# Last Modified: 26/06/22; 11/10/17
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, "SRC_SCLC_tpm-gene-median0_CNA-vs-TPM_q_n70.RData"))
overlaps <- intersect(rownames(src.tpm.gene), gene.sclc.l)
de1 <- src.tpm.gene[overlaps,]
de1$Effect <- de1$RHO
load(file=file.path(wd.de.data, "SRC_SCLC_tpm-gene-median0_SORTING_q_n70.RData"))
de2 <- src.tpm.gene[overlaps,]
de2$Effect <- de2$RHO

overlaps <- intersect(rownames(de1), rownames(de2))
de <- as.data.frame(cbind(de1[overlaps, c("P", "Q", "Effect")], de2[overlaps, c("P", "Q", "Effect")]))
rownames(de) <- overlaps
colnames(de) <- c("P1", "Q1", "Effect1", "P2", "Q2", "Effect2")
de <- cbind(de1[overlaps,]$external_gene_name, de)
colnames(de) <- c("external_gene_name", "P1", "Q1", "Effect1", "P2", "Q2", "Effect2")

de.negative <- subset(subset(de, Effect1 < 0), Effect2 < 0)
de.negative.sig <- subset(subset(de.negative, P1 <= 0.01), P2 <= 0.01)

###
## Cannoli plots
xlab.text <- "Expression vs. CNA"
ylab.text <- expression("Expression vs."~italic('in silico')~"sorting")

plot.de <- file.path(wd.de.plots, "cannoliplot_SRC_SCLC_median0_TPM-CNA-SORTING_P1E02_LATE")
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("SCLC late-repli genes (n=4,009)", "")
plotCannoli(de1, de2, 0.01, genes, file.de, file.main, xlab.text, ylab.text, "topleft", c("", ""), c(red, blue), c(red, blue), fold=0)

# -----------------------------------------------------------------------------
# Cannoli plot (Early)
# Last Modified: 26/06/22; 11/10/17
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, "SRC_SCLC_tpm-gene-median0_CNA-vs-TPM_q_n70.RData"))
overlaps <- intersect(rownames(src.tpm.gene), gene.sclc.e)
de1 <- src.tpm.gene[overlaps,]
de1$Effect <- de1$RHO
load(file=file.path(wd.de.data, "SRC_SCLC_tpm-gene-median0_SORTING_q_n70.RData"))
de2 <- src.tpm.gene[overlaps,]
de2$Effect <- de2$RHO

overlaps <- intersect(rownames(de1), rownames(de2))
de <- as.data.frame(cbind(de1[overlaps, c("P", "Q", "Effect")], de2[overlaps, c("P", "Q", "Effect")]))
rownames(de) <- overlaps
colnames(de) <- c("P1", "Q1", "Effect1", "P2", "Q2", "Effect2")
de <- cbind(de1[overlaps,]$external_gene_name, de)
colnames(de) <- c("external_gene_name", "P1", "Q1", "Effect1", "P2", "Q2", "Effect2")

de.positive <- subset(subset(de, Effect1 > 0), Effect2 > 0)
de.positive.sig <- subset(subset(de.positive, P1 <= 0.01), P2 <= 0.01)

###
## Cannoli plots
xlab.text <- "Expression vs. CNA"
ylab.text <- expression("Expression vs."~italic('in silico')~"sorting")

plot.de <- file.path(wd.de.plots, "cannoliplot_SRC_SCLC_median0_TPM-CNA-SORTING_P1E02_EARLY")
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("SCLC early-repli genes (n=18,011)", "")
plotCannoli(de1, de2, 0.01, genes, file.de, file.main, xlab.text, ylab.text, "topleft", c("", ""), c(red, blue), c(red, blue), fold=0)

# -----------------------------------------------------------------------------
# Cannoli plot (IZ)
# Last Modified: 26/06/22; 11/10/17
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, "SRC_SCLC_tpm-gene-median0_CNA-vs-TPM_q_n70.RData"))
overlaps <- intersect(rownames(src.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.iz))
de1 <- src.tpm.gene[overlaps,]
de1$Effect <- de1$RHO
load(file=file.path(wd.de.data, "SRC_SCLC_tpm-gene-median0_SORTING_q_n70.RData"))
de2 <- src.tpm.gene[overlaps,]
de2$Effect <- de2$RHO

overlaps <- intersect(rownames(de1), rownames(de2))
de <- as.data.frame(cbind(de1[overlaps, c("P", "Q", "Effect")], de2[overlaps, c("P", "Q", "Effect")]))
rownames(de) <- overlaps
colnames(de) <- c("P1", "Q1", "Effect1", "P2", "Q2", "Effect2")
de <- cbind(de1[overlaps,]$external_gene_name, de)
colnames(de) <- c("external_gene_name", "P1", "Q1", "Effect1", "P2", "Q2", "Effect2")

de.positive <- subset(subset(de, Effect1 > 0), Effect2 > 0)
de.positive.sig <- subset(subset(de.positive, P1 <= 0.01), P2 <= 0.01)

de.negative <- subset(subset(de, Effect1 < 0), Effect2 < 0)
de.negative.sig <- subset(subset(de.negative, P1 <= 0.01), P2 <= 0.01)

###
## Cannoli plots
xlab.text <- "Expression vs. CNA"
ylab.text <- expression("Expression vs."~italic('in silico')~"sorting")

plot.de <- file.path(wd.de.plots, "cannoliplot_SRC_SCLC_median0_TPM-CNA-SORTING_P1E02_IZ")
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("SCLC IZ genes (n=3,075)", "")
plotCannoli(de1, de2, 0.01, genes, file.de, file.main, xlab.text, ylab.text, "topleft", c("", ""), c(red, blue), c(red, blue), fold=0)

# -----------------------------------------------------------------------------
# Cannoli plot (TZ)
# Last Modified: 26/06/22; 11/10/17
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, "SRC_SCLC_tpm-gene-median0_CNA-vs-TPM_q_n70.RData"))
overlaps <- intersect(rownames(src.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.tz))
de1 <- src.tpm.gene[overlaps,]
de1$Effect <- de1$RHO
load(file=file.path(wd.de.data, "SRC_SCLC_tpm-gene-median0_SORTING_q_n70.RData"))
de2 <- src.tpm.gene[overlaps,]
de2$Effect <- de2$RHO

overlaps <- intersect(rownames(de1), rownames(de2))
de <- as.data.frame(cbind(de1[overlaps, c("P", "Q", "Effect")], de2[overlaps, c("P", "Q", "Effect")]))
rownames(de) <- overlaps
colnames(de) <- c("P1", "Q1", "Effect1", "P2", "Q2", "Effect2")
de <- cbind(de1[overlaps,]$external_gene_name, de)
colnames(de) <- c("external_gene_name", "P1", "Q1", "Effect1", "P2", "Q2", "Effect2")

de.positive <- subset(subset(de, Effect1 > 0), Effect2 > 0)
de.positive.sig <- subset(subset(de.positive, P1 <= 0.01), P2 <= 0.01)

de.negative <- subset(subset(de, Effect1 < 0), Effect2 < 0)
de.negative.sig <- subset(subset(de.negative, P1 <= 0.01), P2 <= 0.01)

###
## Cannoli plots
xlab.text <- "Expression vs. CNA"
ylab.text <- expression("Expression vs."~italic('in silico')~"sorting")

plot.de <- file.path(wd.de.plots, "cannoliplot_SRC_SCLC_median0_TPM-CNA-SORTING_P1E02_TZ")
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("SCLC TZ genes (n=1,944)", "")
plotCannoli(de1, de2, 0.01, genes, file.de, file.main, xlab.text, ylab.text, "topleft", c("", ""), c(red, blue), c(red, blue), fold=0)

# -----------------------------------------------------------------------------
# Cannoli plot (TTR)
# Last Modified: 26/06/22; 11/10/17
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, "SRC_SCLC_tpm-gene-median0_CNA-vs-TPM_q_n70.RData"))
overlaps <- intersect(rownames(src.tpm.gene), rownames(tpm.gene.log2.m.rfd.ttr))
de1 <- src.tpm.gene[overlaps,]
de1$Effect <- de1$RHO
load(file=file.path(wd.de.data, "SRC_SCLC_tpm-gene-median0_SORTING_q_n70.RData"))
de2 <- src.tpm.gene[overlaps,]
de2$Effect <- de2$RHO

overlaps <- intersect(rownames(de1), rownames(de2))
de <- as.data.frame(cbind(de1[overlaps, c("P", "Q", "Effect")], de2[overlaps, c("P", "Q", "Effect")]))
rownames(de) <- overlaps
colnames(de) <- c("P1", "Q1", "Effect1", "P2", "Q2", "Effect2")
de <- cbind(de1[overlaps,]$external_gene_name, de)
colnames(de) <- c("external_gene_name", "P1", "Q1", "Effect1", "P2", "Q2", "Effect2")

de.positive <- subset(subset(de, Effect1 > 0), Effect2 > 0)
de.positive.sig <- subset(subset(de.positive, P1 <= 0.01), P2 <= 0.01)

de.negative <- subset(subset(de, Effect1 < 0), Effect2 < 0)
de.negative.sig <- subset(subset(de.negative, P1 <= 0.01), P2 <= 0.01)

###
## Cannoli plots
xlab.text <- "Expression vs. CNA"
ylab.text <- expression("Expression vs."~italic('in silico')~"sorting")

plot.de <- file.path(wd.de.plots, "cannoliplot_SRC_SCLC_median0_TPM-CNA-SORTING_P1E02_TTR")
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("SCLC TTR genes (n=17,001)", "")
plotCannoli(de1, de2, 0.01, genes, file.de, file.main, xlab.text, ylab.text, "topleft", c("", ""), c(red, blue), c(red, blue), fold=0)

















# -----------------------------------------------------------------------------
# RFD vs. TPM
# Last Modified: 31/10/18
# -----------------------------------------------------------------------------
nrds.RT.NRFD <- nrds.RT.NRFD.sclc   ## 17/05/20 WHY nrds.RT.NRFD.sclc.nl? ## MUY MUY IMPORTANTE!!

tpm.gene.log2.m <- tpm.gene.log2.m
nrow(tpm.gene.log2.m)
# [1] 34908

tpm.gene.log2.m <- tpm.gene.log2.m[rownames(de.tpm.gene),]
nrow(tpm.gene.log2.m)
# [1] 23572

tpm.gene.log2.m <- tpm.gene.log2.m[rownames(de.tpm.gene),]
nrow(tpm.gene.log2.m)
# [1] 19131

tpm.gene.log2.m.rfd <- getTRC(tpm.gene.log2.m, nrds.RT.NRFD)
tpm.gene.log2.m.rfd$length <- abs(tpm.gene.log2.m.rfd$end_position - tpm.gene.log2.m.rfd$start_position)
nrow(tpm.gene.log2.m.rfd)
# [1] 31616
# [1] 22023
# [1] 18007
nrow(subset(tpm.gene.log2.m.rfd, TRC == 0))
# [1] 3
# [1] 1
# [1] 1

###
## TTR and CTR (IZ +TZ)
save(tpm.gene.log2.m.rfd, file=file.path(wd.de.data, "tpm_gene_log2_m_rfd0.RData"))
#save(tpm.gene.log2.m.rfd, file=file.path(wd.de.data, "tpm_gene_median0_log2_m_rfd0.RData"))
#save(tpm.gene.log2.m.rfd, file=file.path(wd.de.data, "tpm_gene_r5p47_log2_m_rfd0.RData"))

load(file=file.path(wd.de.data, "tpm_gene_log2_m_rfd0.RData"))
#load(file=file.path(wd.de.data, "tpm_gene_median0_log2_m_rfd0.RData"))
#load(file=file.path(wd.de.data, "tpm_gene_r5p47_log2_m_rfd0.RData"))
file.name <- file.path(wd.de.data, "tpm_gene_log2_m_rfd.RData")
#file.name <- file.path(wd.de.data, "tpm_gene_median0_log2_m_rfd.RData")
#file.name <- file.path(wd.de.data, "tpm_gene_r5p47_log2_m_rfd.RData")
setTRC(tpm.gene.log2.m.rfd, rfd=0.9, file.name)
load(file.name)

## CTR
length(which(tpm.gene.log2.m.rfd.ctr.iz$GENE_NRFD < 0))
# [1] 83
# [1] 70
# [1] 60
length(which(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD >= 0))
# [1] 95
# [1] 85
# [1] 80

# -----------------------------------------------------------------------------
# Shared genes (All)
# Last Modified: 10/09/20
# -----------------------------------------------------------------------------
load("/Users/tpyang/Work/uni-koeln/tyang2/SCLC/analysis/expression/kallisto/sclc-tpm-de/data/tpm_gene_log2_m_rfd.RData")

gene.sclc   <- list(rownames(tpm.gene.log2.m.rfd.ctr.iz),   rownames(tpm.gene.log2.m.rfd.ctr.tz),   rownames(tpm.gene.log2.m.rfd.ttr))
gene.sclc.e <- list(rownames(tpm.gene.log2.m.rfd.ctr.iz.e), rownames(tpm.gene.log2.m.rfd.ctr.tz.e), rownames(tpm.gene.log2.m.rfd.ttr.e))
gene.sclc.l <- list(rownames(tpm.gene.log2.m.rfd.ctr.iz.l), rownames(tpm.gene.log2.m.rfd.ctr.tz.l), rownames(tpm.gene.log2.m.rfd.ttr.l))

length(intersect(intersect(gene.sclc[[1]], gene.nbl[[1]]), gene.cll[[1]]))/length(unique(unique(gene.sclc[[1]], gene.nbl[[1]]), gene.cll[[1]]))
# [1] 0.2006778
length(intersect(intersect(gene.sclc[[2]], gene.nbl[[2]]), gene.cll[[2]]))/length(unique(unique(gene.sclc[[2]], gene.nbl[[2]]), gene.cll[[2]]))
# [1] 0.1838505
length(intersect(intersect(gene.sclc[[3]], gene.nbl[[3]]), gene.cll[[3]]))/length(unique(unique(gene.sclc[[3]], gene.nbl[[3]]), gene.cll[[3]]))
# [1] 0.6462588

length(intersect(intersect(gene.sclc.e[[1]], gene.nbl.e[[1]]), gene.cll.e[[1]]))/length(unique(unique(gene.sclc.e[[1]], gene.nbl.e[[1]]), gene.cll.e[[1]]))
# [1] 0.1994245
length(intersect(intersect(gene.sclc.e[[2]], gene.nbl.e[[2]]), gene.cll.e[[2]]))/length(unique(unique(gene.sclc.e[[2]], gene.nbl.e[[2]]), gene.cll.e[[2]]))
# [1] 0.1519174
length(intersect(intersect(gene.sclc.e[[3]], gene.nbl.e[[3]]), gene.cll.e[[3]]))/length(unique(unique(gene.sclc.e[[3]], gene.nbl.e[[3]]), gene.cll.e[[3]]))
# [1] 0.5306519

length(intersect(intersect(gene.sclc.l[[1]], gene.nbl.l[[1]]), gene.cll.l[[1]]))/length(unique(unique(gene.sclc.l[[1]], gene.nbl.l[[1]]), gene.cll.l[[1]]))
# [1] 0.132622
length(intersect(intersect(gene.sclc.l[[2]], gene.nbl.l[[2]]), gene.cll.l[[2]]))/length(unique(unique(gene.sclc.l[[2]], gene.nbl.l[[2]]), gene.cll.l[[2]]))
# [1] 0.1619938
length(intersect(intersect(gene.sclc.l[[3]], gene.nbl.l[[3]]), gene.cll.l[[3]]))/length(unique(unique(gene.sclc.l[[3]], gene.nbl.l[[3]]), gene.cll.l[[3]]))
# [1] 0.6378696

###
## Early
1 - (length(intersect(intersect(gene.sclc.e[[1]], gene.nbl.e[[1]]), gene.cll.e[[1]]))/length(unique(unique(gene.sclc.e[[1]], gene.nbl.e[[1]]), gene.cll.e[[1]])))
# [1] 0.8005755
1 - (length(intersect(intersect(gene.sclc.e[[2]], gene.nbl.e[[2]]), gene.cll.e[[2]]))/length(unique(unique(gene.sclc.e[[2]], gene.nbl.e[[2]]), gene.cll.e[[2]])))
# [1] 0.8480826
1 - (length(intersect(intersect(gene.sclc.e[[3]], gene.nbl.e[[3]]), gene.cll.e[[3]]))/length(unique(unique(gene.sclc.e[[3]], gene.nbl.e[[3]]), gene.cll.e[[3]])))
# [1] 0.4693481

###
## Late
1 - (length(intersect(intersect(gene.sclc.l[[1]], gene.nbl.l[[1]]), gene.cll.l[[1]]))/length(unique(unique(gene.sclc.l[[1]], gene.nbl.l[[1]]), gene.cll.l[[1]])))
# [1] 0.867378
1 - (length(intersect(intersect(gene.sclc.l[[2]], gene.nbl.l[[2]]), gene.cll.l[[2]]))/length(unique(unique(gene.sclc.l[[2]], gene.nbl.l[[2]]), gene.cll.l[[2]])))
# [1] 0.8380062
1 - (length(intersect(intersect(gene.sclc.l[[3]], gene.nbl.l[[3]]), gene.cll.l[[3]]))/length(unique(unique(gene.sclc.l[[3]], gene.nbl.l[[3]]), gene.cll.l[[3]])))
# [1] 0.3621304

###
## Fisher's
test <- toTable(0, 2, 2, c("IZ", "TZ"))
test[1,] <- c(20, 15)
test[2,] <- c(13, 16)
 
fisher.test(test)[[1]]



# -----------------------------------------------------------------------------
# Density (All)
# Last Modified: 18/03/21
# -----------------------------------------------------------------------------
nrow(tpm.gene.log2.m.rfd.ctr.iz)/nrow(nrds.RT.NRFD.sclc.ctr.iz)*1000
# [1] 16.3995
nrow(tpm.gene.log2.m.rfd.ctr.tz)/nrow(nrds.RT.NRFD.sclc.ctr.tz)*1000
# [1] 11.61529
nrow(tpm.gene.log2.m.rfd.ttr)/nrow(nrds.RT.NRFD.sclc.ttr)*1000
# [1] 11.44295

nrow(tpm.gene.log2.m.rfd.ctr.iz.e)/nrow(nrds.RT.NRFD.sclc.ctr.iz.e)*1000
# [1] 21.0098
nrow(tpm.gene.log2.m.rfd.ctr.iz.l)/nrow(nrds.RT.NRFD.sclc.ctr.iz.l)*1000
# [1] 7.583903

nrow(tpm.gene.log2.m.rfd.ctr.tz.e)/nrow(nrds.RT.NRFD.sclc.ctr.tz.e)*1000
# [1] 16.20317
nrow(tpm.gene.log2.m.rfd.ctr.tz.l)/nrow(nrds.RT.NRFD.sclc.ctr.tz.l)*1000
# [1] 7.268418

nrow(tpm.gene.log2.m.rfd.ttr.e)/nrow(nrds.RT.NRFD.sclc.ttr.e)*1000
# [1] 15.26942
nrow(tpm.gene.log2.m.rfd.ttr.l)/nrow(nrds.RT.NRFD.sclc.ttr.l)*1000
# [1] 6.733779

# -----------------------------------------------------------------------------
# E vs. L (All)
# Last Modified: 30/08/20
# -----------------------------------------------------------------------------
nrow(tpm.gene.log2.m.rfd.ttr)/31616
# [1] 0.774418
nrow(tpm.gene.log2.m.rfd.ctr.iz)/31616
# [1] 0.1306617
nrow(tpm.gene.log2.m.rfd.ctr.tz)/31616
# [1] 0.09479378

nrow(tpm.gene.log2.m.rfd.ttr.e)/31616
# [1] 0.5701227
nrow(tpm.gene.log2.m.rfd.ttr.l)/31616
# [1] 0.2042953

nrow(tpm.gene.log2.m.rfd.ctr.iz.e)/31616
# [1] 0.1099127
nrow(tpm.gene.log2.m.rfd.ctr.iz.l)/31616
# [1] 0.02074899

nrow(tpm.gene.log2.m.rfd.ctr.tz.e)/31616
# [1] 0.06433451
nrow(tpm.gene.log2.m.rfd.ctr.tz.l)/31616
# [1] 0.03045926

###
## median0
nrow(tpm.gene.log2.m.rfd.ttr)/22023
# [1] 0.7719657
nrow(tpm.gene.log2.m.rfd.ctr.iz)/22023
# [1] 0.1396268
nrow(tpm.gene.log2.m.rfd.ctr.tz)/22023
# [1] 0.08827135

nrow(tpm.gene.log2.m.rfd.ctr.iz.e)/22023
# [1] 0.1248695
nrow(tpm.gene.log2.m.rfd.ctr.iz.l)/22023
# [1] 0.0147573
 
nrow(tpm.gene.log2.m.rfd.ctr.tz.e)/22023
# [1] 0.07038097
nrow(tpm.gene.log2.m.rfd.ctr.tz.l)/22023
# [1] 0.01789039

###
## r5p47
nrow(tpm.gene.log2.m.rfd.ttr)/18007
# [1] 0.7698117
nrow(tpm.gene.log2.m.rfd.ctr.iz)/18007
# [1] 0.1424446
nrow(tpm.gene.log2.m.rfd.ctr.tz)/18007
# [1] 0.08757705

nrow(tpm.gene.log2.m.rfd.ctr.iz.e)/18007
# [1] 
nrow(tpm.gene.log2.m.rfd.ctr.iz.l)/18007
# [1] 

nrow(tpm.gene.log2.m.rfd.ctr.tz.e)/18007
# [1] 
nrow(tpm.gene.log2.m.rfd.ctr.tz.l)/18007
# [1] 

# -----------------------------------------------------------------------------
# RFD vs. TPM (All)
# Last Modified: 27/11/20; 01/09/20; 29/05/20
# -----------------------------------------------------------------------------
file.name <- file.path(wd.de.data, "tpm_gene_log2_m_rfd+1.RData")
load(file.name)

ylim <- c(min(tpm.gene.log2.m.rfd$MEDIAN), 11)

file.name <- paste0("boxplot3_sclc_tpm.gene+1_RFD_")
plotBox3(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr, tpm.gene.log2.m.rfd.ctr.tz, tpm.gene.log2.m.rfd.ctr.iz, main="SCLC expression", names=c("TTR", "TZ", "IZ"), cols=c("black", blue, red), ylim, text="")

file.name <- paste0("boxplot4_sclc_tpm.gene+1_RFD")
plotBox4(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.l, tpm.gene.log2.m.rfd.ctr.tz.e, tpm.gene.log2.m.rfd.ctr.iz.l, tpm.gene.log2.m.rfd.ctr.iz.e, main="SCLC expression (CTR)", names=c("L", "E", "L", "E"), cols=c(blue, blue, red, red), ylim)

file.name <- paste0("boxplot6_sclc_tpm.gene+1_RFD")
plotBox6(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr.l, tpm.gene.log2.m.rfd.ttr.e, tpm.gene.log2.m.rfd.ctr.tz.l, tpm.gene.log2.m.rfd.ctr.tz.e, tpm.gene.log2.m.rfd.ctr.iz.l, tpm.gene.log2.m.rfd.ctr.iz.e, main="SCLC expression", names=c("L", "E", "L", "E", "L", "E"), cols=c("black", "black", blue, blue, red, red), ylim)

# -----------------------------------------------------------------------------
# Significant mutated genes (All)
# Last Modified: 01/12/20
# -----------------------------------------------------------------------------
muts <- unique(readTable(file.path(wd.meta, "nature14664-s1_ST6_Figure1.txt"), header=F, rownames=F, sep=""))
genes <- getGenes(muts)

colnames <- c("TTR", "TZ", "IZ")
report <- toTable("", length(colnames), 2, colnames)
rownames(report) <- c("E", "L")

getGeneRFD <- function(genes, tpm.gene.log2.m.rfd.ttr.e) {
   ids <- intersect(rownames(genes), rownames(tpm.gene.log2.m.rfd.ttr.e))
   
   return(tpm.gene.log2.m.rfd.ttr.e[ids,]$external_gene_name)
}

report["E", "TTR"] <- toString(getGeneRFD(genes, tpm.gene.log2.m.rfd.ttr.e))
report["L", "TTR"] <- toString(getGeneRFD(genes, tpm.gene.log2.m.rfd.ttr.l))
report["E", "TZ"]  <- toString(getGeneRFD(genes, tpm.gene.log2.m.rfd.ctr.tz.e))
report["L", "TZ"]  <- toString(getGeneRFD(genes, tpm.gene.log2.m.rfd.ctr.tz.l))
report["E", "IZ"]  <- toString(getGeneRFD(genes, tpm.gene.log2.m.rfd.ctr.iz.e))
report["L", "IZ"]  <- toString(getGeneRFD(genes, tpm.gene.log2.m.rfd.ctr.iz.l))

# -----------------------------------------------------------------------------
# RFD vs. TPM (All)
# Last Modified: 01/09/20; 29/05/20
# -----------------------------------------------------------------------------
ylim <- c(min(tpm.gene.log2.m.rfd$MEDIAN), 14.5)

file.name <- paste0("boxplot0_sclc_tpm.gene_rfd_TTR+IZ")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr, tpm.gene.log2.m.rfd.ctr.iz, main="SCLC genes", names=c("TTR", "IZ"), cols=c("black", red), ylim)

file.name <- paste0("boxplot0_sclc_tpm.gene_rfd_TTR+TZ")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr, tpm.gene.log2.m.rfd.ctr.tz, main="SCLC genes", names=c("TTR", "TZ"), cols=c("black", blue), ylim)

file.name <- paste0("boxplot0_sclc_tpm.gene_rfd_TZ+IZ")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz, tpm.gene.log2.m.rfd.ctr.iz, main="SCLC genes (CTR)", names=c("TZ", "IZ"), cols=c(blue, red), ylim)

##
file.name <- paste0("boxplot0_sclc_tpm.gene_rfd_IZ_E+L")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.l, tpm.gene.log2.m.rfd.ctr.iz.e, main="SCLC IZ genes", names=c("Late", "Early"), cols=c(red, red), ylim)

file.name <- paste0("boxplot0_sclc_tpm.gene_rfd_TZ_E+L")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.l, tpm.gene.log2.m.rfd.ctr.tz.e, main="SCLC TZ genes", names=c("Late", "Early"), cols=c(blue, blue), ylim)

file.name <- paste0("boxplot0_sclc_tpm.gene_rfd_TTR_E+L")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr.l, tpm.gene.log2.m.rfd.ttr.e, main="SCLC TTR genes", names=c("Late", "Early"), cols=c("black", "black"), ylim)

# -----------------------------------------------------------------------------
# RFD vs. TPM (r5p47)
# Last Modified: 22/02/21; 02/09/20; 25/08/20; 29/05/20
# -----------------------------------------------------------------------------
ylim <- c(min(tpm.gene.log2.m.rfd$MEDIAN), 14.5)

file.name <- paste0("boxplot0_sclc_tpm.gene.r5p47_rfd_TTR+IZ")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr, tpm.gene.log2.m.rfd.ctr.iz, main="SCLC genes", names=c("TTR", "IZ"), cols=c("black", red), ylim)
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_rfd_TTR+IZ")
plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr, tpm.gene.log2.m.rfd.ctr.iz, main="SCLC genes", names=c("TTR", "IZ"), cols=c("black", red), ylim)

file.name <- paste0("boxplot0_sclc_tpm.gene.r5p47_rfd_TTR+TZ")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr, tpm.gene.log2.m.rfd.ctr.tz, main="SCLC genes", names=c("TTR", "TZ"), cols=c("black", blue), ylim)

#file.name <- paste0("boxplot0_sclc_tpm.gene.median0_rfd_IZ+TZ")
#plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz, tpm.gene.log2.m.rfd.ctr.tz, main="SCLC genes (CTR)", names=c("IZ", "TZ"), cols=c(red, blue), ylim)

file.name <- paste0("boxplot0_sclc_tpm.gene.r5p47_rfd_TZ+IZ")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz, tpm.gene.log2.m.rfd.ctr.iz, main="SCLC genes (CTR)", names=c("TZ", "IZ"), cols=c(blue, red), ylim)
#file.name <- paste0("boxplot_sclc_tpm.gene.median0_rfd_TZ+IZ_long")
#plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz, tpm.gene.log2.m.rfd.ctr.iz, main="SCLC genes", names=c("TZ", "IZ"), cols=c(blue, red), ylim)

##
file.name <- paste0("boxplot0_sclc_tpm.gene.r5p47_rfd_IZ_E+L")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.l, tpm.gene.log2.m.rfd.ctr.iz.e, main="SCLC IZ genes", names=c("Late", "Early"), cols=c(red, red), ylim)
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_rfd_IZ_E+L")
plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.l, tpm.gene.log2.m.rfd.ctr.iz.e, main="IZ", names=c("Late", "Early"), cols=c(red, red), ylim)

file.name <- paste0("boxplot0_sclc_tpm.gene.r5p47_rfd_IZ_E_HO+CD")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.e.ho, tpm.gene.log2.m.rfd.ctr.iz.e.cd, main="SCLC early IZ genes", names=c("HO", "CD"), cols=c(red, red), ylim)
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_rfd_IZ_E_HO+CD")
plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.e.ho, tpm.gene.log2.m.rfd.ctr.iz.e.cd, main="Early IZ", names=c("HO", "CD"), cols=c(red, red), ylim)

##
file.name <- paste0("boxplot0_sclc_tpm.gene.r5p47_rfd_TZ_E+L")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.l, tpm.gene.log2.m.rfd.ctr.tz.e, main="SCLC TZ genes", names=c("Late", "Early"), cols=c(blue, blue), ylim)
#file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_rfd_TZ_E+L_long")
#plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.l, tpm.gene.log2.m.rfd.ctr.tz.e, main="SCLC TZ genes", names=c("Late", "Early"), cols=c(blue, blue), ylim)

file.name <- paste0("boxplot0_sclc_tpm.gene.r5p47_rfd_TZ_E_HO+CD")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.e.ho, tpm.gene.log2.m.rfd.ctr.tz.e.cd, main="SCLC early TZ genes", names=c("HO", "CD"), cols=c(blue, blue), ylim)
#file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_rfd_TZ_E_HO+CD_long")
#plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.e.ho, tpm.gene.log2.m.rfd.ctr.tz.e.cd, main="SCLC early TZ genes", names=c("HO", "CD"), cols=c(blue, blue), ylim)

##
file.name <- paste0("boxplot0_sclc_tpm.gene.r5p47_rfd_TTR_E+L")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr.l, tpm.gene.log2.m.rfd.ttr.e, main="SCLC TTR genes", names=c("Late", "Early"), cols=c("black", "black"), ylim)
#file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_rfd_TTR_E+L_long")
#plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr.l, tpm.gene.log2.m.rfd.ttr.e, main="SCLC TTR genes", names=c("Late", "Early"), cols=c("black", "black"), ylim)

file.name <- paste0("boxplot0_sclc_tpm.gene.r5p47_rfd_TTR_E_HO+CD")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr.e.ho, tpm.gene.log2.m.rfd.ttr.e.cd, main="SCLC early TTR genes", names=c("HO", "CD"), cols=c("black", "black"), ylim)
#file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_rfd_TTR_E_HO+CD_long")
#plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr.e.ho, tpm.gene.log2.m.rfd.ttr.e.cd, main="SCLC early TTR genes", names=c("HO", "CD"), cols=c("black", "black"), ylim)

# -----------------------------------------------------------------------------
# GENE_NRFD plotBoxTSSNRFD (median0)
# Last Modified: 31/08/20
# -----------------------------------------------------------------------------
ylim <- c(-0.01, 0.033)
plotBoxNRFD(base, BASE, ylim, ylim2=c(-36, 35), tpm.gene.log2.m.rfd.ctr.iz.e.cd, tpm.gene.log2.m.rfd.ctr.iz.e.ho, tpm.gene.log2.m.rfd.ctr.tz.e.cd, tpm.gene.log2.m.rfd.ctr.tz.e.ho)

# -----------------------------------------------------------------------------
# RFD vs. TPM (median0)
# Last Modified: 11/02/21; 27/11/20; 01/09/20; 29/05/20
# -----------------------------------------------------------------------------
file.name <- file.path(wd.de.data, "tpm_gene_log2_m_rfd.RData")
load(file.name)

ylim <- c(min(tpm.gene.log2.m.rfd$MEDIAN), 14.5)

file.name <- paste0("boxplot3_sclc_tpm.gene.median0_RFD_3.5_FigureS1")
plotBox3(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr, tpm.gene.log2.m.rfd.ctr.tz, tpm.gene.log2.m.rfd.ctr.iz, main="SCLC expression", names=c("TTR", "TZ", "IZ"), cols=c("black", blue, red), ylim)

file.name <- paste0("boxplot4_sclc_tpm.gene.median0_RFD_3.2")
plotBox4(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.l, tpm.gene.log2.m.rfd.ctr.tz.e, tpm.gene.log2.m.rfd.ctr.iz.l, tpm.gene.log2.m.rfd.ctr.iz.e, main="SCLC expression (CTR)", names=c("L", "E", "L", "E"), cols=c(blue, blue, red, red), ylim)

file.name <- paste0("boxplot6_sclc_tpm.gene.median0_RFD_3.5")
plotBox6(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr.l, tpm.gene.log2.m.rfd.ttr.e, tpm.gene.log2.m.rfd.ctr.tz.l, tpm.gene.log2.m.rfd.ctr.tz.e, tpm.gene.log2.m.rfd.ctr.iz.l, tpm.gene.log2.m.rfd.ctr.iz.e, main="SCLC expression", names=c("L", "E", "L", "E", "L", "E"), cols=c("black", "black", blue, blue, red, red), ylim)








## TTR
main.txt <- paste0(BASE, " early TTR genes")
cols=c("white", "white")
ylim <- c(-1E-16, 1E-16)

file.name <- paste0("boxplot0_", base, "_tpm.gene_r5p47_TSS-NRFD_TTR_E_HO+CD")
plotBoxTSSNRFD(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr.e.ho, tpm.gene.log2.m.rfd.ttr.e.cd, main=main.txt, names=c("HO", "CD"), "TTR efficiency [TSS]", cols, ylim, isFlip=T, height=5, width=3.2)
file.name <- paste0("boxplot1_", base, "_tpm.gene_r5p47_TSS-NRFD_TTR_E_HO+CD")
plotBoxTSSNRFD(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr.e.ho, tpm.gene.log2.m.rfd.ttr.e.cd, main="TSS", names=c("HO", "CD"), "TTR efficiency", cols, ylim, isFlip=T, height=5, width=3.2)

#file.name <- paste0("boxplot_", base, "_tpm.gene_median0_GENE-NRFD_TTR_E_HO+CD")
#plotBoxGENENRFD(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr.e.ho, tpm.gene.log2.m.rfd.ttr.e.cd, main=main.txt, names=c("HO", "CD"), "TTR efficiency [Gene body]", cols, ylim, isFlip=T)
#file.name <- paste0("boxplot0_", base, "_tpm.gene_median0_GENE-NRFD_TTR_E_HO+CD")
#plotBoxGENENRFD(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr.e.ho, tpm.gene.log2.m.rfd.ttr.e.cd, main=main.txt, names=c("HO", "CD"), "TTR efficiency [Gene body]", cols, ylim, isFlip=T, height=5, width=3.2)
file.name <- paste0("boxplot1_", base, "_tpm.gene_r5p47_GENE-NRFD_TTR_E_HO+CD")
plotBoxGENENRFD(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr.e.ho, tpm.gene.log2.m.rfd.ttr.e.cd, main="Gene body", names=c("HO", "CD"), "TTR efficiency", cols, ylim, isFlip=T, height=5, width=3.2)

# -----------------------------------------------------------------------------
# 
# Last Modified: 29/05/20
# -----------------------------------------------------------------------------
genes <- c("PIF1", "MARS", "KIF18B", "BRCA2", "RAD9A", "FOXH1", "TONSL", "AL807752.1", "TOR1AIP1", "TOR1A", "TOR1B", "TP53", "RB1")
genes <- c("GTF3C2", "SUPT7L", "RAD9A", "E2F3")
genes <- c("MARS", "KIF18B", "BRCA2")
genes <- c("BLM", "POLE")
genes <- c("PIF1", "MARS", "KIF18B", "BRCA2")
genes <- c("RP3-407E4.3", "BRD9")
genes <- c("RP11-730B22.1")
genes <- c("RAD9A", "PIF1", "AL049840.1", "KIF18B", "GTPBP3", "BRCA2")
#genes <- c("RSAD2", "IRF2", "IFI35", "IFI6")
genes <- c("BRD9")
for (g in 1:length(genes)) {
   id <- subset(ensGene, external_gene_name == genes[g])$ensembl_gene_id
   plotSRC(genes[g], as.numeric(tpm.gene.log2[id,]), samples$COR, 1, yellow, "bottomright")
}

##
plotNLSCLC <- function(id, gene) {
   median(as.numeric(tpm.gene.lung.log2[id,]))
   median(as.numeric(tpm.gene.log2[id,]))
   testU(as.numeric(tpm.gene.lung.log2[id,]), as.numeric(tpm.gene.log2[id,]))
   file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_", gene, "_")
   plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2[id,], tpm.gene.log2[id,], main=gene, names=c("Lung", "SCLC"))
}

genes <- c("GTF3C2", "SUPT7L", "RAD9A", "E2F3", "ALK")
genes <- c("MARS", "KIF18B", "BRCA2")
genes <- c("BLM", "POLE")
for (g in 1:length(genes)) {
   id <- subset(ensGene, external_gene_name == genes[g])$ensembl_gene_id
   plotNLSCLC(id, genes[g])
}

# -----------------------------------------------------------------------------
# SCLC early IZ genes (Eff vs. Length)
# Last Modified: 16/02/21
# -----------------------------------------------------------------------------
plotEffvsLength <- function(n3, snr3, file.name, main.text, xlab.text, ylab.text, col, col2, pos) {
   pdf(paste0(file.name, ".pdf"), height=6, width=6)
   plot(n3 ~ snr3, ylab="", xlab=xlab.text, main=main.text[1], pch=15, col=col2, lwd=0, cex=2, cex.axis=1.5, cex.lab=1.6, cex.main=1.7)
   lm.fit <- lm(n3 ~ snr3)
   abline(lm.fit, col=col, lwd=4)
 
   cor3 <- cor.test(n3, snr3, method="spearman", exact=F)
   #cor3 <- round0(cor3[[4]], digits=2)
   #legend(pos, paste("rho = ", cor3), text.col=col, pch=15, col=col2, pt.cex=2.5, cex=1.5, pt.lwd=0, text.font=1)
   legend("bottomright", c(paste0("rho = ", round0(cor3[[4]], digits=2)), paste0("p-value = ", scientific(cor3[[3]]))), text.col=cols, text.font=2, bty="n", cex=1.5)
 
   #axis(side=2, at=seq(-0.2, 0.2, by=0.2), labels=c(-0.2, 0, 0.2), cex.axis=1.5)
   #axis(side=2, at=seq(-0.3, 0.1, by=0.2), labels=c("", "", ""), cex.axis=1.5)
   mtext(ylab.text, side=2, line=2.75, cex=1.6)
   mtext(main.text[2], line=0.3, cex=1.6)
   dev.off()
}

##
file.name <- file.path(wd.rt.plots, "SCLC_r5p47_IZ_E-vs-Length")
main.text <- c(paste("IZ efficiency vs. Gene length"), "")
xlab.text <- "Gene length [log10]"
ylab.text <- "IZ efficiency"                                                                         ## "#619CFF", "#F8766D", "#00BA38"      "skyblue3", "lightcoral", "#59a523"
cols <- "black"
cols2 <- "lightgray"
plotEffvsLength(tpm.gene.log2.m.rfd.ctr.iz.e$GENE_NRFD, log10(tpm.gene.log2.m.rfd.ctr.iz.e$length), file.name, main.text, xlab.text, ylab.text, cols, cols2, "topright")

##
file.name <- file.path(wd.rt.plots, "SCLC_r5p47_TZ_E-vs-Length")
main.text <- c(paste("TZ efficiency vs. Gene length"), "")
xlab.text <- "Gene length [log10]"
ylab.text <- "TZ efficiency"                                                                         ## "#619CFF", "#F8766D", "#00BA38"      "skyblue3", "lightcoral", "#59a523"
cols <- "black"
cols2 <- "lightgray"
plotEffvsLength(tpm.gene.log2.m.rfd.ctr.tz.e$GENE_NRFD, log10(tpm.gene.log2.m.rfd.ctr.tz.e$length), file.name, main.text, xlab.text, ylab.text, cols, cols2, "topright")

##
plotBoxLength <- function(wd.de.plots, file.name, tpm.1, tpm.2, main.txt, names, ylab.txt, cols, height=6, width=3) {
   trait <- rep(0, nrow(tpm.1))
   trait <- c(trait, rep(1, nrow(tpm.2)))
   trait <- as.factor(trait)
 
   expr <- as.numeric(c(log10(tpm.1$length), log10(tpm.2$length)))
 
   pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=height, width=width)
   boxplot(expr ~ trait, xaxt="n",  ylab=ylab.txt, main=main.txt, col=cols, cex.axis=1.2, cex.lab=1.25, cex.main=1.3)
 
   p <- testU(tpm.1$length, tpm.2$length)
   text <- ""
   if (p < 1E-3)  text <- "*"
   if (p < 1E-6)  text <- "**"
   if (p < 1E-9) text <- "***"
   text(1.5, ylim[2], text, col="black", cex=2.5)

   axis(side=1, at=seq(1, 2, by=1), labels=names, font=2, cex.axis=1.25)
   #axis(side=1, at=2, labels="Total n=30,978", line=1.3, col=NA, cex.axis=1.2)
   axis(side=1, at=1, labels=paste0("n=", format(nrow(tpm.1), big.mark=",", scientific=F)), line=1.3, col=NA, cex.axis=1.2)
   axis(side=1, at=2, labels=paste0("n=", format(nrow(tpm.2), big.mark=",", scientific=F)), line=1.3, col=NA, cex.axis=1.2)
 
   mtext(paste0("p-value = ", scientific(wilcox.test(expr ~ trait, exact=F)$p.value)), cex=1.25, line=0.3)
   dev.off()
}

names <- c("HO", "CD")
file.name <- paste0("boxplot1_", base, "_tpm.gene_r5p47_LENGTH_IZ_E_HO+CD")
plotBoxLength(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.e.ho, tpm.gene.log2.m.rfd.ctr.iz.e.cd, main="Early IZ genes", names, "Gene length [log10]", "lightgray", height=5, width=3.2)

file.name <- paste0("boxplot1_", base, "_tpm.gene_r5p47_LENGTH_TZ_E_HO+CD")
plotBoxLength(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.e.ho, tpm.gene.log2.m.rfd.ctr.tz.e.cd, main="Early TZ genes", names, "Gene length [log10]", "lightgray", height=5, width=3.2)

# -----------------------------------------------------------------------------
# SCLC early IZ genes
# Last Modified: 31/08/20
# -----------------------------------------------------------------------------
## Expressed IZ-E, CD genes
overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.iz.e.cd))
de.tpm.gene.rfd.iz.e.cd <- cbind(de.tpm.gene[overlaps,], tpm.gene.log2.m.rfd.ctr.iz.e.cd[overlaps, c("GENE_NRFD", "length")])
de.tpm.gene.rfd.iz.e.cd$Q <- qvalue(de.tpm.gene.rfd.iz.e.cd$P)$qvalue
writeTable(de.tpm.gene.rfd.iz.e.cd, file.path(wd.de.data, "de_sclc_tpm-gene-median0_src_q_rfd_iz_e_cd_n70.txt"), colnames=T, rownames=F, sep="\t")

overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.iz.e.ho))
de.tpm.gene.rfd.iz.e.ho <- cbind(de.tpm.gene[overlaps,], tpm.gene.log2.m.rfd.ctr.iz.e.ho[overlaps, c("GENE_NRFD", "length")])
de.tpm.gene.rfd.iz.e.ho$Q <- qvalue(de.tpm.gene.rfd.iz.e.ho$P)$qvalue
writeTable(de.tpm.gene.rfd.iz.e.ho, file.path(wd.de.data, "de_sclc_tpm-gene-median0_src_q_rfd_iz_e_ho_n70.txt"), colnames=T, rownames=F, sep="\t")

## Expressed TZ-E, CD genes
overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.tz.e.cd))
de.tpm.gene.rfd.tz.e.cd <- cbind(de.tpm.gene[overlaps,], tpm.gene.log2.m.rfd.ctr.tz.e.cd[overlaps, c("GENE_NRFD", "length")])
de.tpm.gene.rfd.tz.e.cd$Q <- qvalue(de.tpm.gene.rfd.tz.e.cd$P)$qvalue
writeTable(de.tpm.gene.rfd.tz.e.cd, file.path(wd.de.data, "de_sclc_tpm-gene-median0_src_q_rfd_tz_e_cd_n70.txt"), colnames=T, rownames=F, sep="\t")

overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.tz.e.ho))
de.tpm.gene.rfd.tz.e.ho <- cbind(de.tpm.gene[overlaps,], tpm.gene.log2.m.rfd.ctr.tz.e.ho[overlaps, c("GENE_NRFD", "length")])
de.tpm.gene.rfd.tz.e.ho$Q <- qvalue(de.tpm.gene.rfd.tz.e.ho$P)$qvalue
writeTable(de.tpm.gene.rfd.tz.e.ho, file.path(wd.de.data, "de_sclc_tpm-gene-median0_src_q_rfd_tz_e_ho_n70.txt"), colnames=T, rownames=F, sep="\t")

###
## Expressed TTR-E genes
overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd.ttr.e.cd))
de.tpm.gene.rfd.ttr.e.cd <- cbind(de.tpm.gene[overlaps,], tpm.gene.log2.m.rfd.ttr.e.cd[overlaps, c("GENE_NRFD", "length")])
de.tpm.gene.rfd.ttr.e.cd$Q <- qvalue(de.tpm.gene.rfd.ttr.e.cd$P)$qvalue
writeTable(de.tpm.gene.rfd.ttr.e.cd, file.path(wd.de.data, "de_sclc_tpm-gene-median0_src_q_rfd_ttr_e_cd_n70.txt"), colnames=T, rownames=F, sep="\t")

overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd.ttr.e.ho))
de.tpm.gene.rfd.ttr.e.ho <- cbind(de.tpm.gene[overlaps,], tpm.gene.log2.m.rfd.ttr.e.ho[overlaps, c("GENE_NRFD", "length")])
de.tpm.gene.rfd.ttr.e.ho$Q <- qvalue(de.tpm.gene.rfd.ttr.e.ho$P)$qvalue
writeTable(de.tpm.gene.rfd.ttr.e.ho, file.path(wd.de.data, "de_sclc_tpm-gene-median0_src_q_rfd_ttr_e_ho_n70.txt"), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Volcano plots
# Last Modified: 29/05/20
# -----------------------------------------------------------------------------
fdrToP <- function(fdr, de) {
   de.sig <- subset(de, Q <= fdr)
   de.sig$log10P <- -log10(de.sig$P)
 
   return(max(de.sig$P))
}

pvalueToFDR <- function(pvalue, de) {
   de.sig <- subset(de, P <= pvalue)
 
   return(round0(max(de.sig$Q)*100, digits=0))
}

plotVolcano <- function(de, pvalue, genes, file.de, file.main, xlab.text, ymax=0, cols, fc) {
   #pvalue <- fdrToP(fdr, de)
   #fdr <- pvalueToFDR(pvalue, de)
   de.sig <- subset(de, P <= pvalue)
   de.sig$log10P <- -log10(de.sig$P)
 
   de$log10P <- -log10(de$P)
   xlim <- c(min(de$LOG2_FC), max(de$LOG2_FC))
   if (ymax ==0) ymax <- max(de$log10P)
   
   pdf(file.de, height=6, width=6)
   plot(de$LOG2_FC, de$log10P, pch=16, xlim=xlim, ylim=c(0, ymax), xaxt="n", xlab=xlab.text, ylab="P-value significance [-log10]", col="gray88", main=file.main[1], cex=1.4, cex.axis=1.2, cex.lab=1.25, cex.main=1.3)

   #text(xmax*-1 + 2*xmax/15, -log10(pvalue) - ymax/30, paste0("FDR=", fdr, "%"), cex=1.1)    ## SCLC (IZ)
   #text(xmax*-1 + 2*xmax/13, -log10(pvalue) - ymax/30, paste0("FDR=", fdr*100, "%"), cex=1.1)   ## SCLC (AA)
   #text(xmax*-1 + 2*xmax/9.5, -log10(pvalue) - ymax/30, paste0("BH=1.00E-16"), cex=1.1)
   
   de.up   <- subset(de.sig, LOG2_FC > log2(fc))
   points(de.up$LOG2_FC, de.up$log10P, pch=16, col=cols[1], cex=1.4)
   de.down <- subset(de.sig, LOG2_FC < -log2(fc))
   points(de.down$LOG2_FC, de.down$log10P, pch=16, col=cols[2], cex=1.4)

   abline(v=log2(fc), lty=5, lwd=1.8, col=cols[1]) 
   abline(v=-log2(fc), lty=5, lwd=1.8, col=cols[2])
   #abline(v=0, col="gray")
   abline(h=c(-log10(pvalue)), lty=5, lwd=1.5)
   
   if (nrow(genes) != 0) {
      for (g in 1:nrow(genes)) {
         gene <- subset(de, external_gene_name == genes[g,]$GENE)
         gene <- cbind(gene, genes[g,])
      
         if (nrow(gene) > 0) {
            points(gene$LOG2_FC, gene$log10P, pch=1, col="black", cex=1.4)
         
            if (!is.na(gene$ADJ_1))
               if (is.na(gene$ADJ_2))
                  text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=gene$ADJ_1, cex=1.25)
               else
                  text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=c(gene$ADJ_1, gene$ADJ_2), cex=1.25)
            else
               if (gene$LOG2_FC > 0)
                  text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=c(0, -0.6), cex=1.25)
               else
                  text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=c(1, -0.6), cex=1.25)
         } else
            print(genes[g])
      }
   }
   
   #axis(side=1, at=seq(-3, 6, by=1), labels=c(-3, -2, -1, 0, 1, 2, 3, 4, 5, 6), cex.axis=1.2)
   axis(side=1, at=seq(-12, 12, by=4), labels=c(-12, -8, -4, 0, 4, 8, 12), cex.axis=1.2)
   axis(side=1, at=seq(-10, 10, by=4), labels=c(-10, -6, -2, 2, 6, 10), cex.axis=1.2)
   mtext(file.main[2], cex=1.25, line=0.3)
   legend("topleft", legend=c("Positively", "Negatively"), col=cols, pch=19, cex=1.25)
   dev.off()
}

###
## SCLC r5p47 genes
xlab.text <- "S to G1 fold change [log2]"
file.main <- c("SCLC expressed genes (n=19,131)", "Expression vs. In-silico sorting")
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_r5p47_rfd_p1e-3_all_PIF1")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene, 0.001, genes, file.de, file.main, xlab.text, ymax=5, c(yellow, lightblue), 1.3)

plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_r5p47_rfd_p1e-3_all_BRD9")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene, 0.001, genes, file.de, file.main, xlab.text, ymax=5, c(yellow, lightblue), 1.3)

plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_r5p47_rfd_p1e-3_all_TONSL")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene, 0.001, genes, file.de, file.main, xlab.text, ymax=5, c(yellow, lightblue), 1.3)

plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_r5p47_George_RFD")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene, 0.001, genes, file.de, file.main, xlab.text, ymax=5, c(yellow, lightblue), 1.3)







## SCLC IZ-E, CD genes
xlab.text <- "S to G1 fold change [log2]"
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_r5p47_rfd_p1e-3_iz_e_cd")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c("SCLC expression vs. In-silico sorting", "Early IZ, CD genes (n=1,159)")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene.rfd.iz.e.cd, 0.001, genes, file.de, file.main, xlab.text, ymax=4.3, c(yellow, lightblue), 1.3)







###
## SCLC ALL genes
xlab.text <- "SCLC S/G1 fold change [log2]"
file.main <- c("SCLC expressed genes (n=23,572)", "Expression vs. In-silico sorting")
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_median0_rfd_p1e-3_all_PIF1")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene, 0.001, genes, file.de, file.main, xlab.text, ymax=5.5)

xlab.text <- "SCLC S/G1 fold change [log2]"
file.main <- c("SCLC expressed genes (n=23,572)", "Expression vs. In-silico sorting")
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_median0_rfd_p1e-3_all_Helicases")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene, 0.001, genes, file.de, file.main, xlab.text, ymax=5.5)

xlab.text <- "SCLC S/G1 fold change [log2]"
file.main <- c("SCLC expressed genes (n=23,572)", "Expression vs. In-silico sorting")
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_median0_rfd_p1e-3_all_TFBS")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene, 0.001, genes, file.de, file.main, xlab.text, ymax=5.5)

xlab.text <- "SCLC S/G1 fold change [log2]"
file.main <- c("SCLC expressed genes (n=23,572)", "Expression vs. In-silico sorting")
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_median0_rfd_p1e-3_all_BRD9")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene, 0.001, genes, file.de, file.main, xlab.text, ymax=5.5)




xlab.text <- "SCLC S/G1 [log2 fold change]"
file.main <- c("SCLC expressed genes (n=23,572)", "Correlation between expression and in-silico sorting")
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_median0_rfd_p1e-2_all")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene, 0.01, genes, file.de, file.main, xlab.text, ymax=5.5)

plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_median0_rfd_p1e-2_all_up")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene, 0.01, genes, file.de, file.main, xlab.text, ymax=5.5)

plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_median0_rfd_p1e-2_all_down")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene, 0.01, genes, file.de, file.main, xlab.text, ymax=5.5)


## SCLC IZ-E, CD genes
xlab.text <- "SCLC S/G1 fold change [log2]"
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_median0_rfd_p1e-3_iz_e_cd")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c("SCLC expression vs. In-silico sorting", "Early IZ, CD genes (n=1,377)")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene.rfd.iz.e.cd, 0.001, genes, file.de, file.main, xlab.text, ymax=4.3)






## SCLC IZ-E, CD genes
xlab.text <- "SCLC S/G1 [log2 fold change]"
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_median0_rfd_p1e-2_iz_e_cd")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c("SCLC early IZ, CD genes (n=1,377)", "Correlation between expression and in-silico sorting")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene.rfd.iz.e.cd, 0.01, genes, file.de, file.main, xlab.text, ymax=4.3)

## SCLC TZ-E, HO genes
xlab.text <- "SCLC S/G1 [log2 fold change]"
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_median0_rfd_p1e-2_tz_e_ho")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c("SCLC early TZ, HO genes (n=759)", "Correlation between expression and in-silico sorting")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene.rfd.tz.e.ho, 0.01, genes, file.de, file.main, xlab.text, ymax=4.3)

###
## SCLC IZ-E, HO genes
xlab.text <- "SCLC S/G1 [log2 fold change]"
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_median0_rfd_p1e-2_iz_e_ho")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c("SCLC early IZ, HO genes (n=1,372)", "Correlation between expression and in-silico sorting")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene.rfd.iz.e.ho, 0.01, genes, file.de, file.main, xlab.text, ymax=4.3)

## SCLC TZ-E, CD genes
xlab.text <- "SCLC S/G1 [log2 fold change]"
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_median0_rfd_p1e-2_tz_e_cd")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c("SCLC early TZ, CD genes (n=791)", "Correlation between expression and in-silico sorting")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene.rfd.tz.e.cd, 0.01, genes, file.de, file.main, xlab.text, ymax=4.3)

###
## SCLC TTR-E, HO genes
xlab.text <- "SCLC S/G1 [log2 fold change]"
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_median0_rfd_p1e-2_ttr_e_ho")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c("SCLC early TTR, HO genes (n=6,936)", "Correlation between expression and in-silico sorting")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene.rfd.ttr.e.ho, 0.01, genes, file.de, file.main, xlab.text, ymax=4.3)










# -----------------------------------------------------------------------------
# Replace Ensembl Gene IDs to gene name in Reactome results (Up- and Down-regulation)
# -----------------------------------------------------------------------------
wd.de.pathway <- file.path(wd.de, "pathway_median0")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_IZ-E-CD_p1e-2_up")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_IZ-E-CD_p1e-2_down")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_TZ-E-HO_p1e-2_up")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_TZ-E-HO_p1e-2_down")

setTable <- function(wd.de.reactome) {
   list <- ensGene[,c("ensembl_gene_id",	"external_gene_name")]

   reactome <- read.csv(file.path(wd.de.reactome, "result.csv"))
   colnames(reactome) <- gsub("X.", "", colnames(reactome))
   reactome$Submitted.entities.found <- as.vector(reactome$Submitted.entities.found)
   for (r in 1:nrow(reactome)) {
      ids <- as.vector(reactome$Submitted.entities.found[r])
      ids <- unlist(strsplit(ids, ";"))
 
      for (i in 1:length(ids))
         if (nrow(list[ids[i],]) != 0)
            ids[i] <- list[ids[i],]$external_gene_name
 
      reactome$Submitted.entities.found[r] <- paste(ids, collapse=";")
   }
   writeTable(reactome, file.path(wd.de.reactome, "result.tsv"), colnames=T, rownames=F, sep="\t")
}

# -----------------------------------------------------------------------------
# r5p47 genes (21/02/21)
# -----------------------------------------------------------------------------
wd.de.pathway <- file.path(wd.de, "pathway_r5p47_p1e-3_fc1.3")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_ALL_p1e-3_up")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_ALL_p1e-3_down")

###
## Link: https://www.statmethods.net/graphs/bar.html
wd.de.reactome <- file.path(wd.de.pathway, "reactome_ALL_p1e-3_up")
setTable(wd.de.reactome)
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")
reactome[2, 2] <- "HDR through HRR or SSA"

reactome.up <- subset(reactome, Entities.pValue <= 1e-2)[, 2:7]
reactome.up$log10P <- -log10(reactome.up$Entities.pValue)
reactome.up <- reactome.up[order(reactome.up$log10P),]

main.text <- c("Positively-correlated", "41 genes")
file.de <- file.path(wd.de.reactome, "genes_ALL_p1e-2_n41_up.pdf")

pdf(file.de, height=4.2, width=7.5)
par(mar=c(4,27.5,4,1))   # increase y-axis margin.
barplot(reactome.up$log10P, main=main.text[1], las=1, horiz=T, xlim=c(0, 4), xaxt="n", names.arg=reactome.up$Pathway.name, col=yellow, xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
abline(v=2, lty=5)

axis(side=1, at=seq(0, 4, by=1))
mtext(main.text[2], line=0.3)
dev.off()

###
## Link: https://www.statmethods.net/graphs/bar.html
wd.de.reactome <- file.path(wd.de.pathway, "reactome_ALL_p1e-3_down")
setTable(wd.de.reactome)
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")

reactome.down <- subset(reactome, Entities.pValue <= 1e-2)[, 2:7]
reactome.down$log10P <- -log10(reactome.down$Entities.pValue)
reactome.down <- reactome.down[order(reactome.down$log10P),]
reactome.down$log10P <- reactome.down$log10P * -1

main.text <- c("Negatively-correlated", "17 genes")
file.de <- file.path(wd.de.reactome, "genes_ALL_p1e-3_n17_down.pdf")

pdf(file.de, height=2.9, width=7.5)
par(mar=c(4,1,4,27.5))   # increase y-axis margin.
posbar <- barplot(reactome.down$log10P, main=main.text[1], las=1, horiz=T, xlim=c(-4, 0), xaxt="n", names.arg="", col="steelblue1", xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=reactome.down$Pathway.name, xpd=NA)
abline(v=-2, lty=5)

axis(side=1, at=seq(-4, 0, by=1))
mtext(main.text[2], line=0.3)
dev.off()





# -----------------------------------------------------------------------------
# Heatmap
# Links: https://davetang.org/muse/2010/12/06/making-a-heatmap-with-r/
# Last Modified: 21/04/21; 13/01/20
# -----------------------------------------------------------------------------
#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq")
#install.packages("gplots")

library("DESeq")
library(gplots)

genes <- subset(subset(de.tpm.gene, P < 1e-3), LOG2_FC >= log2(1.4))
genes.ensembl  <- genes$ensembl_gene_id
genes.rownames <- genes$external_gene_name
 
b <- tpm.gene.log2[genes.ensembl,]
colnames(b) <- samples$SAMPLE_ID
rownames(b) <- genes.rownames
b <- data.matrix(b)

heatmap(b)

## Returning the values used for the heatmap
test <- heatmap.2(b, scale="row")
b[rev(test$rowInd), test$colInd]

##
hr <- hclust(as.dist(1-cor(t(b), method="pearson")),  method="complete")
hc <- hclust(as.dist(1-cor(b,    method="spearman")), method="complete")

###
##
colfunc <- colorRampPalette(c("green", "black", "red"))
heatmap.2(b, col=colfunc(7), scale="row", trace="none")














# -----------------------------------------------------------------------------
# TTR
# -----------------------------------------------------------------------------
wd.de.pathway <- file.path(wd.de, "pathway")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_TTR-E-HO_p1e-2_up")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_TTR-E-HO_p1e-2_down")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_TTR-E-CD_p1e-2_up")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_TTR-E-CD_p1e-2_down")

###
## Link: https://www.statmethods.net/graphs/bar.html
wd.de.reactome <- file.path(wd.de.pathway, "reactome_TTR-E-HO_p1e-2_up")
setTable(wd.de.reactome)
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")
reactome[2, 2] <- "Resolution of D-loop Structures through Holliday Junction Intermediate"
reactome[5, 2] <- "HDR through HRR or Single Strand Annealing (SSA)"

reactome.up <- subset(reactome, Entities.pValue <= 1e-2)[, 2:7]
reactome.up$log10P <- -log10(reactome.up$Entities.pValue)
reactome.up <- reactome.up[order(reactome.up$log10P),]

main.text <- c("Positively-correlated", "130 genes")
file.de <- file.path(wd.de.reactome, "genes_TTR-E-HO_p1e-2_n130_up.pdf")

pdf(file.de, height=3.6, width=7.5)
par(mar=c(4,27.5,4,1))   # increase y-axis margin.
barplot(reactome.up$log10P, main=main.text[1], las=1, horiz=T, xlim=c(0, 4), xaxt="n", names.arg=reactome.up$Pathway.name, col="gold", xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
abline(v=3, lty=5)

axis(side=1, at=seq(0, 4, by=1))
mtext(main.text[2], line=0.3)
dev.off()

###
## Link: https://www.statmethods.net/graphs/bar.html
wd.de.reactome <- file.path(wd.de.pathway, "reactome_TTR-E-HO_p1e-2_down")
setTable(wd.de.reactome)
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")

reactome.down <- subset(reactome, Entities.pValue <= 1e-2)[, 2:7]
reactome.down$log10P <- -log10(reactome.down$Entities.pValue)
reactome.down <- reactome.down[order(reactome.down$log10P),]
reactome.down$log10P <- reactome.down$log10P * -1

main.text <- c("Negatively-correlated", "72 genes")
file.de <- file.path(wd.de.reactome, "genes_TTR-E-HO_p1e-2_n72_down.pdf")

pdf(file.de, height=2.2, width=7.5)
par(mar=c(4,1,4,27.5))   # increase y-axis margin.
posbar <- barplot(reactome.down$log10P, main=main.text[1], las=1, horiz=T, xlim=c(-4, 0), xaxt="n", names.arg="", col="steelblue1", xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=reactome.down$Pathway.name, xpd=NA)
abline(v=-3, lty=5)

axis(side=1, at=seq(-4, 0, by=1))
mtext(main.text[2], line=0.3)
dev.off()

###
## Link: https://www.statmethods.net/graphs/bar.html
wd.de.reactome <- file.path(wd.de.pathway, "reactome_TZ-E-CD_p1e-2_up")
setTable(wd.de.reactome)
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")

reactome.up <- subset(reactome, Entities.pValue <= 2e-2)[, 2:7]
reactome.up$log10P <- -log10(reactome.up$Entities.pValue)
reactome.up <- reactome.up[order(reactome.up$log10P),]

main.text <- c("Positively-correlated", "6 genes")
file.de <- file.path(wd.de.reactome, "genes_TZ-E-CD_p1e-2_n6_up.pdf")

pdf(file.de, height=2.2, width=7.5)
par(mar=c(4,27.5,4,1))   # increase y-axis margin.
barplot(reactome.up$log10P, main=main.text[1], las=1, horiz=T, xlim=c(0, 4), xaxt="n", names.arg=reactome.up$Pathway.name, col="gold", xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
abline(v=2, lty=5)

axis(side=1, at=seq(0, 4, by=1))
mtext(main.text[2], line=0.3)
dev.off()

##
wd.de.reactome <- file.path(wd.de.pathway, "reactome_TZ-E-CD_p1e-2_down")
setTable(wd.de.reactome)
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")

reactome.down <- subset(reactome, Entities.pValue <= 1e-2)[, 2:7]
reactome.down$log10P <- -log10(reactome.down$Entities.pValue)
reactome.down <- reactome.down[order(reactome.down$log10P),]
reactome.down$log10P <- reactome.down$log10P * -1

main.text <- c("Negatively-correlated", "6 genes")
file.de <- file.path(wd.de.reactome, "genes_TS-E-CD_p1e-2_n6_down.pdf")

pdf(file.de, height=3.2, width=7.5)
par(mar=c(4,1,4,27.5))   # increase y-axis margin.
posbar <- barplot(reactome.down$log10P, main=main.text[1], las=1, horiz=T, xlim=c(-4, 0), xaxt="n", names.arg="", col="steelblue1", xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=reactome.down$Pathway.name, xpd=NA)

axis(side=1, at=seq(-4, 0, by=1))
mtext(main.text[2], line=0.3)
dev.off()


















## All
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_r5p47_p3e04")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c("SCLC expressed genes", paste0("n=18,004"))
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene, 0.0003, genes, file.de, file.main, xlab.text, ymax=6.8)

##
xlab.text <- "SCLC M2/M1 [log2FC]"
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_r5p47_rfd_fdr0.1")

## IZ-S
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_r5p47_rfd_p3e04_iz_sp_cd")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c("SCLC expressed, IZ-S, CD genes", paste0("n=642"))
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene.rfd.iz.sp.cd, 0.0003, genes, file.de, file.main, xlab.text, ymax=4.5)

plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_r5p47_rfd_fdr0.08_iz_sp_ho")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c("SCLC IZ-S, HO genes", paste0("n=657"))
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene.rfd.iz.sp.ho, 0.08, genes, file.de, file.main, xlab.text)

## TOR1
de.tpm.gene.rfd.ctr.sp <- rbind(de.tpm.gene.rfd.iz.sp, de.tpm.gene.rfd.tz.sp)
de.tpm.gene.rfd.ctr.sp$Q <- qvalue(de.tpm.gene.rfd.ctr.sp$P)$qvalue

plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_r5p47_rfd_fdr0.08_ctr_sp")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c("SCLC CTR-S genes", paste0("n=2,173"))
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene.rfd.ctr.sp, 0.08, genes, file.de, file.main, xlab.text)






















## IZ (CD vs. HO)
# > median(de.tpm.gene.rfd.iz.sp.cd$MEDIAN)
# [1] 3.47189
# > median(de.tpm.gene.rfd.iz.sp.ho$MEDIAN)
# [1] 3.397869
# > testU(de.tpm.gene.rfd.iz.sp.cd$MEDIAN, de.tpm.gene.rfd.iz.sp.ho$MEDIAN)
# [1] 0.5321021
file.name <- paste0("boxplot_sclc_tpm.gene.rfd_IZ_SP_HO-vs-CD")
plotBox(wd.de.plots, file.name, de.tpm.gene.rfd.iz.sp.ho, de.tpm.gene.rfd.iz.sp.cd, main="IZ-S", names=c("HO", "CD"), cols=c("red", "red"), ylim)

## TZ-SP (HO vs CD)  ## 2020/05/17
de.tpm.gene.rfd.tz.sp <- readTable(file.path(wd.de.data, "de_sclc_tpm-gene-r5p47-rfd_src_q_tz_sp_n70.txt"), header=T, rownames=T, sep="\t")

overlaps <- intersect(rownames(de.tpm.gene.rfd.tz.sp), rownames(subset(tpm.gene.log2.m.rfd.ctr.tz, TRC == 1)))
de.tpm.gene.rfd.tz.sp.cd <- cbind(de.tpm.gene.rfd.tz.sp[overlaps,], tpm.gene.log2.m.rfd.ctr.tz[overlaps,])
de.tpm.gene.rfd.tz.sp.cd$Q <- qvalue(de.tpm.gene.rfd.tz.sp.cd$P)$qvalue
writeTable(de.tpm.gene.rfd.tz.sp.cd, file.path(wd.de.data, "de_sclc_tpm-gene-r5p47-rfd_src_q_tz_sp_cd_n70.txt"), colnames=T, rownames=F, sep="\t")

overlaps <- intersect(rownames(de.tpm.gene.rfd.tz.sp), rownames(subset(tpm.gene.log2.m.rfd.ctr.tz, TRC == -1)))
de.tpm.gene.rfd.tz.sp.ho <- cbind(de.tpm.gene.rfd.tz.sp[overlaps,], tpm.gene.log2.m.rfd.ctr.tz[overlaps,])
de.tpm.gene.rfd.tz.sp.ho$Q <- qvalue(de.tpm.gene.rfd.tz.sp.ho$P)$qvalue
writeTable(de.tpm.gene.rfd.tz.sp.ho, file.path(wd.de.data, "de_sclc_tpm-gene-r5p47-rfd_src_q_tz_sp_ho_n70.txt"), colnames=T, rownames=F, sep="\t")

## TZ (CD vs. HO)
# > median(de.tpm.gene.rfd.tz.sp.cd$MEDIAN)
# [1] 3.125119
# > median(de.tpm.gene.rfd.tz.sp.ho$MEDIAN)
# [1] 3.504559
# > testU(de.tpm.gene.rfd.tz.sp.cd$MEDIAN, de.tpm.gene.rfd.tz.sp.ho$MEDIAN)
# [1] 0.0518692
file.name <- paste0("boxplot_sclc_tpm.gene.rfd_TZ_SP_HO-vs-CD")
plotBox(wd.de.plots, file.name, de.tpm.gene.rfd.tz.sp.ho, de.tpm.gene.rfd.tz.sp.cd, main="TZ-S", names=c("HO", "CD"), cols=c("blue", "blue"), ylim)





###
##
overlaps <- intersect(rownames(tpm.gene.log2.m.rfd.ctr.tz), rownames(tpm.gene.log2.m.rfd.ctr.tz.nl))
tpm.gene.log2.m.rfd.ctr.tz.sp <- tpm.gene.log2.m.rfd.ctr.tz[setdiff(rownames(tpm.gene.log2.m.rfd.ctr.tz), overlaps),]

tpm.gene.log2.m.rfd.ctr.tz.sp <- tpm.gene.log2.m.rfd.ctr.tz.sp[intersect(rownames(tpm.gene.log2.m.rfd.ctr.tz.sp), rownames(de.tpm.gene)),]

## IZ
overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.tz.sp))
de.tpm.gene.rfd.tz.sp <- de.tpm.gene[overlaps,]
de.tpm.gene.rfd.tz.sp$Q <- qvalue(de.tpm.gene.rfd.tz.sp$P)$qvalue
writeTable(de.tpm.gene.rfd.tz.sp, file.path(wd.de.data, "de_sclc_tpm-gene-r5p47-rfd_src_q_tz_sp_n70.txt"), colnames=T, rownames=F, sep="\t")

## IZ (Overlapping SCLC and NBL)
overlaps <- intersect(rownames(tpm.gene.log2.m.rfd.ctr.iz.sp), rownames(de.tpm.gene.rfd.iz.s))









###
##
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47.rfd_IZ-NSP-vs-IZ-SP_3")
plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.nsp, tpm.gene.log2.m.rfd.ctr.iz.sp, main="SCLC-specific IZ (IZ-S)", names=c("IZ-NS", "IZ-S"), cols=c("red", "red"), ylim)

###
## IZ vs TZ
# > median(tpm.gene.log2.m.rfd.ctr.iz$MEDIAN)
# [1] 3.455006
# > median(tpm.gene.log2.m.rfd.ctr.tz$MEDIAN)
# [1] 3.236536
# > testU(tpm.gene.log2.m.rfd.ctr.iz$MEDIAN, tpm.gene.log2.m.rfd.ctr.tz$MEDIAN)
# [1] 0.0180636
file.name <- paste0("boxplot_sclc_tpm.gene.rfd_IZ-vs-TZ")
plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz, tpm.gene.log2.m.rfd.ctr.tz, main="CTR", names=c("IZ", "TZ"), cols=c("red", "blue"), ylim)

## IZ (CD vs. HO)
# > median(tpm.gene.log2.m.rfd.ctr.iz.cd$MEDIAN)
# [1] 3.44142
# > median(tpm.gene.log2.m.rfd.ctr.iz.ho$MEDIAN)
# [1] 3.470568
# > testU(tpm.gene.log2.m.rfd.ctr.iz.cd$MEDIAN, tpm.gene.log2.m.rfd.ctr.iz.ho$MEDIAN)
# [1] 0.6792462
file.name <- paste0("boxplot_sclc_tpm.gene.rfd_IZ_HO-vs-CD")
plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.ho, tpm.gene.log2.m.rfd.ctr.iz.cd, main="IZ", names=c("HO", "CD"), cols=c("red", "red"), ylim)

## TZ (CD vs. HO)
# > median(tpm.gene.log2.m.rfd.ctr.tz.cd$MEDIAN)
# [1] 0.6792462
# > median(tpm.gene.log2.m.rfd.ctr.tz.ho$MEDIAN)
# [1] 2.994693
# > testU(tpm.gene.log2.m.rfd.ctr.tz.cd$MEDIAN, tpm.gene.log2.m.rfd.ctr.tz.ho$MEDIAN)
# [1] 0.008081627
file.name <- paste0("boxplot_sclc_tpm.gene.rfd_TZ_HO-vs-CD")
plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.ho, tpm.gene.log2.m.rfd.ctr.tz.cd, main="TZ", names=c("HO", "CD"), cols=c("blue", "blue"), ylim)











# -----------------------------------------------------------------------------
# Gene length vs. RFD slopes
# Last Modified: 29/05/20
# -----------------------------------------------------------------------------
tpm.gene.log2.m.rfd.ctr.iz <- cbind(de.tpm.gene.rfd.iz.sp, tpm.gene.log2.m.rfd.ctr.iz[rownames(de.tpm.gene.rfd.iz.sp),])
tpm.gene.log2.m.rfd.ctr.iz$length <- abs(tpm.gene.log2.m.rfd.ctr.iz$end_position - tpm.gene.log2.m.rfd.ctr.iz$start_position) / 1000
tpm.gene.log2.m.rfd.ctr.iz.cd <- subset(tpm.gene.log2.m.rfd.ctr.iz, TRC > 0)
tpm.gene.log2.m.rfd.ctr.iz.ho <- subset(tpm.gene.log2.m.rfd.ctr.iz, TRC < 0)

##
xlim <- c(min(log10(tpm.gene.log2.m.rfd.ctr.iz$length)), max(log10(tpm.gene.log2.m.rfd.ctr.iz$length)))
ylim <- c(min(tpm.gene.log2.m.rfd.ctr.iz$GENE_NRFD), min(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD) * -1)
ylab.text <- "IZ efficiency"

file.name <- file.path(wd.de.plots, "SCLC-SP_TRC_LENGTH-vs-GENE-NRFD_IZ")
plotTRC(tpm.gene.log2.m.rfd.ctr.iz$GENE_NRFD, log10(tpm.gene.log2.m.rfd.ctr.iz$length), c("SCLC, IZ-S genes", "n=1,299"), file.name, xlim, ylim, col=c("red", "black"), "topright", ylab.text)

file.name <- file.path(wd.de.plots, "SCLC-SP_TRC_LENGTH-vs-GENE-NRFD_IZ_CD")
plotTRC(tpm.gene.log2.m.rfd.ctr.iz.cd$GENE_NRFD, log10(tpm.gene.log2.m.rfd.ctr.iz.cd$length), c("SCLC, IZ-S and CD genes", "n=642"), file.name, xlim, ylim, col=c("red", "black"), "topright", ylab.text)

file.name <- file.path(wd.de.plots, "SCLC-SP_TRC_LENGTH-vs-GENE-NRFD_IZ_HO")
plotTRC(tpm.gene.log2.m.rfd.ctr.iz.ho$GENE_NRFD, log10(tpm.gene.log2.m.rfd.ctr.iz.ho$length), c("SCLC, IZ-S and HO genes", "n=657"), file.name, xlim, ylim, col=c("red", "black"), "topright", ylab.text)

##
tpm.gene.log2.m.rfd.ctr.tz <- cbind(de.tpm.gene.rfd.tz.sp, tpm.gene.log2.m.rfd.ctr.tz[rownames(de.tpm.gene.rfd.tz.sp),])
tpm.gene.log2.m.rfd.ctr.tz$length <- abs(tpm.gene.log2.m.rfd.ctr.tz$end_position - tpm.gene.log2.m.rfd.ctr.tz$start_position) / 1000
tpm.gene.log2.m.rfd.ctr.tz.cd <- subset(tpm.gene.log2.m.rfd.ctr.tz, TRC > 0)
tpm.gene.log2.m.rfd.ctr.tz.ho <- subset(tpm.gene.log2.m.rfd.ctr.tz, TRC < 0)

##
xlim <- c(min(log10(tpm.gene.log2.m.rfd.ctr.iz$length)), max(log10(tpm.gene.log2.m.rfd.ctr.tz$length)))
ylim <- c(max(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD), min(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD)) * -1
ylab.text <- "TZ efficiency"

file.name <- file.path(wd.de.plots, "SCLC-SP_TRC_LENGTH-vs-GENE-NRFD_TZ")
plotTRC(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD * -1, log10(tpm.gene.log2.m.rfd.ctr.tz$length), c("SCLC, TZ-S genes", "n=874"), file.name, xlim, ylim, col=c("blue", "black"), "topright", ylab.text, isFlip=T)

file.name <- file.path(wd.de.plots, "SCLC-SP_TRC_LENGTH-vs-GENE-NRFD_TZ_CD")
plotTRC(tpm.gene.log2.m.rfd.ctr.tz.cd$GENE_NRFD * -1, log10(tpm.gene.log2.m.rfd.ctr.tz.cd$length), c("SCLC, TZ-S and CD genes", "n=446"), file.name, xlim, ylim, col=c("blue", "black"), "topright", ylab.text, isFlip=T)

file.name <- file.path(wd.de.plots, "SCLC-SP_TRC_LENGTH-vs-GENE-NRFD_TZ_HO")
plotTRC(tpm.gene.log2.m.rfd.ctr.tz.ho$GENE_NRFD * -1, log10(tpm.gene.log2.m.rfd.ctr.tz.ho$length), c("SCLC, TZ-S and HO genes", "n=428"), file.name, xlim, ylim, col=c("blue", "black"), "topright", ylab.text, isFlip=T)




###
##
tpm.gene.log2.m.rfd.ctr.iz <- cbind(tpm.gene.log2.m.rfd.ctr.iz, de.tpm.gene[rownames(tpm.gene.log2.m.rfd.ctr.iz),])
tpm.gene.log2.m.rfd.ctr.iz$length <- abs(tpm.gene.log2.m.rfd.ctr.iz$end_position - tpm.gene.log2.m.rfd.ctr.iz$start_position) / 1000
tpm.gene.log2.m.rfd.ctr.iz.cd <- subset(tpm.gene.log2.m.rfd.ctr.iz, TRC > 0)
tpm.gene.log2.m.rfd.ctr.iz.ho <- subset(tpm.gene.log2.m.rfd.ctr.iz, TRC < 0)

##
xlim <- c(min(log10(tpm.gene.log2.m.rfd.ctr.iz$length)), max(log10(tpm.gene.log2.m.rfd.ctr.iz$length)))
ylim <- c(min(tpm.gene.log2.m.rfd.ctr.iz$GENE_NRFD), min(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD) * -1)
ylab.text <- "Initiation efficiency"

file.name <- file.path(wd.de.plots, "TRC_LENGTH-vs-GENE-NRFD_IZ")
plotTRC(tpm.gene.log2.m.rfd.ctr.iz$GENE_NRFD, log10(tpm.gene.log2.m.rfd.ctr.iz$length), c("SCLC initiation zone (IZ)", "(n=2,565)"), file.name, xlim, ylim, col=c("red", "black"), "topright", ylab.text)

file.name <- file.path(wd.de.plots, "TRC_LENGTH-vs-GENE-NRFD_IZ_CD")
plotTRC(tpm.gene.log2.m.rfd.ctr.iz.cd$GENE_NRFD, log10(tpm.gene.log2.m.rfd.ctr.iz.cd$length), c("SCLC IZ (CD) genes", "(n=1,275)"), file.name, xlim, ylim, col=c("red", "black"), "topright", ylab.text)

file.name <- file.path(wd.de.plots, "TRC_LENGTH-vs-GENE-NRFD_IZ_HO")
plotTRC(tpm.gene.log2.m.rfd.ctr.iz.ho$GENE_NRFD, log10(tpm.gene.log2.m.rfd.ctr.iz.ho$length), c("SCLC IZ (HO) genes", "(n=1,289)"), file.name, xlim, ylim, col=c("red", "black"), "topright", ylab.text)

##
tpm.gene.log2.m.rfd.ctr.tz <- cbind(tpm.gene.log2.m.rfd.ctr.tz, de.tpm.gene[rownames(tpm.gene.log2.m.rfd.ctr.tz),])
tpm.gene.log2.m.rfd.ctr.tz$length <- abs(tpm.gene.log2.m.rfd.ctr.tz$end_position - tpm.gene.log2.m.rfd.ctr.tz$start_position) / 1000
tpm.gene.log2.m.rfd.ctr.tz.cd <- subset(tpm.gene.log2.m.rfd.ctr.tz, TRC > 0)
tpm.gene.log2.m.rfd.ctr.tz.ho <- subset(tpm.gene.log2.m.rfd.ctr.tz, TRC < 0)

##
xlim <- c(min(log10(tpm.gene.log2.m.rfd.ctr.iz$length)), max(log10(tpm.gene.log2.m.rfd.ctr.tz$length)))
ylim <- c(max(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD), min(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD)) * -1
ylab.text <- "Termination efficiency"

file.name <- file.path(wd.de.plots, "TRC_LENGTH-vs-GENE-NRFD_TZ")
plotTRC(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD * -1, log10(tpm.gene.log2.m.rfd.ctr.tz$length), c("SCLC termination zone (TZ) genes", "(n=1,577)"), file.name, xlim, ylim, col=c("blue", "black"), "topright", ylab.text, isFlip=T)

file.name <- file.path(wd.de.plots, "TRC_LENGTH-vs-GENE-NRFD_TZ_CD")
plotTRC(tpm.gene.log2.m.rfd.ctr.tz.cd$GENE_NRFD * -1, log10(tpm.gene.log2.m.rfd.ctr.tz.cd$length), c("SCLC TZ (CD) genes", "(n=802)"), file.name, xlim, ylim, col=c("blue", "black"), "topright", ylab.text, isFlip=T)

file.name <- file.path(wd.de.plots, "TRC_LENGTH-vs-GENE-NRFD_TZ_HO")
plotTRC(tpm.gene.log2.m.rfd.ctr.tz.ho$GENE_NRFD * -1, log10(tpm.gene.log2.m.rfd.ctr.tz.ho$length), c("SCLC TZ (HO) genes", "(n=775)"), file.name, xlim, ylim, col=c("blue", "black"), "topright", ylab.text, isFlip=T)
















###
## RHO
xlim <- c(min(tpm.gene.log2.m.rfd.ctr.iz$RHO), max(tpm.gene.log2.m.rfd.ctr.iz$RHO))
ylim <- c(min(tpm.gene.log2.m.rfd.ctr.iz$GENE_NRFD), max(tpm.gene.log2.m.rfd.ctr.iz$GENE_NRFD))

file.name <- file.path(wd.de.plots, "TRC_RHO-vs-GENE-NRFD_IZ")
plotTRC(tpm.gene.log2.m.rfd.ctr.iz$GENE_NRFD, tpm.gene.log2.m.rfd.ctr.iz$RHO, c("Initiation zone (IZ)", "(n=2,565)"), file.name, xlim, ylim, col=c("red", "black"), "topright")

file.name <- file.path(wd.de.plots, "TRC_RHO-vs-GENE-NRFD_IZ_CD")
plotTRC(tpm.gene.log2.m.rfd.ctr.iz.cd$GENE_NRFD, tpm.gene.log2.m.rfd.ctr.iz.cd$RHO, c("IZ (CD)", "(n=1,275)"), file.name, xlim, ylim, col=c("red", "black"), "topright")

file.name <- file.path(wd.de.plots, "TRC_RHO-vs-GENE-NRFD_IZ_HO")
plotTRC(tpm.gene.log2.m.rfd.ctr.iz.ho$GENE_NRFD, tpm.gene.log2.m.rfd.ctr.iz.ho$RHO, c("IZ (HO)", "(n=1,289)"), file.name, xlim, ylim, col=c("red", "black"), "topright")

##
xlim <- c(min(tpm.gene.log2.m.rfd.ctr.tz$RHO), max(tpm.gene.log2.m.rfd.ctr.tz$RHO))
ylim <- c(min(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD), max(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD))

file.name <- file.path(wd.de.plots, "TRC_RHO-vs-GENE-NRFD_TZ")
plotTRC(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD, tpm.gene.log2.m.rfd.ctr.tz$RHO, c("Termination zone (TZ)", "(n=1,577)"), file.name, xlim, ylim, col=c("blue", "black"), "bottomright")

file.name <- file.path(wd.de.plots, "TRC_RHO-vs-GENE-NRFD_TZ_CD")
plotTRC(tpm.gene.log2.m.rfd.ctr.tz.cd$GENE_NRFD, tpm.gene.log2.m.rfd.ctr.tz.cd$RHO, c("TZ (CD)", "(n=802)"), file.name, xlim, ylim, col=c("blue", "black"), "bottomright")

file.name <- file.path(wd.de.plots, "TRC_RHO-vs-GENE-NRFD_TZ_HO")
plotTRC(tpm.gene.log2.m.rfd.ctr.tz.ho$GENE_NRFD, tpm.gene.log2.m.rfd.ctr.tz.ho$RHO, c("TZ (HO)", "(n=775)"), file.name, xlim, ylim, col=c("blue", "black"), "bottomright")











ylim <- c(min(tpm.gene.log2.m.rfd.ctr.iz$GENE_NRFD), max(tpm.gene.log2.m.rfd.ctr.iz$GENE_NRFD))
xlim <- c(min(-log10(tpm.gene.log2.m.rfd.ctr.iz$length)), max(-log10(tpm.gene.log2.m.rfd.ctr.iz$)))

## 
file.name <- file.path(wd.de.plots, "TRC_P-vs-G-NRFD_IZ")
plotTRC(-log10(de.tpm.gene.iz$P), tpm.gene.log2.m.rfd.ctr.iz$GENE_NRFD, c("Initiation (IZ)", "(n=2,381)"), file.name, xlim, ylim, col=c("red", "black"), "topright")

file.name <- file.path(wd.de.plots, "TRC_P-vs-G-NRFD_IZ_HO")
plotTRC(-log10(de.tpm.gene.iz.ho$P), tpm.gene.log2.m.rfd.ctr.iz.ho$GENE_NRFD, c("IZ (HO)", "(n=1,198)"), file.name, xlim, ylim, col=c("red", "black"), "topright")

file.name <- file.path(wd.de.plots, "TRC_P-vs-G-NRFD_IZ_CD")
plotTRC(-log10(de.tpm.gene.iz.cd$P), tpm.gene.log2.m.rfd.ctr.iz.cd$GENE_NRFD, c("IZ (CD)", "(n=1,182)"), file.name, xlim, ylim, col=c("red", "black"), "topright")





# -----------------------------------------------------------------------------
# D.E.
# Last Modified: 08/01/20
# -----------------------------------------------------------------------------
#tpm.gene.log2.0      <- tpm.gene.log2[rownames(de.tpm.gene.iz),]
#tpm.gene.lung.log2.0 <- tpm.gene.lung.log2[rownames(de.tpm.gene.iz),]
overlaps <- intersect(intersect(rownames(tpm.gene.log2.m.rfd.ctr.iz.cd), rownames(tpm.gene.log2)), rownames(tpm.gene.lung.log2))
length(overlaps)
tpm.gene.log2.0      <- tpm.gene.log2[overlaps,]
tpm.gene.lung.log2.0 <- tpm.gene.lung.log2[overlaps,]

colnames <- c("P", "FDR", "N", "T", "LOG2_FC")
de2 <- toTable(0, length(colnames), nrow(tpm.gene.log2.0), colnames)
rownames(de2) <- rownames(tpm.gene.log2.0)

## SRC
de2$P <- mapply(x = 1:nrow(tpm.gene.log2.0), function(x) testU(as.numeric(tpm.gene.log2.0[x,]), as.numeric(tpm.gene.lung.log2.0[x,])))

## Log2 fold change
de2$N <- mapply(x = 1:nrow(tpm.gene.log2.0), function(x) median(as.numeric(tpm.gene.lung.log2.0[x,])))
de2$T <- mapply(x = 1:nrow(tpm.gene.log2.0), function(x) median(as.numeric(tpm.gene.log2.0[x,])))
de2$LOG2_FC <- de2$T - de2$N

de2.iz.sp <- de2[intersect(rownames(de2), rownames(tpm.gene.log2.m.rfd.ctr.iz.cd)),]
de2.iz.sp$FDR <- p.adjust(de2.iz.sp$P, method="BH", n=length(de2.iz.sp$P))
#de2.iz.sp$Q <- qvalue(de2.iz.sp$P)$qvalue
de2.iz.sp <- de2.iz.sp[order(de2.iz.sp$P),]

## FDR
#library(qvalue)
#de2$Q   <- qvalue(de2$P)$qvalue
de2$FDR <- p.adjust(de2$P, method="BH", n=length(de2$P))
#de$ANOVA_Q <- qvalue(de$ANOVA_P)$qvalue
de2 <- de2[order(de2$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de2.tpm.gene <- cbind(annot[rownames(de2.iz.sp),], de2.iz.sp)   ## BE EXTRA CAREFUL!!

save(de2.tpm.gene, samples, file=file.path(wd.de.data, "de2_sclc_tpm-gene-r5p47_lung_U_q_n70+41_iz_sp_cd.RData"))
writeTable(de2.tpm.gene, file.path(wd.de.data, "de2_sclc_tpm-gene-r5p47_lung_U_q_n70+41_iz_sp_cd.txt"), colnames=T, rownames=F, sep="\t")
# nrow(de2.tpm.gene)
# [1] 616

## Volcano
xlab.text <- "SCLC/Lung [log2FC]"
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_r5p47_iz_sp_cd_bh1e-16_lung")

genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c("SCLC, IZ-S and CD genes", paste0("n=616"))
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de2.tpm.gene, 1E-16, genes, file.de, file.main, xlab.text)









overlaps <- intersect(rownames(de2.tpm.gene), rownames(tpm.gene))
length(overlaps)
# [1] 2199
de2.tpm.gene <- de2.tpm.gene[overlaps,]
save(de2.tpm.gene, samples, file=file.path(wd.de.data, "de2_sclc_tpm-gene-r5p47_lung_U_q_n70+41_iz_sp_expressed-in-lung.RData"))
writeTable(de2.tpm.gene, file.path(wd.de.data, "de2_sclc_tpm-gene-r5p47_lung_U_q_n70+41_iz_sp_expressed-in-lung.txt"), colnames=T, rownames=F, sep="\t")

##
genes <- readTable(paste0(plot.de, "_lung.tab"), header=T, rownames=F, sep="\t")
file.main <- c("TTR + CTR", paste0("(n=17,311)"))
file.de <- paste0(plot.de, "_lung.pdf")
plotVolcano(de2.tpm.gene, 1E-16, genes, file.de, file.main, xlab.text)

###
## Volcano (All genes)
xlab.text <- "SCLC/Lung [log2FC]"
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_r5p47_iz_sp_bh1e-16")

genes <- readTable(paste0(plot.de, "_lung2.tab"), header=T, rownames=F, sep="\t")
file.main <- c("SCLC-specificly initiated (IZ-S), expressed genes", paste0("(n=1,245)"))
file.de <- paste0(plot.de, "_lung2.pdf")
plotVolcano(de2.tpm.gene, 1E-16, genes, file.de, file.main, xlab.text)

genes <- readTable(paste0(plot.de, "_lung_FOXH1.tab"), header=T, rownames=F, sep="\t")
file.main <- c("Initiation (IZ)", paste0("(n=1,299)"))
file.de <- paste0(plot.de, "_lung_FOXH1.pdf")
plotVolcano(de2.tpm.gene, 1E-16, genes, file.de, file.main, xlab.text)

de2.tpm.gene <- de2.tpm.gene[rownames(de.tpm.gene),]
genes <- readTable(paste0(plot.de, "_iz_lung_all.tab"), header=T, rownames=F, sep="\t")
file.main <- c("TTR + CTR", paste0("(n=31,207)"))
file.de <- paste0(plot.de, "_iz_lung_all.pdf")
plotVolcano(de2.tpm.gene, 1E-16, genes, file.de, file.main, xlab.text)

# -----------------------------------------------------------------------------
# 
# Last Modified: 31/10/18
# -----------------------------------------------------------------------------
## B2M
median(as.numeric(tpm.gene.lung.log2["ENSG00000166710",]))
median(as.numeric(tpm.gene.log2["ENSG00000166710",]))
testU(as.numeric(tpm.gene.lung.log2["ENSG00000166710",]), as.numeric(tpm.gene.log2["ENSG00000166710",]))
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_B2M_")
plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2["ENSG00000166710",], tpm.gene.log2["ENSG00000166710",], main="B2M", names=c("Lung", "SCLC"))

## MARS
median(as.numeric(tpm.gene.lung.log2["ENSG00000166986",]))
median(as.numeric(tpm.gene.log2["ENSG00000166986",]))
testU(as.numeric(tpm.gene.lung.log2["ENSG00000166986",]), as.numeric(tpm.gene.log2["ENSG00000166986",]))
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_MARS_")
plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2["ENSG00000166986",], tpm.gene.log2["ENSG00000166986",], main="MARS", names=c("Lung", "SCLC"))

## KIF18B
median(as.numeric(tpm.gene.lung.log2["ENSG00000186185",]))
median(as.numeric(tpm.gene.log2["ENSG00000186185",]))
testU(as.numeric(tpm.gene.lung.log2["ENSG00000186185",]), as.numeric(tpm.gene.log2["ENSG00000186185",]))
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_KIF18B_")
plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2["ENSG00000186185",], tpm.gene.log2["ENSG00000186185",], main="KIF18B", names=c("Lung", "SCLC"))

## BRCA2
median(as.numeric(tpm.gene.lung.log2["ENSG00000139618",]))
median(as.numeric(tpm.gene.log2["ENSG00000139618",]))
testU(as.numeric(tpm.gene.lung.log2["ENSG00000139618",]), as.numeric(tpm.gene.log2["ENSG00000139618",]))
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_BRCA2_")
plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2["ENSG00000139618",], tpm.gene.log2["ENSG00000139618",], main="BRCA2", names=c("Lung", "SCLC"))

## TP53
median(as.numeric(tpm.gene.lung.log2["ENSG00000141510",]))
median(as.numeric(tpm.gene.log2["ENSG00000141510",]))
testU(as.numeric(tpm.gene.lung.log2["ENSG00000141510",]), as.numeric(tpm.gene.log2["ENSG00000141510",]))
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_TP53_")
plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2["ENSG00000141510",], tpm.gene.log2["ENSG00000141510",], main="TP53", names=c("Lung", "SCLC"))

## RB1
median(as.numeric(tpm.gene.lung.log2["ENSG00000139687",]))
median(as.numeric(tpm.gene.log2["ENSG00000139687",]))
testU(as.numeric(tpm.gene.lung.log2["ENSG00000139687",]), as.numeric(tpm.gene.log2["ENSG00000139687",]))
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_RB1_")
plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2["ENSG00000139687",], tpm.gene.log2["ENSG00000139687",], main="RB1", names=c("Lung", "SCLC"))



## TONSL
median(as.numeric(tpm.gene.lung.log2["ENSG00000160949",]))
median(as.numeric(tpm.gene.log2["ENSG00000160949",]))
testU(as.numeric(tpm.gene.lung.log2["ENSG00000160949",]), as.numeric(tpm.gene.log2["ENSG00000160949",]))
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_TONSL_")
plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2["ENSG00000160949",], tpm.gene.log2["ENSG00000160949",], main="TONSL", names=c("Lung", "SCLC"))

## RECQL4
median(as.numeric(tpm.gene.lung.log2["ENSG00000160957",]))
median(as.numeric(tpm.gene.log2["ENSG00000160957",]))
testU(as.numeric(tpm.gene.lung.log2["ENSG00000160957",]), as.numeric(tpm.gene.log2["ENSG00000160957",]))
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_RECQL4_")
plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2["ENSG00000160957",], tpm.gene.log2["ENSG00000160957",], main="RECQL4", names=c("Lung", "SCLC"))

## CPSF1
median(as.numeric(tpm.gene.lung.log2["ENSG00000071894",]))
median(as.numeric(tpm.gene.log2["ENSG00000071894",]))
testU(as.numeric(tpm.gene.lung.log2["ENSG00000071894",]), as.numeric(tpm.gene.log2["ENSG00000071894",]))
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_CPSF1_")
plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2["ENSG00000071894",], tpm.gene.log2["ENSG00000071894",], main="CPSF1", names=c("Lung", "SCLC"))

## FOXH1
median(as.numeric(tpm.gene.lung.log2["ENSG00000160973",]))
median(as.numeric(tpm.gene.log2["ENSG00000160973",]))
testU(as.numeric(tpm.gene.lung.log2["ENSG00000160973",]), as.numeric(tpm.gene.log2["ENSG00000160973",]))
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_FOXH1_")
plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2["ENSG00000160973",], tpm.gene.log2["ENSG00000160973",], main="FOXH1", names=c("Lung", "SCLC"))

## PIF1
median(as.numeric(tpm.gene.lung.log2["ENSG00000140451",]))
median(as.numeric(tpm.gene.log2["ENSG00000140451",]))
testU(as.numeric(tpm.gene.lung.log2["ENSG00000140451",]), as.numeric(tpm.gene.log2["ENSG00000140451",]))
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_PIF1_")
plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2["ENSG00000140451",], tpm.gene.log2["ENSG00000140451",], main="PIF1", names=c("Lung", "SCLC"))

## TOR1AIP1
median(as.numeric(tpm.gene.lung.log2["ENSG00000143337",]))
median(as.numeric(tpm.gene.log2["ENSG00000143337",]))
testU(as.numeric(tpm.gene.lung.log2["ENSG00000143337",]), as.numeric(tpm.gene.log2["ENSG00000143337",]))
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_TOR1AIP1_")
plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2["ENSG00000143337",], tpm.gene.log2["ENSG00000143337",], main="TOR1AIP1", names=c("Lung", "SCLC"))

## BRCA2
median(as.numeric(tpm.gene.lung.log2["ENSG00000139618",]))
median(as.numeric(tpm.gene.log2["ENSG00000139618",]))
testU(as.numeric(tpm.gene.lung.log2["ENSG00000139618",]), as.numeric(tpm.gene.log2["ENSG00000139618",]))
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_BRCA2_")
plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2["ENSG00000139618",], tpm.gene.log2["ENSG00000139618",], main="BRCA2", names=c("Lung", "SCLC"))

## BRD9
median(as.numeric(tpm.gene.lung.log2["ENSG00000028310",]))
median(as.numeric(tpm.gene.log2["ENSG00000028310",]))
testU(as.numeric(tpm.gene.lung.log2["ENSG00000028310",]), as.numeric(tpm.gene.log2["ENSG00000028310",]))
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_BRD9")
plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2["ENSG00000028310",], tpm.gene.log2["ENSG00000028310",], main="BRD9 (ENSG00000028310)", names=c("Lung", "SCLC"))















## IZ / CD
genes <- readTable(paste0(plot.de, "_iz_cd.tab"), header=T, rownames=F, sep="\t")
file.main <- c("Initiation (IZ)", "Co-directional (CD)")
file.de <- paste0(plot.de, "_iz_cd.pdf")
plotVolcano(de.tpm.gene.iz.cd, 0.1, genes, file.de, file.main, xlab.text)

## IZ / HO
genes <- readTable(paste0(plot.de, "_iz_ho.tab"), header=T, rownames=F, sep="\t")
file.main <- c("Initiation (IZ)", "Head-on (HO)")
file.de <- paste0(plot.de, "_iz_ho.pdf")
plotVolcano(de.tpm.gene.iz.ho, 0.1, genes, file.de, file.main, xlab.text)

## TZ
genes <- readTable(paste0(plot.de, "_tz.tab"), header=T, rownames=F, sep="\t")
file.main <- c("Termination", "")
file.de <- paste0(plot.de, "_tz.pdf")
plotVolcano(de.tpm.gene.tz, 0.1, genes, file.de, file.main, xlab.text)

## TZ / CD
genes <- readTable(paste0(plot.de, "_tz_cd.tab"), header=T, rownames=F, sep="\t")
file.main <- c("Termination (TZ)", "Co-directional (CD)")
file.de <- paste0(plot.de, "_tz_cd.pdf")
plotVolcano(de.tpm.gene.tz.cd, 0.1, genes, file.de, file.main, xlab.text)

## TZ / HO
genes <- readTable(paste0(plot.de, "_tz_ho.tab"), header=T, rownames=F, sep="\t")
file.main <- c("Termination (TZ)", "Head-on (HO)")
file.de <- paste0(plot.de, "_tz_ho.pdf")
plotVolcano(de.tpm.gene.tz.ho, 0.1, genes, file.de, file.main, xlab.text)







##
plot.main <- "Differential expression between M2 and M1 in SCLC"
xlab.text <- "log2FC(SCLC M2/M1)"
plot.de <- file.path(wd.de.plots, "volcanoplot-r5p47-wgs_sclc_p1e-4_m2_")

## E2F3
genes <- readTable(paste0(plot.de, "E2F3.tab"), header=T, rownames=F, sep="\t")
rownames(genes) <- genes$GENE
genes <- genes[intersect(genes$GENE, de.tpm.gene$external_gene_name),]

file.main <- c(plot.main, "")
#file.de <- paste0(plot.de, "_chr2_log2FC1.pdf")
file.de <- paste0(plot.de, "E2F3.pdf")
plotVolcano(de.tpm.gene, 1.00E-04, genes, file.de, file.main, xlab.text)
















# -----------------------------------------------------------------------------
# 
# Last Modified: 
# -----------------------------------------------------------------------------
tpm.gene.log2.rl <- subset(tpm.gene.log2, TSS_RFD > 0)
tpm.gene.log2.ll <- subset(tpm.gene.log2, TSS_RFD < 0)

tpm.gene.log2.rl.cd <- subset(tpm.gene.log2.rl, TRC > 0)
tpm.gene.log2.rl.ho <- subset(tpm.gene.log2.rl, TRC < 0)
tpm.gene.log2.ll.cd <- subset(tpm.gene.log2.ll, TRC > 0)
tpm.gene.log2.ll.ho <- subset(tpm.gene.log2.ll, TRC < 0)

tpm.gene.log2.cd <- subset(tpm.gene.log2, TRC > 0)
tpm.gene.log2.ho <- subset(tpm.gene.log2, TRC < 0)

##
plotTRC <- function(tpm, rfd, main.text, file.name, xlim, ylim, col, pos) {
 xlab.text <- "RFD"
 ylab.text <- "log2(TPM + 0.01)"
 
 pdf(paste0(file.name, ".pdf"), height=6, width=6)
 plot(tpm ~ rfd, ylim=ylim, xlim=xlim, ylab="", xlab=xlab.text, main=main.text[1], col=col[1], pch=1, cex=1, cex.axis=1.7, cex.lab=1.7, cex.main=1.8)
 
 lm.fit <- lm(tpm ~ rfd)
 abline(lm.fit, col=col[2], lwd=3)
 
 cor <- cor.test(tpm, rfd, method="spearman", exact=F)
 legend(pos, c(paste0("rho = ", round0(cor[[4]], digits=2)), paste0("p-value = ", scientific(cor[[3]], digits=2))), text.col=col[2], bty="n", cex=1.75)
 
 mtext(ylab.text, side=2, line=2.85, cex=1.7)
 mtext(main.text[2], cex=1.2, line=0.3)
 dev.off()
}

xlim <- c(min(tpm.gene.log2$TSS_RFD), max(tpm.gene.log2$TSS_RFD))
ylim <- c(min(tpm.gene.log2$MEDIAN), max(tpm.gene.log2$MEDIAN))

## Right-leading, Co-directional
file.name <- file.path(wd.rt.plots, "TRC_RFD_RL_CD")
plotTRC(tpm.gene.log2.rl.cd$MEDIAN, tpm.gene.log2.rl.cd$TSS_RFD, c("Right-leading", "Co-directional"), file.name, c(0, 1), ylim, col=c("sandybrown", "blue"), "topright")

file.name <- file.path(wd.rt.plots, "TRC_RFD_RL_HO")
plotTRC(tpm.gene.log2.rl.ho$MEDIAN, tpm.gene.log2.rl.ho$TSS_RFD, c("Right-leading", "Head-on"), file.name, c(0, 1), ylim, col=c("sandybrown", "red"), "topright")

## Left-leading, Co-directional
file.name <- file.path(wd.rt.plots, "TRC_RFD_LL_CD")
plotTRC(tpm.gene.log2.ll.cd$MEDIAN, tpm.gene.log2.ll.cd$TSS_RFD, c("Left-leading", "Co-directional"), file.name, c(-1, 0), ylim, col=c("steelblue1", "blue"), "topright")

file.name <- file.path(wd.rt.plots, "TRC_RFD_LL_HO")
plotTRC(tpm.gene.log2.ll.ho$MEDIAN, tpm.gene.log2.ll.ho$TSS_RFD, c("Left-leading", "Head-on"), file.name, c(-1, 0), ylim, col=c("steelblue1", "red"), "topright")











# -----------------------------------------------------------------------------
# Principal component analysis (PCA)
# Last Modified: 23/11/17
# -----------------------------------------------------------------------------
#load(file.path(wd, base, "analysis/expression/kallisto", paste0(base.lcnec, "-tpm-de/data/", "de_", base, "_tpm-gene-r5p47_rb1_wilcox_q_n54.RData")))
genes.rb1.q0.1  <- rownames(subset(de.tpm.gene, FDR <= 0.1))
genes.rb1.q0.05 <- rownames(subset(de.tpm.gene, FDR <= 0.05))
## > length(genes.rb1.q0.1)
## [1] 639
## > length(genes.rb1.q0.05)
## [1] 145

## RB1 status on D.E genes
test <- tpm.gene[genes.rb1.q0.05, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
#test <- tpm.gene[genes.rb1.q0.1, rownames(samples)]
pca.de <- getPCA(t(test))

trait <- as.numeric(samples[,"RB1_MUT"])
trait[which(trait == 0)] <- "WT"
trait[which(trait == 1)] <- "RB1"

##
file.main <- "LCNEC RB1 status on 145 D.E. (Q < 0.05) genes"
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_RB1_Q0.05_145DE", size=6.5, file.main, "topleft", c("gray", "red", "dodgerblue"), NULL, flip.x=1, flip.y=-1)

file.main <- "LCNEC RB1 status on 639 D.E. (Q < 0.1) genes"
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_RB1_Q0.1_639DE", size=6.5, file.main, "topleft", c("gray", "red", "dodgerblue"), NULL, flip.x=1, flip.y=-1)
