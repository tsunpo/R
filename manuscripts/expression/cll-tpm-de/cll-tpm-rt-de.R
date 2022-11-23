# =============================================================================
# Manuscript   : 
# Chapter      : 
# Name         : manuscripts/expression/cll-tpm-rt-de.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 29/05/20
# =============================================================================
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Commons.R", "Transcription.R", "TranscriptionReplicationConflict.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.ensGene.bed.1kb.RData"))

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE <- "CLL"
base <- tolower(BASE)

wd.wgs   <- file.path(wd, BASE, "ngs/WGS")
wd.rna   <- file.path(wd, BASE, "ngs/RNA")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.de       <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

samples.wgs <- readTable(file.path(wd.wgs, "cll_wgs_n96.txt"), header=T, rownames=T, sep="\t")
samples.rna <- readTable(file.path(wd.rna, "cll_rna_n71.txt"), header=T, rownames=T, sep="\t")
overlaps <- intersect(samples.wgs$SAMPLE_ID, samples.rna$ID_WGS)

samples.tpm.cll <- cbind(samples.wgs[overlaps,], samples.rna[overlaps,])
#samples <- subset(samples, COR < 0)
samples$M2 <- as.factor(samples$M2)
# > length(which(samples$M2 == 1))
# [1] 36
# > length(which(samples$M2 == 0))
# [1] 35

#load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.RData")))
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
#load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.r5p47.RData")))
tpm.gene <- tpm.gene[, samples$ID2_RNA]   ## VERY VERY VERY IMPORTANT!!!
tpm.gene.log2   <- log2(tpm.gene + 1)
tpm.gene.log2.m <- getLog2andMedian(tpm.gene, 1)
nrow(tpm.gene.log2.m)
# [1] 34908
# [1] 22807
# [1] 18502

# -----------------------------------------------------------------------------
# Purities
# Last Modified: 10/07/22
# -----------------------------------------------------------------------------
samples.wgs$purity <- NA
samples.wgs$ploidy <- NA
#samples.wgs$purity2 <- NA
#samples.wgs$ploidy2 <- NA
samples.wgs$purity3 <- NA
#samples.wgs$ploidy3 <- NA

for (s in 1:nrow(samples.wgs)) {
   sample <- samples.wgs$SAMPLE_ID[s]
   purities <- readTable(file.path(wd.wgs, "peiflyne", sample, paste0(sample, "_ANALYSIS"), paste0(sample, "_pupl.txt")), header=F, rownames=F, sep="\t")
   #purities2 <- readTable(file.path(wd.wgs, "peiflyne/2015", sample, paste0(sample, "_ANALYSIS"), paste0(sample, "_pupl.txt")), header=F, rownames=F, sep="\t")
 
   samples.wgs$purity[s]  <- purities$V2
   samples.wgs$ploidy[s]  <- purities$V3
   #samples.wgs$purity2[s] <- purities2$V2
   #samples.wgs$ploidy2[s] <- purities2$V3
   
   samples.wgs$purity3[s] <- samples.purity[sample, "Purity"]
}

#file.name <- file.path(wd.driver.plots, "2015", paste0("Correlation_NBL_purity_n57-1"))
#x <- samples.wgs$purity
#y <- samples.wgs$purity2
#plotCorrelation(file.name, "NB purity", "BWA-MEM", "Peifer 2015", x, y, pos="bottomright", cols=c("dimgray", "black"), size=5)

#file.name <- file.path(wd.driver.plots, "2015", paste0("Correlation_NBL_ploidy_n57-1"))
#x <- samples.wgs$ploidy
#y <- samples.wgs$ploidy2
#plotCorrelation(file.name, "NB ploidy", "BWA-MEM", "Peifer 2015", x, y, pos="bottomright", cols=c("dimgray", "black"), size=5)

file.name <- file.path(wd.de.plots, paste0("Correlation_CLL_Purity3"))
samples.wgs.nona <- samples.wgs[!is.na(samples.wgs$purity3),]
x <- samples.wgs.nona$purity
y <- samples.wgs.nona$purity3
plotCorrelation(file.name, "CLL purities (n=86)", "PeifLyne", "PCAWG", x, y, pos="bottomright", cols=c("dimgray", "black"), size=5)

for (q in 1:4) {
   samples.tpm.cll.Q4 <- subset(samples.tpm.cll, Q4 == q)
 
   file.name <- file.path(wd.de.plots, paste0("Correlation_CLL_Purity3_Q", q))
   samples.tpm.cll.Q4.nona <- samples.tpm.cll.Q4[!is.na(samples.tpm.cll.Q4$purity3),]
   x <- samples.tpm.cll.Q4.nona$purity
   y <- samples.tpm.cll.Q4.nona$purity3
   plotCorrelation(file.name, paste0("CLL purities (Q", q, ")"), "PeifLyne", "PCAWG", x, y, pos="bottomright", cols=c("dimgray", "black"), size=5)
}

# -----------------------------------------------------------------------------
# Correlation bwteen TPM and in-silico sorting
# Last Modified: 12/12/19; 08/01/19; 17/08/17
# -----------------------------------------------------------------------------
samples.tpm.cll$Purity <- samples.tpm.cll$purity
samples.tpm.cll$M2 <- samples.tpm.cll$M2 + 1
 
xlab.text <- "Expression vs. SCF index [rho]"
ylab.text <- "Expression vs. Purity [rho]"
pvalue <- 0.01

colnames <- c("GENE", "ADJ_1", "ADJ_2")
genes0 <- c("MKI67", "BCL2")
genes <- toTable(NA, length(colnames), length(genes0), colnames)
genes$GENE <- genes0
#genes[2, 2] <- -0.15
#genes[2, 3] <- 1.2

for (q in 1:4) {
   samples.tpm.cll.Q4 <- subset(samples.tpm.cll, Q4 == q)
   n <- nrow(samples.tpm.cll.Q4)
 
   load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
   tpm.gene <- tpm.gene[, samples.tpm.cll.Q4$ID2_RNA]   ## VERY VERY VERY IMPORTANT!!!
   tpm.gene.log2   <- log2(tpm.gene + 1)
 
   #getSRC(wd.de.data, BASE, tpm.gene.log2, samples.tpm.cll.Q4, "COR", n, "Q", q)
   #getSRC(wd.de.data, BASE, tpm.gene.log2, samples.tpm.cll.Q4, "Purity", n, "Q", q)
 
   ##
   expressed <- rownames(tpm.gene)
   de <- getCannoli(wd.de.data, BASE, n, expressed, TEST="COR", TEST2="Purity", M2=paste0("_Q", q))
   plot.de <- file.path(wd.de.plots, paste0("cannoliplot_SRC_", BASE, "_TPM-COR-Purity_P1E02_MEDIAN0_Q", q))
   genes <- getCannoliGenes(de, pvalue, genes0)
   #genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
   file.de <- paste0(plot.de, ".pdf")
   file.main <- c(paste0("Lymph-CLL Q", q, " (n=", n, ")"), "")
   plotCannoli(de, pvalue, genes, file.de, file.main, xlab.text, ylab.text, "bottomleft", c("", ""), c(red.lighter, blue.lighter), c(red, blue), fold=0, pos="bottomright")
}

for (m in 1:2) {
   samples.tpm.cll.M2 <- subset(samples.tpm.cll, M2 == m)
   n <- nrow(samples.tpm.cll.M2)
 
   load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
   tpm.gene <- tpm.gene[, samples.tpm.cll.M2$ID2_RNA]   ## VERY VERY VERY IMPORTANT!!!
   tpm.gene.log2   <- log2(tpm.gene + 1)
 
   #getSRC(wd.de.data, BASE, tpm.gene.log2, samples.tpm.cll.M2, "COR", n, "M", m)
   #getSRC(wd.de.data, BASE, tpm.gene.log2, samples.tpm.cll.M2, "Purity", n, "M", m)
 
   ##
   expressed <- rownames(tpm.gene)
   de <- getCannoli(wd.de.data, BASE, n, expressed, TEST="COR", TEST2="Purity", M2=paste0("_M", m))
   plot.de <- file.path(wd.de.plots, paste0("cannoliplot_SRC_", BASE, "_TPM-COR-Purity_P1E02_MEDIAN0_M", m))
   genes <- getCannoliGenes(de, pvalue, genes0)
   #genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
   file.de <- paste0(plot.de, ".pdf")
   file.main <- c(paste0("Lymph-CLL M", m, " (n=", n, ")"), "")
   plotCannoli(de, pvalue, genes, file.de, file.main, xlab.text, ylab.text, "bottomleft", c("", ""), c(red.lighter, blue.lighter), c(red, blue), fold=0, pos="bottomright")
}

n <- nrow(samples.tpm.cll)

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
tpm.gene <- tpm.gene[, samples.tpm.cll$ID2_RNA]   ## VERY VERY VERY IMPORTANT!!!
tpm.gene.log2   <- log2(tpm.gene + 1)

#getSRC0(wd.de.data, BASE, tpm.gene.log2, samples.tpm.cll, "COR", n)
#getSRC0(wd.de.data, BASE, tpm.gene.log2, samples.tpm.cll, "Purity", n)

##
expressed <- rownames(tpm.gene)
de <- getCannoli(wd.de.data, BASE, n, expressed, TEST="COR", TEST2="Purity", M2="")
plot.de <- file.path(wd.de.plots, paste0("cannoliplot_SRC_", BASE, "_TPM-COR-Purity_P1E02_MEDIAN0"))
genes <- getCannoliGenes(de, pvalue, genes0)
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c(paste0("Lymph-CLL (n=", n, ")"), "")
plotCannoli(de, pvalue, genes, file.de, file.main, xlab.text, ylab.text, "bottomleft", c("", ""), c(red.lighter, blue.lighter), c(red, blue), fold=0, pos="bottomright")

# -----------------------------------------------------------------------------
# 
# Last Modified: 19/10/22
# -----------------------------------------------------------------------------
genes.G1S <- readTable(file.path(wd.meta, "Dominguez_G1-S.list"), header=F, rownames=F, sep="")
genes.G2M <- readTable(file.path(wd.meta, "Dominguez_G2-M.list"), header=F, rownames=F, sep="")

for (q in 1:4) {
   n <- nrow(subset(samples.tpm.cll, Q4 == q))
   expressed <- rownames(tpm.gene)
   de <- getCannoli(wd.de.data, BASE, n, expressed, TEST="COR", TEST2="Purity", M2=paste0("_Q", q))
 
   plot.de <- file.path(wd.de.plots, paste0("cannoliplot_SRC_", BASE, "_TPM-COR-Purity_MEDIAN0_G1-S_Q", q))
   file.de <- paste0(plot.de, ".pdf")
   file.main <- c(paste0("Lymph-CLL Q", q, " (n=", n, ")"), "")
   plotCannoliGenes(de, genes.G1S, file.de, file.main, xlab.text, ylab.text, col=red, pos="bottomright")
 
   plot.de <- file.path(wd.de.plots, paste0("cannoliplot_SRC_", BASE, "_TPM-COR-Purity_MEDIAN0_G2-M_Q", q))
   file.de <- paste0(plot.de, ".pdf")
   file.main <- c(paste0("Lymph-CLL Q", q, " (n=", n, ")"), "")
   plotCannoliGenes(de, genes.G2M, file.de, file.main, xlab.text, ylab.text, col=green, pos="bottomright")
}

##
file.name <- file.path(wd.de.plots, paste0("Correlation_CLL_purity_SCF"))
x <- samples.tpm.cll$COR
y <- samples.tpm.cll$purity
plotCorrelation(file.name, "Lymph-CLL", "SCF index", "Purity", x, y, pos="bottomright", cols=c("dimgray", "black"), size=5)

for (q in 1:4) {
   samples.tpm.cll.Q4 <- subset(samples.tpm.cll, Q4 == q)
 
   file.name <- file.path(wd.de.plots, paste0("Correlation_CLL_purity_SCF_Q", q))
   x <- samples.tpm.cll.Q4$COR
   y <- samples.tpm.cll.Q4$purity
   plotCorrelation(file.name, paste0("Lymph-CLL (Q", q, ")"), "SCF index", "Purity", x, y, pos="bottomright", cols=c("dimgray", "black"), size=5)
}

##
file.name <- file.path(wd.de.plots, paste0("Correlation_CLL_ploidy_SCF"))
x <- samples.tpm.cll$COR
y <- samples.tpm.cll$ploidy
plotCorrelation(file.name, "Lymph-CLL", "SCF index", "Ploidy", x, y, pos="bottomright", cols=c("dimgray", "black"), size=5)

for (q in 1:4) {
   samples.tpm.cll.Q4 <- subset(samples.tpm.cll, Q4 == q)
 
   file.name <- file.path(wd.de.plots, paste0("Correlation_CLL_ploidy_SCF_Q", q))
   x <- samples.tpm.cll.Q4$COR
   y <- samples.tpm.cll.Q4$ploidy
   plotCorrelation(file.name, paste0("Lymph-CLL (Q", q, ")"), "SCF index", "Ploidy", x, y, pos="bottomright", cols=c("dimgray", "black"), size=5)
}




# -----------------------------------------------------------------------------
# Relationship between expression and in-silico sorting
# Last Modified: 29/05/20
# -----------------------------------------------------------------------------
colnames <- c("RHO", "P", "Q", "M1", "M2", "LOG2_FC")   ##"ANOVA_P", "ANOVA_Q", 
de <- toTable(0, length(colnames), nrow(tpm.gene.log2), colnames)
rownames(de) <- rownames(tpm.gene.log2)

## SRC
de$RHO <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples$COR, method="spearman", exact=F)[[4]])
de$P   <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples$COR, method="spearman", exact=F)[[3]])
de <- de[!is.na(de$P),]

## Log2 fold change
de$M1 <- median00(tpm.gene.log2, subset(samples, M2 == 0)$ID2_RNA)
de$M2 <- median00(tpm.gene.log2, subset(samples, M2 == 1)$ID2_RNA)
de$LOG2_FC <- de$M2 - de$M1

## FDR
library(qvalue)
de$Q <- qvalue(de$P)$qvalue
de <- de[order(de$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)   ## BE EXTRA CAREFUL!!

writeTable(de.tpm.gene, file.path(wd.de.data, "de_cll_tpm-gene-median0_src_q_n71.txt"), colnames=T, rownames=F, sep="\t")
save(de.tpm.gene, samples, file=file.path(wd.de.data, "de_cll_tpm-gene-median0_src_q_n71.RData"))
nrow(de.tpm.gene)
# [1] 22807

# -----------------------------------------------------------------------------
# RFD vs. TPM
# Last Modified: 29/05/20
# -----------------------------------------------------------------------------
nrds.RT.NRFD <- nrds.RT.NRFD.cll   ## 17/05/20 WHY nrds.RT.NRFD.sclc.nl? ## MUY MUY IMPORTANTE!!
tpm.gene.log2.m <- tpm.gene.log2.m
nrow(tpm.gene.log2.m)
# [1] 34908

nrds.RT.NRFD <- nrds.RT.NRFD.cll   ## MUY MUY IMPORTANTE!!
tpm.gene.log2.m <- tpm.gene.log2.m[rownames(de.tpm.gene),]
nrow(tpm.gene.log2.m)
# [1] 22807
# [1] 

tpm.gene.log2.m.rfd <- getTRC(tpm.gene.log2.m, nrds.RT.NRFD)
tpm.gene.log2.m.rfd$length <- abs(tpm.gene.log2.m.rfd$end_position - tpm.gene.log2.m.rfd$start_position)
nrow(tpm.gene.log2.m.rfd)
# [1] 30347
# [1] 20281
# [1] 
nrow(subset(tpm.gene.log2.m.rfd, TRC == 0))
# [1] 9
# [1] 6
# [1] 

###
## TTR and CTR (IZ +TZ)
save(tpm.gene.log2.m.rfd, file=file.path(wd.de.data, "tpm_gene_log2+1_m_rfd0.RData"))

load(file=file.path(wd.de.data, "tpm_gene_log2+1_m_rfd0.RData"))
file.name <- file.path(wd.de.data, "tpm_gene_log2+1_m_rfd.RData")
#file.name <- file.path(wd.de.data, "tpm-gene-median0_log2_m_de_rfd.RData")
setTRC(tpm.gene.log2.m.rfd, rfd=0.9, file.name)
load(file.name)

length(which(tpm.gene.log2.m.rfd.ctr.iz$GENE_NRFD < 0))
# [1] 106
# [1] 88
# [1] 
length(which(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD > 0))
# [1] 92
# [1] 77
# [1] 

# -----------------------------------------------------------------------------
# RFD vs. TPM (Tables)
# Last Modified: 29/05/20
# -----------------------------------------------------------------------------
#overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd))
#de.tpm.gene.rfd <- de.tpm.gene[overlaps,]
#de.tpm.gene.rfd$Q <- qvalue(de.tpm.gene.rfd$P)$qvalue
#writeTable(de.tpm.gene.rfd, file.path(wd.de.data, "de_cll_tpm-gene-median0_src_q_rfd_n70.txt"), colnames=T, rownames=F, sep="\t")

## IZ
overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.iz.e.cd))
de.tpm.gene.rfd.iz.e.cd <- de.tpm.gene[overlaps,]
de.tpm.gene.rfd.iz.e.cd$Q <- qvalue(de.tpm.gene.rfd.iz.e.cd$P)$qvalue
writeTable(de.tpm.gene.rfd.iz.e.cd, file.path(wd.de.data, "de_cll_tpm-gene-median0_src_q_rfd_iz_e_cd_n71.txt"), colnames=T, rownames=F, sep="\t")

overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.iz.e.ho))
de.tpm.gene.rfd.iz.e.ho <- de.tpm.gene[overlaps,]
de.tpm.gene.rfd.iz.e.ho$Q <- qvalue(de.tpm.gene.rfd.iz.e.ho$P)$qvalue
writeTable(de.tpm.gene.rfd.iz.e.ho, file.path(wd.de.data, "de_cll_tpm-gene-median0_src_q_rfd_iz_e_ho_n71.txt"), colnames=T, rownames=F, sep="\t")

## TZ
overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.tz.e.cd))
de.tpm.gene.rfd.tz.e.cd <- de.tpm.gene[overlaps,]
de.tpm.gene.rfd.tz.e.cd$Q <- qvalue(de.tpm.gene.rfd.tz.e.cd$P)$qvalue
writeTable(de.tpm.gene.rfd.tz.e.cd, file.path(wd.de.data, "de_cll_tpm-gene-median0_src_q_rfd_tz_e_cd_n71.txt"), colnames=T, rownames=F, sep="\t")

overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.tz.e.ho))
de.tpm.gene.rfd.tz.e.ho <- de.tpm.gene[overlaps,]
de.tpm.gene.rfd.tz.e.ho$Q <- qvalue(de.tpm.gene.rfd.tz.e.ho$P)$qvalue
writeTable(de.tpm.gene.rfd.tz.e.ho, file.path(wd.de.data, "de_cll_tpm-gene-median0_src_q_rfd_tz_e_ho_n71.txt"), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Shared genes (All)
# Last Modified: 10/09/20
# -----------------------------------------------------------------------------
load("/Users/tpyang/Work/uni-koeln/tyang2/CLL/analysis/expression/kallisto/cll-tpm-de/data/tpm_gene_log2_m_rfd.RData")

gene.cll   <- list(rownames(tpm.gene.log2.m.rfd.ctr.iz),   rownames(tpm.gene.log2.m.rfd.ctr.tz),   rownames(tpm.gene.log2.m.rfd.ttr))
gene.cll.e <- list(rownames(tpm.gene.log2.m.rfd.ctr.iz.e), rownames(tpm.gene.log2.m.rfd.ctr.tz.e), rownames(tpm.gene.log2.m.rfd.ttr.e))
gene.cll.l <- list(rownames(tpm.gene.log2.m.rfd.ctr.iz.l), rownames(tpm.gene.log2.m.rfd.ctr.tz.l), rownames(tpm.gene.log2.m.rfd.ttr.l))

# -----------------------------------------------------------------------------
# Density (All)
# Last Modified: 17/03/21
# -----------------------------------------------------------------------------
nrow(tpm.gene.log2.m.rfd.ctr.iz)/nrow(nrds.RT.NRFD.cll.ctr.iz)*1000
# [1] 15.30022
nrow(tpm.gene.log2.m.rfd.ctr.tz)/nrow(nrds.RT.NRFD.cll.ctr.tz)*1000
# [1] 10.78863
nrow(tpm.gene.log2.m.rfd.ttr)/nrow(nrds.RT.NRFD.cll.ttr)*1000
# [1] 10.97982

nrow(tpm.gene.log2.m.rfd.ctr.iz.e)/nrow(nrds.RT.NRFD.cll.ctr.iz.e)*1000
# [1] 22.1724
nrow(tpm.gene.log2.m.rfd.ctr.iz.l)/nrow(nrds.RT.NRFD.cll.ctr.iz.l)*1000
# [1] 6.715485

nrow(tpm.gene.log2.m.rfd.ctr.tz.e)/nrow(nrds.RT.NRFD.cll.ctr.tz.e)*1000
# [1] 17.90273
nrow(tpm.gene.log2.m.rfd.ctr.tz.l)/nrow(nrds.RT.NRFD.cll.ctr.tz.l)*1000
# [1] 6.572995

nrow(tpm.gene.log2.m.rfd.ttr.e)/nrow(nrds.RT.NRFD.cll.ttr.e)*1000
# [1] 17.58857
nrow(tpm.gene.log2.m.rfd.ttr.l)/nrow(nrds.RT.NRFD.cll.ttr.l)*1000
# [1] 6.164826

# -----------------------------------------------------------------------------
# E vs. L (All)
# Last Modified: 30/08/20
# -----------------------------------------------------------------------------
nrow(tpm.gene.log2.m.rfd.ttr)/30347
# [1] 0.7221142
nrow(tpm.gene.log2.m.rfd.ctr.iz)/30347
# [1] 0.1587307
nrow(tpm.gene.log2.m.rfd.ctr.tz)/30347
# [1] 0.1182654

nrow(tpm.gene.log2.m.rfd.ttr.e)/30347
# [1] 0.4875605
nrow(tpm.gene.log2.m.rfd.ttr.l)/30347
# [1] 0.2345537

nrow(tpm.gene.log2.m.rfd.ctr.iz.e)/30347
# [1] 0.1277556
nrow(tpm.gene.log2.m.rfd.ctr.iz.l)/30347
# [1] 0.03097506

nrow(tpm.gene.log2.m.rfd.ctr.tz.e)/30347
# [1] 0.07302205
nrow(tpm.gene.log2.m.rfd.ctr.tz.l)/30347
# [1] 0.04524335

###
## median0
nrow(tpm.gene.log2.m.rfd.ttr)/20281
# [1] 0.7200828
nrow(tpm.gene.log2.m.rfd.ctr.iz)/20281
# [1] 0.1711454
nrow(tpm.gene.log2.m.rfd.ctr.tz)/20281
# [1] 0.1081308

nrow(tpm.gene.log2.m.rfd.ctr.iz.e)/20281
# [1] 
nrow(tpm.gene.log2.m.rfd.ctr.iz.l)/20281
# [1] 

nrow(tpm.gene.log2.m.rfd.ctr.tz.e)/20281
# [1] 
nrow(tpm.gene.log2.m.rfd.ctr.tz.l)/20281
# [1] 

# -----------------------------------------------------------------------------
# RFD vs. TPM (All)
# Last Modified: 27/11/20; 01/09/20; 29/05/20
# -----------------------------------------------------------------------------
file.name <- file.path(wd.de.data, "tpm_gene_log2+1_m_rfd.RData")
load(file.name)

ylim <- c(min(tpm.gene.log2.m.rfd$MEDIAN), 12)
# > min(tpm.gene.log2.m.rfd$MEDIAN)
# [1] 0
# > max(tpm.gene.log2.m.rfd$MEDIAN)
# [1] 13.94057

file.name <- paste0("boxplot3_cll_tpm.gene+1_RFD")
plotBox3(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr, tpm.gene.log2.m.rfd.ctr.tz, tpm.gene.log2.m.rfd.ctr.iz, main="CLL expression", names=c("TTR", "TZ", "IZ"), cols=c("black", blue, red), ylim)

file.name <- paste0("boxplot4_cll_tpm.gene+1_RFD")
plotBox4(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.l, tpm.gene.log2.m.rfd.ctr.tz.e, tpm.gene.log2.m.rfd.ctr.iz.l, tpm.gene.log2.m.rfd.ctr.iz.e, main="CLL expression (CTR)", names=c("L", "E", "L", "E"), cols=c(blue, blue, red, red), ylim)

file.name <- paste0("boxplot6_cll_tpm.gene+1_RFD")
plotBox6(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr.l, tpm.gene.log2.m.rfd.ttr.e, tpm.gene.log2.m.rfd.ctr.tz.l, tpm.gene.log2.m.rfd.ctr.tz.e, tpm.gene.log2.m.rfd.ctr.iz.l, tpm.gene.log2.m.rfd.ctr.iz.e, main="CLL expression", names=c("L", "E", "L", "E", "L", "E"), cols=c("black", "black", blue, blue, red, red), ylim)




# -----------------------------------------------------------------------------
# RFD vs. TPM (All)
# Last Modified: 01/09/20; 29/05/20
# -----------------------------------------------------------------------------
ylim <- c(min(tpm.gene.log2.m.rfd$MEDIAN), 14.5)

file.name <- paste0("boxplot_cll_tpm.gene_rfd_TTR+IZ_nasa.blue")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr, tpm.gene.log2.m.rfd.ctr.iz, main="CLL genes", names=c("TTR", "IZ"), cols=c("black", red), ylim)

file.name <- paste0("boxplot_cll_tpm.gene_rfd_TTR+TZ_nasa.blue")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr, tpm.gene.log2.m.rfd.ctr.tz, main="CLL genes", names=c("TTR", "TZ"), cols=c("black", blue), ylim)

file.name <- paste0("boxplot_cll_tpm.gene_rfd_TZ+IZ_nasa.blue")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz, tpm.gene.log2.m.rfd.ctr.iz, main="CLL genes (CTR)", names=c("TZ", "IZ"), cols=c(blue, red), ylim)

##
file.name <- paste0("boxplot_cll_tpm.gene_rfd_IZ_E+L_nasa.blue")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.l, tpm.gene.log2.m.rfd.ctr.iz.e, main="CLL IZ genes", names=c("Late", "Early"), cols=c(red, red), ylim)

file.name <- paste0("boxplot_cll_tpm.gene_rfd_TZ_E+L_nasa.blue")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.l, tpm.gene.log2.m.rfd.ctr.tz.e, main="CLL TZ genes", names=c("Late", "Early"), cols=c(blue, blue), ylim)

file.name <- paste0("boxplot0_cll_tpm.gene_rfd_TTR_E+L")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr.l, tpm.gene.log2.m.rfd.ttr.e, main="CLL TTR genes", names=c("Late", "Early"), cols=c("black", "black"), ylim)

# -----------------------------------------------------------------------------
# RFD vs. TPM (median0)
# Last Modified: 02/09/20; 25/08/20; 29/05/20
# -----------------------------------------------------------------------------
ylim <- c(min(tpm.gene.log2.m.rfd$MEDIAN), 14.5)

file.name <- paste0("boxplot0_cll_tpm.gene.median0_rfd_TTR+IZ")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr, tpm.gene.log2.m.rfd.ctr.iz, main="CLC genes", names=c("TTR", "IZ"), cols=c("black", "red"), ylim)
file.name <- paste0("boxplot_cll_tpm.gene.median0_rfd_TTR+IZ")
plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr, tpm.gene.log2.m.rfd.ctr.iz, main="CLL genes", names=c("TTR", "IZ"), cols=c("black", "red"), ylim)

file.name <- paste0("boxplot0_cll_tpm.gene.median0_rfd_TTR+TZ")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr, tpm.gene.log2.m.rfd.ctr.tz, main="CLL genes", names=c("TTR", "TZ"), cols=c("black", "blue"), ylim)

#file.name <- paste0("boxplot0_cll_tpm.gene.median0_rfd_IZ+TZ")
#plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz, tpm.gene.log2.m.rfd.ctr.tz, main="CLL genes (CTR)", names=c("IZ", "TZ"), cols=c("red", "blue"), ylim)

file.name <- paste0("boxplot0_cll_tpm.gene.median0_rfd_TZ+IZ")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz, tpm.gene.log2.m.rfd.ctr.iz, main="CLL genes (CTR)", names=c("TZ", "IZ"), cols=c("blue", "red"), ylim)
#file.name <- paste0("boxplot_cll_tpm.gene.median0_rfd_TZ+IZ_long")
#plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz, tpm.gene.log2.m.rfd.ctr.iz, main="CLL genes", names=c("TZ", "IZ"), cols=c("blue", "red"), ylim)

##
file.name <- paste0("boxplot0_cll_tpm.gene.median0_rfd_IZ_E+L")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.l, tpm.gene.log2.m.rfd.ctr.iz.e, main="CLL IZ genes", names=c("Late", "Early"), cols=c("red", "red"), ylim)
file.name <- paste0("boxplot_cll_tpm.gene.median0_rfd_IZ_E+L")
plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.l, tpm.gene.log2.m.rfd.ctr.iz.e, main="CLL IZ genes", names=c("Late", "Early"), cols=c("red", "red"), ylim)

file.name <- paste0("boxplot0_cll_tpm.gene.median0_rfd_IZ_E_HO+CD")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.e.ho, tpm.gene.log2.m.rfd.ctr.iz.e.cd, main="CLL early IZ genes", names=c("HO", "CD"), cols=c("red", "red"), ylim)
file.name <- paste0("boxplot_cll_tpm.gene.median0_rfd_IZ_E_HO+CD")
plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.e.ho, tpm.gene.log2.m.rfd.ctr.iz.e.cd, main="CLL early IZ genes", names=c("HO", "CD"), cols=c("red", "red"), ylim)

##
file.name <- paste0("boxplot0_cll_tpm.gene.median0_rfd_TZ_E+L")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.l, tpm.gene.log2.m.rfd.ctr.tz.e, main="CLL TZ genes", names=c("Late", "Early"), cols=c("blue", "blue"), ylim)
#file.name <- paste0("boxplot_cll_tpm.gene.median0_rfd_TZ_E+L_long")
#plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.l, tpm.gene.log2.m.rfd.ctr.tz.e, main="CLL TZ genes", names=c("Late", "Early"), cols=c("blue", "blue"), ylim)

file.name <- paste0("boxplot0_cll_tpm.gene.median0_rfd_TZ_E_HO+CD")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.e.ho, tpm.gene.log2.m.rfd.ctr.tz.e.cd, main="CLL early TZ genes", names=c("HO", "CD"), cols=c("blue", "blue"), ylim)
#file.name <- paste0("boxplot_cll_tpm.gene.median0_rfd_TZ_E_HO+CD_long")
#plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.e.ho, tpm.gene.log2.m.rfd.ctr.tz.e.cd, main="CLL early TZ genes", names=c("HO", "CD"), cols=c("blue", "blue"), ylim)

# -----------------------------------------------------------------------------
# GENE_NRFD plotBoxTSSNRFD (median0)
# Last Modified: 31/08/20
# -----------------------------------------------------------------------------
ylim <- c(-0.01, 0.032)
plotBoxNRFD(base, BASE, ylim, ylim2=c(-36, 35), tpm.gene.log2.m.rfd.ctr.iz.e.cd, tpm.gene.log2.m.rfd.ctr.iz.e.ho, tpm.gene.log2.m.rfd.ctr.tz.e.cd, tpm.gene.log2.m.rfd.ctr.tz.e.ho)

# -----------------------------------------------------------------------------
# 
# Last Modified: 29/05/20
# -----------------------------------------------------------------------------
genes <- c("CCR8", "PXN", "TAP1", "MAX", "USP41")
for (g in 1:length(genes)) {
   id <- subset(ensGene, external_gene_name == genes[g])$ensembl_gene_id
   plotCYS(genes[g], as.numeric(tpm.gene.log2[id,]), samples$COR, 1, "black", "bottomright")
}

# -----------------------------------------------------------------------------
# 
# Last Modified: 18/09/20
# -----------------------------------------------------------------------------
## SCLC's top genes   ## ADD: 04/11/20
xlab.text <- "CLL S/G1 [log2 fold change]"
plot.de <- file.path(wd.de.plots, "volcanoplot_cll_median0_rfd_p1e-2_all_PIF1")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("CLL expressed genes (n=22,807)", "Correlation between expression and in-silico sorting")
plotVolcano(de.tpm.gene, 0.01, genes, file.de, file.main, xlab.text, ymax=7)




## CLL IZ-E, CD genes
xlab.text <- "CLL S/G1 [log2 fold change]"
plot.de <- file.path(wd.de.plots, "volcanoplot_cll_median0_rfd_p1e-2_iz_e_cd")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c("CLL early IZ, CD genes (n=32)", "Correlation between expression and in-silico sorting")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene.rfd.iz.e.cd, 0.01, genes, file.de, file.main, xlab.text, ymax=5.3)

## CLL TZ-E, HO genes
xlab.text <- "CLL S/G1 [log2 fold change]"
plot.de <- file.path(wd.de.plots, "volcanoplot_cll_median0_rfd_p1e-2_tz_e_ho")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c("CLL early TZ, HO genes (n=23)", "Correlation between expression and in-silico sorting")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene.rfd.tz.e.ho, 0.01, genes, file.de, file.main, xlab.text, ymax=5.3)

# -----------------------------------------------------------------------------
# Replace Ensembl Gene IDs to gene name in Reactome results (Up- and Down-regulation)
# -----------------------------------------------------------------------------
wd.de.pathway <- file.path(wd.de, "pathway")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_IZ-E-CD_p1e-2_up")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_IZ-E-CD_p1e-2_down")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_TZ-E-HO_p1e-2_up")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_TZ-E-HO_p1e-2_down")

###
## Link: https://www.statmethods.net/graphs/bar.html
wd.de.reactome <- file.path(wd.de.pathway, "reactome_IZ-E-CD_p1e-2_up")
setTable(wd.de.reactome)
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")
reactome[2, 2] <- "Regulation of cytoskeletal remodeling and cell spreading by IPP"

reactome.up <- subset(reactome, Entities.pValue <= 2e-2)[, 2:7]
reactome.up$log10P <- -log10(reactome.up$Entities.pValue)
reactome.up <- reactome.up[order(reactome.up$log10P),]

main.text <- c("Positively-correlated", "16 genes")
file.de <- file.path(wd.de.reactome, "genes_IZ-E-CD_p1e-2_n16_up.pdf")

pdf(file.de, height=2.2, width=7.5)
par(mar=c(4,27.5,4,1))   # increase y-axis margin.
barplot(reactome.up$log10P, main=main.text[1], las=1, horiz=T, xlim=c(0, 8), xaxt="n", names.arg=reactome.up$Pathway.name, col="gold", xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
abline(v=2, lty=5)

axis(side=1, at=seq(0, 8, by=1))
mtext(main.text[2], line=0.3)
dev.off()

###
## Link: https://www.statmethods.net/graphs/bar.html
wd.de.reactome <- file.path(wd.de.pathway, "reactome_IZ-E-CD_p1e-2_down")
setTable(wd.de.reactome)
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")

reactome.down <- subset(reactome, Entities.pValue <= 2e-2)[, 2:7]
reactome.down$log10P <- -log10(reactome.down$Entities.pValue)
reactome.down <- reactome.down[order(reactome.down$log10P),]
reactome.down$log10P <- reactome.down$log10P * -1

main.text <- c("Negatively-correlated", "16 genes")
file.de <- file.path(wd.de.reactome, "genes_IS-E-CD_p1e-3_n16_down.pdf")

pdf(file.de, height=1.87, width=7.5)
par(mar=c(4,1,4,27.5))   # increase y-axis margin.
posbar <- barplot(reactome.down$log10P, main=main.text[1], las=1, horiz=T, xlim=c(-8, 0), xaxt="n", names.arg="", col="steelblue1", xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=reactome.down$Pathway.name, xpd=NA)
abline(v=-2, lty=5)

axis(side=1, at=seq(-8, 0, by=1))
mtext(main.text[2], line=0.3)
dev.off()

###
## Link: https://www.statmethods.net/graphs/bar.html
wd.de.reactome <- file.path(wd.de.pathway, "reactome_TZ-E-HO_p1e-2_up")
setTable(wd.de.reactome)
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")

reactome.up <- subset(reactome, Entities.pValue <= 2e-2)[, 2:7]
reactome.up$log10P <- -log10(reactome.up$Entities.pValue)
reactome.up <- reactome.up[order(reactome.up$log10P),]

main.text <- c("Positively-correlated", "12 genes")
file.de <- file.path(wd.de.reactome, "genes_TZ-E-HO_p1e-2_n12_up.pdf")

pdf(file.de, height=2.2, width=7.5)
par(mar=c(4,27.5,4,1))   # increase y-axis margin.
barplot(reactome.up$log10P, main=main.text[1], las=1, horiz=T, xlim=c(0, 8), xaxt="n", names.arg=reactome.up$Pathway.name, col="gold", xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
abline(v=2, lty=5)

axis(side=1, at=seq(0, 8, by=1))
mtext(main.text[2], line=0.3)
dev.off()

##
wd.de.reactome <- file.path(wd.de.pathway, "reactome_TZ-E-HO_p1e-2_down")
setTable(wd.de.reactome)
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")
reactome[3, 2] <- "Antigen Presentation: Folding, assembly and peptide loading of MHC I"

reactome.down <- subset(reactome, Entities.pValue <= 2e-2)[, 2:7]
reactome.down$log10P <- -log10(reactome.down$Entities.pValue)
reactome.down <- reactome.down[order(reactome.down$log10P),]
reactome.down$log10P <- reactome.down$log10P * -1

main.text <- c("Negatively-correlated", "11 genes")
file.de <- file.path(wd.de.reactome, "genes_TS-E-HO_p1e-2_n11_down.pdf")

pdf(file.de, height=3.6, width=7.5)
par(mar=c(4,1,4,27.5))   # increase y-axis margin.
posbar <- barplot(reactome.down$log10P, main=main.text[1], las=1, horiz=T, xlim=c(-8, 0), xaxt="n", names.arg="", col="steelblue1", xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=reactome.down$Pathway.name, xpd=NA)
abline(v=-2, lty=5)

axis(side=1, at=seq(-8, 0, by=1))
mtext(main.text[2], line=0.3)
dev.off()














# -----------------------------------------------------------------------------
# Figures
# Last Modified: 06/01/20
# -----------------------------------------------------------------------------
plotBox3 <- function(wd.de.plots, file.name, tpm.1, tpm.2, tpm.3, main, names, cols, ylim) {
   trait <- rep(0, nrow(tpm.1))
   trait <- c(trait, rep(1, nrow(tpm.2)))
   trait <- c(trait, rep(2, nrow(tpm.3)))
   trait <- as.factor(trait)
 
   expr <- as.numeric(c(tpm.1$MEDIAN, tpm.2$MEDIAN, tpm.3$MEDIAN))
 
   pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=5, width=3.2)
   boxplot(expr ~ trait, outline=T, xaxt="n", ylab=paste0("log2(TPM + 0.01)"), main=main, boxcol=cols, whiskcol=cols, outcol=cols, medcol=cols, staplecol=cols, ylim=ylim, cex.axis=1.1, cex.lab=1.2, cex.main=1.25)
   text(2, 15, "***", col="black", cex=2.5)
 
   ##
   trait <- rep(0, nrow(tpm.1))
   trait <- c(trait, rep(1, nrow(tpm.2)))
   trait <- as.factor(trait)
   expr <- as.numeric(c(tpm.1$MEDIAN, tpm.2$MEDIAN))
 
   axis(side=1, at=seq(1, 3, by=1), labels=names, font=2, cex.axis=1.2)
   axis(side=1, at=2, labels="Total n=29,417", line=1.3, col=NA, cex.axis=1.2)
 
   mtext(paste0("p-value = ", scientific(wilcox.test(expr ~ trait, exact=F)$p.value)), cex=1.25, line=0.3)
   dev.off()
}

ylim <- c(min(tpm.gene.log2.m.rfd$MEDIAN), 15)
file.name <- paste0("boxplot_cll_tpm.gene.rfd_TTR+IZ+TZ_3.2*5_n70")
plotBox3(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr, tpm.gene.log2.m.rfd.ctr.iz, tpm.gene.log2.m.rfd.ctr.tz, main="CLL expression", names=c("TTR", "IZ", "TZ"), cols=c("black", "red", "blue"), ylim)






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
# 
# Last Modified: 31/10/18
# -----------------------------------------------------------------------------
plotCYS <- function(gene, cn, snr, pch, col, pos) {
   unit <- (max(snr) - min(snr))/10
   xlim <- c(min(snr) - unit, max(snr) + unit)
   unit <- (max(cn) - min(cn))/10
   ylim <- c(min(cn) - unit, max(cn) + unit)
 
   xlab.text <- "In-silico sorting [rho]"
   ylab.text <- "log2(TPM + 0.01)"
   id <- subset(ensGene, external_gene_name == gene)$ensembl_gene_id
   file.name <- file.path(wd.de.plots, paste0("TPM-vs-SORTING_n70_", genes[g]))
   
   pdf(paste0(file.name, ".pdf"), height=6, width=6)
   plot(cn ~ snr, ylim=ylim, xlim=xlim, ylab="", xlab=xlab.text, main=paste0(gene, " (", id, ")"), col=col, pch=pch, cex=2, cex.axis=1.7, cex.lab=1.7, cex.main=1.8)
 
   lm.fit <- lm(cn ~ snr)
   abline(lm.fit, col=col, lwd=3)
 
   cor <- cor.test(cn, snr, method="spearman", exact=F)
   legend(pos, c(paste0("rho = ", round0(cor[[4]], digits=2)), paste0("p-value = ", scientific(cor[[3]], digits=2))), text.col=col, bty="n", cex=1.8)

   mtext(ylab.text, side=2, line=2.85, cex=1.7)
   #mtext(main.text[2], cex=1.2, line=0.3)
   dev.off()
}

##
genes <- c("TECRL", "AL355490.1", "DHRS7B", "BIRC5", "RDM1")
for (g in 1:length(genes)) {
   id <- subset(ensGene, external_gene_name == genes[g])$ensembl_gene_id
   plotCYS(genes[g], as.numeric(tpm.gene.log2[id,]), samples$COR, 1, "black", "bottomright")
}

# -----------------------------------------------------------------------------
# Volcano plots
# Last Modified: 15/08/18
# -----------------------------------------------------------------------------
fdrToP <- function(fdr, de) {
   de.sig <- subset(de, Q <= fdr)
   de.sig$log10P <- -log10(de.sig$P)
 
   return(max(de.sig$P))
}

plotVolcano <- function(de, pvalue, genes, file.de, file.main, xlab.text, ymax=0) {
   #pvalue <- fdrToP(fdr, de)
   #fdr <- pvalueToFDR(pvalue, de)
   de.sig <- subset(de, P <= pvalue)
   de.sig$log10P <- -log10(de.sig$P)
 
   de$log10P <- -log10(de$P)
   xmax <- max(de$LOG2_FC)
   xmin <- min(de$LOG2_FC)
   if (ymax ==0) ymax <- max(de$log10P)
   
   pdf(file.de, height=6, width=6)
   plot(de$LOG2_FC, de$log10P, pch=16, xlim=c(-xmax, xmax), ylim=c(0, ymax), xlab=xlab.text, ylab="P-value significance [-log10]", col="lightgray", main=file.main[1], cex=1.4, cex.axis=1.2, cex.lab=1.25, cex.main=1.3)

   abline(h=c(-log10(pvalue)), lty=5)
   #text(xmax*-1 + 2*xmax/15, -log10(pvalue) - ymax/30, paste0("FDR=", fdr*100, "%"), cex=1.1)    ## SCLC (IZ)
   #text(xmax*-1 + 2*xmax/13, -log10(pvalue) - ymax/30, paste0("FDR=", fdr*100, "%"), cex=1.1)   ## SCLC (AA)

   de.up   <- subset(de.sig, LOG2_FC > 0)
   points(de.up$LOG2_FC, de.up$log10P, pch=16, col=yellow, cex=1.4)
   de.down <- subset(de.sig, LOG2_FC < 0)
   points(de.down$LOG2_FC, de.down$log10P, pch=16, col=lightblue, cex=1.4)
 
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
                  text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=c(0, -0.5), cex=1.25)
               else
                  text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=c(1, -0.5), cex=1.25)
         } else
            print(genes[g])
      }
   }
   
   mtext(file.main[2], cex=1.25, line=0.3)
   legend("topleft", legend=c("Positively", "Negatively"), col=c(yellow, lightblue), pch=19, cex=1.25)
   dev.off()
}

## CLL ALL genes
xlab.text <- "CLL S/G1 fold change [log2]"
plot.de <- file.path(wd.de.plots, "volcanoplot_cll_median0_rfd_p1e-3_all_Helicases")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("CLL expressed genes (n=22,807)", "Expression vs. In-silico sorting")
plotVolcano(de.tpm.gene, 0.001, genes, file.de, file.main, xlab.text)

#overlaps <- intersect(helicases, de.tpm.gene$external_gene_name)
#setdiff(helicases, overlaps)

xlab.text <- "CLL S/G1 fold change [log2]"
plot.de <- file.path(wd.de.plots, "volcanoplot_cll_median0_rfd_p1e-3_all_TFBS")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("CLL expressed genes (n=22,807)", "Expression vs. In-silico sorting")
plotVolcano(de.tpm.gene, 0.001, genes, file.de, file.main, xlab.text)

#overlaps <- intersect(tfbs, de.tpm.gene$external_gene_name)
#setdiff(tfbs, overlaps)




##
xlab.text <- "CLL M2/M1 [log2FC]"
plot.de <- file.path(wd.de.plots, "volcanoplot_cll_r5p47_rfd_fdr0.1")

## IZ
genes <- readTable(paste0(plot.de, "_iz.tab"), header=T, rownames=F, sep="\t")
genes$GENE <- rownames(genes)
file.main <- c("CLL initiated (IZ), expressed genes", paste0("(n=2,893)"))
file.de <- paste0(plot.de, "_iz.pdf")
plotVolcano(de.tpm.gene.rfd.iz, 0.1, genes, file.de, file.main, xlab.text)

# -----------------------------------------------------------------------------
# Gene length vs. RFD slopes
# Last Modified: 10/01/20
# -----------------------------------------------------------------------------
plotTRC <- function(rfd, tpm, main.text, file.name, xlim, ylim, col, pos, ylab.text, isFlip=F) {
   xlab.text <- "Gene length [log10]"

   pdf(paste0(file.name, ".pdf"), height=6, width=6)
   plot(rfd ~ tpm, ylim=ylim, xlim=xlim, ylab="", xlab=xlab.text, main="", col=col[1], pch=1, cex=1.2, cex.axis=1.5, cex.lab=1.55, yaxt="n")

   lm.fit <- lm(rfd ~ tpm)
   abline(lm.fit, col=col[2], lwd=3)
 
   cor <- cor.test(rfd, tpm, method="spearman", exact=F)
   legend(pos, c(paste0("rho = ", round0(cor[[4]], digits=2)), paste0("p-value = ", scientific(cor[[3]], digits=2))), text.col=col[2], bty="n", cex=1.55)
 
   if (isFlip)
      axis(side=2, at=seq(0, 0.2, by=0.1), labels=c(0, -0.1, -0.2), cex.axis=1.5)
   else
      axis(side=2, at=seq(0, 0.2, by=0.1), labels=c(0, 0.1, 0.2), cex.axis=1.5)
   mtext(ylab.text, side=2, line=2.85, cex=1.5)
   mtext(main.text[1], line=1.8, cex=1.55, font=2)
   mtext(main.text[2], line=0.3, cex=1.55)
   dev.off()
}

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
overlaps <- intersect(intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2)), rownames(tpm.gene.lung.log2))
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

## FDR
#library(qvalue)
#de2$Q   <- qvalue(de2$P)$qvalue
de2$FDR <- p.adjust(de2$P, method="BH", n=length(de2$P))
#de$ANOVA_Q <- qvalue(de$ANOVA_P)$qvalue
de2 <- de2[order(de2$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de2.tpm.gene <- cbind(annot[rownames(de2),], de2)   ## BE EXTRA CAREFUL!!

save(de2.tpm.gene, samples, file=file.path(wd.de.data, "de2_sclc_tpm-gene-r5p47_lung_U_q_n70+41.RData"))
writeTable(de2.tpm.gene, file.path(wd.de.data, "de2_sclc_tpm-gene-r5p47_lung_U_q_n70+41.txt"), colnames=T, rownames=F, sep="\t")

## Volcano
xlab.text <- "SCLC/Lung [log2FC]"
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_r5p47_bh1e-16")

genes <- readTable(paste0(plot.de, "_lung.tab"), header=T, rownames=F, sep="\t")
file.main <- c("TTR + CTR", paste0("(n=17,311)"))
file.de <- paste0(plot.de, "_lung.pdf")
plotVolcano(de2.tpm.gene, 1E-16, genes, file.de, file.main, xlab.text)

###
## Volcano (All genes)
xlab.text <- "SCLC/Lung [log2FC]"
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_q1e-06")

genes <- readTable(paste0(plot.de, "_iz_lung.tab"), header=T, rownames=F, sep="\t")
file.main <- c("Initiation (IZ)", paste0("(n=4,108)"))
file.de <- paste0(plot.de, "_iz_lung.pdf")
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
## PIF1
median(as.numeric(tpm.gene.lung.log2["ENSG00000140451",]))
median(as.numeric(tpm.gene.log2["ENSG00000140451",]))
testU(as.numeric(tpm.gene.lung.log2["ENSG00000140451",]), as.numeric(tpm.gene.log2["ENSG00000140451",]))
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_PIF1")
plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2["ENSG00000140451",], tpm.gene.log2["ENSG00000140451",], main="PIF1 (ENSG00000140451)", names=c("Lung", "SCLC"))

## TOR1AIP1
median(as.numeric(tpm.gene.lung.log2["ENSG00000143337",]))
median(as.numeric(tpm.gene.log2["ENSG00000143337",]))
testU(as.numeric(tpm.gene.lung.log2["ENSG00000143337",]), as.numeric(tpm.gene.log2["ENSG00000143337",]))
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_TOR1AIP1")
plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2["ENSG00000143337",], tpm.gene.log2["ENSG00000143337",], main="TOR1AIP1 (ENSG00000143337)", names=c("Lung", "SCLC"))

## BRCA2
median(as.numeric(tpm.gene.lung.log2["ENSG00000139618",]))
median(as.numeric(tpm.gene.log2["ENSG00000139618",]))
testU(as.numeric(tpm.gene.lung.log2["ENSG00000139618",]), as.numeric(tpm.gene.log2["ENSG00000139618",]))
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_BRCA2")
plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2["ENSG00000139618",], tpm.gene.log2["ENSG00000139618",], main="BRCA2 (ENSG00000139618)", names=c("Lung", "SCLC"))

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
