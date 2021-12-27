# =============================================================================
# Manuscript   : 
# Name         : manuscripts/expression/esad-tpm-de.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 21/08/20
# =============================================================================
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Commons.R", "Transcription.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
BASE <- "ESAD"
base <- tolower(BASE)

#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
wd.rna   <- file.path(wd, BASE, "ngs/3RNA")
wd.anlys <- file.path(wd, BASE, "analysis")
wd.meta  <- file.path(wd, BASE, "metadata")

wd.de    <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

samples <- readTable(file.path(wd.rna, "esad_3rna_n68.txt2"), header=T, rownames=T, sep="\t")
samples <- samples[!is.na(samples[, "tumorcontent_per_sample"]),]
samples$tumorcontent_per_sample <- as.numeric(sub("%", "", samples$tumorcontent_per_sample))/100

b <- rownames(subset(samples, GROUP_ID == 0))
x <- rownames(subset(samples, GROUP_ID == 1))
samples$GROUP_ID2 <- 0
samples[b,]$GROUP_ID2 <- 0
samples[x,]$GROUP_ID2 <- 1
samples$GROUP_ID2 <- as.factor(samples$GROUP_ID2)

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
tpm.gene.log2 <- log2(tpm.gene + 1)   ## Use pseudocount=0.01

tpm.gene.log2 <- tpm.gene.log2[, rownames(samples)]
nrow(tpm.gene.log2)
# [1] 15658

# -----------------------------------------------------------------------------
# Spearman's rank correlation
# Last Modified: 03/-6/21; 22/08/20
# -----------------------------------------------------------------------------
testANOVA <- function(expr, group, content) {
   fit1 <- lm(expr ~ group + content)
   fit2 <- lm(expr ~ group)
   a1 <- anova(fit1, fit2)
 
   return(a1$Pr[2])
}

colnames <- c("P", "FDR", "B", "X", "FC_X_B")
de <- toTable(0, length(colnames), nrow(tpm.gene.log2), colnames)
rownames(de) <- rownames(tpm.gene.log2)

## ANNOVA
de$P   <- mapply(x = 1:nrow(tpm.gene.log2), function(x) testANOVA(as.numeric(tpm.gene.log2[x,]), samples$GROUP_ID2, samples$tumorcontent_per_sample))

## Log2 fold change
de$B <- median00(tpm.gene.log2, rownames(subset(samples, GROUP_ID2 == 0)))
de$X <- median00(tpm.gene.log2, rownames(subset(samples, GROUP_ID2 == 1)))
de$FC_X_B <- de$X - de$B

## FDR
library(qvalue)
de$FDR   <- qvalue(de$P)$qvalue
de <- de[order(de$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)   ## BE EXTRA CAREFUL!!

save(de.tpm.gene, samples, file=file.path(wd.de.data, "DE_EAC_tpm-gene-median0_B-vs-X_anova_content_q_n44.RData"))
writeTable(de.tpm.gene, file.path(wd.de.data, "DE_EAC_tpm-gene-median0_B-vs-X_anova_content_q_n44.txt"), colnames=T, rownames=F, sep="\t")

###
##
plot.main <- "45 differentially expressed genes in EAC"
plot.de <- file.path(wd.de.plots, "volcanoplot_DE_EAC_median0_content_anova_B-vs-X_p1e-3_CNA")

genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c(plot.main, "Before (B) vs. After treatment (X)")
file.de <- paste0(plot.de, ".pdf")
de.tpm.gene$LOG2_FC <- de.tpm.gene$FC_X_B
plotVolcano(de.tpm.gene, 1.00E-03, genes, file.de, file.main, ymax=5)

# -----------------------------------------------------------------------------
# Replace Ensembl Gene IDs to gene name in Reactome results (Up- and Down-regulation)
# -----------------------------------------------------------------------------
wd.de.pathway <- file.path(wd.de, "pathway")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_B-vs-X_content-anova_p1e-3_up")
#wd.de.reactome <- file.path(wd.de.pathway, "reactome_B-vs-X_content-anova_p1e-3_down")

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

###
## Link: https://www.statmethods.net/graphs/bar.html
wd.de.reactome <- file.path(wd.de.pathway, "reactome_B-vs-X_content-anova_p1e-3_up")
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")
#reactome[9, 2] <- "Amp. from unattached kinetochores via MAD2"

reactome.up <- subset(reactome, Entities.pValue <= 1e-2)[, 2:7]
reactome.up$log10P <- -log10(reactome.up$Entities.pValue)
reactome.up <- reactome.up[order(reactome.up$log10P),]

main.text <- c("Up-regulated pathways in EAC", "15 genes")
file.de <- file.path(wd.de.reactome, "genes_B-vs-X_p1e-3_n15_up.pdf")

pdf(file.de, height=4.5, width=7.5)
par(mar=c(4,18,4,3.1))   # increase y-axis margin.
barplot(reactome.up$log10P, main=main.text[1], las=1, horiz=T, xlim=c(0, 3), xaxt="n", names.arg=reactome.up$Pathway.name, col=purple, xlab="-log10(p-value)", width=1)   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
abline(v=2, lty=5)

axis(side=1, at=seq(0, 3, by=1))
mtext(main.text[2], line=0.3)
dev.off()

##
wd.de.reactome <- file.path(wd.de.pathway, "reactome_B-vs-X_content-anova_p1e-3_down")
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")

reactome.down <- subset(reactome, Entities.pValue <= 1e-2)[, 2:7]
reactome.down$log10P <- -log10(reactome.down$Entities.pValue)
reactome.down <- reactome.down[order(reactome.down$log10P),]
reactome.down$log10P <- reactome.down$log10P * -1

main.text <- c("Down-regulated pathways in EAC", "30 genes")
file.de <- file.path(wd.de.reactome, "genes_B-vs-X_p1e-3_n30_down.pdf")

pdf(file.de, height=2.3, width=7.5)
par(mar=c(4,3,4,18))   # increase y-axis margin.
posbar <- barplot(reactome.down$log10P, main=main.text[1], las=1, horiz=T, xlim=c(-3, 0), xaxt="n", names.arg="", col=red, xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=reactome.down$Pathway.name, xpd=NA)
abline(v=-2, lty=5)

axis(side=1, at=seq(-3, 0, by=1), labels=c(-3, -2, -1, 0))
mtext(main.text[2], line=0.3)
dev.off()

# -----------------------------------------------------------------------------
# PCA on G1-S/G2-M genes
# Last Modified: 09/12/21; 18/05/19
# -----------------------------------------------------------------------------
plotCellCyclePCA <- function(wd.de.plots, BASE, tpm.gene, samples, trait, genes, isID, cols, flip.x, flip.y) {
   genes.G1S <- intersect(genes.G1S, rownames(tpm.gene))
   genes <- intersect(genes, rownames(tpm.gene))
   traits <- samples[,trait]
   ids <- NA
   if (isID) ids <- rownames(samples)
 
   ##
   test <- tpm.gene[genes.G1S, rownames(samples)]
   pca.de <- getPCA(t(test))
 
   file.main <- paste0(BASE, " on ", length(genes.G1S), "/304 G1-S genes")
   plotPCA(1, 2, pca.de, traits, wd.de.plots, paste0("pca_", tolower(BASE), "_G1-S"), size=6.5, file.main, "topleft", cols, ids, flip.x, flip.y)
   plotPCA(1, 3, pca.de, traits, wd.de.plots, paste0("pca_", tolower(BASE), "_G1-S"), size=6.5, file.main, "topleft", cols, ids, flip.x, flip.y)
   plotPCA(2, 3, pca.de, traits, wd.de.plots, paste0("pca_", tolower(BASE), "_G1-S"), size=6.5, file.main, "topleft", cols, ids, flip.x, flip.y)
 
   ##
   test <- tpm.gene[genes.G2M, rownames(samples)]
   pca.de <- getPCA(t(test))
 
   file.main <- paste0(BASE, " on ", length(genes.G2M), "/876 G2-M genes")
   plotPCA(1, 2, pca.de, traits, wd.de.plots, paste0("pca_", tolower(BASE), "_G2-M"), size=6.5, file.main, "topleft", cols, ids, flip.x, flip.y)
   plotPCA(1, 3, pca.de, traits, wd.de.plots, paste0("pca_", tolower(BASE), "_G2-M"), size=6.5, file.main, "topleft", cols, ids, flip.x, flip.y)
   plotPCA(2, 3, pca.de, traits, wd.de.plots, paste0("pca_", tolower(BASE), "_G2-M"), size=6.5, file.main, "topleft", cols, ids, flip.x, flip.y)
 
   ##
   test <- tpm.gene[genes, rownames(samples)]
   pca.de <- getPCA(t(test))
 
   file.main <- paste0(BASE, " on ", length(genes), "/151 read depth and expression correlated genes")
   plotPCA(1, 2, pca.de, traits, wd.de.plots, paste0("pca_", tolower(BASE), "_DRE"), size=6.5, file.main, "bottomleft", cols, ids, flip.x, flip.y)
   plotPCA(1, 3, pca.de, traits, wd.de.plots, paste0("pca_", tolower(BASE), "_DRE"), size=6.5, file.main, "topleft", cols, ids, flip.x, flip.y)
   plotPCA(2, 3, pca.de, traits, wd.de.plots, paste0("pca_", tolower(BASE), "_DEE"), size=6.5, file.main, "topleft", cols, ids, flip.x, flip.y)
}

##
plotCellCyclePCA(wd.de.plots, "EAC", tpm.gene, samples, "PCA_Tumour-Content-Genes", genes, isID=T, c("purple3", "forestgreen", "gold"), 1, 1)


