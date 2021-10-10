# =============================================================================
# Manuscript   : 
# Name         : manuscripts/expression/esad-tpm-de.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 26/08/20
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

wd.de    <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

samples <- readTable(file.path(wd.rna, "esad_3rna_n68.txt"), header=T, rownames=T, sep="\t")

#load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.median0.RData")))
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
tpm.gene.log2 <- log2(tpm.gene + 0.01)   ## Use pseudocount=0.01

b <- rownames(subset(samples, GROUP_ID == 0))
x <- rownames(subset(samples, GROUP_ID == 1))
n <- rownames(subset(samples, GROUP_ID == 2))
samples$GROUP_ID3 <- 0
samples[n,]$GROUP_ID3 <- 0
samples[x,]$GROUP_ID3 <- 1
samples[b,]$GROUP_ID3 <- 2
samples$GROUP_ID3 <- as.numeric(samples$GROUP_ID3)

tpm.gene.log2 <- tpm.gene.log2[, rownames(samples)]
nrow(tpm.gene.log2)
# [1] 15658

# -----------------------------------------------------------------------------
# Spearman's rank correlation
# Last Modified: 22/08/20
# -----------------------------------------------------------------------------
colnames <- c("RHO", "P", "FDR", "N", "X", "B", "FC_X_N", "FC_B_N")
de <- toTable(0, length(colnames), nrow(tpm.gene.log2), colnames)
rownames(de) <- rownames(tpm.gene.log2)

## SRC
de$RHO <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples$GROUP_ID3, method="spearman", exact=F)[[4]])
de$P   <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples$GROUP_ID3, method="spearman", exact=F)[[3]])
## ANOVA
#samples$Cancer_Type <- as.factor(pheno.all$Cancer_Type)
#de$ANOVA_P <- mapply(x = 1:nrow(expr.pheno.log2), function(x) testANOVA(x, expr.pheno.log2, pheno.all))

## Log2 fold change
de$N <- median00(tpm.gene.log2, rownames(subset(samples, GROUP_ID3 == 0)))
de$X <- median00(tpm.gene.log2, rownames(subset(samples, GROUP_ID3 == 1)))
de$B <- median00(tpm.gene.log2, rownames(subset(samples, GROUP_ID3 == 2)))
de$FC_X_N <- de$X - de$N
de$FC_B_N <- de$B - de$N

## FDR
library(qvalue)
de$FDR   <- qvalue(de$P)$qvalue
#de$ANOVA_Q <- qvalue(de$ANOVA_P)$qvalue
de <- de[order(de$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)   ## BE EXTRA CAREFUL!!

save(de.tpm.gene, samples, file=file.path(wd.de.data, "src_esad_tpm-gene-r5p47_N-X-B_src_q_n68.RData"))
writeTable(de.tpm.gene, file.path(wd.de.data, "src_esad_tpm-gene-r5p47_N-X-B_src_q_n68.txt"), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Volcano plots of RB1-loss DE genes in LCNEC
# Figure(s)    : Figure S1 (A)
# Last Modified: 07/01/19
# -----------------------------------------------------------------------------
plot.main <- "325 additive expression after treatment"
plot.de <- file.path(wd.de.plots, "volcanoplot_EAC_median0_N-X-B_p1e-12_CNA")
de.tpm.gene$LOG2_FC <- de.tpm.gene$FC_B_N

genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c(plot.main, "Between N, X and B")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene, 1.00E-12, genes, file.de, file.main)

# -----------------------------------------------------------------------------
# Replace Ensembl Gene IDs to gene name in Reactome results (Up- and Down-regulation)
# -----------------------------------------------------------------------------
wd.de.pathway <- file.path(wd.de, "pathway")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_TR-T-N_p1e-6_up")
#wd.de.reactome <- file.path(wd.de.pathway, "reactome_TR-T-N_p1e-6_down")

list <- ensGene[,c("ensembl_gene_id",	"external_gene_name")]

reactome <- read.csv(file.path(wd.de.reactome, "result.csv"))
colnames(reactome) <- gsub("X.", "", colnames(reactome))
reactome$Submitted.entities.found <- as.vector(reactome$Submitted.entities.found)
reactome$Pathway.name             <- as.vector(reactome$Pathway.name)
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
wd.de.reactome <- file.path(wd.de.pathway, "reactome_N-X-B_p1e-12_up")
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")

reactome.up <- subset(reactome, Entities.pValue <= 1e-4)[, 2:7]
reactome.up$log10P <- -log10(reactome.up$Entities.pValue)
reactome.up <- reactome.up[order(reactome.up$log10P),]

main.text <- c("Up-regulated in cancer proliferation", "176 genes")
file.de <- file.path(wd.de.reactome, "genes_N-X-B_p1e-12_n176_up.pdf")

pdf(file.de, height=4, width=7.5)
par(mar=c(4,20.7,4,1.3))   # increase y-axis margin.
barplot(reactome.up$log10P, main=main.text[1], las=1, horiz=T, xlim=c(0, 12), xaxt="n", names.arg=reactome.up$Pathway.name, col=red, xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
abline(v=4, lty=5)

axis(side=1, at=seq(0, 12, by=2))
mtext(main.text[2], line=0.3)
dev.off()

##
wd.de.reactome <- file.path(wd.de.pathway, "reactome_N-X-B_p1e-12_down")
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")
#reactome[3, 2] <- "TCA cycle and respiratory electron transport"

reactome.down <- subset(reactome, Entities.pValue <= 1e-4)[, 2:7]
reactome.down$log10P <- -log10(reactome.down$Entities.pValue)
reactome.down <- reactome.down[order(reactome.down$log10P),]
reactome.down$log10P <- reactome.down$log10P * -1
#reactome.down[2, 1] <- "Formyl peptide receptors bind formyl peptides and ligands"

main.text <- c("Down-regulated in cancer proliferation", "149 genes")
file.de <- file.path(wd.de.reactome, "genes_N-X-B_p1e-12_n149_down.pdf")

pdf(file.de, height=2.8, width=7.5)
par(mar=c(4,1.3,4,20.7))   # increase y-axis margin.
posbar <- barplot(reactome.down$log10P, main=main.text[1], las=1, horiz=T, xlim=c(-12, 0), xaxt="n", names.arg="", col=green, xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=reactome.down$Pathway.name, xpd=NA)
abline(v=-4, lty=5)

axis(side=1, at=seq(-12, 0, by=2), labels=c(12, 10, 8, 6, 4, 2, 0))
mtext(main.text[2], line=0.3)
dev.off()
