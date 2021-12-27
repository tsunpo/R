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

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
tpm.gene.log2 <- log2(tpm.gene + 0.01)   ## Use pseudocount=0.01

b <- rownames(subset(samples, GROUP_ID == 0))
x <- rownames(subset(samples, GROUP_ID == 1))
n <- rownames(subset(samples, GROUP_ID == 2))
samples$GROUP_ID3 <- 0
samples[x,]$GROUP_ID3 <- 0
samples[n,]$GROUP_ID3 <- 1
samples[b,]$GROUP_ID3 <- 2
samples$GROUP_ID3 <- as.numeric(samples$GROUP_ID3)

tpm.gene.log2 <- tpm.gene.log2[, rownames(samples)]

# -----------------------------------------------------------------------------
# Residuals of expression
# Last Modified: 26/03/21; 28/08/20
# -----------------------------------------------------------------------------
tpm.gene.log2 <- tpm.gene.log2[, rownames(samples)]
scores <- pcaScores(pca.de)
pc1 <- scores[rownames(samples), "PC1"]
#as.numeric(resid(lm(as.numeric(tpm.gene.log2[1,]) ~ pc1)))

tpm.gene.res <- t(mapply(x = 1:nrow(tpm.gene), function(x) as.numeric(resid(lm(as.numeric(tpm.gene[x,]) ~ pc1)))))
colnames(tpm.gene.res) <- rownames(samples)
rownames(tpm.gene.res) <- rownames(tpm.gene)
save(tpm.gene.res, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.meidan0.res.RData")))

tpm.gene.log2.res <- t(mapply(x = 1:nrow(tpm.gene.log2), function(x) as.numeric(resid(lm(as.numeric(tpm.gene.log2[x,]) ~ pc1)))))
colnames(tpm.gene.log2.res) <- rownames(samples)
rownames(tpm.gene.log2.res) <- rownames(tpm.gene.log2)
save(tpm.gene.log2.res, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.median0.log2.res.RData")))

#tpm.gene.log2 <- tpm.gene.log2.res

# -----------------------------------------------------------------------------
# Spearman's rank correlation
# Last Modified: 22/08/20
# -----------------------------------------------------------------------------
colnames <- c("RHO", "P", "FDR", "X", "N", "B", "FC_N_X", "FC_B_X")
de <- toTable(0, length(colnames), nrow(tpm.gene.log2), colnames)
rownames(de) <- rownames(tpm.gene.log2)

## SRC
de$RHO <- mapply(x = 1:nrow(tpm.gene.log2.res), function(x) cor.test(as.numeric(tpm.gene.log2.res[x,]), samples$GROUP_ID3, method="spearman", exact=F)[[4]])
de$P   <- mapply(x = 1:nrow(tpm.gene.log2.res), function(x) cor.test(as.numeric(tpm.gene.log2.res[x,]), samples$GROUP_ID3, method="spearman", exact=F)[[3]])
## ANOVA
#samples$Cancer_Type <- as.factor(pheno.all$Cancer_Type)
#de$ANOVA_P <- mapply(x = 1:nrow(expr.pheno.log2), function(x) testANOVA(x, expr.pheno.log2, pheno.all))

## Log2 fold change
de$X <- median00(tpm.gene.log2, rownames(subset(samples, GROUP_ID3 == 0)))
de$N <- median00(tpm.gene.log2, rownames(subset(samples, GROUP_ID3 == 1)))
de$B <- median00(tpm.gene.log2, rownames(subset(samples, GROUP_ID3 == 2)))
de$FC_N_X <- de$N - de$X
de$FC_B_X <- de$B - de$X

## FDR
library(qvalue)
de$FDR   <- qvalue(de$P)$qvalue
#de$ANOVA_Q <- qvalue(de$ANOVA_P)$qvalue
de <- de[order(de$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)   ## BE EXTRA CAREFUL!!

save(de.tpm.gene, samples, file=file.path(wd.de.data, "src_esad_tpm-gene-median0-residual-pc1_X-N-B_src_q_n68.RData"))
writeTable(de.tpm.gene, file.path(wd.de.data, "src_esad_tpm-gene-median0-residual-pc1_X-N-B_src_q_n68.txt"), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Volcano plots of RB1-loss DE genes in LCNEC
# Figure(s)    : Figure S1 (A)
# Last Modified: 07/01/19
# -----------------------------------------------------------------------------
plotVolcano <- function(de, pvalue, genes, file.de, file.main) {
   #pvalue <- fdrToP(fdr, de)
   fdr <- pvalueToFDR(pvalue, de)
   de.sig <- subset(de, P <= pvalue)
   de.sig$log10P <- -log10(de.sig$P)
 
   de$log10P <- -log10(de$P)
   xmax <- max(de$LOG2_FC)
   #ymax <- max(de$log10P)
   ymax <- 8   ## 28/08/20
 
   pdf(file.de, height=6, width=6)
   plot(de$LOG2_FC, de$log10P, pch=16, xlim=c(-xmax, xmax), ylim=c(0, ymax), xlab="EAC B versus N [log2 fold change]", ylab="Significance [-log10(p-value)]", col="lightgray", main=file.main[1])

   de.up   <- subset(de.sig, LOG2_FC > 0)
   points(de.up$LOG2_FC, de.up$log10P, pch=16, col=red)
   de.down <- subset(de.sig, LOG2_FC < 0)
   points(de.down$LOG2_FC, de.down$log10P, pch=16, col=green)
 
   abline(v=c(-log2(2), log2(2)), lty=5, col="darkgray")
   abline(h=c(-log10(pvalue)), lty=5, col="black")
 
   for (g in 1:nrow(genes)) {
      gene <- subset(de, external_gene_name == genes[g,]$GENE)
      gene <- cbind(gene, genes[g,])
  
      if (nrow(gene) > 0) {
         points(gene$LOG2_FC, gene$log10P, pch=1, col="black")
   
         if (!is.na(gene$ADJ_1))
            if (is.na(gene$ADJ_2))
               text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=gene$ADJ_1, cex=0.85)
            else
               text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=c(gene$ADJ_1, gene$ADJ_2), cex=0.85)
         else
            if (gene$LOG2_FC > 0)
               text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=c(0, -0.5), cex=0.85)
            else
               text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=c(1, -0.5), cex=0.85)
      } else
         print(genes[g])
   }
 
   mtext(file.main[2], cex=1.2, line=0.3)
   legend("topleft", legend=c("Up-regulated (N < X < B)", "Down-regulated (N > X > B)"), col=c(red, green), pch=19)
   dev.off()
}

##
plot.main <- "80 additive expression after treatment"
plot.de <- file.path(wd.de.plots, "volcanoplot_EAC_r5p47_residual_pc1_N-X-B_p1e-3")
de.tpm.gene$LOG2_FC <- de.tpm.gene$FC_B_N

genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c(plot.main, "Between N, X and B (controlled for cell cycle)")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene, 1.00E-03, genes, file.de, file.main)

# -----------------------------------------------------------------------------
# Replace Ensembl Gene IDs to gene name in Reactome results (Up- and Down-regulation)
# -----------------------------------------------------------------------------
wd.de.pathway <- file.path(wd.de, "pathway")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_X-N-B_PC1_p1e-6_up")
#wd.de.reactome <- file.path(wd.de.pathway, "reactome_X-N-B_PC1_p1e-6_down")

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
#wd.de.reactome <- file.path(wd.de.pathway, "reactome_p1e-6_up")
#reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")
#reactome[4, 2] <- "Assembly of collagen fibrils and multimeric structures"

reactome.up <- subset(reactome, Entities.pValue <= 1e-2)[, 2:7]
reactome.up$log10P <- -log10(reactome.up$Entities.pValue)
reactome.up <- reactome.up[order(reactome.up$log10P),]
reactome.up[5, 1] <- "TP53 regulates transcription of additional cell death genes"
reactome.up[2, 1] <- "Defective SLC34A2 causes pulmonary alveolar microlithiasis"

main.text <- c("Up-regulated pathways", "64 genes")
file.de <- file.path(wd.de.reactome, "genes_N-X-B_res_p1e-3_n64_up.pdf")

pdf(file.de, height=3.2, width=7.5)
par(mar=c(4,23.5,4,0))   # increase y-axis margin.
barplot(reactome.up$log10P, main=main.text[1], las=1, horiz=T, xlim=c(0, 3), xaxt="n", names.arg=reactome.up$Pathway.name, col=red, xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
abline(v=2, lty=5)

axis(side=1, at=seq(0, 4, by=2))
mtext(main.text[2], line=0.3)
dev.off()

##
#wd.de.reactome <- file.path(wd.de.pathway, "reactome_p1e-6_down")
#reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")
#reactome[3, 2] <- "TCA cycle and respiratory electron transport"

reactome.down <- subset(reactome, Entities.pValue <= 1.9e-2)[, 2:7]
reactome.down$log10P <- -log10(reactome.down$Entities.pValue)
reactome.down <- reactome.down[order(reactome.down$log10P),]
reactome.down$log10P <- reactome.down$log10P * -1
reactome.down[2, 1] <- "Formyl peptide receptors bind formyl peptides and ligands"

main.text <- c("Down-regulated pathways", "23 genes")
file.de <- file.path(wd.de.reactome, "genes_N-X-B_res_p1e-3_n23_down.pdf")

pdf(file.de, height=3.2, width=7.5)
par(mar=c(4,0,4,23.5))   # increase y-axis margin.
posbar <- barplot(reactome.down$log10P, main=main.text[1], las=1, horiz=T, xlim=c(-3, 0), xaxt="n", names.arg="", col=green, xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=reactome.down$Pathway.name, xpd=NA)
abline(v=-2, lty=5)

axis(side=1, at=seq(-12, 0, by=2), labels=c(12, 10, 8, 6, 4, 2, 0))
mtext(main.text[2], line=0.3)
dev.off()
