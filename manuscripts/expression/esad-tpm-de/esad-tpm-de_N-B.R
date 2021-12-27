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
samples <- samples[rownames(subset(samples, GROUP_ID != 1)),]
n <- rownames(subset(samples, GROUP_ID == 2))
b <- rownames(subset(samples, GROUP_ID == 0))
samples$GROUP_ID2 <- 0
samples[n,]$GROUP_ID2 <- 0
samples[b,]$GROUP_ID2 <- 1
samples$GROUP_ID2 <- as.factor(samples$GROUP_ID2)

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
tpm.gene.log2 <- log2(tpm.gene + 1)   ## Use pseudocount=0.01

tpm.gene.log2 <- tpm.gene.log2[, rownames(samples)]
nrow(tpm.gene.log2)
# [1] 15658

# -----------------------------------------------------------------------------
# Wilcoxon rank sum test (non-parametric; n=45, 22 TR vs 23 UN)
# Last Modified: 22/08/20
# -----------------------------------------------------------------------------
## Test: Wilcoxon/Mann–Whitney/U/wilcox.test
##       Student's/t.test
## FDR : Q/BH
## DE  : N (0) vs B (1) as factor
argv      <- data.frame(predictor="GROUP_ID2", predictor.wt=0, test="Wilcoxon", test.fdr="Q", stringsAsFactors=F)
file.name <- paste0("DE_EAC_tpm-gene-median0_N-vs-B_wilcox_q_n46")
file.main <- paste0("N (n=23) vs. B (n=23) in ", BASE)

de <- differentialAnalysis(tpm.gene.log2, samples, argv$predictor, argv$predictor.wt, argv$test, argv$test.fdr)

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)

save(de.tpm.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

de.tpm.gene.NB.up   <- subset(de.tpm.gene, LOG2_FC >= 0)
de.tpm.gene.NB.down <- subset(de.tpm.gene, LOG2_FC < 0)

# -----------------------------------------------------------------------------
# N, B, X
# Last Modified: 20/12/21
# -----------------------------------------------------------------------------
de.tpm.gene.NB.NBX <- de.tpm.gene.NB[c(genes.NBX.up, genes.NBX.down),]

genes.NB.NBX.up   <- rownames(subset(subset(de.tpm.gene.NB.NBX, P < 1E-6), LOG2_FC >= 0))
genes.NB.NBX.down <- rownames(subset(subset(de.tpm.gene.NB.NBX, P < 1E-6), LOG2_FC < 0))
genes.NB.NBX <- c(genes.NB.NBX.up, genes.NB.NBX.down)

##
genes <- genes.NB.NBX
test <- tpm.gene[genes, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))
file.main <- c(paste0("N vs. B genes (n=", length(genes), ")"), "")
plotPCA(1, 2, pca.de, trait, wd.de.plots, "PCA_EAC_N-B-X_n68_N-vs-B", size=6, file.main, "topleft", trait.v, cols, NULL, flip.x=-1, flip.y=-1, legend.title=NA)

# -----------------------------------------------------------------------------
# N, X, B
# Last Modified: 21/12/21
# -----------------------------------------------------------------------------
de.tpm.gene.NB.NXB <- de.tpm.gene.NB[c(genes.NXB.up, genes.NXB.down),]

genes.NB.NXB.up   <- rownames(subset(subset(de.tpm.gene.NB.NXB, P < 1E-6), LOG2_FC >= 0))
genes.NB.NXB.down <- rownames(subset(subset(de.tpm.gene.NB.NXB, P < 1E-6), LOG2_FC < 0))
genes.NB.NXB <- c(genes.NB.NXB.up, genes.NB.NXB.down)

##
genes <- genes.NB.NXB
test <- tpm.gene[genes, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))
file.main <- c(paste0("N vs. B genes (n=", length(genes), ")"), "")
plotPCA(1, 2, pca.de, trait, wd.de.plots, "PCA_EAC_N-X-B_n68_N-vs-B", size=6, file.main, "bottomright", trait.v, cols, NULL, flip.x=-1, flip.y=-1, legend.title=NA)










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
   ymax <- 9
 
   pdf(file.de, height=6, width=6)
   plot(de$LOG2_FC, de$log10P, pch=16, xlim=c(-xmax, xmax), ylim=c(0, ymax), xlab="B to N fold change [log2]", ylab="P-value significance [-log10]", col="lightgray", main=file.main[1])

   #text(xmax*-1 + 2*xmax/28, -log10(pvalue) + ymax/42, paste0("FDR=", fdr, "%"), cex=0.85)
   #text(xmax*-1 + 2*xmax/35, -log10(pvalue) + ymax/42, "FDR=0.05", cex=0.85)
   #abline(h=c(-log10(fdrToP(0.1, de))), lty=5, col="darkgray")
   #text(xmax*-1 + 2*xmax/50, -log10(fdrToP(0.1, de)) + ymax/42, "FDR=0.1", col="darkgray", cex=0.85)

   de.up   <- subset(de.sig, LOG2_FC > 0)
   points(de.up$LOG2_FC, de.up$log10P, pch=16, col=red)
   de.down <- subset(de.sig, LOG2_FC < 0)
   points(de.down$LOG2_FC, de.down$log10P, pch=16, col=blue)
 
   abline(v=c(-log2(2), log2(2)), lty=5, col="darkgray")
   abline(h=c(-log10(pvalue)), lty=5, col="black")
   
   for (g in 1:nrow(genes)) {
      gene <- subset(de, external_gene_name == genes[g,]$GENE)
      gene <- cbind(gene, genes[g,])
      
      if (nrow(gene) > 0) {
         points(gene$LOG2_FC, gene$log10P, pch=1, col="black")
         
         if (!is.na(gene$ADJ_1))
            if (is.na(gene$ADJ_2))
               text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=gene$ADJ_1, cex=0.75)
            else
               text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=c(gene$ADJ_1, gene$ADJ_2), cex=0.75)
         else
            if (gene$LOG2_FC > 0)
               text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=c(0, -0.5), cex=0.75)
            else
               text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=c(1, -0.5), cex=0.75)
      } else
         print(genes[g])
   }
   
   mtext(file.main[2], cex=1.2, line=0.3)
   legend("topleft", legend=c("Up-regulated (N < B < X)", "Down-regulated (N > B > X)"), col=c(red, blue), pch=19)
   dev.off()
}

##
plot.main <- "1,240 differentially expressed genes in EAC"
plot.de <- file.path(wd.de.plots, "volcanoplot_DE_EAC_median0_N-vs-B_p1e-6_CNA")

genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c(plot.main, "Normal (N) vs. Tumour (B)")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene, 1.00E-06, genes, file.de, file.main)








de.tpm.gene.NB.NX.BX <- de.tpm.gene.NB[c(overlaps.up2, overlaps.down2),]
dim(subset(subset(de.tpm.gene.NB.NX.BX, P < 1E-3), LOG2_FC >= 0))
dim(subset(subset(de.tpm.gene.NB.NX.BX, P < 1E-3), LOG2_FC < 0))
genes.NB.NBX.up   <- rownames(subset(subset(de.tpm.gene.NB.NX.BX, FDR < 0.1), LOG2_FC >= 0))
genes.NB.NBX.down <- rownames((subset(subset(de.tpm.gene.NB.NX.BX, FDR < 0.1), LOG2_FC < 0)))
genes.NB.NBX <- c(genes.NB.NBX.up, genes.NB.NBX.down)

file.name <- paste0("DE_EAC_tpm-gene-median0-NBX_N-vs-B_wilcox_q_n46")
de <- de.tpm.gene.NB.NX.BX
de <- de[order(de$P),]
de$FDR <- qvalue(de$P)$qvalue
de.tpm.gene.NB.NX.BX <- de
save(de.tpm.gene.NB.NX.BX, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene.NB.NX.BX, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

##
file.name <- paste0("DE_EAC_tpm-gene-median0-NBX_N-vs-X_wilcox_q_n46")
de <- de.tpm.gene.NX[c(overlaps.up2, overlaps.down2),]
de <- de[order(de$P),]
de$FDR <- qvalue(de$P)$qvalue
de.tpm.gene.NX.NBX <- de
save(de.tpm.gene.NX.NBX, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene.NX.NBX, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")
dim(subset(subset(de.tpm.gene.NX.NBX, P < 1E-3), LOG2_FC >= 0))
dim(subset(subset(de.tpm.gene.NX.NBX, P < 1E-3), LOG2_FC < 0))
genes.NX.NBX.up   <- rownames(subset(subset(de.tpm.gene.NX.NBX, FDR < 0.1), LOG2_FC >= 0))
genes.NX.NBX.down <- rownames(subset(subset(de.tpm.gene.NX.NBX, FDR < 0.1), LOG2_FC < 0))
genes.NX.NBX <- c(genes.NX.NBX.up, genes.NX.NBX.down)



##
file.name <- paste0("DE_EAC_tpm-gene-median0-NBX_B-vs-X_wilcox_q_n46")
de <- de.tpm.gene.BX[c(overlaps.up2, overlaps.down2),]
de <- de[order(de$P),]
de$FDR <- qvalue(de$P)$qvalue
de.tpm.gene.BX.NBX <- de
save(de.tpm.gene.BX.NBX, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene.BX.NBX, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")
#dim(subset(subset(de.tpm.gene.BX.NBX, FDR < 0.001), LOG2_FC >= 0))
#dim(subset(subset(de.tpm.gene.BX.NBX, FDR < 0.001), LOG2_FC < 0))
dim(subset(subset(de.tpm.gene.BX.NBX, P < 1E-6), LOG2_FC >= 0))
dim(subset(subset(de.tpm.gene.BX.NBX, P < 1E-6), LOG2_FC < 0))
genes.BX.NBX.up   <- rownames(subset(subset(de.tpm.gene.BX.NBX, P < 1E-6), LOG2_FC >= 0))
genes.BX.NBX.down <- rownames(subset(subset(de.tpm.gene.BX.NBX, P < 1E-6), LOG2_FC < 0))
genes.BX.NBX <- c(genes.BX.NBX.up, genes.BX.NBX.down)



##
#genes.NBX.199 <- intersect(intersect(genes.NX.NBX, genes.NB.NBX), genes.BX.NBX)
#genes.NBX.25  <- intersect(intersect(genes.NX.NBX, genes.NB.NBX), genes.BX.NBX)

file.name <- paste0("DE_EAC_tpm-gene-median0-NBX_25_wilcox_q_n46")
de <- de.tpm.gene.NB[genes.NBX.25,]
de <- de[order(de$P),]
de$FDR <- qvalue(de$P)$qvalue
#de.tpm.gene.BX.NBX <- de
#save(de.tpm.gene.BX.NBX, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

## Volcano plots
plot.main <- "30 differentially expressed genes in EAC"
plot.de <- file.path(wd.de.plots, "volcanoplot_DE_EAC_median0_N-vs-B_p1e-6_N-B-X_CNA")

genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
rownames(genes) <- rownames(getGenes(genes$GENE))
genes <- genes[intersect(rownames(de.tpm.gene.NB.NX.BX), rownames(getGenes(genes$GENE))), ]
file.main <- c(plot.main, "Normal (N) vs. Tumour (B)")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene.NB.NX.BX, 1.00E-06, genes, file.de, file.main)

###
##
genes.NXB.up   <- intersect(overlaps.up,   rownames(de.tpm.gene.BX.down))
genes.NXB.down <- intersect(overlaps.down, rownames(de.tpm.gene.BX.up))

de.tpm.gene.NB.NX.XB <- de.tpm.gene.NB[c(overlaps.up3, overlaps.down3),]
dim(subset(subset(de.tpm.gene.NB.NX.XB, P < 1E-6), LOG2_FC >= 0))
dim(subset(subset(de.tpm.gene.NB.NX.XB, P < 1E-6), LOG2_FC < 0))

file.name <- paste0("DE_EAC_tpm-gene-median0_N-vs-B_wilcox_q_n46_N-X-B")
save(de.tpm.gene.NB.NX.XB, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene.NB.NX.XB, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

## Volcano plots
plot.main <- "1,087 differentially expressed genes in EAC"
plot.de <- file.path(wd.de.plots, "volcanoplot_DE_EAC_median0_N-vs-B_p1e-6_N-X-B_CNA")

genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
rownames(genes) <- rownames(getGenes(genes$GENE))
genes <- genes[intersect(rownames(de.tpm.gene.NB.NX.XB), rownames(getGenes(genes$GENE))), ]
file.main <- c(plot.main, "Normal (N) vs. Tumour (B)")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene.NB.NX.XB, 1.00E-06, genes, file.de, file.main)

# -----------------------------------------------------------------------------
# 
# Last Modified: 11/12/21
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Replace Ensembl Gene IDs to gene name in Reactome results (Up- and Down-regulation)
# -----------------------------------------------------------------------------
wd.de.pathway <- file.path(wd.de, "pathway")
#wd.de.reactome <- file.path(wd.de.pathway, "reactome_p1e-6_up")
#wd.de.reactome <- file.path(wd.de.pathway, "reactome_p1e-6_down")

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
wd.de.reactome <- file.path(wd.de.pathway, "reactome_N-vs-B_p1e-6_up")
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")
reactome[9, 2] <- "Amp. from unattached kinetochores via MAD2"

reactome.up <- subset(reactome, Entities.pValue <= 1e-6)[, 2:7]
reactome.up$log10P <- -log10(reactome.up$Entities.pValue)
reactome.up <- reactome.up[order(reactome.up$log10P),]

main.text <- c("Up-regulated pathways in EAC", "763 genes")
file.de <- file.path(wd.de.reactome, "genes_B-vs-N_p1e-6_n763_up.pdf")

pdf(file.de, height=4.5, width=7.5)
par(mar=c(4,18,4,3.1))   # increase y-axis margin.
barplot(reactome.up$log10P, main=main.text[1], las=1, horiz=T, xlim=c(0, 16), xaxt="n", names.arg=reactome.up$Pathway.name, col=red, xlab="-log10(p-value)", width=1)   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
abline(v=6, lty=5)

axis(side=1, at=seq(0, 16, by=2))
mtext(main.text[2], line=0.3)
dev.off()

##
wd.de.reactome <- file.path(wd.de.pathway, "reactome_N-vs-B_p1e-6_down")
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")

reactome.down <- subset(reactome, Entities.pValue <= 1e-6)[, 2:7]
reactome.down$log10P <- -log10(reactome.down$Entities.pValue)
reactome.down <- reactome.down[order(reactome.down$log10P),]
reactome.down$log10P <- reactome.down$log10P * -1

main.text <- c("Down-regulated pathways in EAC", "477 genes")
file.de <- file.path(wd.de.reactome, "genes_N-vs-B_p1e-6_n477_down.pdf")

pdf(file.de, height=2.1, width=7.5)
par(mar=c(4,3,4,18))   # increase y-axis margin.
posbar <- barplot(reactome.down$log10P, main=main.text[1], las=1, horiz=T, xlim=c(-16, 0), xaxt="n", names.arg="", col=blue, xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=reactome.down$Pathway.name, xpd=NA)
abline(v=-6, lty=5)

axis(side=1, at=seq(-16, 0, by=2), labels=c(16, 14, 12, 10, 8, 6, 4, 2, 0))
mtext(main.text[2], line=0.3)
dev.off()











# -----------------------------------------------------------------------------
# PCA vs. Purites
# Last Modified: 02/06/21
# -----------------------------------------------------------------------------
rownames(samples) <- samples$PATIENT_ID2
overlaps <- intersect(samples$PATIENT_ID2, purities$sample_name)
samples <- samples[overlaps,]
samples$purity <- purities[overlaps,]$purity

