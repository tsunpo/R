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
load(file.path(wd.src.ref, "mm10.RData"))

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
BASE <- "SFB"
base <- tolower(BASE)

#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2/"   ## tpyang@localhost
wd.rna   <- file.path(wd, BASE, "ngs/3RNA")
wd.anlys <- file.path(wd, BASE, "analysis")
wd.meta  <- file.path(wd, BASE, "metadata")

wd.de    <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

samples <- readTable(file.path(wd.rna, "sfb_3rna_n41.txt"), header=T, rownames=T, sep="\t")
samples <- samples[rownames(subset(samples, GROUP_NAME %in% c("NEMO", "WT"))),]

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
tpm.gene.log2 <- log2(tpm.gene + 1)   ## Use pseudocount=1
tpm.gene.log2 <- tpm.gene.log2[, rownames(samples)]
dim(tpm.gene.log2)
# [1] 18383     11

# -----------------------------------------------------------------------------
# PCA
# Last Modified: 10/10/21
# -----------------------------------------------------------------------------
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
test <- tpm.gene[, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))

trait <- samples$GROUP_NAME
file.main <- c("NEMO + WT", "")
plotPCA(1, 2, pca.de, trait, wd.de.plots, "PCA_SFB_median0_NEMO+WT_n11", size=6, file.main, "topright", c("NEMO", "WT"), c(green, grey), flip.x=-1, flip.y=1)

# -----------------------------------------------------------------------------
# Wilcoxon rank sum test (non-parametric; n=45, 22 TR vs 23 UN)
# Last Modified: 22/08/20
# -----------------------------------------------------------------------------
## Test: Wilcoxon/Mann–Whitney/U/wilcox.test
##       Student's/t.test
## FDR : Q/BH
## DE  : C (0) vs MA (1) as factor
argv      <- data.frame(predictor="GROUP_ID", predictor.wt=0, test="Wilcoxon", test.fdr="Q", stringsAsFactors=F)
file.name <- paste0("DE_SFB_tpm-gene-median0_NEMO-vs-WT_wilcox_q_n11")
file.main <- paste0("NEMO (n=5) vs. WT (n=6) in ", BASE)

de <- differentialAnalysis(tpm.gene.log2, samples, argv$predictor, argv$predictor.wt, argv$test, argv$test.fdr)

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)

save(de.tpm.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

# > de.tpm.gene[1:10,]
# ensembl_gene_id external_gene_name chromosome_name strand start_position end_position           gene_biotype           P       FDR GROUP_ID_WT
# ENSMUSG00000021680 ENSMUSG00000021680              Crhbp           chr13     -1       95431371     95444924         protein_coding 0.003890487 0.9215105   0.0000000
# ENSMUSG00000018595 ENSMUSG00000018595              Glra4            chrX     -1      136757674    136780141         protein_coding 0.005494111 0.9215105   0.0000000
# ENSMUSG00000034467 ENSMUSG00000034467            Dynlrb2            chr8      1      116504974    116515915         protein_coding 0.005494111 0.9215105   1.7690385
# ENSMUSG00000102175 ENSMUSG00000102175             Gm6119            chr1     -1        4692219      4693424   processed_pseudogene 0.005494111 0.9215105   0.0000000
# ENSMUSG00000109750 ENSMUSG00000109750         Hmgb1-rs17            chr8      1       33484196     33485167 unprocessed_pseudogene 0.005494111 0.9215105   0.0000000
# ENSMUSG00000113053 ENSMUSG00000113053            Gm18262           chr13      1       12864124     12866069   processed_pseudogene 0.005494111 0.9215105   0.0000000
# ENSMUSG00000115799 ENSMUSG00000115799            Vmn1r19            chr6      1       57403780     57406034         protein_coding 0.005494111 0.9215105   0.3814651
# ENSMUSG00000036198 ENSMUSG00000036198           Arhgap36            chrX      1       49463945     49500244         protein_coding 0.006735947 0.9215105   1.2876382
# ENSMUSG00000019194 ENSMUSG00000019194              Scn1b            chr7     -1       31116524     31127003         protein_coding 0.007546234 0.9215105   0.7049706
# ENSMUSG00000015337 ENSMUSG00000015337              Endog            chr2      1       30171493     30174069         protein_coding 0.007969413 0.9215105   1.5530695

# -----------------------------------------------------------------------------
# Volcano plots of RB1-loss DE genes in LCNEC
# Figure(s)    : Figure S1 (A)
# Last Modified: 07/01/19
# -----------------------------------------------------------------------------
load("/Users/tpyang/Work/uni-koeln/tyang2/SFB/analysis/expression/kallisto/sfb-tpm-de/data/DE_SFB_tpm-gene-median0_NEMO-vs-WT_wilcox_q_n11.RData")

xlab.text <- expression("Fold change " * "[log" * ""[2] * "]")
ylab.text <- expression("Significance " * "[-log" * ""[10] * "(" * italic("P") * ")]")

plot.de <- file.path(wd.de.plots, "volcanoplot_DE_SFB_median0_NEMO-vs-WT_p0.01_fc1")
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("NEMO vs. WT", "")
plotVolcano(de.tpm.gene, 0.01, genes, file.de, file.main, xlab.text, ylab.text, "topright", c("NEMO > WT", "NEMO < WT"), c(red, blue), c(red, blue), fold=1)

pvalue <- 0.01
fold <- 1
de.sig <- subset(de.tpm.gene, P <= pvalue)
de.down <- subset(de.sig, LOG2_FC < -fold)
de.up   <- subset(de.sig, LOG2_FC > fold)
nrow(de.up)
nrow(de.down)

# -----------------------------------------------------------------------------
# 
# Last Modified: 11/12/21
# -----------------------------------------------------------------------------
overlaps.up   <- intersect(rownames(de.tpm.gene.NB.up),   rownames(de.tpm.gene.NX.up))
overlaps.down <- intersect(rownames(de.tpm.gene.NB.down), rownames(de.tpm.gene.NX.down))

de.tpm.gene.NB.NX <- de.tpm.gene.NB[c(overlaps.up2, overlaps.down2),]
dim(subset(subset(de.tpm.gene.NB.NX, P < 1E-6), LOG2_FC >= 0))
dim(subset(subset(de.tpm.gene.NB.NX, P < 1E-6), LOG2_FC < 0))

###
##
overlaps.up2   <- intersect(overlaps.up,   rownames(de.tpm.gene.BX.up))
overlaps.down2 <- intersect(overlaps.down, rownames(de.tpm.gene.BX.down))

de.tpm.gene.NB.NX.BX <- de.tpm.gene.NB[c(overlaps.up2, overlaps.down2),]
dim(subset(subset(de.tpm.gene.NB.NX.BX, P < 1E-6), LOG2_FC >= 0))
dim(subset(subset(de.tpm.gene.NB.NX.BX, P < 1E-6), LOG2_FC < 0))
genes.NBX <- c(rownames(subset(subset(de.tpm.gene.NB.NX.BX, P < 1E-6), LOG2_FC >= 0)), rownames(subset(subset(de.tpm.gene.NB.NX.BX, P < 1E-6), LOG2_FC < 0)))

file.name <- paste0("DE_EAC_tpm-gene-median0_N-vs-B_wilcox_q_n46_N-B-X")
save(de.tpm.gene.NB.NX.BX, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene.NB.NX.BX, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

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
overlaps.up3   <- intersect(overlaps.up,   rownames(de.tpm.gene.BX.down))
overlaps.down3 <- intersect(overlaps.down, rownames(de.tpm.gene.BX.up))

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

