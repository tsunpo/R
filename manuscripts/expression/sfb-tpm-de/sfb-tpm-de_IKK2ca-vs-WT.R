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
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
wd.rna   <- file.path(wd, BASE, "ngs/3RNA")
wd.anlys <- file.path(wd, BASE, "analysis")
wd.meta  <- file.path(wd, BASE, "metadata")

wd.de    <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.gsea  <- file.path(wd.de, "gsea", "GSEA_Ranked Gene Lists_IKK2ca-vs-WT", "MH_mouse-ortholog hallmark gene sets")
wd.de.plots <- file.path(wd.de, "plots")

samples <- readTable(file.path(wd.rna, "sfb_3rna_n41.txt"), header=T, rownames=T, sep="\t")
samples <- samples[rownames(subset(samples, GROUP_NAME %in% c("IKK2ca", "WT"))),]

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
tpm.gene.log2 <- log2(tpm.gene + 1)   ## Use pseudocount=1
tpm.gene.log2 <- tpm.gene.log2[, rownames(samples)]
dim(tpm.gene.log2)
# [1] 15381

# -----------------------------------------------------------------------------
# PCA
# Last Modified: 10/10/21
# -----------------------------------------------------------------------------
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
test <- tpm.gene[, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))

trait <- samples$GROUP_NAME
file.main <- c("IKK2ca + WT", "")
plotPCA(1, 2, pca.de, trait, wd.de.plots, "PCA_SFB_median0_IKK2ca+WT_n12", size=6, file.main, "topright", c("IKK2ca", "WT"), c(orange, grey), flip.x=1, flip.y=1)

# -----------------------------------------------------------------------------
# Wilcoxon rank sum test (non-parametric; n=45, 22 TR vs 23 UN)
# Last Modified: 22/08/20
# -----------------------------------------------------------------------------
## Test: Wilcoxon/Mann–Whitney/U/wilcox.test
##       Student's/t.test
## FDR : Q/BH
## DE  : C (0) vs MA (1) as factor
argv      <- data.frame(predictor="GROUP_ID", predictor.wt=0, test="Wilcoxon", test.fdr="Q", stringsAsFactors=F)
file.name <- paste0("DE_SFB_tpm-gene-median0_IKK2ca-vs-WT_wilcox_q_n12")
file.main <- paste0("IKK2ca (n=6) vs. WT (n=6) in ", BASE)

de <- differentialAnalysis(tpm.gene.log2, samples, argv$predictor, argv$predictor.wt, argv$test, argv$test.fdr)

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)

save(de.tpm.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

# > de.tpm.gene[1:10,]
# ensembl_gene_id external_gene_name chromosome_name strand start_position end_position           gene_biotype           P       FDR GROUP_ID_WT
# ENSMUSG00000112640 ENSMUSG00000112640            Gm32687           chr10      1       81863934     81881081         protein_coding 0.003664548 0.2796019   0.0000000
# ENSMUSG00000071893 ENSMUSG00000071893             Vmn1r4            chr6      1       56924015     56958100         protein_coding 0.004336913 0.2796019   0.0000000
# ENSMUSG00000033227 ENSMUSG00000033227               Wnt6            chr1      1       74771892     74785322         protein_coding 0.004771822 0.2796019   0.1055003
# ENSMUSG00000044029 ENSMUSG00000044029            Olfr178           chr16     -1       58888308     58892063         protein_coding 0.004771822 0.2796019   0.7162802
# ENSMUSG00000055413 ENSMUSG00000055413              H2-Q5           chr17      1       35394126     35397800 polymorphic_pseudogene 0.004771822 0.2796019   0.1911959
# ENSMUSG00000033966 ENSMUSG00000033966              Cdkl4           chr17     -1       80523550     80577813         protein_coding 0.004998125 0.2796019   0.2051306
# ENSMUSG00000036198 ENSMUSG00000036198           Arhgap36            chrX      1       49463945     49500244         protein_coding 0.004998125 0.2796019   1.2876382
# ENSMUSG00000043931 ENSMUSG00000043931             Gimap7            chr6      1       48718621     48724636         protein_coding 0.004998125 0.2796019   0.3933584
# ENSMUSG00000044014 ENSMUSG00000044014              Npy5r            chr8     -1       66679965     66688128         protein_coding 0.004998125 0.2796019   1.9342948
# ENSMUSG00000069515 ENSMUSG00000069515               Lyz1           chr10     -1      117287797    117292868         protein_coding 0.004998125 0.2796019   1.2967185

file.name <- paste0("DE_SFB_tpm-gene-median0_IKK2ca-vs-WT_wilcox_q_n12")
load(file=file.path(wd.de.data, paste0(file.name, ".RData")))
nrow(subset(de.tpm.gene, P < 0.01))
nrow(subset(de.tpm.gene, P < 0.05))

de.tpm.gene[536, "FDR"]
de.tpm.gene[2284, "FDR"]

# -----------------------------------------------------------------------------
# Volcano plots of RB1-loss DE genes in LCNEC
# Figure(s)    : Figure S1 (A)
# Last Modified: 07/01/19
# -----------------------------------------------------------------------------
load("/Users/tpyang/Work/uni-koeln/tyang2/SFB/analysis/expression/kallisto/sfb-tpm-de/data/DE_SFB_tpm-gene-median0_IKK2ca-vs-WT_wilcox_q_n12.RData")

xlab.text <- expression("Fold change " * "[log" * ""[2] * "]")
ylab.text <- expression("Significance " * "[-log" * ""[10] * "(" * italic("P") * ")]")

plot.de <- file.path(wd.de.plots, "volcanoplot_DE_SFB_median0_IKK2ca-vs-WT_p0.01_fc1_log2")
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("IKK2ca vs. WT", "")
plotVolcano(de.tpm.gene, 0.01, genes, file.de, file.main, xlab.text, ylab.text, "topright", c("IKK2ca > WT", "IKK2ca < WT"), c(red, blue), c(red, blue), fold=1)

pvalue <- 0.01
fold <- 1
de.sig <- subset(de.tpm.gene, P <= pvalue)
de.down <- subset(de.sig, LOG2_FC < -fold)
de.up   <- subset(de.sig, LOG2_FC > fold)
nrow(de.up)
nrow(de.down)

# -----------------------------------------------------------------------------
# Volcano plots of RB1-loss DE genes in LCNEC
# Figure(s)    : Figure S1 (A)
# Last Modified: 07/01/19
# -----------------------------------------------------------------------------
hallmark1 <- readTable(file.path(wd.de.gsea, "HALLMARK_INFLAMMATORY_RESPONSE.tsv"), header=T, rownames=T, sep="\t")
hallmark2 <- readTable(file.path(wd.de.gsea, "HALLMARK_TNFA_SIGNALING_VIA_NFKB.tsv"), header=T, rownames=T, sep="\t")
hallmarks <- rbind(hallmark1, hallmark2)
hallmarks.core <- subset(hallmarks, CORE.ENRICHMENT == "Yes")
nrow(hallmarks.core)
# [1] 229
hallmarks.core.gene <- unique(hallmarks.core$SYMBOL)
length(hallmarks.core.gene)
# [1] 194

plot.de <- file.path(wd.de.plots, "volcanoplot_DE_SFB_median0_IKK2ca-vs-WT_p0.05")
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("IKK2ca vs. WT", "")
plotVolcano(de.tpm.gene, 0.01, genes, file.de, file.main, xlab.text, ylab.text, "topright", c("IKK2ca > WT", "IKK2ca < WT"), c(red, blue), c(red, blue), fold=1)

###
##
overlaps <- intersect(de.up$external_gene_name, hallmarks.core.gene)
colnames <- c("GENE", "ADJ_1", "ADJ_2")
genes <- toTable(NA, length(colnames), length(overlaps), colnames)
genes$GENE <- overlaps
writeTable(genes, file.path(wd.de.plots, "volcanoplot_DE_SFB_median0_IKK2ca-vs-WT_p0.01_fc1_log2_GSEA.tab"), colnames=T, rownames=F, sep="\t")

file.name <- "DE_SFB_tpm-gene-median0_IKK2ca-vs-WT_wilcox_q_n12_GSEA"
writeTable(subset(de.tpm.gene, external_gene_name %in% genes$GENE), file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

plot.de <- file.path(wd.de.plots, "volcanoplot_DE_SFB_median0_IKK2ca-vs-WT_p0.01_fc1_log2_GSEA_O")
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("IKK2ca vs. WT", "")
plotVolcano(de.tpm.gene, 0.01, genes, file.de, file.main, xlab.text, ylab.text, "topright", c("IKK2ca > WT", "IKK2ca < WT"), c(red, blue), c(red, blue), fold=1)






# -----------------------------------------------------------------------------
# Replace Ensembl Gene IDs to gene name in Reactome results (Up- and Down-regulation)
# -----------------------------------------------------------------------------
wd.de.pathway <- file.path(wd.de, "pathway", "")
wd.de.reactome <- file.path(wd.de.pathway, "up")
#wd.de.reactome <- file.path(wd.de.pathway, "down")

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
wd.de.reactome <- file.path(wd.de.pathway, "up")
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")
#reactome[9, 2] <- "Amp. from unattached kinetochores via MAD2"

reactome.up <- subset(reactome, Entities.pValue <= 1e-3)[, 2:7]
reactome.up$log10P <- -log10(reactome.up$Entities.pValue)
reactome.up <- reactome.up[order(reactome.up$log10P),]

main.text <- c("IKK2ca > WT", "731 genes")
file.de <- file.path(wd.de.reactome, "IKK2ca-vs-WT_p0.05_fc1.5_up_n731.pdf")

pdf(file.de, height=2.9, width=7.5)
par(mar=c(4,27,4,2))   # increase y-axis margin.
barplot(reactome.up$log10P, main=main.text[1], las=1, horiz=T, xlim=c(0, 7), xaxt="n", names.arg=reactome.up$Pathway.name, col=red, xlab="-log10(p-value)", width=1)   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
abline(v=4, lty=5)

axis(side=1, at=seq(0, 6, by=2))
mtext(main.text[2], line=0.3)
dev.off()

##
wd.de.reactome <- file.path(wd.de.pathway, "down")
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")
reactome[2, 2] <- "RUNX1 regulates genes involved in megakaryocyte differentiation"

reactome.down <- subset(reactome, Entities.pValue <= 1.00e-3)[, 2:7]
reactome.down$log10P <- -log10(reactome.down$Entities.pValue)
reactome.down <- reactome.down[order(reactome.down$log10P),]
reactome.down$log10P <- reactome.down$log10P * -1

main.text <- c("IKK2ca < WT", "423 genes")
file.de <- file.path(wd.de.reactome, "IKK2ca-vs-WT_p0.05_fc1.5_down_n423.pdf")

pdf(file.de, height=4.5, width=7.5)
par(mar=c(4,2,4,27))   # increase y-axis margin.
posbar <- barplot(reactome.down$log10P, main=main.text[1], las=1, horiz=T, xlim=c(-7, 0), xaxt="n", names.arg="", col=blue, xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=reactome.down$Pathway.name, xpd=NA)
abline(v=-4, lty=5)

axis(side=1, at=seq(-6, 0, by=2), labels=c(6, 4, 2, 0))
mtext(main.text[2], line=0.3)
dev.off()

###
## Converted to Human

###
## Link: https://www.statmethods.net/graphs/bar.html
wd.de.reactome <- file.path(wd.de.pathway, "up")
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")
#reactome[9, 2] <- "Amp. from unattached kinetochores via MAD2"

reactome.up <- subset(reactome, Entities.pValue <= 1e-2)[, 2:7]
reactome.up$log10P <- -log10(reactome.up$Entities.pValue)
reactome.up <- reactome.up[order(reactome.up$log10P),]

main.text <- c("IKK2ca > WT", "731 genes")
file.de <- file.path(wd.de.reactome, "IKK2ca-vs-WT_p0.05_fc1.5_up_n731.pdf")

pdf(file.de, height=3.2, width=7.5)
par(mar=c(4,27,4,2))   # increase y-axis margin.
barplot(reactome.up$log10P, main=main.text[1], las=1, horiz=T, xlim=c(0, 4), xaxt="n", names.arg=reactome.up$Pathway.name, col=red, xlab="-log10(p-value)", width=1)   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
abline(v=3, lty=5)

axis(side=1, at=seq(0, 4, by=2))
mtext(main.text[2], line=0.3)
dev.off()

##
wd.de.reactome <- file.path(wd.de.pathway, "down")
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")
reactome[1, 2] <- "Respiratory electron transport, ATP synthesis by chemiosmotic coupling"

reactome.down <- subset(reactome, Entities.pValue <= 1e-2)[, 2:7]
reactome.down$log10P <- -log10(reactome.down$Entities.pValue)
reactome.down <- reactome.down[order(reactome.down$log10P),]
reactome.down$log10P <- reactome.down$log10P * -1

main.text <- c("IKK2ca < WT", "423 genes")
file.de <- file.path(wd.de.reactome, "IKK2ca-vs-WT_p0.05_fc1.5_down_n423.pdf")

pdf(file.de, height=3, width=7.5)
par(mar=c(4,1,4,28))   # increase y-axis margin.
posbar <- barplot(reactome.down$log10P, main=main.text[1], las=1, horiz=T, xlim=c(-4, 0), xaxt="n", names.arg="", col=blue, xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=reactome.down$Pathway.name, xpd=NA)
abline(v=-3, lty=5)

axis(side=1, at=seq(-4, 0, by=2), labels=c(4, 2, 0))
mtext(main.text[2], line=0.3)
dev.off()

# -----------------------------------------------------------------------------
# GSVA
# Last Modified: 24/06/22
# -----------------------------------------------------------------------------
install.packages("BiocManager")
BiocManager::install("GSVA")

library(GSVA)


# -----------------------------------------------------------------------------
# Gene set enrichment analysis (GSEA) on LCNEC RB1/WT ranked gene lists
# Figure(s)    : Figure S1 (A and B)
# Last Modified: 08/01/19
# -----------------------------------------------------------------------------
file.name <- paste0("GSEA_SFB_tpm-gene-median0_IKK2ca-vs-WT_wilcox_q_n12")
writeRNKformat(de.tpm.gene, wd.de.gsea, file.name)










# -----------------------------------------------------------------------------
# PCA vs. Purites
# Last Modified: 02/06/21
# -----------------------------------------------------------------------------
rownames(samples) <- samples$PATIENT_ID2
overlaps <- intersect(samples$PATIENT_ID2, purities$sample_name)
samples <- samples[overlaps,]
samples$purity <- purities[overlaps,]$purity

