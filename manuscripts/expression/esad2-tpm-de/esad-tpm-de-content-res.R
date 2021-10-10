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

#load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.median0.RData")))
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
tpm.gene.log2 <- log2(tpm.gene + 0.01)   ## Use pseudocount=0.01

tpm.gene      <- tpm.gene[, rownames(samples)]
tpm.gene.log2 <- tpm.gene.log2[, rownames(samples)]

# -----------------------------------------------------------------------------
# Residuals of expression
# Last Modified: 26/03/21; 28/08/20
# -----------------------------------------------------------------------------
tpm.gene.res <- t(mapply(x = 1:nrow(tpm.gene), function(x) as.numeric(resid(lm(as.numeric(tpm.gene[x,]) ~ samples$tumorcontent_per_sample)))))
colnames(tpm.gene.res) <- rownames(samples)
rownames(tpm.gene.res) <- rownames(tpm.gene)
save(tpm.gene.res, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.meidan0.content.res.RData")))

tpm.gene.log2.res <- t(mapply(x = 1:nrow(tpm.gene.log2), function(x) as.numeric(resid(lm(as.numeric(tpm.gene.log2[x,]) ~ samples$tumorcontent_per_sample)))))
colnames(tpm.gene.log2.res) <- rownames(samples)
rownames(tpm.gene.log2.res) <- rownames(tpm.gene.log2)
save(tpm.gene.log2.res, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.median0.log2.content.res.RData")))

# -----------------------------------------------------------------------------
# Wilcoxon rank sum test (non-parametric; n=45, 22 X vs 23 B) RES
# Last Modified: 30/03/21; 22/08/20
# -----------------------------------------------------------------------------
## Test: Wilcoxon/Mann–Whitney/U/wilcox.test
##       Student's/t.test
## FDR : Q/BH
## DE  : TR (1) vs UN (0) as factor
argv      <- data.frame(predictor="GROUP_ID2", predictor.wt=0, test="Wilcoxon", test.fdr="Q", stringsAsFactors=F)
file.name <- paste0("de_", base, "_tpm-gene-median0-content-res_B-vs-X_wilcox_q_n44")
file.main <- paste0("B (n=23) vs X (n=21) in ", BASE)

de <- differentialAnalysis(tpm.gene.log2.res, samples, argv$predictor, argv$predictor.wt, argv$test, argv$test.fdr)

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)

save(de.tpm.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Volcano plots of RB1-loss DE genes in LCNEC
# Figure(s)    : Figure S1 (A)
# Last Modified: 07/01/19
# -----------------------------------------------------------------------------
plot.main <- "Correction for tumour content (Linear regression)"
plot.de <- file.path(wd.de.plots, "volcanoplot_DE_EAC_median0_content_residual_B-vs-X_p1e-3")

genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c(plot.main, "Before (B) vs. After treatment (X)")
file.de <- paste0(plot.de, ".pdf")
de.tpm.gene$LOG2_FC <- de.tpm.gene$GROUP_ID2 - de.tpm.gene$GROUP_ID2_WT
plotVolcano(de.tpm.gene, 1.00E-03, genes, file.de, file.main, ymax=5.5)

# -----------------------------------------------------------------------------
# Replace Ensembl Gene IDs to gene name in Reactome results (Up- and Down-regulation)
# -----------------------------------------------------------------------------
wd.de.pathway <- file.path(wd.de, "pathway")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_B-vs-X_content-res_p1e-3_up")
#wd.de.reactome <- file.path(wd.de.pathway, "reactome_B-vs-X_content-res_p1e-3_down")

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
wd.de.reactome <- file.path(wd.de.pathway, "reactome_B-vs-X_content-res_p1e-3_up")
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")
#reactome[5, 2] <- "FOXO-mediated oxidative stress, metabolic and neuronal genes"

reactome.up <- subset(reactome, Entities.pValue <= 1e-2)[, 2:7]
reactome.up$log10P <- -log10(reactome.up$Entities.pValue)
reactome.up <- reactome.up[order(reactome.up$log10P),]

main.text <- c("Up-regulated after treatment", "12 genes")
file.de <- file.path(wd.de.reactome, "genes_B-vs-X_content-res_p1e-3_n12_up.pdf")

pdf(file.de, height=2.5, width=7.5)
par(mar=c(4,25,4,2))   # increase y-axis margin.
barplot(reactome.up$log10P, main=main.text[1], las=1, horiz=T, xlim=c(0, 4), xaxt="n", names.arg=reactome.up$Pathway.name, col=purple, xlab="-log10(p-value)", width=1)   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
abline(v=2, lty=5)

axis(side=1, at=seq(0, 8, by=2))
mtext(main.text[2], line=0.3)
dev.off()

##
wd.de.reactome <- file.path(wd.de.pathway, "reactome_B-vs-X_content-res_p1e-3_down")
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")
#reactome[2, 2] <- "Transcription of E2F targets under negative control by DREAM"
#reactome[3, 2] <- "HDR through Homologous Recombination or Single Strand Annealing"

reactome.down <- subset(reactome, Entities.pValue <= 1e-2)[, 2:7]
reactome.down$log10P <- -log10(reactome.down$Entities.pValue)
reactome.down <- reactome.down[order(reactome.down$log10P),]
reactome.down$log10P <- reactome.down$log10P * -1

main.text <- c("Down-regulated after treatment", "8 genes")
file.de <- file.path(wd.de.reactome, "genes_B-vs-X_content-res_p1e-3_n8_down.pdf")

pdf(file.de, height=3.5, width=7.5)
par(mar=c(4,2,4,25))   # increase y-axis margin.
posbar <- barplot(reactome.down$log10P, main=main.text[1], las=1, horiz=T, xlim=c(-4, 0), xaxt="n", names.arg="", col=red, xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=reactome.down$Pathway.name, xpd=NA)
abline(v=-2, lty=5)

axis(side=1, at=seq(-8, 0, by=2), labels=c(8, 6, 4, 2, 0))
mtext(main.text[2], line=0.3)
dev.off()
