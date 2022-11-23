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
wd.de.plots <- file.path(wd.de, "plots")

samples <- readTable(file.path(wd.rna, "sfb_3rna_n41.txt"), header=T, rownames=T, sep="\t")
#samples <- samples[rownames(subset(samples, GROUP_NAME %in% c("TNFR1KO", "WT"))),]

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
tpm.gene.log2 <- log2(tpm.gene + 1)   ## Use pseudocount=1
tpm.gene.log2 <- tpm.gene.log2[, rownames(samples)]
dim(tpm.gene.log2)
# [1] 15381

# -----------------------------------------------------------------------------
# PCA
# Last Modified: 10/10/21
# -----------------------------------------------------------------------------
samples <- samples[rownames(subset(samples, GROUP_NAME %in% c("TNFR1KO", "TNFR1LEC-KO", "WT"))),]

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
test <- tpm.gene[, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))

trait <- samples$GROUP_NAME
trait.v <- c("TNFR1KO", "TNFR1LEC-KO", "WT")
cols    <- c(blue, blue.lighter, grey)
file.main <- c("TNFR1KO + TNFR1LEC-KO + WT", "")
plotPCA(1, 2, pca.de, trait, wd.de.plots, "PCA_SFB_TNFR1KO+TNFR1LEC-KO+WT_n12_median0", size=6, file.main, "topright", trait.v, cols, flip.x=1, flip.y=1)

# -----------------------------------------------------------------------------
# Wilcoxon rank sum test (non-parametric; n=45, 22 TR vs 23 UN)
# Last Modified: 22/08/20
# -----------------------------------------------------------------------------
## Test: Wilcoxon/Mann–Whitney/U/wilcox.test
##       Student's/t.test
## FDR : Q/BH
## DE  : C (0) vs MA (1) as factor
argv      <- data.frame(predictor="GROUP_ID", predictor.wt=0, test="Wilcoxon", test.fdr="Q", stringsAsFactors=F)
file.name <- paste0("DE_SFB_tpm-gene-median0_TNFR1KO-vs-WT_wilcox_q_n11")
file.main <- paste0("TNFR1KO (n=5) vs. WT (n=6) in ", BASE)

de <- differentialAnalysis(tpm.gene.log2, samples, argv$predictor, argv$predictor.wt, argv$test, argv$test.fdr)

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)

save(de.tpm.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

# > de.tpm.gene[1:10,]
# ensembl_gene_id external_gene_name chromosome_name strand start_position end_position                     gene_biotype           P       FDR
# ENSMUSG00000069609 ENSMUSG00000069609           Cd300ld4           chr11     -1      115020728    115027012                   protein_coding 0.006735947 0.8566251
# ENSMUSG00000071893 ENSMUSG00000071893             Vmn1r4            chr6      1       56924015     56958100                   protein_coding 0.006735947 0.8566251
# ENSMUSG00000089722 ENSMUSG00000089722           Cd300ld5           chr11     -1      115031486    115037769                   protein_coding 0.006735947 0.8566251
# ENSMUSG00000115512 ENSMUSG00000115512          Rpl17-ps7           chr14      1      103010246    103010786             processed_pseudogene 0.006735947 0.8566251
# ENSMUSG00000074704 ENSMUSG00000074704             Rad21l            chr2     -1      151645404    151668533                   protein_coding 0.007546234 0.8566251
# ENSMUSG00000020234 ENSMUSG00000020234      4930404N11Rik           chr10     -1       81363067     81365816                   protein_coding 0.007969413 0.8566251
# ENSMUSG00000039095 ENSMUSG00000039095                En2            chr5      1       28165694     28172166                   protein_coding 0.007969413 0.8566251
# ENSMUSG00000047237 ENSMUSG00000047237             Fbxw21            chr9     -1      109139447    109162041                   protein_coding 0.007969413 0.8566251
# ENSMUSG00000117492 ENSMUSG00000117492            Gm19052           chr18      1       10217043     10238921 transcribed_processed_pseudogene 0.007969413 0.8566251
# ENSMUSG00000000378 ENSMUSG00000000378               Ccm2           chr11      1        6546887      6596744                   protein_coding 0.008113117 0.8566251

file.name <- paste0("DE_SFB_tpm-gene-median0_TNFR1KO-vs-WT_wilcox_q_n11")
load(file=file.path(wd.de.data, paste0(file.name, ".RData")))
nrow(subset(de.tpm.gene, P < 0.01))
nrow(subset(de.tpm.gene, P < 0.05))

de.tpm.gene[77, "FDR"]
de.tpm.gene[577, "FDR"]

# -----------------------------------------------------------------------------
# Volcano plots of RB1-loss DE genes in LCNEC
# Figure(s)    : Figure S1 (A)
# Last Modified: 07/01/19
# -----------------------------------------------------------------------------
load("/Users/tpyang/Work/uni-koeln/tyang2/SFB/analysis/expression/kallisto/sfb-tpm-de/data/DE_SFB_tpm-gene-median0_TNFR1KO-vs-WT_wilcox_q_n11.RData")

de.sig.pos <- subset(subset(de.tpm.gene, P <= 0.01), LOG2_FC >= log2(2))
de.sig.neg <- subset(subset(de.tpm.gene, P <= 0.01), LOG2_FC <= -log2(2))
de.sig <- rbind(de.sig.pos, de.sig.neg)

colnames <- c("GENE", "ADJ_1", "ADJ_2")
genes0 <- de.sig$external_gene_name
genes <- toTable(NA, length(colnames), length(genes0), colnames)
genes$GENE <- genes0

xlab.text <- expression("Fold change " * "[log" * ""[2] * "]")
ylab.text <- expression("Significance " * "[-log" * ""[10] * "(" * italic("P") * ")]")

plot.de <- file.path(wd.de.plots, "volcanoplot_DE_SFB_median0_TNFR1KO-vs-WT_p0.01_fc1_log2_GENE")
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("TNFR1KO vs. WT", "")
plotVolcano(de.tpm.gene, 0.01, genes, file.de, file.main, xlab.text, ylab.text, "topright", c("TNFR1KO > WT", "TNFR1KO < WT"), c(red, blue), c(red, blue), fold=1)

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

