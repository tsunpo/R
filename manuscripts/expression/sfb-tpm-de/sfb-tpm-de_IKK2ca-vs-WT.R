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
samples <- samples[rownames(subset(samples, GROUP_NAME %in% c("IKK2ca", "WT"))),]

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
tpm.gene.log2 <- log2(tpm.gene + 1)   ## Use pseudocount=1
tpm.gene.log2 <- tpm.gene.log2[, rownames(samples)]
dim(tpm.gene.log2)
# [1] 15381

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
file.name <- paste0("DE_SFB_tpm-gene-median0_IKK2ca-vs-WT_wilcox_q_n12")
load(file=file.path(wd.de.data, paste0(file.name, ".RData")))

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
   legend("topleft", legend=c("Up-regulated (WT < NEMO)", "Down-regulated (WT > NEMO)"), col=c(red, blue), pch=19)
   dev.off()
}

##
plot.main <- "1,240 differentially expressed genes in EAC"
plot.de <- file.path(wd.de.plots, "volcanoplot_DE_EAC_median0_N-vs-B_p1e-6_CNA")

genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c(plot.main, "Normal (N) vs. Tumour (B)")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene, 1.00E-06, genes, file.de, file.main)

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

