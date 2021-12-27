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
samples[x,]$GROUP_ID3 <- 0
samples[n,]$GROUP_ID3 <- 1
samples[b,]$GROUP_ID3 <- 2
samples$GROUP_ID3 <- as.numeric(samples$GROUP_ID3)

tpm.gene.log2 <- tpm.gene.log2[, rownames(samples)]
nrow(tpm.gene.log2)
# [1] 15658

# -----------------------------------------------------------------------------
# Spearman's rank correlation
# Last Modified: 22/08/20
# -----------------------------------------------------------------------------
colnames <- c("RHO", "P", "FDR", "X", "N", "B", "FC_N_X", "FC_B_X")
de <- toTable(0, length(colnames), nrow(tpm.gene.log2), colnames)
rownames(de) <- rownames(tpm.gene.log2)

## SRC
de$RHO <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples$GROUP_ID3, method="spearman", exact=F)[[4]])
de$P   <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples$GROUP_ID3, method="spearman", exact=F)[[3]])
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

save(de.tpm.gene, samples, file=file.path(wd.de.data, "src_esad_tpm-gene-median0_X-N-B_src_q_n68.RData"))
writeTable(de.tpm.gene, file.path(wd.de.data, "src_esad_tpm-gene-median0_X-N-B_src_q_n68.txt"), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
plotBox3 <- function(wd.de.plots, file.name, tpm.1, tpm.2, tpm.3, main, names, cols, ylim, height=5, width=3.5, ylab.txt="log2(TPM + 0.01)") {
   trait <- rep(0, length(tpm.1))
   trait <- c(trait, rep(1, length(tpm.2)))
   trait <- c(trait, rep(2, length(tpm.3)))
   trait <- as.factor(trait)
   expr <- as.numeric(c(tpm.1, tpm.2, tpm.3))
 
   pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=height, width=width)
   boxplot(expr ~ trait, outline=T, xaxt="n", ylab=ylab.txt, main=main, boxcol=cols, whiskcol=cols, outcol=cols, medcol=cols, staplecol=cols, ylim=ylim, cex.axis=1.2, cex.lab=1.25, cex.main=1.3)

   axis(side=1, at=seq(1, 3, by=1), labels=names, font=2, cex.axis=1.25)
   #axis(side=1, at=2, labels="Total n=30,978", line=1.3, col=NA, cex.axis=1.2)
   axis(side=1, at=1, labels=paste0("n=", format(length(tpm.1), big.mark=",", scientific=F)), line=1.2, col=NA, cex.axis=1.25)
   axis(side=1, at=2, labels=paste0("n=", format(length(tpm.2), big.mark=",", scientific=F)), line=1.2, col=NA, cex.axis=1.25)
   axis(side=1, at=3, labels=paste0("n=", format(length(tpm.3), big.mark=",", scientific=F)), line=1.2, col=NA, cex.axis=1.25)

   #mtext(paste0("p-value = ", scientific(wilcox.test(expr ~ trait, exact=F)$p.value)), cex=1.2, line=0.3)
   dev.off()
}

###
##
n <- tpm.gene.log2["ENSG00000087245", rownames(subset(samples, GROUP_ID2 == 0))]
b <- tpm.gene.log2["ENSG00000087245", rownames(subset(samples, GROUP_ID2 == 1))]
x <- tpm.gene.log2["ENSG00000087245", rownames(subset(samples, GROUP_ID2 == 2))]
ylim <- c(min(tpm.gene.log2["ENSG00000087245",]), max(tpm.gene.log2["ENSG00000087245",]))

file.name <- paste0("boxplot_EAC_tpm.gene.median0_MMP2")
plotBox3(wd.de.plots, file.name, n, b, x, main="MMP2 (ENSG00000087245)", names=c("N", "B", "X"), cols=c(green, red, purple), ylim=ylim)

###
##
n <- tpm.gene.res["ENSG00000087245", rownames(subset(samples, GROUP_ID2 == 0))]
b <- tpm.gene.res["ENSG00000087245", rownames(subset(samples, GROUP_ID2 == 1))]
x <- tpm.gene.res["ENSG00000087245", rownames(subset(samples, GROUP_ID2 == 2))]
ylim <- c(min(tpm.gene.res["ENSG00000087245",]), max(tpm.gene.res["ENSG00000087245",]))

file.name <- paste0("boxplot_EAC_tpm.gene.median0_MMP2_res")
plotBox3(wd.de.plots, file.name, n, b, x, main="MMP2 (ENSG00000087245)", names=c("N", "B", "X"), cols=c(green, red, purple), ylim, ylab.txt="Residual of expression (controlled for PC1)")

###
##
n <- tpm.gene.log2.res["ENSG00000087245", rownames(subset(samples, GROUP_ID2 == 0))]
b <- tpm.gene.log2.res["ENSG00000087245", rownames(subset(samples, GROUP_ID2 == 1))]
x <- tpm.gene.log2.res["ENSG00000087245", rownames(subset(samples, GROUP_ID2 == 2))]
ylim <- c(min(tpm.gene.log2.res["ENSG00000087245",]), max(tpm.gene.log2.res["ENSG00000087245",]))

file.name <- paste0("boxplot_EAC_tpm.gene.median0_MMP2_log2.res")
plotBox3(wd.de.plots, file.name, n, b, x, main="MMP2 (ENSG00000087245)", names=c("N", "B", "X"), cols=c(green, red, purple), ylim, ylab.txt="Residual of expression (controlled for PC1)")


















##
n  <- tpm.gene.log2["ENSG00000133110", rownames(subset(samples, GROUP_ID == 0))]
t  <- tpm.gene.log2["ENSG00000133110", rownames(subset(samples, GROUP_ID == 1))]
tr <- tpm.gene.log2["ENSG00000133110", rownames(subset(samples, GROUP_ID == 2))]

file.name <- paste0("boxplot_esad_tpm.gene.median_POSTN")
plotBox3(wd.de.plots, file.name, n, t, tr, main="POSTN (ENSG00000133110)", names=c("N", "T", "TR"), cols=c("dodgerblue", "red", "#01DF01"), ylab.txt="log2(TPM + 0.01)")

n  <- tpm.gene.res["ENSG00000133110", rownames(subset(samples, GROUP_ID == 0))]
t  <- tpm.gene.res["ENSG00000133110", rownames(subset(samples, GROUP_ID == 1))]
tr <- tpm.gene.res["ENSG00000133110", rownames(subset(samples, GROUP_ID == 2))]

file.name <- paste0("boxplot_esad_tpm.gene.median_POSTN_RES")
plotBox3(wd.de.plots, file.name, n, t, tr, main="POSTN (ENSG00000133110)", names=c("N", "T", "TR"), cols=c("dodgerblue", "red", "#01DF01"), ylab.txt="Residual of expression (controlled for PC1)")

##
n  <- tpm.gene.log2["ENSG00000154096", rownames(subset(samples, GROUP_ID == 0))]
t  <- tpm.gene.log2["ENSG00000154096", rownames(subset(samples, GROUP_ID == 1))]
tr <- tpm.gene.log2["ENSG00000154096", rownames(subset(samples, GROUP_ID == 2))]

file.name <- paste0("boxplot_esad_tpm.gene.median_THY1")
plotBox3(wd.de.plots, file.name, n, t, tr, main="THY1 (ENSG00000154096)", names=c("N", "T", "TR"), cols=c("dodgerblue", "red", "#01DF01"), ylab.txt="log2(TPM + 0.01)")

n  <- tpm.gene.res["ENSG00000154096", rownames(subset(samples, GROUP_ID == 0))]
t  <- tpm.gene.res["ENSG00000154096", rownames(subset(samples, GROUP_ID == 1))]
tr <- tpm.gene.res["ENSG00000154096", rownames(subset(samples, GROUP_ID == 2))]

file.name <- paste0("boxplot_esad_tpm.gene.median_THY1_RES")
plotBox3(wd.de.plots, file.name, n, t, tr, main="THY1 (ENSG00000154096)", names=c("N", "T", "TR"), cols=c("dodgerblue", "red", "#01DF01"), ylab.txt="Residual of expression (controlled for PC1)")



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
   ymax <- max(de$log10P)
 
   pdf(file.de, height=7, width=7)
   plot(de$LOG2_FC, de$log10P, pch=16, xlim=c(-xmax, xmax), ylim=c(0, ymax), xlab="log2FC(TR/N)", ylab="-log10(p-value)", col="darkgray", main=file.main[1])

   abline(h=c(-log10(pvalue)), lty=5)
   #text(xmax*-1 + 2*xmax/28, -log10(pvalue) + ymax/42, paste0("FDR=", fdr, "%"), cex=0.85)
   #text(xmax*-1 + 2*xmax/35, -log10(pvalue) + ymax/42, "FDR=0.05", cex=0.85)
   #abline(h=c(-log10(fdrToP(0.1, de))), lty=5, col="darkgray")
   #text(xmax*-1 + 2*xmax/50, -log10(fdrToP(0.1, de)) + ymax/42, "FDR=0.1", col="darkgray", cex=0.85)

   de.up   <- subset(de.sig, LOG2_FC > 0)
   points(de.up$LOG2_FC, de.up$log10P, pch=16, col="#01DF01")
   de.down <- subset(de.sig, LOG2_FC < 0)
   points(de.down$LOG2_FC, de.down$log10P, pch=16, col="dodgerblue")
 
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
   legend("topleft", legend=c("Upregulated in TR", "Downregulated in TR"), col=c("#01DF01", "dodgerblue"), pch=19)
   dev.off()
}

##
plot.main <- "374 additively expressed genes after treatment"
plot.de <- file.path(wd.de.plots, "volcanoplot_esad_median0_median0_TR-T-N_p1e-6")
de.tpm.gene$LOG2_FC <- de.tpm.gene$FC_TR_N

genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c(plot.main, "Between TR, T and N")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene, 1.00E-06, genes, file.de, file.main)

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
wd.de.reactome <- file.path(wd.de.pathway, "reactome_TR-T-N_p1e-6_up")
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")
reactome[3, 2] <- "Assembly of collagen fibrils and multimeric structures"

reactome.up <- subset(reactome, Entities.pValue <= 1e-6)[, 2:7]
reactome.up$log10P <- -log10(reactome.up$Entities.pValue)
reactome.up <- reactome.up[order(reactome.up$log10P),]

main.text <- c("Up-regulated pathways after treatment", "158 genes")
file.de <- file.path(wd.de.reactome, "genes_TR-T-N_p1e-6_n158_up.pdf")

pdf(file.de, height=1.85, width=7.5)
par(mar=c(4,20.7,4,1.3))   # increase y-axis margin.
barplot(reactome.up$log10P, main=main.text[1], las=1, horiz=T, xlim=c(0, 12), xaxt="n", names.arg=reactome.up$Pathway.name, col="#01DF01", xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3

axis(side=1, at=seq(0, 12, by=2))
mtext(main.text[2], line=0.3)
dev.off()

##
wd.de.reactome <- file.path(wd.de.pathway, "reactome_TR-T-N_p1e-6_down")
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")
reactome[3, 2] <- "TCA cycle and respiratory electron transport"

reactome.down <- subset(reactome, Entities.pValue <= 1e-6)[, 2:7]
reactome.down$log10P <- -log10(reactome.down$Entities.pValue)
reactome.down <- reactome.down[order(reactome.down$log10P),]
reactome.down$log10P <- reactome.down$log10P * -1

main.text <- c("Down-regulated pathways after treatment", "216 genes")
file.de <- file.path(wd.de.reactome, "genes_TR-T-N_p1e-4_n216_down.pdf")

pdf(file.de, height=2.1, width=7.5)
par(mar=c(4,2,4,20))   # increase y-axis margin.
posbar <- barplot(reactome.down$log10P, main=main.text[1], las=1, horiz=T, xlim=c(-12, 0), xaxt="n", names.arg="", col="dodgerblue", xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=reactome.down$Pathway.name, xpd=NA)

axis(side=1, at=seq(-12, 0, by=2), labels=c(12, 10, 8, 6, 4, 2, 0))
mtext(main.text[2], line=0.3)
dev.off()
