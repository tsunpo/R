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

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.median0.RData")))
tpm.gene.log2 <- log2(tpm.gene + 0.01)   ## Use pseudocount=0.01

t  <- rownames(subset(samples, GROUP_ID == 0))
tr <- rownames(subset(samples, GROUP_ID == 1))
n  <- rownames(subset(samples, GROUP_ID == 2))
samples[n,]$GROUP_ID  <- 0
samples[tr,]$GROUP_ID <- 1
samples[t,]$GROUP_ID  <- 2
samples$GROUP_ID <- as.numeric(samples$GROUP_ID)

tpm.gene.log2 <- tpm.gene.log2[, rownames(samples)]

# -----------------------------------------------------------------------------
# Spearman's rank correlation
# Last Modified: 22/08/20
# -----------------------------------------------------------------------------
colnames <- c("RHO", "P", "FDR", "N", "TR", "T", "FC_TR_N", "FC_T_N")
de <- toTable(0, length(colnames), nrow(tpm.gene.log2), colnames)
rownames(de) <- rownames(tpm.gene.log2)

## SRC
de$RHO <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples$GROUP_ID, method="spearman", exact=F)[[4]])
de$P   <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples$GROUP_ID, method="spearman", exact=F)[[3]])
## ANOVA
#samples$Cancer_Type <- as.factor(pheno.all$Cancer_Type)
#de$ANOVA_P <- mapply(x = 1:nrow(expr.pheno.log2), function(x) testANOVA(x, expr.pheno.log2, pheno.all))

## Log2 fold change
de$N  <- median00(tpm.gene.log2, rownames(subset(samples, GROUP_ID == 0)))
de$TR <- median00(tpm.gene.log2, rownames(subset(samples, GROUP_ID == 1)))
de$T  <- median00(tpm.gene.log2, rownames(subset(samples, GROUP_ID == 2)))
de$FC_TR_N <- de$TR - de$N
de$FC_T_N  <- de$T - de$N

## FDR
library(qvalue)
de$FDR   <- qvalue(de$P)$qvalue
#de$ANOVA_Q <- qvalue(de$ANOVA_P)$qvalue
de <- de[order(de$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)   ## BE EXTRA CAREFUL!!

save(de.tpm.gene, samples, file=file.path(wd.de.data, "src_esad_tpm-gene-median0-median0_T-TR-N_src_q_n68.RData"))
writeTable(de.tpm.gene, file.path(wd.de.data, "src_esad_tpm-gene-median0-median0_T-TR-N_src_q_n68.txt"), colnames=T, rownames=F, sep="\t")

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
