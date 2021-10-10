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

samples  <- readTable(file.path(wd.rna, "esad_3rna_n68.txt2"), header=T, rownames=T, sep="\t")
purities <- readTable(file.path(wd.meta, "EAD-pupl.txt"), header=T, rownames=T, sep="")

#load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
#tpm.gene.log2 <- log2(tpm.gene + 0.01)   ## Use pseudocount=0.01

#samples <- samples[rownames(subset(samples, GROUP_ID != 1)),]
#samples$GROUP_ID <- samples$GROUP_ID - 2
#tpm.gene.log2 <- tpm.gene.log2[, rownames(samples)]

# -----------------------------------------------------------------------------
# PCA vs. Purites
# Last Modified: 02/06/21
# -----------------------------------------------------------------------------
rownames(samples) <- samples$PATIENT_ID2
overlaps <- intersect(samples$PATIENT_ID2, purities$sample_name)
samples <- samples[overlaps,]
samples$purity <- purities[overlaps,]$purity
rownames(samples) <- samples$SAMPLE_ID

scores <- pcaScores(pca.de)
scores <- scores[rownames(samples),]

# -----------------------------------------------------------------------------
# 
# Last Modified: 03/06/21
# -----------------------------------------------------------------------------
plotFACS3 <- function(p, pc, file.name, main.text, xlab.text, ylab.text, col, col2) {
   pdf(paste0(file.name, ".pdf"), height=6, width=6)
   plot(p ~ pc, ylab="", xlab=xlab.text, main=main.text[1], col="white", pch=19, cex=2, lwd=0, cex.axis=1.5, cex.lab=1.6, cex.main=1.7)
   points(pc[ 1:10], p[ 1:10], col=col2[1], pch=19, cex=1.7)
   points(pc[11:20], p[11:20], col=col2[2], pch=19, cex=1.7)
   
   lm.fit <- lm(p ~ pc)
   abline(lm.fit, col=col, lwd=4)
 
   cor3 <- cor.test(p, pc, method="spearman", exact=F)
   legend("bottomright", c(paste0("rho = ", round0(cor3[[4]], digits=2)), paste0("p-value = ", scientific(cor3[[3]]))), text.col=cols, text.font=2, bty="n", cex=1.5)
 
   #axis(side=2, at=seq(-0.2, 0.2, by=0.2), labels=c(-0.2, 0, 0.2), cex.axis=1.5)
   #axis(side=2, at=seq(-0.3, 0.1, by=0.2), labels=c(-0.3, -0.1, 0.1), cex.axis=1.5)
   mtext(ylab.text, side=2, line=2.75, cex=1.6)
   #mtext(main.text[2], line=0.3, cex=1.6)
   dev.off()
}

file.name <- file.path(wd.de.plots, "ESAD_Purity-vs-PC4_n20")
main.text <- c(paste("Purity vs. PC4 (n=20)"), "")
xlab.text <- paste0("PC", 4, " (", pcaProportionofVariance(pca.de, 4), "%)")
ylab.text <- "Purity"                                                                         ## "#619CFF", "#F8766D", "#00BA38"      "skyblue3", "lightcoral", "#59a523"
cols <- "black"
cols2 <- c(red, purple)
plotFACS3(samples$purity, scores$PC4, file.name, main.text, xlab.text, ylab.text, cols, cols2)

file.name <- file.path(wd.de.plots, "ESAD_Purity-vs-tumorcontent_n20")
main.text <- c(paste("Purity vs. Tumour content (n=20)"), "")
xlab.text <- "Tumour content"
ylab.text <- "Purity"                                                                         ## "#619CFF", "#F8766D", "#00BA38"      "skyblue3", "lightcoral", "#59a523"
cols <- "black"
cols2 <- c(red, purple)
plotFACS3(samples$purity, as.numeric(sub("%", "", samples$tumorcontent_per_sample))/100, file.name, main.text, xlab.text, ylab.text, cols, cols2)





# -----------------------------------------------------------------------------
# Wilcoxon rank sum test (non-parametric; n=45, 22 TR vs 23 UN)
# Last Modified: 22/08/20
# -----------------------------------------------------------------------------
## Test: Wilcoxon/Mann–Whitney/U/wilcox.test
##       Student's/t.test
## FDR : Q/BH
## DE  : TR (1) vs UN (0) as factor
argv      <- data.frame(predictor="GROUP_ID", predictor.wt=0, test="Wilcoxon", test.fdr="Q", stringsAsFactors=F)
file.name <- paste0("de_", base, "_tpm-gene-median0-median0_B-vs-N_wilcox_q_n46")
file.main <- paste0("TR (n=22) vs UN (n=23) in ", BASE)

de <- differentialAnalysis(tpm.gene.log2, samples, argv$predictor, argv$predictor.wt, argv$test, argv$test.fdr)

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)

save(de.tpm.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# PCA
# Last Modified: 26/03/21; 23/08/20
# -----------------------------------------------------------------------------
test <- tpm.gene[, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))
save(pca.de, file=file.path(wd.de.data, paste0("PCA_EAC_N-B-X.RData")))

trait <- as.numeric(samples[, "GROUP_ID2"])
trait[which(trait == 0)] <- "N"
trait[which(trait == 1)] <- "B"
trait[which(trait == 2)] <- "X"

file.main <- c("EAC samples (n=68)", "")
plotPCA(1, 2, pca.de, trait, wd.de.plots, "PCA_EAC_N-B-X_n68_EAC", size=6.5, file.main, "topright", c(green, red, purple), NULL, flip.x=-1, flip.y=1, legend.title=NA)

# -----------------------------------------------------------------------------
# PCA (Controlled for PC1)
# Last Modified: 26/03/21; 23/01/21
# -----------------------------------------------------------------------------
test <- tpm.gene.res[, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de.res <- getPCA(t(test))
save(pca.de.res, file=file.path(wd.de.data, paste0("PCA_EAC_N-B-X_RES.RData")))

trait <- as.numeric(samples[, "GROUP_ID2"])
trait[which(trait == 0)] <- "N"
trait[which(trait == 1)] <- "B"
trait[which(trait == 2)] <- "X"

file.main <- c("EAC samples (Controlled for PC1)", "")
plotPCA(1, 2, pca.de.res, trait, wd.de.plots, "PCA_EAC_B-X-N_n68_RES", size=6.5, file.main, "topright", c(green, red, purple), NULL, flip.x=-1, flip.y=1, legend.title=NA)

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
   ymax <- 9.5
 
   pdf(file.de, height=6, width=6)
   plot(de$LOG2_FC, de$log10P, pch=16, xlim=c(-xmax, xmax), ylim=c(0, ymax), xlab="EAC B versus N [log2 fold change]", ylab="Significance [-log10(p-value)]", col="lightgray", main=file.main[1])

   #text(xmax*-1 + 2*xmax/28, -log10(pvalue) + ymax/42, paste0("FDR=", fdr, "%"), cex=0.85)
   #text(xmax*-1 + 2*xmax/35, -log10(pvalue) + ymax/42, "FDR=0.05", cex=0.85)
   #abline(h=c(-log10(fdrToP(0.1, de))), lty=5, col="darkgray")
   #text(xmax*-1 + 2*xmax/50, -log10(fdrToP(0.1, de)) + ymax/42, "FDR=0.1", col="darkgray", cex=0.85)

   de.up   <- subset(de.sig, LOG2_FC > 0)
   points(de.up$LOG2_FC, de.up$log10P, pch=16, col=red)
   de.down <- subset(de.sig, LOG2_FC < 0)
   points(de.down$LOG2_FC, de.down$log10P, pch=16, col=green)
 
   abline(v=c(-log2(2), log2(2)), lty=5, col="darkgray")
   abline(h=c(-log10(pvalue)), lty=5, col="black")
   
   if(nrow(genes) != 0) {
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
   }
   
   mtext(file.main[2], cex=1.2, line=0.3)
   legend("topleft", legend=c("Up-regulated (N < B)", "Down-regulated (N > B)"), col=c(red, green), pch=19)
   dev.off()
}

##
plot.main <- "1,240 differentially expressed genes in EAC"
plot.de <- file.path(wd.de.plots, "volcanoplot_esad_median0_N-vs-B_p1e-6_CNA")

genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c(plot.main, "Normal (N) vs. Tumour (B)")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene, 1.00E-06, genes, file.de, file.main)

##
plot.main <- "1,240 differentially expressed genes in EAC"
plot.de <- file.path(wd.de.plots, "volcanoplot_esad_median0_N-vs-B_p1e-6")

genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c(plot.main, "Normal (N) vs. Tumour (B)")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene, 1.00E-06, genes, file.de, file.main)

# -----------------------------------------------------------------------------
# Replace Ensembl Gene IDs to gene name in Reactome results (Up- and Down-regulation)
# -----------------------------------------------------------------------------
wd.de.pathway <- file.path(wd.de, "pathway")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_p1e-6_up")
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
wd.de.reactome <- file.path(wd.de.pathway, "reactome_B-vs-N_p1e-6_up")
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
wd.de.reactome <- file.path(wd.de.pathway, "reactome_B-vs-N_p1e-6_down")
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")

reactome.down <- subset(reactome, Entities.pValue <= 1e-6)[, 2:7]
reactome.down$log10P <- -log10(reactome.down$Entities.pValue)
reactome.down <- reactome.down[order(reactome.down$log10P),]
reactome.down$log10P <- reactome.down$log10P * -1

main.text <- c("Down-regulated pathways in EAC", "477 genes")
file.de <- file.path(wd.de.reactome, "genes_T-vs-N_p1e-6_n477_down.pdf")

pdf(file.de, height=2.1, width=7.5)
par(mar=c(4,3,4,18))   # increase y-axis margin.
posbar <- barplot(reactome.down$log10P, main=main.text[1], las=1, horiz=T, xlim=c(-16, 0), xaxt="n", names.arg="", col=green, xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=reactome.down$Pathway.name, xpd=NA)
abline(v=-6, lty=5)

axis(side=1, at=seq(-16, 0, by=2), labels=c(16, 14, 12, 10, 8, 6, 4, 2, 0))
mtext(main.text[2], line=0.3)
dev.off()
