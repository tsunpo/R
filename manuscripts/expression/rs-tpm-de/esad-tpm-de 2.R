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

wd.de    <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

samples <- readTable(file.path(wd.rna, "esad_3rna_n68.txt"), header=T, rownames=T, sep="\t")
#for (s in 1:nrow(samples)) {
#   if (as.numeric(gsub("S", "", samples[s,]$SAMPLE_ID)) != samples[s,]$CCG_ID) print(paste0("HERE: ", s))
#   if (as.numeric(gsub("P1", "", samples[s,]$PATIENT_ID)) != samples[s,]$Patient) print(paste0("HERE: ", s))
#   if (as.numeric(unlist(strsplit(samples[s,]$FILE_NAME, "_"))[2]) != samples[s,]$CCG_ID) print(paste0("HERE: ", s))
#}

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
tpm.gene.log2 <- log2(tpm.gene + 0.01)   ## Use pseudocount=0.01

samples <- samples[rownames(subset(samples, GROUP_ID != 2)),]
tpm.gene.log2 <- tpm.gene.log2[, rownames(samples)]

###
##
tpm.gene.un <- tpm.gene[, rownames(subset(samples, GROUP_ID == 0))]
tpm.gene.un <- removeMedian0(tpm.gene.un)
nrow(tpm.gene.un)
# [1] 156

tpm.gene.tr <- tpm.gene[, rownames(subset(samples, GROUP_ID == 1))]
tpm.gene.tr <- removeMedian0(tpm.gene.tr)
nrow(tpm.gene.tr)
# [1] 412

tpm.gene.n <- tpm.gene[, rownames(subset(samples, GROUP_ID == 2))]
tpm.gene.n <- removeMedian0(tpm.gene.n)
nrow(tpm.gene.n)
# [1] 984

overlaps <- intersect(intersect(rownames(tpm.gene.un), rownames(tpm.gene.tr)), rownames(tpm.gene.n))
length(overlaps)
tpm.gene <- tpm.gene[overlaps,]
save(tpm.gene, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.median0.median0.RData")))
nrow(tpm.gene)
# [1] 14334

###
##
#load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.median0.RData")))
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
tpm.gene.log2 <- log2(tpm.gene + 0.01)   ## Use pseudocount=0.01

samples <- samples[rownames(subset(samples, GROUP_ID != 1)),]
samples$GROUP_ID <- samples$GROUP_ID - 2
tpm.gene.log2 <- tpm.gene.log2[, rownames(samples)]

# -----------------------------------------------------------------------------
# Wilcoxon rank sum test (non-parametric; n=45, 22 TR vs 23 UN)
# Last Modified: 22/08/20
# -----------------------------------------------------------------------------
## Test: Wilcoxon/Mann–Whitney/U/wilcox.test
##       Student's/t.test
## FDR : Q/BH
## DE  : TR (1) vs UN (0) as factor
argv      <- data.frame(predictor="GROUP_ID", predictor.wt=0, test="Wilcoxon", test.fdr="Q", stringsAsFactors=F)
file.name <- paste0("de_", base, "_tpm-gene-median0-median0_N-vs-UN_wilcox_q_n45")
file.main <- paste0("TR (n=22) vs UN (n=23) in ", BASE)

de <- differentialAnalysis(tpm.gene.log2, samples, argv$predictor, argv$predictor.wt, argv$test, argv$test.fdr)

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)

save(de.tpm.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# PCA
# Last Modified: 23/08/20
# -----------------------------------------------------------------------------
test <- tpm.gene[, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))

trait <- as.numeric(samples[, "GROUP_ID"])
trait[which(trait == 0)] <- "B"
trait[which(trait == 1)] <- "X"
trait[which(trait == 2)] <- "N"

file.main <- c("ESAD samples (n=68)", "")
plotPCA(1, 2, pca.de, trait, wd.de.plots, "PCA_ESAD_N-B-X_median0_n68", size=6.5, file.main, "topright", c("red", "#01DF01", "dodgerblue"), NA, flip.x=-1, flip.y=1, legend.title=NA)

##
test <- tpm.gene.res[, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.res <- getPCA(t(test))

trait <- as.numeric(samples[, "GROUP_ID"])
trait[which(trait == 0)] <- "B"
trait[which(trait == 1)] <- "X"
trait[which(trait == 2)] <- "N"

file.main <- c("ESAD samples (n=68)", "After controlling for cell cycle")
plotPCA(1, 2, pca.res, trait, wd.de.plots, "PCA_ESAD_TR-T-N_median0_RES_n68", size=6.5, file.main, "topleft", c("red", "#01DF01", "dodgerblue"), NA, flip.x=1, flip.y=1, legend.title=NA)

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
   plot(de$LOG2_FC, de$log10P, pch=16, xlim=c(-xmax, xmax), ylim=c(0, ymax), xlab="log2FC(TR/T)", ylab="-log10(p-value)", col="darkgray", main=file.main[1])

   abline(h=c(-log10(pvalue)), lty=5)
   #text(xmax*-1 + 2*xmax/28, -log10(pvalue) + ymax/42, paste0("FDR=", fdr, "%"), cex=0.85)
   #text(xmax*-1 + 2*xmax/35, -log10(pvalue) + ymax/42, "FDR=0.05", cex=0.85)
   #abline(h=c(-log10(fdrToP(0.1, de))), lty=5, col="darkgray")
   #text(xmax*-1 + 2*xmax/50, -log10(fdrToP(0.1, de)) + ymax/42, "FDR=0.1", col="darkgray", cex=0.85)

   de.up   <- subset(de.sig, LOG2_FC > 0)
   points(de.up$LOG2_FC, de.up$log10P, pch=16, col="#01DF01")
   de.down <- subset(de.sig, LOG2_FC < 0)
   points(de.down$LOG2_FC, de.down$log10P, pch=16, col="red")
 
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
   legend("topleft", legend=c("Upregulated in TR", "Downregulated in TR"), col=c("#01DF01", "red"), pch=19)
   dev.off()
}

##
plot.main <- "368 differentially expressed genes before and after treatment"
plot.de <- file.path(wd.de.plots, "volcanoplot_esad_median0_median0_TR-vs-T_p1e-6")

genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c(plot.main, "TR vs. T")
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
reactome.up <- subset(reactome, Entities.pValue <= 1e-6)[, 2:7]
reactome.up$log10P <- -log10(reactome.up$Entities.pValue)
reactome.up <- reactome.up[order(reactome.up$log10P),]

main.text <- c("Up-regulated pathways after treatment", "224 genes")
file.de <- file.path(wd.de.reactome, "genes_tr_p1e-6_n224_up.pdf")

pdf(file.de, height=3, width=7.5)
par(mar=c(4,18,4,3.1))   # increase y-axis margin.
barplot(reactome.up$log10P, main=main.text[1], las=1, horiz=T, xlim=c(0, 15), xaxt="n", names.arg=reactome.up$Pathway.name, col="#01DF01", xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3

axis(side=1, at=seq(0, 14, by=2))
mtext(main.text[2], line=0.3)
dev.off()

##
wd.de.reactome <- file.path(wd.de.pathway, "reactome_TR-vs-T_p1e-6_down")
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")

reactome.down <- subset(reactome, Entities.pValue <= 1e-6)[, 2:7]
reactome.down$log10P <- -log10(reactome.down$Entities.pValue)
reactome.down <- reactome.down[order(reactome.down$log10P),]
reactome.down$log10P <- reactome.down$log10P * -1

main.text <- c("Down-regulated pathways after treatment", "144 genes")
file.de <- file.path(wd.de.reactome, "genes_TR-vs-T_p1e-6_n144_down.pdf")

pdf(file.de, height=5.8, width=7.5)
par(mar=c(4,2,4,18))   # increase y-axis margin.
posbar <- barplot(reactome.down$log10P, main=main.text[1], las=1, horiz=T, xlim=c(-16, 0), xaxt="n", names.arg="", col="red", xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=reactome.down$Pathway.name, xpd=NA)

axis(side=1, at=seq(-14, 0, by=2), labels=c(14, 12, 10, 8, 6, 4, 2, 0))
mtext(main.text[2], line=0.3)
dev.off()
