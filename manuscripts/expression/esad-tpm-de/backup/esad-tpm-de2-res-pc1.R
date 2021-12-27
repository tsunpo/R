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

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
tpm.gene.log2 <- log2(tpm.gene + 0.01)   ## Use pseudocount=0.01

samples <- samples[rownames(subset(samples, GROUP_ID != 2)),]
#tpm.gene.log2 <- tpm.gene.log2[, rownames(samples)]
tpm.gene.log2.res <- tpm.gene.log2.res[, rownames(samples)]

# -----------------------------------------------------------------------------
# Wilcoxon rank sum test (non-parametric; n=45, 22 TR vs 23 UN)
# Last Modified: 22/08/20
# -----------------------------------------------------------------------------
## Test: Wilcoxon/Mann–Whitney/U/wilcox.test
##       Student's/t.test
## FDR : Q/BH
## DE  : TR (1) vs UN (0) as factor
argv      <- data.frame(predictor="GROUP_ID", predictor.wt=0, test="Wilcoxon", test.fdr="Q", stringsAsFactors=F)
file.name <- paste0("de_", base, "_tpm-gene-median0-res-pc1_B-vs-X_wilcox_q_n45")
file.main <- paste0("TR (n=22) vs UN (n=23) in ", BASE)

de <- differentialAnalysis(tpm.gene.log2.res, samples, argv$predictor, argv$predictor.wt, argv$test, argv$test.fdr)

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)

save(de.tpm.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Wilcoxon rank sum test (non-parametric; n=45, 22 X vs 23 B) RES
# Last Modified: 30/03/21; 22/08/20
# -----------------------------------------------------------------------------
## Test: Wilcoxon/Mann–Whitney/U/wilcox.test
##       Student's/t.test
## FDR : Q/BH
## DE  : TR (1) vs UN (0) as factor
argv      <- data.frame(predictor="GROUP_ID", predictor.wt=0, test="Wilcoxon", test.fdr="Q", stringsAsFactors=F)
file.name <- paste0("de_", base, "_tpm-gene-median0-res_X-vs-B_wilcox_q_n45")
file.main <- paste0("X (n=22) vs B (n=23) in ", BASE)

de <- differentialAnalysis(tpm.gene.log2.res, samples, argv$predictor, argv$predictor.wt, argv$test, argv$test.fdr)

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

trait <- as.numeric(samples[, "GROUP_ID2"])
trait[which(trait == 0)] <- "N"
trait[which(trait == 1)] <- "B"
trait[which(trait == 2)] <- "X"

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
   #ymax <- max(de$log10P)
   ymax <- 9.5
 
   pdf(file.de, height=6, width=6)
   plot(de$LOG2_FC, de$log10P, pch=16, xlim=c(-xmax, xmax), ylim=c(0, ymax), xlab="EAC X versus B [log2 fold change]", ylab="Significance [-log10(p-value)]", col="lightgray", main=file.main[1])

   #text(xmax*-1 + 2*xmax/28, -log10(pvalue) + ymax/42, paste0("FDR=", fdr, "%"), cex=0.85)
   #text(xmax*-1 + 2*xmax/35, -log10(pvalue) + ymax/42, "FDR=0.05", cex=0.85)
   #abline(h=c(-log10(fdrToP(0.1, de))), lty=5, col="darkgray")
   #text(xmax*-1 + 2*xmax/50, -log10(fdrToP(0.1, de)) + ymax/42, "FDR=0.1", col="darkgray", cex=0.85)

   de.up   <- subset(de.sig, LOG2_FC > 0)
   points(de.up$LOG2_FC, de.up$log10P, pch=16, col=purple)
   de.down <- subset(de.sig, LOG2_FC < 0)
   points(de.down$LOG2_FC, de.down$log10P, pch=16, col=red)
 
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
   legend("topleft", legend=c("Up-regulated (B < X)", "Down-regulated (B > X)"), col=c(purple, red), pch=19)
   dev.off()
}

##
plot.main <- "387 differentially expressed genes in EAC"
plot.de <- file.path(wd.de.plots, "volcanoplot_esad_median0_X-vs-B_p1e-6_CNA")

genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c(plot.main, "Before (B) vs. After treatment (X)")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene, 1.00E-06, genes, file.de, file.main)

##
plot.main <- "387 differentially expressed genes in EAC"
plot.de <- file.path(wd.de.plots, "volcanoplot_esad_median0_B-vs-X_p1e-6")

genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c(plot.main, "Before (B) vs. After treatment (X)")
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
wd.de.reactome <- file.path(wd.de.pathway, "reactome_X-vs-B_p1e-6_up")
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")

reactome.up <- subset(reactome, Entities.pValue <= 1e-6)[, 2:7]
reactome.up$log10P <- -log10(reactome.up$Entities.pValue)
reactome.up <- reactome.up[order(reactome.up$log10P),]

main.text <- c("Up-regulated pathways after treatment", "228 genes")
file.de <- file.path(wd.de.reactome, "genes_X-vs-B_p1e-6_n228_up.pdf")

pdf(file.de, height=3, width=7.5)
par(mar=c(4,18,4,3.1))   # increase y-axis margin.
barplot(reactome.up$log10P, main=main.text[1], las=1, horiz=T, xlim=c(0, 15), xaxt="n", names.arg=reactome.up$Pathway.name, col=purple, xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
abline(v=6, lty=5)

axis(side=1, at=seq(0, 14, by=2))
mtext(main.text[2], line=0.3)
dev.off()

##
wd.de.reactome <- file.path(wd.de.pathway, "reactome_X-vs-B_p1e-6_down")
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")

reactome.down <- subset(reactome, Entities.pValue <= 1e-6)[, 2:7]
reactome.down$log10P <- -log10(reactome.down$Entities.pValue)
reactome.down <- reactome.down[order(reactome.down$log10P),]
reactome.down$log10P <- reactome.down$log10P * -1

main.text <- c("Down-regulated pathways after treatment", "159 genes")
file.de <- file.path(wd.de.reactome, "genes_X-vs-B_p1e-6_n159_down.pdf")

pdf(file.de, height=5.8, width=7.5)
par(mar=c(4,2,4,18))   # increase y-axis margin.
posbar <- barplot(reactome.down$log10P, main=main.text[1], las=1, horiz=T, xlim=c(-16, 0), xaxt="n", names.arg="", col=red, xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=reactome.down$Pathway.name, xpd=NA)
abline(v=-6, lty=5)

axis(side=1, at=seq(-14, 0, by=2), labels=c(14, 12, 10, 8, 6, 4, 2, 0))
mtext(main.text[2], line=0.3)
dev.off()






# -----------------------------------------------------------------------------
# Heatmap
# Last Modified: 15/04/21; 11/01/20
# -----------------------------------------------------------------------------
#load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.RData")))
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
#load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.r5p47.RData")))
tpm.gene <- tpm.gene[, rownames(samples)]   ## VERY VERY VERY IMPORTANT!!!
tpm.gene.log2   <- log2(tpm.gene + 0.01)
tpm.gene.log2.m <- getLog2andMedian(tpm.gene, 0.01)

#genes <- c("ENSG00000137699", "ENSG00000158825", "ENSG00000196754", "ENSG00000121552", "ENSG00000175906", "ENSG00000108839", "ENSG00000204421", "ENSG00000167768", "ENSG00000205076", "ENSG00000205420", "ENSG00000172005", "ENSG00000170477", "ENSG00000126233", "ENSG00000057149", "ENSG00000229035", "ENSG00000198807", "ENSG00000088002", "ENSG00000170322", "ENSG00000108828")
genes <- c("ENSG00000092853", "ENSG00000117724", "ENSG00000143870", "ENSG00000088325", "ENSG00000110917", "ENSG00000086232", "ENSG00000138182", "ENSG00000148773", "ENSG00000119969", "ENSG00000169045", "ENSG00000112308", "ENSG00000122565", "ENSG00000122643", "ENSG00000066279", "ENSG00000139734", "ENSG00000118193", "ENSG00000198901", "ENSG00000163781", "ENSG00000131747", "ENSG00000198826", "ENSG00000163840", "ENSG00000143321", "ENSG00000040275", "ENSG00000096433", "ENSG00000106012", "ENSG00000138160", "ENSG00000087586", "ENSG00000046604", "ENSG00000144959", "ENSG00000111816", "ENSG00000072571", "ENSG00000092199", "ENSG00000111602")
genes <- c("ENSG00000092853", "ENSG00000117724", "ENSG00000143870", "ENSG00000088325", "ENSG00000110917", "ENSG00000086232", "ENSG00000138182", "ENSG00000148773", "ENSG00000119969", "ENSG00000169045", "ENSG00000112308", "ENSG00000122565", "ENSG00000122643", "ENSG00000066279", "ENSG00000139734", "ENSG00000118193", "ENSG00000198901", "ENSG00000163781", "ENSG00000131747", "ENSG00000198826", "ENSG00000163840", "ENSG00000143321", "ENSG00000040275", "ENSG00000096433", "ENSG00000106012", "ENSG00000138160", "ENSG00000087586", "ENSG00000046604", "ENSG00000144959", "ENSG00000111816", "ENSG00000072571", "ENSG00000092199", "ENSG00000111602", "ENSG00000090889", "ENSG00000155640", "ENSG00000184445", "ENSG00000104154", "ENSG00000076554", "ENSG00000148219", "ENSG00000169679", "ENSG00000187951", "ENSG00000137804", "ENSG00000197301", "ENSG00000101361", "ENSG00000177565", "ENSG00000107897", "ENSG00000125885", "ENSG00000171320", "ENSG00000149136", "ENSG00000170312", "ENSG00000142687", "ENSG00000058085", "ENSG00000175538", "ENSG00000204909", "ENSG00000108846", "ENSG00000127328", "ENSG00000143375", "ENSG00000184661")

genes <- c("ENSG00000143870", "ENSG00000092853", "ENSG00000110917", "ENSG00000112308", "ENSG00000169045", "ENSG00000088325", "ENSG00000122565", "ENSG00000117724", "ENSG00000163781", "ENSG00000086232", "ENSG00000096433", "ENSG00000118193", "ENSG00000163840", "ENSG00000119969", "ENSG00000144959", "ENSG00000040275", "ENSG00000155640", "ENSG00000121621", "ENSG00000111816", "ENSG00000092199", "ENSG00000106012", "ENSG00000072571", "ENSG00000148773", "ENSG00000148219", "ENSG00000204909", "ENSG00000122643", "ENSG00000138160", "ENSG00000090889", "ENSG00000184445", "ENSG00000198901", "ENSG00000197301", "ENSG00000184661", "ENSG00000125885", "ENSG00000107897", "ENSG00000198826", "ENSG00000139734", "ENSG00000138182", "ENSG00000177565", "ENSG00000131747", "ENSG00000111602")
genes <- c("ENSG00000143870", "ENSG00000092853", "ENSG00000110917", "ENSG00000112308", "ENSG00000169045", "ENSG00000088325", "ENSG00000122565", "ENSG00000117724", "ENSG00000163781", "ENSG00000086232", "ENSG00000096433", "ENSG00000118193", "ENSG00000163840", "ENSG00000119969", "ENSG00000144959", "ENSG00000040275", "ENSG00000155640", "ENSG00000121621", "ENSG00000111816", "ENSG00000092199", "ENSG00000106012", "ENSG00000072571", "ENSG00000148773", "ENSG00000148219", "ENSG00000204909", "ENSG00000122643", "ENSG00000138160", "ENSG00000090889", "ENSG00000184445", "ENSG00000198901", "ENSG00000197301", "ENSG00000184661", "ENSG00000125885", "ENSG00000107897", "ENSG00000198826", "ENSG00000139734", "ENSG00000138182", "ENSG00000177565", "ENSG00000131747", "ENSG00000111602", "ENSG00000134690", "ENSG00000106245", "ENSG00000034510", "ENSG00000175538", "ENSG00000105438", "ENSG00000149136", "ENSG00000138496", "ENSG00000143375", "ENSG00000094804", "ENSG00000142687", "ENSG00000143924", "ENSG00000112312", "ENSG00000136205", "ENSG00000137804", "ENSG00000101868", "ENSG00000102699", "ENSG00000109685", "ENSG00000120802", "ENSG00000171320", "ENSG00000046604", "ENSG00000187951", "ENSG00000116406", "ENSG00000146918", "ENSG00000087586", "ENSG00000112029", "ENSG00000138434", "ENSG00000155100", "ENSG00000196187", "ENSG00000170312", "ENSG00000154370", "ENSG00000165118", "ENSG00000108511", "ENSG00000115267", "ENSG00000165105", "ENSG00000111247", "ENSG00000147689", "ENSG00000037965", "ENSG00000127328")
#tpm.gene.log2.m[genes.clspn,]$MEDAIN

## B subgroup
samples.b <- samples[rownames(subset(samples, GROUP_ID == 0)),]
tpm.gene.b <- tpm.gene[, rownames(samples.b)]   ## VERY VERY VERY IMPORTANT!!!
tpm.gene.b.log2   <- log2(tpm.gene.b + 0.01)
tpm.gene.b.log2.m <- getLog2andMedian(tpm.gene.b, 0.01)

## X subgroup
samples.x <- samples[rownames(subset(samples, GROUP_ID == 1)),]
tpm.gene.x <- tpm.gene[, rownames(samples.x)]   ## VERY VERY VERY IMPORTANT!!!
tpm.gene.x.log2   <- log2(tpm.gene.x + 0.01)
tpm.gene.x.log2.m <- getLog2andMedian(tpm.gene.x, 0.01)

## N subgroup
samples.n <- samples[rownames(subset(samples, GROUP_ID == 2)),]
tpm.gene.n <- tpm.gene[, rownames(samples.n)]   ## VERY VERY VERY IMPORTANT!!!
tpm.gene.n.log2   <- log2(tpm.gene.n + 0.01)
tpm.gene.n.log2.m <- getLog2andMedian(tpm.gene.n, 0.01)

# -----------------------------------------------------------------------------
# Heatmap
# Links: https://davetang.org/muse/2010/12/06/making-a-heatmap-with-r/
# Last Modified: 13/01/20
# -----------------------------------------------------------------------------
#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq")
#install.packages("gplots")

library("DESeq")
library(gplots)

#genes.rownames <- c("TRIM29", "CDA", "S100A2", "CSTA", "ARL4D", "ALOX12", "LY6G6C", "KRT1", "LGALS7", "KRT6A", "MAL", "KRT4", "SLURP1", "SERPINB3", "SPRR2C", "PAX9", "SULT2B1", "NFRKB", "VAT1")
genes.rownames <- c("CLSPN", "CENPF", "PDIA6", "TPX2", "MLEC", "EIF2AK1", "KIF20B", "MKI67", "HELLS", "HNRNPH1", "C6orf62", "CBX3", "NT5C3A", "ASPM", "DIAPH3", "KIF14", "PRC1", "TOPBP1", "TOP2A", "ARHGAP11A", "DTX3L", "HDGF", "SPDL1", "ITPR3", "IQCE", "KIF11", "AURKA", "DSG2", "NCEH1", "FRK", "HMMR", "HNRNPC", "TIMELESS")
genes.rownames <- c("CLSPN", "CENPF", "PDIA6", "TPX2", "MLEC", "EIF2AK1", "KIF20B", "MKI67", "HELLS", "HNRNPH1", "C6orf62", "CBX3", "NT5C3A", "ASPM", "DIAPH3", "KIF14", "PRC1", "TOPBP1", "TOP2A", "ARHGAP11A", "DTX3L", "HDGF", "SPDL1", "ITPR3", "IQCE", "KIF11", "AURKA", "DSG2", "NCEH1", "FRK", "HMMR", "HNRNPC", "TIMELESS", "KIF4A", "C10orf12", "KNTC1", "SLC30A4", "TPD52", "ASTN2", "BUB1", "ARHGAP11B", "NUSAP1", "RP11-366L20.2", "NOP56", "TBL1XR1", "ACBD5", "MCM8", "ESCO2", "SSRP1", "CDK1", "KIAA0319L", "LAMC2", "KCNE3", "SPINK9", "ABCC3", "RAB3IP", "CGN", "CDCA2")

genes.rownames <- c("PDIA6", "CLSPN", "MLEC", "C6orf62", "HNRNPH1", "TPX2", "CBX3", "CENPF", "TOPBP1", "EIF2AK1", "ITPR3", "KIF14", "DTX3L", "HELLS", "NCEH1", "SPDL1", "C10orf12", "KIF18A", "FRK", "HNRNPC", "IQCE", "HMMR", "MKI67", "ASTN2", "SPINK9", "NT5C3A", "KIF11", "KIF4A", "KNTC1", "PRC1", "RP11-366L20.2", "CDCA2", "MCM8", "ACBD5", "ARHGAP11A", "DIAPH3", "KIF20B", "TBL1XR1", "TOP2A", "TIMELESS")
genes.rownames <- c("PDIA6", "CLSPN", "MLEC", "C6orf62", "HNRNPH1", "TPX2", "CBX3", "CENPF", "TOPBP1", "EIF2AK1", "ITPR3", "KIF14", "DTX3L", "HELLS", "NCEH1", "SPDL1", "C10orf12", "KIF18A", "FRK", "HNRNPC", "IQCE", "HMMR", "MKI67", "ASTN2", "SPINK9", "NT5C3A", "KIF11", "KIF4A", "KNTC1", "PRC1", "RP11-366L20.2", "CDCA2", "MCM8", "ACBD5", "ARHGAP11A", "DIAPH3", "KIF20B", "TBL1XR1", "TOP2A", "TIMELESS", "CDCA8", "BUD31", "TMSB10", "KCNE3", "KDELR1", "SSRP1", "PARP9", "CGN", "CDC6", "KIAA0319L", "EML4", "GMNN", "TNS3", "NUSAP1", "POLA1", "PARP4", "WHSC1", "TMPO", "ESCO2", "DSG2", "ARHGAP11B", "EDEM3", "NCAPG2", "AURKA", "FBXO5", "SSFA2", "OTUD6B", "TMEM63A", "CDK1", "TRIM11", "C9orf64", "HOXB6", "IFIH1", "RASEF", "RAD51AP1", "FAM83A", "HOXC8", "RAB3IP")

b <- tpm.gene.b.log2[genes,]
colnames(b) <- samples.b$NR_intern
rownames(b) <- genes.rownames
b <- data.matrix(b)

x <- tpm.gene.x.log2[genes,]
colnames(x) <- samples.x$NR_intern
rownames(x) <- genes.rownames
x <- data.matrix(x)

n <- tpm.gene.n.log2[genes,]
colnames(n) <- samples.n$NR_intern
rownames(n) <- genes.rownames
n <- data.matrix(n)

xb <- cbind(x, b)
nb <- cbind(n, b)
xn <- cbind(x, n)

heatmap(b)

## Returning the values used for the heatmap
test <- heatmap.2(b, scale="row")
b[rev(test$rowInd), test$colInd]

##
hr <- hclust(as.dist(1-cor(t(b), method="pearson")),  method="complete")
hc <- hclust(as.dist(1-cor(b,    method="spearman")), method="complete")

###
##
colfunc <- colorRampPalette(c("green", "black", "red"))
heatmap.2(b, col=colfunc(7), scale="row", trace="none")
