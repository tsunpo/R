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

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
tpm.gene.log2 <- log2(tpm.gene + 0.01)   ## Use pseudocount=0.01

b <- rownames(subset(samples, GROUP_ID == 0))
x <- rownames(subset(samples, GROUP_ID == 1))
n <- rownames(subset(samples, GROUP_ID == 2))
samples$GROUP_ID3 <- 0
samples[n,]$GROUP_ID3 <- 0
samples[x,]$GROUP_ID3 <- 1
samples[b,]$GROUP_ID3 <- 2
samples$GROUP_ID3 <- as.numeric(samples$GROUP_ID3)

tpm.gene.log2 <- tpm.gene.log2[, rownames(samples)]

# -----------------------------------------------------------------------------
# Residuals of expression
# Last Modified: 26/03/21; 28/08/20
# -----------------------------------------------------------------------------
#tpm.gene.log2 <- tpm.gene.log2[, rownames(samples)]
scores <- pcaScores(pca.de)
pc2 <- scores[rownames(samples), "PC2"]
#as.numeric(resid(lm(as.numeric(tpm.gene.log2[1,]) ~ pc1)))

tpm.gene.res <- t(mapply(x = 1:nrow(tpm.gene), function(x) as.numeric(resid(lm(as.numeric(tpm.gene[x,]) ~ pc2)))))
colnames(tpm.gene.res) <- rownames(samples)
rownames(tpm.gene.res) <- rownames(tpm.gene)
save(tpm.gene.res, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.median0.res.pc2.RData")))

tpm.gene.log2.res <- t(mapply(x = 1:nrow(tpm.gene.log2), function(x) as.numeric(resid(lm(as.numeric(tpm.gene.log2[x,]) ~ pc2)))))
colnames(tpm.gene.log2.res) <- rownames(samples)
rownames(tpm.gene.log2.res) <- rownames(tpm.gene.log2)
save(tpm.gene.log2.res, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.median0.log2.res.pc2.RData")))

#tpm.gene.log2 <- tpm.gene.log2.res

# -----------------------------------------------------------------------------
# Spearman's rank correlation
# Last Modified: 22/08/20
# -----------------------------------------------------------------------------
colnames <- c("RHO", "P", "FDR", "N", "X", "B", "FC_X_N", "FC_B_N")
de <- toTable(0, length(colnames), nrow(tpm.gene.log2), colnames)
rownames(de) <- rownames(tpm.gene.log2)

## SRC
de$RHO <- mapply(x = 1:nrow(tpm.gene.log2.res), function(x) cor.test(as.numeric(tpm.gene.log2.res[x,]), samples$GROUP_ID3, method="spearman", exact=F)[[4]])
de$P   <- mapply(x = 1:nrow(tpm.gene.log2.res), function(x) cor.test(as.numeric(tpm.gene.log2.res[x,]), samples$GROUP_ID3, method="spearman", exact=F)[[3]])
## ANOVA
#samples$Cancer_Type <- as.factor(pheno.all$Cancer_Type)
#de$ANOVA_P <- mapply(x = 1:nrow(expr.pheno.log2), function(x) testANOVA(x, expr.pheno.log2, pheno.all))

## Log2 fold change
de$N <- median00(tpm.gene.log2, rownames(subset(samples, GROUP_ID3 == 0)))
de$X <- median00(tpm.gene.log2, rownames(subset(samples, GROUP_ID3 == 1)))
de$B <- median00(tpm.gene.log2, rownames(subset(samples, GROUP_ID3 == 2)))
de$FC_X_N <- de$X - de$N
de$FC_B_N <- de$B - de$N

## FDR
library(qvalue)
de$FDR   <- qvalue(de$P)$qvalue
#de$ANOVA_Q <- qvalue(de$ANOVA_P)$qvalue
de <- de[order(de$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)   ## BE EXTRA CAREFUL!!

save(de.tpm.gene, samples, file=file.path(wd.de.data, "src_esad_tpm-gene-median0-residual-pc2_N-X-B_src_q_n68.RData"))
writeTable(de.tpm.gene, file.path(wd.de.data, "src_esad_tpm-gene-median0-residual-pc2_N-X-B_src_q_n68.txt"), colnames=T, rownames=F, sep="\t")

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
   ymax <- 25   ## 28/08/20
 
   pdf(file.de, height=6, width=6)
   plot(de$LOG2_FC, de$log10P, pch=16, xlim=c(-xmax, xmax), ylim=c(0, ymax), xlab="EAC B versus N [log2 fold change]", ylab="Significance [-log10(p-value)]", col="lightgray", main=file.main[1])

   de.up   <- subset(de.sig, LOG2_FC > 0)
   points(de.up$LOG2_FC, de.up$log10P, pch=16, col=red)
   de.down <- subset(de.sig, LOG2_FC < 0)
   points(de.down$LOG2_FC, de.down$log10P, pch=16, col=green)
 
   abline(v=c(-log2(2), log2(2)), lty=5, col="darkgray")
   abline(h=c(-log10(pvalue)), lty=5, col="black")
 
   for (g in 1:nrow(genes)) {
      gene <- subset(de, external_gene_name == genes[g,]$GENE)
      gene <- cbind(gene, genes[g,])
  
      if (nrow(gene) > 0) {
         points(gene$LOG2_FC, gene$log10P, pch=1, col="black")
   
         if (!is.na(gene$ADJ_1))
            if (is.na(gene$ADJ_2))
               text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=gene$ADJ_1, cex=0.85)
            else
               text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=c(gene$ADJ_1, gene$ADJ_2), cex=0.85)
         else
            if (gene$LOG2_FC > 0)
               text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=c(0, -0.5), cex=0.85)
            else
               text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=c(1, -0.5), cex=0.85)
      } else
         print(genes[g])
   }
 
   mtext(file.main[2], cex=1.2, line=0.3)
   legend("topleft", legend=c("Up-regulated (N < X < B)", "Down-regulated (N > X > B)"), col=c(red, green), pch=19)
   dev.off()
}

##
plot.main <- "458 additive expression after treatment"
plot.de <- file.path(wd.de.plots, "volcanoplot_EAC_median0_residual_pc2_N-X-B_p1e-12_CNA")
de.tpm.gene$LOG2_FC <- de.tpm.gene$FC_B_N

genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c(plot.main, "Between N, X and B (controlled for fibroblast)")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene, 1.00E-12, genes, file.de, file.main)

##
plot.main <- "458 additive expression after treatment"
plot.de <- file.path(wd.de.plots, "volcanoplot_EAC_median0_residual_pc2_N-X-B_p1e-12_CNA")
de.tpm.gene$LOG2_FC <- de.tpm.gene$FC_B_N

genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c(plot.main, "Between N, X and B (controlled for fibroblast)")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene, 1.00E-12, genes, file.de, file.main)

# -----------------------------------------------------------------------------
# plotBox3 in esad-tpm-de3.R
# -----------------------------------------------------------------------------
n <- tpm.gene.log2["ENSG00000092853", rownames(subset(samples, GROUP_ID3 == 0))]
x <- tpm.gene.log2["ENSG00000092853", rownames(subset(samples, GROUP_ID3 == 1))]
b <- tpm.gene.log2["ENSG00000092853", rownames(subset(samples, GROUP_ID3 == 2))]
ylim <- c(min(tpm.gene.log2["ENSG00000092853",]), max(tpm.gene.log2["ENSG00000092853",]))

file.name <- paste0("boxplot_EAC_tpm.gene.median0_CLSPN_log2")
plotBox3(wd.de.plots, file.name, n, x, b, main="CLSPN (ENSG00000092853)", names=c("N", "X", "B"), cols=c(green, purple, red), ylim)

n <- tpm.gene.log2.res["ENSG00000092853", rownames(subset(samples, GROUP_ID3 == 0))]
x <- tpm.gene.log2.res["ENSG00000092853", rownames(subset(samples, GROUP_ID3 == 1))]
b <- tpm.gene.log2.res["ENSG00000092853", rownames(subset(samples, GROUP_ID3 == 2))]
ylim <- c(min(tpm.gene.log2.res["ENSG00000092853",]), max(tpm.gene.log2.res["ENSG00000092853",]))

file.name <- paste0("boxplot_EAC_tpm.gene.median0_CLSPN_log2_res")
plotBox3(wd.de.plots, file.name, n, x, b, main="CLSPN (ENSG00000092853)", names=c("N", "X", "B"), cols=c(green, purple, red), ylim, ylab.txt="Residual of expression (controlled for PC2)")

# -----------------------------------------------------------------------------
# Replace Ensembl Gene IDs to gene name in Reactome results (Up- and Down-regulation)
# -----------------------------------------------------------------------------
wd.de.pathway <- file.path(wd.de, "pathway")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_N-X-B_PC1_p1e-3_up")
#wd.de.reactome <- file.path(wd.de.pathway, "reactome_N-X-B_PC1_p1e-3_down")

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
wd.de.reactome <- file.path(wd.de.pathway, "reactome_N-X-B_PC2_p1e-12_up")
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")

reactome.up <- subset(reactome, Entities.pValue <= 1e-4)[, 2:7]
reactome.up$log10P <- -log10(reactome.up$Entities.pValue)
reactome.up <- reactome.up[order(reactome.up$log10P),]

main.text <- c("Up-regulated in cancer proliferation", "203 genes")
file.de <- file.path(wd.de.reactome, "genes_N-X-B_res_pc2_p1e-12_n203_up.pdf")

pdf(file.de, height=4, width=7.5)
par(mar=c(4,20.7,4,1.3))   # increase y-axis margin.
barplot(reactome.up$log10P, main=main.text[1], las=1, horiz=T, xlim=c(0, 12), xaxt="n", names.arg=reactome.up$Pathway.name, col=red, xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
abline(v=4, lty=5)

axis(side=1, at=seq(0, 12, by=2))
mtext(main.text[2], line=0.3)
dev.off()

##
wd.de.reactome <- file.path(wd.de.pathway, "reactome_N-X-B_PC2_p1e-12_down")
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")
#reactome[3, 2] <- "TCA cycle and respiratory electron transport"

reactome.down <- subset(reactome, Entities.pValue <= 1e-4)[, 2:7]
reactome.down$log10P <- -log10(reactome.down$Entities.pValue)
reactome.down <- reactome.down[order(reactome.down$log10P),]
reactome.down$log10P <- reactome.down$log10P * -1
#reactome.down[2, 1] <- "Formyl peptide receptors bind formyl peptides and ligands"

main.text <- c("Down-regulated in cancer proliferation", "255 genes")
file.de <- file.path(wd.de.reactome, "genes_N-X-B_res_pc2_p1e-12_n255_down.pdf")

pdf(file.de, height=2.5, width=7.5)
par(mar=c(4,1.3,4,20.7))   # increase y-axis margin.
posbar <- barplot(reactome.down$log10P, main=main.text[1], las=1, horiz=T, xlim=c(-12, 0), xaxt="n", names.arg="", col=green, xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=reactome.down$Pathway.name, xpd=NA)
abline(v=-4, lty=5)

axis(side=1, at=seq(-12, 0, by=2), labels=c(12, 10, 8, 6, 4, 2, 0))
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

genes <- c("ENSG00000092853", "ENSG00000111602")
#genes <- c("ENSG00000137699", "ENSG00000158825", "ENSG00000196754", "ENSG00000121552", "ENSG00000175906", "ENSG00000108839", "ENSG00000204421", "ENSG00000167768", "ENSG00000205076", "ENSG00000205420", "ENSG00000172005", "ENSG00000170477", "ENSG00000126233", "ENSG00000057149", "ENSG00000229035", "ENSG00000198807", "ENSG00000088002", "ENSG00000170322", "ENSG00000108828")
#genes <- c("ENSG00000092853", "ENSG00000117724", "ENSG00000143870", "ENSG00000088325", "ENSG00000110917", "ENSG00000086232", "ENSG00000138182", "ENSG00000148773", "ENSG00000119969", "ENSG00000169045", "ENSG00000112308", "ENSG00000122565", "ENSG00000122643", "ENSG00000066279", "ENSG00000139734", "ENSG00000118193", "ENSG00000198901", "ENSG00000163781", "ENSG00000131747", "ENSG00000198826", "ENSG00000163840", "ENSG00000143321", "ENSG00000040275", "ENSG00000096433", "ENSG00000106012", "ENSG00000138160", "ENSG00000087586", "ENSG00000046604", "ENSG00000144959", "ENSG00000111816", "ENSG00000072571", "ENSG00000092199", "ENSG00000111602")
genes <- c("ENSG00000092853", "ENSG00000117724", "ENSG00000143870", "ENSG00000088325", "ENSG00000110917", "ENSG00000086232", "ENSG00000138182", "ENSG00000148773", "ENSG00000119969", "ENSG00000169045", "ENSG00000112308", "ENSG00000122565", "ENSG00000122643", "ENSG00000066279", "ENSG00000139734", "ENSG00000118193", "ENSG00000198901", "ENSG00000163781", "ENSG00000131747", "ENSG00000198826", "ENSG00000163840", "ENSG00000143321", "ENSG00000040275", "ENSG00000096433", "ENSG00000106012", "ENSG00000138160", "ENSG00000087586", "ENSG00000046604", "ENSG00000144959", "ENSG00000111816", "ENSG00000072571", "ENSG00000092199", "ENSG00000111602", "ENSG00000090889", "ENSG00000155640", "ENSG00000184445", "ENSG00000104154", "ENSG00000076554", "ENSG00000148219", "ENSG00000169679", "ENSG00000187951", "ENSG00000137804", "ENSG00000197301", "ENSG00000101361", "ENSG00000177565", "ENSG00000107897", "ENSG00000125885", "ENSG00000171320", "ENSG00000149136", "ENSG00000170312", "ENSG00000142687", "ENSG00000058085", "ENSG00000175538", "ENSG00000204909", "ENSG00000108846", "ENSG00000127328", "ENSG00000143375", "ENSG00000184661")
#tpm.gene.log2.m[genes.clspn,]$MEDAIN

## SRC-PC2 (1E-16; n=21)
genes <- c("ENSG00000092853", "ENSG00000117724", "ENSG00000143870", "ENSG00000088325", "ENSG00000110917", "ENSG00000086232", "ENSG00000138182", "ENSG00000148773", "ENSG00000119969", "ENSG00000169045", "ENSG00000112308", "ENSG00000122565", "ENSG00000122643", "ENSG00000066279", "ENSG00000139734", "ENSG00000118193", "ENSG00000198901", "ENSG00000163781", "ENSG00000131747", "ENSG00000198826", "ENSG00000163840")
## SRC (1E-16; n=13)
genes <- c("ENSG00000143870", "ENSG00000092853", "ENSG00000110917", "ENSG00000112308", "ENSG00000169045", "ENSG00000088325", "ENSG00000122565", "ENSG00000117724", "ENSG00000163781", "ENSG00000086232", "ENSG00000096433", "ENSG00000118193", "ENSG00000163840")
## SRC (1E-18; n=7)
genes <- c("ENSG00000143870", "ENSG00000092853", "ENSG00000110917", "ENSG00000112308", "ENSG00000169045", "ENSG00000088325", "ENSG00000122565")

## Fibroblast SRC-PC1 (1E-12; n=21)
genes <- c("ENSG00000087245", "ENSG00000204262", "ENSG00000189334", "ENSG00000164692", "ENSG00000154096", "ENSG00000108821", "ENSG00000168542", "ENSG00000129195", "ENSG00000091490", "ENSG00000186806", "ENSG00000099875", "ENSG00000143320", "ENSG00000144810", "ENSG00000136155", "ENSG00000134317", "ENSG00000140254", "ENSG00000177106", "ENSG00000051620", "ENSG00000133110", "ENSG00000163359", "ENSG00000254709")

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

genes.rownames <- c("CLSPN", "TIMELESS")

#genes.rownames <- c("TRIM29", "CDA", "S100A2", "CSTA", "ARL4D", "ALOX12", "LY6G6C", "KRT1", "LGALS7", "KRT6A", "MAL", "KRT4", "SLURP1", "SERPINB3", "SPRR2C", "PAX9", "SULT2B1", "NFRKB", "VAT1")
genes.rownames <- c("CLSPN", "CENPF", "PDIA6", "TPX2", "MLEC", "EIF2AK1", "KIF20B", "MKI67", "HELLS", "HNRNPH1", "C6orf62", "CBX3", "NT5C3A", "ASPM", "DIAPH3", "KIF14", "PRC1", "TOPBP1", "TOP2A", "ARHGAP11A", "DTX3L", "HDGF", "SPDL1", "ITPR3", "IQCE", "KIF11", "AURKA", "DSG2", "NCEH1", "FRK", "HMMR", "HNRNPC", "TIMELESS")
genes.rownames <- c("CLSPN", "CENPF", "PDIA6", "TPX2", "MLEC", "EIF2AK1", "KIF20B", "MKI67", "HELLS", "HNRNPH1", "C6orf62", "CBX3", "NT5C3A", "ASPM", "DIAPH3", "KIF14", "PRC1", "TOPBP1", "TOP2A", "ARHGAP11A", "DTX3L", "HDGF", "SPDL1", "ITPR3", "IQCE", "KIF11", "AURKA", "DSG2", "NCEH1", "FRK", "HMMR", "HNRNPC", "TIMELESS", "KIF4A", "C10orf12", "KNTC1", "SLC30A4", "TPD52", "ASTN2", "BUB1", "ARHGAP11B", "NUSAP1", "RP11-366L20.2", "NOP56", "TBL1XR1", "ACBD5", "MCM8", "ESCO2", "SSRP1", "CDK1", "KIAA0319L", "LAMC2", "KCNE3", "SPINK9", "ABCC3", "RAB3IP", "CGN", "CDCA2")

## SRC-PC2 (1E-16; n=21)
genes.rownames <- c("CLSPN", "CENPF", "PDIA6", "TPX2", "MLEC", "EIF2AK1", "KIF20B", "MKI67", "HELLS", "HNRNPH1", "C6orf62", "CBX3", "NT5C3A", "ASPM", "DIAPH3", "KIF14", "PRC1", "TOPBP1", "TOP2A", "ARHGAP11A", "DTX3L")
## SRC (1E-16; n=13)
genes.rownames <- c("PDIA6", "CLSPN", "MLEC", "C6orf62", "HNRNPH1", "TPX2", "CBX3", "CENPF", "TOPBP1", "EIF2AK1", "ITPR3", "KIF14", "DTX3L")
## SRC (1E-18; n=7)
genes.rownames <- c("PDIA6", "CLSPN", "MLEC", "C6orf62", "HNRNPH1", "TPX2", "CBX3")

## Fibroblast SRC-PC1 (1E-12; n=21)
genes.rownames <- c("MMP2", "COL5A2", "S100A14", "COL1A2", "THY1", "COL1A1", "COL3A1", "FAM64A", "SEL1L3", "VSIG10L", "MKNK2", "CRABP2", "COL8A1", "SCEL", "GRHL1", "DUOXA1", "EPS8L2", "HEBP2", "POSTN", "COL6A3", "IGLL5")

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
