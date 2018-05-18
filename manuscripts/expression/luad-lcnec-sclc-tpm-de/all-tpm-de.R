# =============================================================================
# Manuscript   : Loss of RB1 escapes cell cycle arrest
# Chapter      : RB-loss additive effect in SCLC
# Name         : manuscripts/expression/all-tpm-de.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 30/11/17
# =============================================================================
#wd.src <- "/projects/cangen/tyang2/dev/R"            ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"                 ## tyang2@gauss
#wd.src <- "/re/home/tyang2/dev/R"                    ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"                  ## tpyang@localhost

wd.src.handbook <- file.path(wd.src, "handbook-of")   ## Required handbooks/libraries for the manuscript
handbooks <- c("Common.R", "DifferentialExpression.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.handbook, x))))

wd.src.guide <- file.path(wd.src, "guide-to-the")     ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.guide, "hg19.RData"))

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
#wd <- "/ngs/cangen/tyang2/"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2/"   ## tyang2@local

wd.all <- paste0(wd, "ALL/")
wd.all.de <- paste0(wd.all, "analysis/expression/kallisto/luad-lcnec-sclc-rnaseq-de/")
setwd(wd.all)

# -----------------------------------------------------------------------------
# Find genes with "additive" effects between SCLC, LUAD and LCNEC
# Last Modified: 17/08/17
# -----------------------------------------------------------------------------
testANOVA <- function(x, expr, pheno) {
   fit1 <- lm(as.numeric(expr[x,]) ~ pheno$Cancer_Type)
   fit2 <- lm(as.numeric(expr[x,]) ~ 1)
   a1 <- anova(fit1, fit2)
 
   return(a1$Pr[2])
}

getMedian <- function(x, expr, pheno, type) {
   return(median(as.numeric(expr[x, rownames(subset(pheno, Cancer_Type == type))])))
}

##
pheno.all <- readTable(paste0(wd.all.de, "data/pheno.all_n198.txt"), header=T, rownames=T, sep="")
pheno.all$RB1[which(is.na(pheno.all$RB1))] <- "NA"
pheno.all$Cancer_Type <- as.factor(pheno.all$Cancer_Type)

load("/Users/tpyang/Work/uni-koeln/tyang2/ALL/analysis/expression/kallisto/luad-lcnec-sclc-rnaseq-de/data/all_kallisto_0.43.1_tpm.gene_r5_p47.RData")
# > dim(tpm.gene)
# [1] 18898   198
tpm.gene.log2 <- log2(tpm.gene + 0.01)
expr.pheno.log2 <- tpm.gene.log2[,rownames(pheno.all)]   ## VERY VERY VERY IMPORTANT!!!

colnames <- c("Ensembl_Gene_ID", "SRC_RHO", "SRC_P", "SRC_FDR", "ANOVA_P", "ANOVA_FDR", "LUAD", "LCNEC", "SCLC", "LCNEC_LUAD", "SCLC_LUAD")
de <- toTable(NA, length(colnames), nrow(expr.pheno.log2), colnames)
rownames(de) <- rownames(expr.pheno.log2)
de$Ensembl_Gene_ID <- rownames(expr.pheno.log2)

## ANOVA
pheno.all$Cancer_Type <- as.factor(pheno.all$Cancer_Type)
de$ANOVA_P <- mapply(x = 1:nrow(expr.pheno.log2), function(x) testANOVA(x, expr.pheno.log2, pheno.all))

## SRC
de$SRC_RHO   <- mapply(x = 1:nrow(expr.pheno.log2), function(x) cor.test(as.numeric(expr.pheno.log2[x,]), pheno.all$RB1_RATE, method="spearman", exact=F)[[4]])
de$SRC_P <- mapply(x = 1:nrow(expr.pheno.log2), function(x) cor.test(as.numeric(expr.pheno.log2[x,]), pheno.all$RB1_RATE, method="spearman", exact=F)[[3]])

## FC
de$LUAD  <- mapply(x = 1:nrow(expr.pheno.log2), function(x) getMedian(x, expr.pheno.log2, pheno.all, 0))
de$LCNEC <- mapply(x = 1:nrow(expr.pheno.log2), function(x) getMedian(x, expr.pheno.log2, pheno.all, 1))
de$SCLC  <- mapply(x = 1:nrow(expr.pheno.log2), function(x) getMedian(x, expr.pheno.log2, pheno.all, 2))
de$LCNEC_LUAD <- de$LCNEC - de$LUAD
de$SCLC_LUAD  <- de$SCLC - de$LUAD

library(qvalue)
de$ANOVA_FDR <- qvalue(de$ANOVA_P)$qvalue
de$SRC_FDR   <- qvalue(de$SRC_P)$qvalue
de <- de[order(de$SRC_P),]

##
annot.gene <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name",  "strand", "start_position", "end_position", "gene_biotype")]
overlaps <- intersect(rownames(de), rownames(annot.gene))   ## BE EXTRA CAREFUL!! NOT intersect(rownames(annot.gene), rownames(de)) 
de.all.tpm.gene <- cbind(annot.gene[overlaps,], de[,-1])    ## BE EXTRA CAREFUL!!

save(de.all.tpm.gene, pheno.all, file=paste0(wd.all.de, "de-all-tpm-gene/de_all_tpm-gene_rb1_src+anova_n198.RData"))
writeTable(de.all.tpm.gene, file=paste0(wd.all.de, "de-all-tpm-gene/de_all_tpm-gene_rb1_src+anova_n198.txt"), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Comparisions between SCLC, LUAD and LCNEC
# Last Modified: 18/08/17
# -----------------------------------------------------------------------------
#install.packages('beeswarm')
library(beeswarm)

plotBeeswarm <- function(gene, wd.de, expr.pheno.log2, pheno.all, position) {
   ensembl_gene_id <- subset(ensGene, external_gene_name == gene)$ensembl_gene_id[1]   ## To avoid ENSG00000269846 (one of RBL1)
   gene.tpms <- cbind(t(expr.pheno.log2)[rownames(pheno.all), ensembl_gene_id], pheno.all)
   colnames(gene.tpms)[1] <- "LOG2_TPM"
 
   pdf(paste0(wd.de, "plots/beeswarm/beeswarm_tpm.gene.log2_", gene, ".pdf"), height=6, width=4)
   ymin <- min(gene.tpms$LOG2_TPM)
   ymax <- max(gene.tpms$LOG2_TPM)
   boxplot(LOG2_TPM ~ Cancer_Type, data=gene.tpms, outline=F, names=c("LUAD", "LCNEC", "SCLC"), ylim=c(ymin, ymax), ylab="log2(TPM + 0.01)", main=paste0(gene, " (", ensembl_gene_id, ")"))
 
   beeswarm(LOG2_TPM ~ Cancer_Type, data=subset(gene.tpms, RB1 == "NA"), col="lightgray", pch=16, add=T)
   beeswarm(LOG2_TPM ~ Cancer_Type, data=subset(gene.tpms, RB1 == "WT"), col="blue", pch=16, add=T)
   beeswarm(LOG2_TPM ~ Cancer_Type, data=subset(gene.tpms, RB1 == "RB1"), col="red", pch=16, add=T)
 
   gene.tpms[intersect(rownames(gene.tpms),  c("S02229", "S01763", "S01093")),]$RB1 <- "LCNEC (WES)"
   gene.tpms["S00396",]$RB1 <- "LCNEC (WES)"
   gene.tpms["S00396",]$LOG2_TPM <- -8
   gene.tpms["S00006",]$RB1 <- "LCNEC (WES)"
   gene.tpms["S00006",]$LOG2_TPM <- -8
   beeswarm(LOG2_TPM ~ Cancer_Type, data=subset(gene.tpms, RB1=="LCNEC (WES)"), col="yellow", pch=16, add=T)
 
   legend(position, legend = c("RB1", "LCNEC (WES)", "WT"), pch=16, col=c("red", "yellow", "blue"))
   dev.off()
}

plotBox <- function(gene, wd.de, expr.pheno.log2, pheno.all) {
   ensembl_gene_id <- subset(ensGene, external_gene_name == gene)$ensembl_gene_id[1]   ## To avoid ENSG00000269846 (one of RBL1)
   gene.tpms <- cbind(t(expr.pheno.log2)[rownames(pheno.all), ensembl_gene_id], pheno.all)
   colnames(gene.tpms)[1] <- "LOG2_TPM"
 
   pdf(paste0(wd.de, "plots/boxplot/boxplot_tpm.gene.log2_", gene, ".pdf"), height=6, width=4)
   ymin <- min(gene.tpms$LOG2_TPM)
   ymax <- max(gene.tpms$LOG2_TPM)
   boxplot(LOG2_TPM ~ Cancer_Type, data=gene.tpms, outline=T, names=c("LUAD", "LCNEC", "SCLC"), ylim=c(ymin, ymax), ylab="log2(TPM + 0.01)", main=paste0(gene, " (", ensembl_gene_id, ")"))
 
   dev.off()
}

genes <- readTable(paste0(wd.all, "plots/beeswarm/gene.list"), header=F, rownames=F, sep="")
for (g in 1:length(genes))
   plotBeeswarm(genes[g], wd.all, tpm.gene.log2, pheno.all, "bottomleft")

genes <- c("BARD1", "CDKN2A", "MCM2", "PSIP1", "RBL1", "TP53", "TP53BP1", "TP53BP2")
genes <- c("POLH", "POLK", "POLQ", "FANCD2")
for (g in 1:length(genes))
   plotBeeswarm(genes[g], wd.all, tpm.gene.log2, pheno.all, "topleft")

genes <- readTable(paste0(wd.all, "plots/boxplot/gene.list"), header=F, rownames=F, sep="")
genes <- c("C7orf49")
genes <- c("XRCC5", "XRCC6", "PPM1D", "MAPK1", "MAPK3", "MAP2K1", "MAP2K2", "MAPK14")
genes <- c("RICTOR", "IRS2")
genes <- c("MMS22L", "FANCG", "FANCC", "FANCI", "FANCE", "FANCA", "FANCL", "FANCB", "FANCM", "FANCF", "FANCD2", "PALB2")
genes <- c("POLK")
genes <- c("POLA1", "POLD4", "POLE", "POLA2", "POLD3", "POLE2", "POLD1", "POLQ", "POLE4", "POLR2G", "POLR3F", "POLR2D", "POLE3",  "POLR3B", "POLRMT", "POLR3K", "POLR2I", "POLR2H", "POLR2E", "POLD2", "POLR1E", "POLH", "POLR3H", "POLDIP2", "POLI", "POLR3A", "POLR2L", "POLDIP3", "POLR3E", "POLK", "POLR2J", "POLR1D", "POLR2M", "POLR2C", "POLL", "POLM",  "POLR2F", "POLR1B", "POLR2A", "POLN", "POLR1C", "POLR2K", "POLG2", "POLR3G", "POLR3C", "POLR1A", "POLR3D", "POLR2J3", "POLR2B", "POLR2J2", "POLG", "POLR2J2", "POLR3GL", "POLRMTP1", "POLR2J4", "POLB")
genes <- c("PRIMPOL")
genes <- c("CNTNAP2", "EZH2")
for (g in 1:length(genes))
   plotBox(genes[g], wd.all.de, tpm.gene.log2, pheno.all)

# -----------------------------------------------------------------------------
# Principal component analysis (PCA)
# Last Modified: 17/09/17
# -----------------------------------------------------------------------------
summaryProportionofVariance <- function(summaries, pc) {
   return(round(summaries[2, pc]*100, 1))
}

plotPCAALL <- function(x, y, scores, summaries, trait, wd.pca, variable, file.main, legend.xy, cols, isID, isFlip) {
   #trait[is.na(trait)] <- "NA"
   trait.v <- sort(unique(trait))
 
   cols <- cols[1:length(trait.v)]
   trait.col <- mapply(x = 1:length(trait), function(x) return(cols[which(trait[x] == trait.v)]))   ## Assign colours to each subtypes
 
   #cols[which(trait.v == "NA")] <- "lightgrey"
   #trait.col[which(trait == "NA")] <- "lightgrey"
 
   xlab <- paste0("Principal component ", x, " (", summaryProportionofVariance(summaries, x), "%)")
   ylab <- paste0("Principal component ", y, " (", summaryProportionofVariance(summaries, y), "%)")
 
   if (isFlip)
      scores[,x] <- -scores[,x]
   scores[,y] <- -scores[,y]
   if (isID)
      pdf(paste0(wd.pca, "pca_", variable, "_", names(scores)[x], "-", names(scores)[y], "_ID.pdf"))
   else
      pdf(paste0(wd.pca, "pca_", variable, "_", names(scores)[x], "-", names(scores)[y], ".pdf"))
   plot(scores[,x], scores[,y], col=trait.col, pch=1, cex=1.5, main=file.main, xlab=xlab, ylab=ylab)
   if (isID)
      text(scores[,x], scores[,y], rownames(scores), cex=0.6, col="black", pos=3)
      
   legend(legend.xy, c("LUAD (n=48)", "LCNEC (n=69)", "SCLC (n=81)"), col=cols, pch=1, cex=1)   ##bty="n")
   dev.off()
}

##
tpm.gene.pcg.exp.log2 <- log2(tpm.gene.pcg.exp + 0.01)   ## Gene-level, protein-coding-gene TPMs (with 0 values)

pca <- getPCA(t(tpm.gene.pcg.exp.log2))
scores <- pcaScores(pca)
summaries <- pcaSummary(pca)

###
##
wd.pca <- paste0(wd, "analysis/expression/kallisto/all-rnaseq-de-kallisto-ensembl/data/")
file.main <- "Cancer Type"
trait <- pheno.all$Cancer_Type

#trait <- gsub(0, "LUAD", trait)
#trait <- gsub(1, "LCNEC", trait)
#trait <- gsub(2, "SCLC", trait)
plotPCAALL(1, 2, scores, summaries, trait, wd.pca, "all_n198_Cancer_Type", file.main, "bottomright", c("blue", "green", "red"), F, F)

###
##
de.tpm.gene.pcg.exp.sig <- rbind(subset(de.tpm.gene.pcg.exp, SRC_RHO >= 0.8), subset(de.tpm.gene.pcg.exp, SRC_RHO <= -0.8))

pca <- getPCA(t(tpm.gene.pcg.exp.log2[rownames(de.tpm.gene.pcg.exp.sig),]))   ## BUG FIX 13/02/17: Perform PCA using normalised data
scores <- pcaScores(pca)
summaries <- pcaSummary(pca)

file.main <- c("Cancer Type on 2 D.E. Genes", "rho >= ±0.8")
trait <- pheno.all$Cancer_Type

plotPCAALL(1, 2, scores, summaries, trait, wd.pca, "all_de_rho0.8_n198_Cancer_Type", file.main, "bottomright", c("blue", "green", "red"), F, F)

# -----------------------------------------------------------------------------
# Multidimensional Scaling (MDS)
# https://www.r-bloggers.com/multidimensional-scaling-mds-with-r/
# Last Modified: 17/09/17
# -----------------------------------------------------------------------------
de.tpm.gene.pcg.exp.sig <- rbind(subset(de.tpm.gene.pcg.exp, SRC_RHO >= 0.6), subset(de.tpm.gene.pcg.exp, SRC_RHO <= -0.6))
mydata <- t(tpm.gene.pcg.exp.log2[rownames(de.tpm.gene.pcg.exp.sig),])

d <- dist(mydata)   ## euclidean distances between the rows
fit <- cmdscale(d, eig=T, k=2)   ## k is the number of dim
fit   ## view results

# plot solution
x <- fit$points[,1]
y <- fit$points[,2]

trait <- pheno.all$Cancer_Type
trait.v <- sort(unique(trait))

cols <- c("blue", "green", "red")
cols <- cols[1:length(trait.v)]
trait.col <- mapply(x = 1:length(trait), function(x) return(cols[which(trait[x] == trait.v)]))   ## Assign colours to each subtypes

pdf(paste0(wd.pca, "mds_de_rho0.6_n198.pdf"))
plot(x, y, col=trait.col, pch=1, cex=1.5, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS")
text(x, y, labels=row.names(mydata), cex=.7) 
dev.off()

# -----------------------------------------------------------------------------
# Volcano plots
# Last Modified: 18/08/17
# -----------------------------------------------------------------------------
rhoToP <- function(rho, de) {
   de.sig.up <- subset(de, SRC_RHO >= rho)
   de.sig.down <- subset(de, SRC_RHO <= rho*-1)

   return(max(max(de.sig.up$SRC_P), max(de.sig.down$SRC_P)))
}

pToRHO <- function(p, de) {
   de.sig <- subset(de, SRC_P <= p)

   return(max(abs(de.sig$SRC_RHO)))
}

plotVolcanoALL <- function(de, rho, genes, file.de, file.main, isIg, up, fc) {
   if (!is.na(up))
      de <- subset(de, SRC_P > rhoToP(up, de))
   if (!is.na(fc))
      de <- subset(subset(de, Effect < fc), Effect > fc*-1)
   
   p <- rhoToP(rho, de)
   de.sig <- subset(de, SRC_P <= p)
   de.sig$log10P <- -log10(de.sig$SRC_P)
   #rho <- pToRHO(p, de)
 
   de$log10P <- -log10(de$SRC_P)
   xmax <- max(de$Effect)
   ymax <- max(de$log10P)
   #p <- max(significantAndVariable(de, effect, fdr)$P)
 
   pdf(file.de, height=7, width=7)
   plot(de$Effect, de$log10P, pch=16, xlim=c(-xmax, xmax), ylim=c(0, ymax), xlab=c("Median fold change (log2 SCLC/LUAD)", "", ""), ylab="Significance (-log10 P-value)", col="darkgray", main=file.main)
   abline(h=c(-log10(rhoToP(0.4, de))), lty=5, col="darkgray")
   text(xmax*-1 + 2*xmax/29, -log10(rhoToP(0.4, de)) + ymax/42, paste0("rho=±", 0.4), col="darkgray")
   
   if (isIg) {   ## ADD 02/11/17: Immunoglobulin (Ig) variable chain and T-cell receptor (TCR) genes
      de.ig <- de[grep("IG", de$gene_biotype),]
      points(de.ig$Effect, de.ig$log10P, pch=16, col="gold")
      de.tr <- de[grep("TR", de$gene_biotype),]
      points(de.tr$Effect, de.tr$log10P, pch=16, col="forestgreen")
   }
   de.up <- subset(de.sig, Effect > 0)
   points(de.up$Effect, de.up$log10P, pch=16, col="red")
   de.down <- subset(de.sig, Effect < 0)
   points(de.down$Effect, de.down$log10P, pch=16, col="dodgerblue")
   
   for (g in 1:length(genes)) {
      gene <- subset(de, external_gene_name == genes[g])
      if (nrow(gene) > 0) {
         points(gene$Effect, gene$log10P, pch=1, col="black")
       
         if (genes[g] == "TP53BP1" || genes[g] == "TP53")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(-0.07, 1.2), cex=0.75)
         else if (genes[g] == "ATM" || genes[g] == "DLL3")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(-0.05, 1.2), cex=0.75)
         else if (genes[g] == "KAT5")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(-0.12, 1.2), cex=0.75)
         else if (genes[g] == "TLR6" || genes[g] == "TLR3" || genes[g] == "CD28")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=1.12, cex=0.75)
         else if (genes[g] == "CDC7" || genes[g] == "XRCC3" || genes[g] == "HMGB1" || genes[g] == "CTLA4" || genes[g] == "CD86" || genes[g] == "NOTCH2" || genes[g] == "NHEJ1" || genes[g] == "SLFN11")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=1.1, cex=0.75)
         else if (genes[g] == "XRCC1")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=1.07, cex=0.75)  
         else if (genes[g] == "SMARCB1")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=1.05, cex=0.75)  
          else if (genes[g] == "PDCD1LG2" || genes[g] == "CDK9")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(1.07, 0), cex=0.75)
         else if (genes[g] == "MYD88")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(1.09, 0), cex=0.75)
         else if (genes[g] == "TLR7")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(1.14, 0.8), cex=0.75)
         else if (genes[g] == "HES1" || genes[g] == "REST")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(0, -0.5), cex=0.75)
         else if (genes[g] == "BRIP1")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(-0.05, -0.45), cex=0.75)
         else if (genes[g] == "BRCA2" || genes[g] == "WRN")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(-0.1, 0), cex=0.75)
         else if (genes[g] == "SMARCA4")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=-0.07, cex=0.75)
         else if (genes[g] == "TOP2A" || genes[g] == "ASF1B")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=-0.09, cex=0.75)
         else if (genes[g] == "CHEK2" || genes[g] == "XRCC2" || genes[g] == "BARD1" || genes[g] == "DNMT1")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=-0.1, cex=0.75)
         else if (genes[g] == "POLQ" || genes[g] == "BRCA1")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=-0.11, cex=0.75)
         else if (genes[g] == "XRCC4" || genes[g] == "NOTCH1" || genes[g] == "HMGB2" || genes[g] == "CDK6" || genes[g] == "CD274" || genes[g] == "ATR" || genes[g] == "RAD51B" || genes[g] == "SETD2" || genes[g] == "TLR8" || genes[g] == "TGFBR3")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(1.05, 1.3), cex=0.75)
         else if (genes[g] == "DNA2" || genes[g] == "RAD51D" || genes[g] == "EZH2" || genes[g] == "MAD2L2" || genes[g] == "ASCL1" || genes[g] == "FBXO5" || genes[g] == "RNF138" || genes[g] == "RNASEH2A" || genes[g] == "MSH6" || genes[g] == "CDK2" || genes[g] == "PCNA" || genes[g] == "PSIP1" || genes[g] == "MSH2" || genes[g] == "SIRT1" || genes[g] == "CHAF1B")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(1, -0.5), cex=0.75)
         else if (genes[g] == "RBBP8")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(1.05, -0.4), cex=0.75)
         else if (genes[g] == "TMPO" || genes[g] == "PPM1D")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(1, -0.45), cex=0.75)
         else if (genes[g] == "PARP2")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(1.06, -0.43), cex=0.75)
         else
            if (gene$Effect > 0)
               text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(0, -0.5), cex=0.75)
            else
               text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(1, -0.5), cex=0.75)
      }
   }
   
   if (isIg)
      legend("topleft", legend=c("Up-regulation", "Down-regulation", "Ig VJC genes (no D)", "TCR VC genes (no DJ)"), col=c("red", "dodgerblue", "gold", "forestgreen"), pch=19)
   else
      legend("topleft", legend=c("Up-regulation", "Down-regulation"), col=c("red", "dodgerblue"), pch=19)
   
   abline(h=c(-log10(p)), lty=5)
   text(xmax*-1 + 2*xmax/29, -log10(p) + ymax/42, paste0("rho=±", 0.6))
   dev.off()
}

###
## Volcano plot
colnames <- c("Ensembl_Gene_ID", "external_gene_name", "chromosome_name",  "strand", "start_position", "end_position", "gene_biotype", "SRC_RHO", "SRC_P", "SRC_FDR", "LUAD", "LCNEC", "SCLC", "LCNEC_LUAD", "SCLC_LUAD")
colnames(de.all.tpm.gene) <- colnames
de.all <- de.all.tpm.gene
#de.all <- de.all[setdiff(rownames(de.all), "ENSG00000269846"),]
de.all$Effect <- de.all$SCLC_LUAD
#de.all <- de.all[which(!is.na(de.all$SRC_P)),]
#de.all <- subset(test, SRC_P != 1)

plot.main <- "RB1mut-correlated differential expression in LUAD, LCNEC and SCLC"
plot.de <- paste0(wd.all.de, "volcanoplot_all_tpm-gene_rb1_rho6")

##
genes <- c("MYCL", "MYCN", "NEUROD1", "ASCL1", "SYP", "NCAM1", "ELAVL4", "CD44", "BMP4", "ATXN1")
file.main <- c(plot.main, "", "Non-/Neuroendocrine marker", "")
file.de   <- paste0(plot.de, "_Neuroendocrine.pdf")
plotVolcanoALL(de.all, 0.6, genes, file.de, file.main, F, NA, NA)

## Cycle   ##"ERCC6L2", "ERCC6L"
#genes <- c("PPM1D", "MAPK14", "CDK6", "CDK4", "MCM2", "RNASEH2A", "CDKN2A", "CDKN2C", "PARP1", "PARP2", "BARD1", "BRIP1", "BRCA1", "BRCA2", "CDC7", "CCNE2", "CCND1", "E2F2", "E2F7", "E2F1", "ATM", "ATR", "CHEK1", "CHEK2", "CDK2", "POLD4", "TOPBP1", "CLSPN", "TOP2A", "TOP1", "RB1", "RBL1", "RBL2", "TP73", "TP53", "TP53BP1", "PCNA")
#genes <- c("STMN1", "FBXO5", "TUBA1A", "TMPO", "MAD2L2", genes)
#genes <- c("CD44", "CHAF1B", "ASF1B", "PPM1D", "TMSB15A", "TMSB15B", "CDK6", "CDK4", "MCM2", "RNASEH2A", "CDKN2A", "CDKN2C", "PARP1", "PARP2", "BARD1", "BRIP1", "BRCA1", "BRCA2", "CDC7", "CCNE2", "CCND1", "E2F2", "E2F7", "E2F1", "ATM", "ATR", "CHEK1", "CHEK2", "CDK2", "POLD4", "TOPBP1", "CLSPN", "TOP2A", "TOP1", "RB1", "RBL1", "RBL2", "TP73", "TP53", "TP53BP1", "PCNA", "STMN1", "FBXO5", "TUBA1A", "TMPO", "MAD2L2")
genes <- c("TP53", "RB1", "E2F1", "E2F7", "E2F8", "E2F4", "E2F5", "CDK6", "CDK4", "MCM2", "RNASEH2A", "CDKN2A", "CDKN1A", "CDKN1B", "BARD1", "CDC7", "CDK2", "CCNE2", "CCND1", "POLD4", "TOPBP1", "CHEK2", "ATR", "ATM")
genes <- c("STMN1", "FBXO5", "TUBA1A", "TMPO", "MAD2L2", genes)
file.main <- c(plot.main, "", "Cell cycle control", "")
file.de   <- paste0(plot.de, "_Figure_2_Cycle.pdf")
plotVolcanoALL(de.all, 0.6, genes, file.de, file.main, F, NA, NA)

## PPM1D   ## MAPK14 (p38), HSPB1 (Hsp27), CBX1/CBX5 (HP1), MAP2K1 (MEK1), MAP2K2 (MEK2), MAPK3 (ERK1), MAPK1 (ERK2)
genes <- c("STMN1", "TUBA1A", "TMPO", "MAD2L2", "TP53", "RB1", "RBL1", "RBL2", "CDKN2A", "CCND1", "CBX1", "PPM1D", "MAPK1", "MAPK3", "MAP2K1", "MAP2K2", "MAPK14", "HSPB1", "CBX5", "EZH2", "SLFN11", "DNMT1", "DNMT3A", "DNMT3B", "BRCA1", "TP53BP1", "ATM")
file.main <- c(plot.main, "", "DNA methylation-mediated silencing", "")
file.de   <- paste0(plot.de, "_Cycle_PPM1D.pdf")
plotVolcanoALL(de.all, 0.6, genes, file.de, file.main, F, NA, NA)

## 53BP1   ##"UBE2N", "RNF8", "DNMT1", "DNMT3A", "DNMT3B", "KAT5", "SIRT1", "HELLS", "MAX", "ARID1A", "C7orf49", "CHAF1B", "ASF1B", "CHEK2",
genes <- c("RAD50", "NBN", "MRE11A", "MSH2", "MSH6", "TOPBP1", "SMARCAD1", "SMARCA4", "SMARCB1",  "NHEJ1", "LIG1", "LIG4", "E2F4", "POLQ", "DNA2", "SETD2", "RNF8", "RNF138", "RNF168", "PSIP1", "HMGB1", "HMGB2", "PARP1", "PARP2", "PARP3", "BARD1", "BRIP1", "BRCA1", "BRCA2", "RB1", "RBL1", "RBL2", "TP53", "TP53BP1", "RIF1", "MUM1", "MDC1", "XRCC1", "XRCC2", "XRCC3", "XRCC4", "RAD51AP1", "RBBP8", "FEN1", "EXO1", "BLM", "WRN")
file.main <- c(plot.main, "", "DNA repair pathway", "")
file.de   <- paste0(plot.de, "_53BP1_rho6_fc3.pdf")
plotVolcanoALL(de.all, 0.6, genes, file.de, file.main, F, 0.45, 3)

## 53BP1   ##"UBE2N", "RNF8", "DNMT1", "DNMT3A", "DNMT3B", "KAT5", "SIRT1", "HELLS", "MAX", "ARID1A", "C7orf49", "CHAF1B", "ASF1B", "CHEK2",
genes <- c("RAD50", "NBN", "MRE11A", "MSH2", "MSH6", "TOPBP1", "SMARCAD1", "SMARCA4", "SMARCB1",  "NHEJ1", "LIG1", "LIG4", "E2F4", "POLQ", "DNA2", "SETD2", "RNF8", "RNF138", "RNF168", "PSIP1", "HMGB1", "HMGB2", "PARP1", "PARP2", "PARP3", "BARD1", "BRIP1", "BRCA1", "BRCA2", "RB1", "RBL1", "RBL2", "TP53", "TP53BP1", "RIF1", "MUM1", "MDC1", "XRCC1", "XRCC2", "XRCC3", "XRCC4", "RAD51AP1", "RBBP8", "FEN1", "EXO1", "BLM", "WRN", "RAD52", "CDK6", "E2F5", "DNA2")
file.main <- c(plot.main, "", "DNA repair pathway", "")
file.de   <- paste0(plot.de, "_53BP1_rho6_fc3_2.pdf")
plotVolcanoALL(de.all, 0.6, genes, file.de, file.main, F, 0.45, 3)

file.de   <- paste0(plot.de, "_53BP1.pdf")
plotVolcanoALL(de.all, 0.6, genes, file.de, file.main, F, NA, NA)

## Immune response
genes <- c("CD274", "PDCD1LG2", "MYD88", "NFKB1", "NFKB2", "TLR2", "BCL3", "IL6", "IL6ST", "IL6R")
file.main <- c("RB1mut-correlated differential expression in LUAD, LCNEC and SCLC", "", "Immune response", "")
file.de   <- paste0(plot.de, "_Immune.pdf")
plotVolcanoALL(de.all, 0.6, genes, file.de, file.main, T, NA, NA)

## MHC I
genes <- c("B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F")
file.main <- c(plot.main, "", "MHC class I", "")
file.de <- paste0(plot.de, "_MHCI.pdf")
plotVolcanoALL(de.all, 0.6, genes, file.de, file.main, T, NA, NA)

## MHC II
genes <- c("CIITA", "HLA-DRA", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1", "HLA-DRB5")
file.main <- c(plot.main, "", "MHC class II", "")
file.de <- paste0(plot.de, "_MHCII.pdf")
plotVolcanoALL(de.all, 0.6, genes, file.de, file.main, T, NA, NA)

#file.main <- c("RB1mut-associated gene expression profiling in SCLC", "(with genes not-expressed in LUAD)", "Neuroendocrine marker", "")
#file.de <- paste0(wd.de, "volcanoplot_all_tpm-gene-pcg_rho6_NOTCH.pdf")
#plotVolcanoALL(de.tpm.gene.pcg[setdiff(rownames(de.tpm.gene.pcg), "ENSG00000269846"),], 0.6, genes, file.de, file.main)

## Immune response   ##"IL6R", "AHRR", "AHR", "IDH1"
#genes <- c("CD28", "CD86", "CD274", "PDCD1LG2", "CTLA4", "PARP3", "MYD88", "TNFSF13", "CYP1A2", "CYP1B1", "TRADD", "TNFRSF1A", "TNFRSF10B", "NFKB1", "NFKB2", "TLR2", "TLR5", "TLR3", "TLR6", "TLR7", "TLR8", "TLR4", "BCL2", "BCL2L11", "BCL3", "RUNX1", "RUNX2")
#file.main <- c("RB1mut-associated gene expression profiling in SCLC", "", "Immune response", "")
#file.de <- paste0(wd.de, "_Immune_TLRX.pdf")
#plotVolcanoALL(de.all, 0.6, genes, file.de, file.main)

###
## Immunoglobulin class switching
genes <- c("IL4R", "IFNGR2", "IFNGR1", "TGFBR2", "TGFB1", "TGFB3", "TGFB2", "TGFBR1", "TGFBR3")
file.main <- c(plot.main, "", "Immunoglobulin class switching", "")
file.de <- paste0(plot.de, "_Ig.pdf")
plotVolcanoALL(de.all, 0.6, genes, file.de, file.main, T)

## Innate immunity
genes <- c("IL6R", "GSK3A", "GSK3B")
file.main <- c(plot.main, "", "Innate immunity", "")
file.de <- paste0(plot.de, "_Innate.pdf")
plotVolcanoALL(de.all, 0.6, genes, file.de, file.main, T)

## XRCC   ##"UBE2N", "RNF8" 
#genes <- c("SMARCA4", "SMARCB1", "NHEJ1", "LIG4", "POLQ", "SETD2", "RNF138", "RNF168", "PSIP1", "HMGB1", "HMGB2", "PARP1", "PARP2", "PARP3", "BARD1", "BRIP1", "BRCA1", "BRCA2", "RB1", "RBL1", "RBL2", "TP53", "TP53BP1", "XRCC1", "EZH2", "XRCC2", "XRCC3", "XRCC4", "RAD51AP1", "RBBP8", "FEN1", "EXO1", "BLM", "WRN")
#file.main <- c(plot.main, "", "DNA repair pathway", "")
#file.de   <- paste0(plot.de, "_53BP1_XRCC.pdf")
#plotVolcanoALL(de.all, 0.6, genes, file.de, file.main, F)

## NOTCH
genes <- c("CDK9", "KRAS", "MYC", "MAX", "MYCL", "MYCN", "NEUROD1", "ASCL1", "DLL3", "PSIP1", "HMGB2", "BRD4", "EZH2", "UCHL1", "HES1", "REST", "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "EGFR", "ERBB2", "ERBB3", "DPYSL5", "CRMP1", "DPYSL3")
file.main <- c(plot.main, "", "Neuroendocrine marker", "")
file.de   <- paste0(plot.de, "_Neuroendocrine.pdf")
plotVolcanoALL(de.all, 0.6, genes, file.de, file.main, F)

# -----------------------------------------------------------------------------
# Volcano plots (Lite)
# Last Modified: 22/10/17
# -----------------------------------------------------------------------------
rhoToP <- function(rho, de) {
   de.sig.up <- subset(de, SRC_RHO >= rho)
   de.sig.down <- subset(de, SRC_RHO <= rho*-1)
 
   return(max(max(de.sig.up$SRC_P), max(de.sig.down$SRC_P)))
}

pToRHO <- function(p, de) {
   de.sig <- subset(de, SRC_P <= p)
 
   return(min(abs(de.sig$SRC_RHO)))
}

plotVolcanoALLLite <- function(de, rho, genes, file.de, file.main) {
   p <- rhoToP(rho, de)
   de.sig <- subset(de, SRC_P <= p)
   de.sig$log10P <- -log10(de.sig$SRC_P)
   #rho <- round(pToRHO(p, de), 2)
 
   de$log10P <- -log10(de$SRC_P)
   xmax <- max(de$Effect)
   ymax <- max(de$log10P)
   #p <- max(significantAndVariable(de, effect, fdr)$P)
 
   pdf(file.de, height=7, width=7)
   plot(de$Effect, de$log10P, pch=16, xlim=c(-xmax, xmax), ylim=c(0, ymax), xlab=c("Median fold change (log2 SCLC/LUAD)", "", ""), ylab="Significance (-log10 P-value)", col="darkgray", main=file.main)
   abline(h=c(-log10(rhoToP(0.4, de))), lty=5, col="darkgray")
   text(xmax*-1 + 2*xmax/29, -log10(rhoToP(0.4, de)) + ymax/42, paste0("rho=±", 0.4), col="darkgray")
 
   de.up <- subset(de.sig, Effect > 0)
   points(de.up$Effect, de.up$log10P, pch=16, col="red")
 
   de.down <- subset(de.sig, Effect < 0)
   points(de.down$Effect, de.down$log10P, pch=16, col="dodgerblue")
 
   for (g in 1:length(genes)) {
      gene <- subset(de, external_gene_name == genes[g])
      if (nrow(gene) > 0) {
         points(gene$Effect, gene$log10P, pch=1, col="black")
   
         if (genes[g] == "RBL2" || genes[g] == "LIG4" || genes[g] == "CDK6")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=1.15, cex=0.75)
         else if (genes[g] == "SMARCA4" || genes[g] == "CREBBP")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=-0.08, cex=0.75)
         else if (genes[g] == "XRCC2" || genes[g] == "XRCC3" || genes[g] == "TP53BP1")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=-0.09, cex=0.75)
         else if (genes[g] == "BARD1" || genes[g] == "BRCA1" || genes[g] == "SLFN11" || genes[g] == "CHEK1" || genes[g] == "RAD51" || genes[g] == "TOPBP1" || genes[g] == "RBL1" || genes[g] == "HMGB1" || genes[g] == "ASF1B" || genes[g] == "CCNE2" || genes[g] == "TUBA1A" || genes[g] == "PSIP1")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=-0.1, cex=0.75)
         else if (genes[g] == "WRN" || genes[g] == "TP53" || genes[g] == "E2F7")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=-0.12, cex=0.75)
         else if (genes[g] == "NHEJ1")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(1.04, 1.4), cex=0.75)
         else if (genes[g] == "ATM" || genes[g] == "RAD52")   ##CBX5
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(-0.05, 1.2), cex=0.75)
         else if (genes[g] == "XRCC1")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(-0.07, 1.2), cex=0.75)
         else if (genes[g] == "")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(-0.08, 1.2), cex=0.75)
         else if (genes[g] == "ATR")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(1.05, 1.3), cex=0.75)
         else if (genes[g] == "MAPK3")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=1.1, cex=0.75)
         else if (genes[g] == "DNA2" || genes[g] == "RNF168" || genes[g] == "RNF8")  
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(1, -0.5), cex=0.75)
         #else if (genes[g] == "PSIP1")
         #   text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(-0.05, 1.3), cex=0.75)
         else
            if (gene$Effect > 0)
               text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(0, -0.5), cex=0.75)
            else
               text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(1, -0.5), cex=0.75)
      } else
         print(genes[g])
   }
 
   abline(h=c(-log10(rhoToP(0.6, de))), lty=5)
   text(xmax*-1 + 2*xmax/29, -log10(rhoToP(0.6, de)) + ymax/42, paste0("rho=±", 0.6))
   legend("topleft", legend=c("Up-regulation", "Down-regulation"), col=c("red", "dodgerblue"), pch=19)
   dev.off()
}

###
## Volcano plot
#de.all <- de.all.tpm.gene.pcg.exp3
de.all <- de.all.tpm.gene
de.all <- de.all[setdiff(rownames(de.all), "ENSG00000269846"),]

colnames <- c("Ensembl_Gene_ID", "external_gene_name", "chromosome_name",  "strand", "start_position", "end_position", "gene_biotype", "SRC_RHO", "SRC_P", "SRC_FDR", "LUAD", "LCNEC", "SCLC", "LCNEC_LUAD", "SCLC_LUAD")
colnames(de.all) <- colnames
de.all$Effect <- de.all$SCLC_LUAD

##
plot.main <- "RB-loss additive effect in SCLC"
plot.de <- paste0(wd.all.de, "volcanoplot_all_tpm-gene_rb1_rho6")

genes <- c("POLA1", "POLD4", "POLE", "POLA2", "POLD3", "POLE2", "POLD1", "POLQ", "POLE4", "POLE3", "POLD2", "POLH", "POLI", "POLDIP3", "POLK", "POLL", "POLM", "POLN", "POLG2", "POLG", "POLB")
file.main <- c(plot.main, "", "DNA polymerases", "")
file.de   <- paste0(plot.de, "_Figure_2_Polymerase.pdf")
plotVolcanoALLLite(de.all, 0.6, genes, file.de, file.main)

genes <- c("MMS22L", "FANCG", "FANCC", "FANCI", "FANCE", "FANCA", "FANCB", "FANCF", "FANCD2", "BRCA1", "PALB2", "BRCA2", "RAD51AP1", "RAD51", "RAD51B", "RAD51C", "RAD51D", "XRCC2")   #"FANCL", "FANCM",
genes <- c("FANCG", "FANCC", "FANCI", "FANCE", "FANCA", "FANCB", "FANCF", "FANCD2", "BRCA1", "PALB2", "RAD51", "XRCC2")   #"FANCL", "FANCM",
file.main <- c(plot.main, "", "Fanconi anaemia repair", "")
file.de   <- paste0(plot.de, "_Figure_2_Fanconi.pdf")
plotVolcanoALLLite(de.all, 0.6, genes, file.de, file.main)

genes <- c("TP53", "RB1", "E2F1", "E2F7", "E2F8", "E2F4", "E2F5", "CDK6", "CDK4", "MCM2", "RNASEH2A", "CDKN2A", "CDKN1B", "CDC7", "CDK2", "CCNE2", "CCND1", "POLD4", "TOPBP1", "CHEK2", "ATR", "ATM")
genes <- c("STMN1", "FBXO5", "TUBA1A", "TMPO", "MAD2L2", genes)
file.main <- c(plot.main, "", "Cell cycle control", "")
file.de   <- paste0(plot.de, "_Figure_2_Cycle.pdf")
plotVolcanoALLLite(de.all, 0.6, genes, file.de, file.main)

genes <- c("RB1", "RBL1", "RBL2", "XRCC1", "XRCC2", "XRCC3", "XRCC4", "NHEJ1", "LIG1", "LIG3", "LIG4", "TP53BP1", "PARP1", "FEN1", "BLM", "EXO1", "DNA2", "WRN", "BRCA1", "BRCA2", "RAD51", "RAD52", "RAD51AP1", "PSIP1", "RBBP8", "RNF138", "RNF168", "BARD1", "TOPBP1")   #"ASF1A", "CHAF1A", "RB1", "RBL1", "RBL2", "E2F4", "TP53BP2", "MAX", "PARP1", "PARP2", "CHEK1", "CHEK2", "MSH2")
file.main <- c(plot.main, "", "DSB repair pathway", "")
file.de <- paste0(plot.de, "_Figure_2_53BP1.pdf")
plotVolcanoALLLite(de.all, 0.6, genes, file.de, file.main)

## Figure 1 (Histone chaperone and DSB repair pathway choice)
genes <- c("TP53BP1", "RB1",  "RBL1", "RBL2", "RNF8", "RNF168", "RNF138", "ASF1B", "CHAF1B", "CHAF1A", "ASF1A", "TOPBP1", "PSIP1", "RBBP8", "RBBP4", "HMGB2", "HMGB1", "SMARCAD1", "SMARCA4", "SMARCB1", "EP300", "CREBBP", "DNMT1")   ## "RAD51", "EXO1", "DNA2", "MDC1", "RIF1", "MAD2L2", "CDH1", "RNF8")
file.main <- c(plot.main, "", "Histone chaperone and chromatin remodeler", "")
file.de <- paste0(plot.de, "_Figure_2_53BP1_ASF1A.pdf")
plotVolcanoALLLite(de.all, 0.6, genes, file.de, file.main)

## PPM1D   ## MAPK14 (p38), HSPB1 (Hsp27), CBX1/CBX5 (HP1), MAP2K1 (MEK1), MAP2K2 (MEK2), MAPK3 (ERK1), MAPK1 (ERK2)
#genes <- c("TP53", "RB1", "RBL1", "RBL2", "CDKN2A", "MAP2K1", "MAPK3", "MAPK1", "PPM1D", "MAPK14", "TP53BP1", "ATR", "ATM", "NBN", "H2AFX", "PRKDC", "CHEK1", "CHEK2")   ## "MAP2K2"
genes <- c("EZH2", "PPM1D", "DNMT1", "DNMT3B", "CBX1", "CBX3", "CBX5", "APOBEC3A", "APOBEC3B", "MAPK14", "MAPK1", "MAPK3", "MAP2K1", "MAP2K2")
file.main <- c(plot.main, "", "DNA methylation-induced silencing", "")
file.de   <- paste0(plot.de, "_Figure_2_PPM1D.pdf")
plotVolcanoALLLite(de.all, 0.6, genes, file.de, file.main)

# -----------------------------------------------------------------------------
# CDKN2A (p16/p14)
# Last Modified: 30/11/17
# -----------------------------------------------------------------------------
transcripts <- rownames(subset(ensGene.transcript, external_gene_name == "CDKN2A"))
transcripts <- intersect(transcripts, rownames(tpm))
tpm.cdkn2a <- log2(tpm[transcripts,] + 0.01)

cdkn2a <- ensGene.transcript[transcripts,]

pdf(paste0(plot.de, ".pdf"), height=7, width=7)
#hist(median0(log2(tpm[cdkn2a,] + 0.01)), xlab="Median log2(TPM + 0.01)", ylab="Frequency")
plot(c(min(cdkn2a$transcript_start), max(cdkn2a$transcript_start) + 15000), c(-1, 5.5), col="white", xlab="CDKN2A (ENSG00000147889)", ylab="Median log2(TPM + 0.01)")
for (t in 1:length(transcripts)) {
   segments(cdkn2a$transcript_start[t], median(tpm.cdkn2a[t,]), cdkn2a$transcript_end[t], median(tpm.cdkn2a[t,]), col="blue")
   text(cdkn2a$transcript_end[t] + 5000, median(tpm.cdkn2a[t,]), paste0(cdkn2a$ensembl_transcript_id[t], " (FC=", round(de.lcnec.tpm.cdkn2a$Effect[t], 2), ")"), cex=0.55)
}
dev.off()

# -----------------------------------------------------------------------------
# Perform differential analysis using non-parametric test (RB1NEW vs WT; n=69-15NA)
# Last Modified: 28/03/17
# -----------------------------------------------------------------------------
pipeDE <- function(expr, pheno.expr, dao, file.de, file.main, annot) {
   ## Differential expression
   de <- differentialAnalysis(expr, pheno.expr, dao$predictor, dao$predictor.wt, dao$test, dao$test.fdr)
   de <- cbind(annot[rownames(de),], de)
 
   writeTable(de, paste0(file.de, ".txt"), colnames=T, rownames=T, sep="\t")
   #writeTable(onlyVariable(de, dao$effect), paste0(file.de, "_Effect>", dao$effect, ".txt"), colnames=T, rownames=T, sep="\t")   ## Only varible genes
 
   ## Volcano plot
   #plotVolcano(de, dao$fdr, dao$effect, file.de, file.main, legend.x=NA)
 
   return(de)
}

filenameDE <- function(wd.de, title, dao) {
   return(paste0(wd.de, title, "_", dao$predictor, "_", dao$test, "_", dao$test.fdr))
}

## Data access object (DAO) for test parameters
## Test: Wilcoxon/Wilcox/U/Students/ttest
## FDR:  Q/BH
## DE:   RB1NEW vs WT(0) as factor
dao <- data.frame(test="Wilcox", test.fdr="Q", fdr=0.05, effect=0, predictor="RB1NEW", predictor.wt=0, stringsAsFactors=F)
expr <- tpm.cdkn2a[,rownames(pheno.expr)]
#pheno.expr.final <- getFinalPhenotype(expr, pheno.expr, dao$predictor)
#expr.pheno <- getFinalExpression(expr, pheno.expr.final)

annot.gene <- ensGene.transcript[,c("ensembl_transcript_id", "external_gene_name", "gene_biotype")]
de.lcnec.tpm.cdkn2a <- pipeDE(expr, pheno.expr, dao, file.de, file.main, annot.gene)










## NHEJ and alternative pathways to DSB repair   ## NHEJ1 (XLF/PAXX)
genes <- c("XRCC1", "XRCC2", "XRCC3", "XRCC4", "NHEJ1", "FEN1", "LIG1", "LIG4", "PARP1")
file.main <- c(plot.main, "", "DSB repair pathway choice", "")
file.de <- paste0(plot.de, "_53BP1_Lite_0.pdf")
plotVolcanoALLLite(de.all, 0.6, genes, file.de, file.main)

genes <- c("XRCC1", "XRCC2", "XRCC3", "XRCC4", "NHEJ1", "LIG1", "LIG3", "LIG4", "PARP1", "FEN1", "BLM", "EXO1", "WRN", "BRCA1", "BRCA2", "RAD51", "RAD51AP1", "POLQ", "PSIP1", "RBBP8", "RNF138", "RNF168", "BARD1", "RAD50", "NBN", "MRE11A")   #"RB1", "RBL1", "RBL2", "E2F4", "TP53BP2", "MAX", "PARP1", "PARP2", "CHEK1", "CHEK2", "MSH2")
file.main <- c(plot.main, "", "DSB repair pathway choice", "")
file.de <- paste0(plot.de, "_53BP1_Lite_1.pdf")
plotVolcanoALLLite(de.all, 0.6, genes, file.de, file.main)

## PPM1D   ## MAPK14 (p38), HSPB1 (Hsp27), CBX1/CBX5 (HP1), MAP2K1 (MEK1), MAP2K2 (MEK2), MAPK3 (ERK1), MAPK1 (ERK2)
genes <- c("TP53", "RB1", "RBL1", "RBL2", "CDKN2A", "MAP2K1", "MAPK3", "MAPK1", "PPM1D", "MAPK14", "TP53BP1", "ATR", "ATM", "NBN", "H2AFX", "PRKDC", "CHEK1", "CHEK2")   ## "MAP2K2"
file.main <- c(plot.main, "", "PPM1D-p53 loop", "")
file.de   <- paste0(plot.de, "_53BP1_PPM1D_Lite.pdf")
plotVolcanoALLLite(de.all, 0.6, genes, file.de, file.main)

## PPM1D   ## MAPK14 (p38), HSPB1 (Hsp27), CBX1/CBX5 (HP1), MAP2K1 (MEK1), MAP2K2 (MEK2), MAPK3 (ERK1), MAPK1 (ERK2)
genes <- c("TP53", "RB1", "RBL1", "RBL2", "CDKN2A", "CBX1", "PPM1D", "EZH2", "SLFN11", "DNMT1", "DNMT3A", "DNMT3B", "BRCA1", "ATM")   ## "MAP2K2"
file.main <- c(plot.main, "", "DNA methylation-mediated silencing", "")
file.de   <- paste0(plot.de, "_53BP1_PPM1D_EZH2_Lite.pdf")
plotVolcanoALLLite(de.all, 0.6, genes, file.de, file.main)






## Proliferation / Apoptosis inhibition
genes <- c("TRADD", "TNFRSF1A", "PARP3", "BCL2", "BCL2L11", "BBC3", "PMAIP1")
file.main <- c("RB1mut-correlated differential expression in LUAD, LCNEC and SCLC", "", "Apoptosis inhibition (Proliferation)", "")
plot.de <- paste0(file.de, "_BCL2_Lite.pdf")
plotVolcanoALLLite(de.all, 0.6, genes, plot.de, file.main)











## Immunoglobulin class switching
genes <- c("IL4R", "IFNGR2", "IFNGR1", "TGFBR2", "TGFB1", "TGFB3", "TGFB2", "TGFBR1", "TGFBR3")
file.main <- c("RB1mut-associated expression profiling (LUAD, LCNEC and SCLC)", "", "Immunoglobulin class switching", "")
file.de <- paste0(wd.de, "volcanoplot_all_tpm-gene-pcg-exp_p6_Ig.pdf")
plotVolcanoALL(de.tpm.gene.pcg.exp[setdiff(rownames(de.tpm.gene.pcg.exp), "ENSG00000269846"),], 0.6, genes, file.de, file.main)

## Immunoglobulin class switching
genes <- c("IL4R", "IFNGR2", "IFNGR1", "TGFBR2", "TGFB1", "TGFB3", "TGFB2", "TGFBR1", "TGFBR3")
file.main <- c("RB1mut-associated gene expression profiling in SCLC", "", "Immunoglobulin class switching", "")
file.de <- paste0(wd.de, "volcanoplot_all_tpm-gene-pcg-exp_rho6_Immunoglobulin.pdf")
plotVolcanoALL(de.tpm.gene.pcg.exp[setdiff(rownames(de.tpm.gene.pcg.exp), "ENSG00000269846"),], 0.6, genes, file.de, file.main)

## Innate immunity
genes <- c("IL6R", "GSK3A", "GSK3B")
file.main <- c("RB1mut-associated gene expression profiling in SCLC", "", "Innate immunity", "")
file.de <- paste0(wd.de, "volcanoplot_all_tpm-gene-pcg-exp_rho6_Innate.pdf")
plotVolcanoALL(de.tpm.gene.pcg.exp[setdiff(rownames(de.tpm.gene.pcg.exp), "ENSG00000269846"),], 0.6, genes, file.de, file.main)

## Volcano plot
## MHC I
genes <- c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F")
file.main <- c("RB1mut-associated expression profiling (LUAD, LCNEC and SCLC)", "", "MHC class I", "")
file.de <- paste0(wd.de, "volcanoplot_all_tpm-gene-pcg-exp_rho6_MHCI.pdf")
plotVolcanoALL(de.tpm.gene.pcg.exp[setdiff(rownames(de.tpm.gene.pcg.exp), "ENSG00000269846"),], 1E-06, genes, file.de, file.main)

## MHC II
genes <- c("CIITA", "HLA-DRA", "HLA-DMA", "HLA-DMB", "HLA-DPA1", "HLA-DPB1", "HLA-DOA", "HLA-DOB", "HLA-DRB1", "HLA-DRB5", "HLA-DQA1", "HLA-DQB1")
file.main <- c("RB1mut-associated expression profiling (LUAD, LCNEC and SCLC)", "", "MHC class II", "")
file.de <- paste0(wd.de, "volcanoplot_all_tpm-gene-pcg-exp_rho6_MHCII.pdf")
plotVolcanoALL(de.tpm.gene.pcg.exp[setdiff(rownames(de.tpm.gene.pcg.exp), "ENSG00000269846"),], 1E-06, genes, file.de, file.main)

# -----------------------------------------------------------------------------
# Replace Ensembl Gene IDs to gene name in Reactome results
# -----------------------------------------------------------------------------
reactome <- read.csv(paste0(wd.all.de, "de-all-tpm-gene/pathway_anova20/result.csv"))
#reactome <- read.csv(paste0(wd.all.de, "de-all-tpm-gene/pathway_rho6/result.csv"))
reactome <- read.csv(paste0(wd.lcnec.de, "pathway_fdr5/result.csv"))
#reactome <- read.csv(paste0(wd.lcnec.de, "de-all-tpm-gene/pathway_fdr10/result.csv"))
#reactome <- read.csv(paste0("/Users/tpyang/Work/uni-koeln/tyang2/ALL/analysis/expression/kallisto/luad-lcnec-sclc-rnaseq-de/de-lung-tpm-gene/", "pathway_p25/result.csv"))

colnames(reactome) <- gsub("X.", "", colnames(reactome))
reactome$Submitted.entities.found <- as.vector(reactome$Submitted.entities.found)
for (r in 1:nrow(reactome)) {
   ids <- as.vector(reactome$Submitted.entities.found[r])
   ids <- unlist(strsplit(ids, ";"))
 
   for (i in 1:length(ids))
      if (nrow(ensGene[ids[i],]) != 0)
         ids[i] <- ensGene[ids[i],]$external_gene_name
   
   reactome$Submitted.entities.found[r] <- paste(ids, collapse=";")
}
writeTable(reactome, paste0(wd.all.de, "de-lcnec-tpm-gene/pathway_fdr5/result.tsv"), colnames=T, rownames=F, sep="\t")

##
file.de <- paste0(wd.reactome, "barplot_all_tpm-gene_rb1_rho6.pdf")
pdf(file.de, height=7, width=6.5)
par(las=1)             ## make label text perpendicular to axis
par(mar=c(5,22,4,2))   ## increase y-axis margin
barplot(-log10(reactome[20:1,]$Entities.pValue), names.arg=reactome[20:1,]$Pathway.name, xlab="Significance (-log10 P-value)", horiz=T)
dev.off()

file.de <- paste0(wd.reactome, "barplot_all_tpm-gene_anova20.pdf")
pdf(file.de, height=7, width=8.7)
par(las=1)             ## make label text perpendicular to axis
par(mar=c(5,31,4,2))   ## increase y-axis margin
barplot(-log10(reactome[20:1,]$Entities.pValue), names.arg=reactome[20:1,]$Pathway.name, xlab="Significance (-log10 P-value)", horiz=T)
dev.off()

file.de <- paste0(wd.reactome, "barplot_lcnec_tpm-gene_rb1_fdr5.pdf")
pdf(file.de, height=7, width=8.9)
par(las=1)             ## make label text perpendicular to axis
par(mar=c(5,34,4,2))   ## increase y-axis margin
barplot(-log10(reactome[20:1,]$Entities.pValue), names.arg=reactome[20:1,]$Pathway.name, xlab="Significance (-log10 P-value)", horiz=T)
dev.off()

file.de <- paste0(wd.reactome, "barplot_lcnec_tpm-gene_rb1_fdr10.pdf")
pdf(file.de, height=7, width=8)
par(las=1)             ## make label text perpendicular to axis
par(mar=c(5,30,4,2))   ## increase y-axis margin
barplot(-log10(reactome[20:1,]$Entities.pValue), names.arg=reactome[20:1,]$Pathway.name, xlab="Significance (-log10 P-value)", horiz=T)
dev.off()

file.de <- paste0(wd.reactome, "barplot_ccle_rpkm-gene_p25.pdf")
pdf(file.de, height=7, width=7.5)
par(las=1)             ## make label text perpendicular to axis
par(mar=c(5,27,4,2))   ## increase y-axis margin
barplot(-log10(reactome[20:1,]$Entities.pValue), names.arg=reactome[20:1,]$Pathway.name, xlab="Significance (-log10 P-value)", horiz=T)
dev.off()

