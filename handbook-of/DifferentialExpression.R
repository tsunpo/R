# =============================================================================
# Library      : Differential Gene Expression
# Name         : handbook-of/DifferentialExpression.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 07/08/18
# =============================================================================
# -----------------------------------------------------------------------------
# Methods: Density plot
# Last Modified: 11/06/18
# -----------------------------------------------------------------------------
getDensityCount <- function(median) {
   d <- density(median)
   d$y <- d$y * d$n
 
   return(d)
}

getMaxDensityCount <- function(median) {
   return(max(getDensityCount(median)$y))
}

plotDensityCount <- function(median, file.main, file.name, ymax) {
   d <- getDensityCount(median)
 
   if (is.null(ymax))
      ymax <- max(d$y)
   
   pdf(file.name, height=6, width=6)
   plot(d, ylab="Frequency", xlab="log2(TPM + 0.01)", main=paste0("Genes (n=", file.main, ")"), ylim=c(0, ymax))
   dev.off()
}

# -----------------------------------------------------------------------------
# Methods: Principal component analysis (PCA)
# Last Modified: 05/02/17
# -----------------------------------------------------------------------------
getPCA <- function(expr) {
   pca <- prcomp(expr, retx=T, center=T, scale=T)   ## Perform PCA using normalised data!!
   #loadings <- pca$rotation
   
   return(pca)
}

pcaScores <- function(pca) {
   return(as.data.frame(pca$x))
}

pcaSummary <- function(pca) {
   return(as.data.frame(summary(pca)[[6]]))
}

pcaProportionofVariance <- function(pca, pc) {
   summary <- pcaSummary(pca)
   
   return(round(summary[2, pc]*100, 1))
}

isNA <- function(input) {
   if (is.na(input) || input == "NA")   ##|| input == "")
      return(T)
   return(F)
}

## http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
plotPCA <- function(x, y, pca, trait, wd.de.data, file.name, file.main, legend.x, legend.y, cols) {
   scores <- pcaScores(pca)
   trait[is.na(trait)] <- "NA"
   trait.v <- sort(unique(trait))
   
   if (isNA(cols))
      cols <- c("red","deepskyblue","forestgreen","purple3","blue","gold","lightsalmon","turquoise1","limegreen")   #,"salmon","tomato","steelblue2","cyan")
   if (length(trait.v) > length(cols))
      cols <- rainbow(length(trait.v))
   else
      cols <- cols[1:length(trait.v)]
   trait.col <- mapply(x = 1:length(trait), function(x) return(cols[which(trait[x] == trait.v)]))   ## Assign colours to each subtypes
   
   cols[which(trait.v == "NA")] <- "lightgrey"
   trait.col[which(trait == "NA")] <- "lightgrey"
   
   xlab <- paste0("Principal component ", x, " (", pcaProportionofVariance(pca, x), "%)")
   ylab <- paste0("Principal component ", y, " (", pcaProportionofVariance(pca, y), "%)")
   
   pdf(file.path(wd.de.data, paste0(file.name, "_", names(scores)[x], "-", names(scores)[y], ".pdf")))
   plot(scores[,x], scores[,y], col=trait.col, pch=16, cex=1.5, main=file.main, xlab=xlab, ylab=ylab)
   
   if (is.na(legend.x))
      legend.x <- min(scores[,x])
   if (is.na(legend.y))
      legend.y <- max(scores[,y])
   legend(legend.x, legend.y, trait.v, col=cols, pch=16, cex=1)   ##bty="n")
   dev.off()
}

plotPCAs <- function(x, y, pca, traits, wd.pca, file.main, legend.x, legend.y, cols) {
   for (t in 1:ncol(traits)) {
      variable <- colnames(traits)[t]
      if (isNA(file.main))
         file.main <- variable
      
      plotPCA(x, y, pca, traits[,t], wd.pca, variable, file.main, legend.x, legend.y, cols)
   }
}

# -----------------------------------------------------------------------------
# Methods: Differential expression analysis
# Last Modified: 06/02/17
# -----------------------------------------------------------------------------
median0 <- function(expr) {
   return(mapply(x = 1:nrow(expr), function(x) median(as.numeric(expr[x,])))) 
}

median00 <- function(expr, samples) {
   return(mapply(x = 1:nrow(expr), function(x) median(as.numeric(expr[x, samples])))) 
}

## Student's t-test
testStudents <- function(expr, samples.mut, samples.wt) {
   expr.mut <- expr[,samples.mut]
   expr.wt  <- expr[,samples.wt]
 
   return(mapply(x = 1:nrow(expr), function(x) t.test(expr.mut[x,], expr.wt[x,])$p.value))
}

## Wilcoxon rank-sum test (Mann-Whitney U test) with continuity correction (logical indicator exact=F)
testWilcoxon <- function(expr, pheno, predictor) {
   trait <- as.factor(pheno[,predictor])
   
   return(mapply(x = 1:nrow(expr), function(x) wilcox.test(as.numeric(expr[x,]) ~ trait, exact=F)$p.value))
}

## As used in all-tpm.de.R
testANOVA <- function(x, expr, pheno) {
   fit1 <- lm(as.numeric(expr[x,]) ~ pheno$Cancer_Type)
   fit2 <- lm(as.numeric(expr[x,]) ~ 1)
   a1 <- anova(fit1, fit2)
 
   return(a1$Pr[2])
}

## source("https://bioconductor.org/biocLite.R")
## biocLite("qvalue")
testFDR <- function(P, test.fdr) {
   if (test.fdr == "Q") {
      library(qvalue)
      return(qvalue(P)$qvalue)
   
   } else if (test.fdr == "BH")
      return(p.adjust(P, "BH"))
}

## Generate final pheno.expr and expr.pheno for D.E. analysis
removeNA <- function(df, colname) {
   return(df[!is.na(df[,colname]),])
}

removeNACovariate <- function(expr, pheno, covariate) {
   pheno.nona <- removeNA(pheno, covariate)
   
   cat(paste0("Original sample size from expression data: ", ncol(expr)), sep="\n")
   cat(paste0("Original sample size from phenotypic data: ", nrow(pheno)), sep="\n")
   cat(paste0("After removed samples with NA covariate ", covariate, ": ", nrow(pheno.nona)), sep="\n")
   return(pheno.nona)
}

findOverlapingSamples <- function(expr, pheno) {
   samples.expr <- intersect(colnames(expr), rownames(pheno))
   
   cat(paste0("Overlapping samples between new phenotypic and expression data: ", length(samples.expr)), sep="\n")   
   return(samples.expr)
}

getFinalPhenotype <- function(expr, pheno, covariate) {
   pheno.nona <- removeNACovariate(expr, pheno, covariate)   ## Even if pheno.expr.nona is inputed (with NAs removed)
   samples.expr <- findOverlapingSamples(expr, pheno.nona)   ## Users can always enter expr.d.e.log2 and let this function to remove samples with NA covariate
   
   return(pheno.nona[samples.expr,])
}

getFinalExpression <- function(expr, pheno.expr) {
   return(expr[,rownames(pheno.expr)])
}

## Main D.E. method
differentialAnalysis <- function(expr, pheno, predictor, predictor.wt, test, test.fdr) {
   pheno.expr <- getFinalPhenotype(expr, pheno, predictor)
   expr.pheno <- getFinalExpression(expr, pheno.expr)

   samples.expr.wt <- rownames(subset0(pheno.expr, predictor, predictor.wt))
   samples.expr.mut <- setdiff(colnames(expr.pheno), samples.expr.wt)
   cat(paste0("Samples with MUT ", predictor, ": ", length(samples.expr.mut)), sep="\n")
   cat(paste0("Samples with  WT ", predictor, ": ", length(samples.expr.wt)), sep="\n")
   
   ## DE
   de <- toTable(0, 4, nrow(expr.pheno), c("P", "FDR", paste0(predictor, "_WT"), predictor))
   rownames(de) <- rownames(expr.pheno)
   if (test == "Wilcoxon" || test == "Wilcox" || test == "U") {
      de$P <- testWilcoxon(expr.pheno, pheno.expr, predictor)
   } else if (test == "Students" || test == "ttest") {
      de$P <- testStudents(expr, samples.expr.mut, samples.expr.wt)
   }
   
   ## FDR
   de$FDR <- testFDR(de$P, test.fdr)
   
   ## Log fold change
   de[,3] <- median00(expr.pheno, samples.expr.wt)
   de[,4] <- median00(expr.pheno, samples.expr.mut)
   de$LOG_FC <- de[,4] - de[,3]
 
   ## NOTE: Must sort AFTER fold change and BEFORE annotation!!
   de <- de[order(de$P),]
   return(de)
}

# -----------------------------------------------------------------------------
# Pipeline: Differential expression analysis
# Last Modified: 28/03/17
# -----------------------------------------------------------------------------
pipeDE <- function(expr, pheno, argv, ensGene) {
   annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
 
   de <- differentialAnalysis(expr, pheno, argv$predictor, argv$predictor.wt, argv$test, argv$test.fdr)
   de <- cbind(annot[rownames(de),], de)
 
   return(de)
}

# -----------------------------------------------------------------------------
# Method: Volcano plots
# Last Modified: 13/02/17
# -----------------------------------------------------------------------------
plotVolcano <- function(de, p, fdr, effect, file.name, file.main, legend.x) {
   de$log10P <- -log10(de$P)
   xmax <- max(de$Effect)
   ymax <- max(de$log10P)
   if (is.null(p))
      p <- max(significantAndVariable(de, effect, fdr)$P)
 
   pdf(paste0(file.de, "_Effect>", effect, ".pdf"))
   plot(de$Effect, de$log10P, xlim=c(-xmax, xmax), ylim=c(0, ymax), xlab="Effect size (log2 FC)", ylab="Significance (-log10 P)", col="darkgray", main=file.main)
 
   de.up <- subset(subset(de, FDR <= fdr), Effect >= effect)   ## UPDATE 13/02/17: Not using P
   points(de.up$Effect, de.up$log10P, col="red")
 
   de.down <- subset(subset(de, FDR <= fdr), Effect <= -effect)   ## UPDATE 13/02/17: Not using P
   points(de.down$Effect, de.down$log10P, col="blue")

   abline(h=c(-log10(p)), lty=5)
   abline(v=c(effect), col="red", lty=5)
   abline(v=c(-effect), col="blue", lty=5)

   if (is.na(legend.x))
      legend.x <- -xmax
   legend(x=legend.x, y=ymax, legend=c("Overexpressed", "Underexpressed"), col=c("red", "blue"), pch=21)
   dev.off()
}

# =============================================================================
# Inner Class  : kallisto/sleuth File Reader (for manuscripts/expression/*-tpm.R)
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 10/06/18
# =============================================================================
# -----------------------------------------------------------------------------
# Methods: Custom filters for using kallisto
# Last Modified: 10/06/18
# -----------------------------------------------------------------------------
getFiltered <- function(reads, min_reads, min_prop) { 
   numbers <- mapply(x = 1:nrow(reads), function(x) length(which(reads[x,] >= min_reads)))
 
  return(which(numbers/ncol(reads) >= min_prop))
}

kallisto_table_to_matrix <- function(kallisto.table, min_reads, min_prop) {
   reads <- list2Matrix(kallisto.table$scaled_reads_per_base, kallisto.table)
   tpm <- list2Matrix(kallisto.table$tpm, kallisto.table)
 
   return(tpm[getFiltered(reads, min_reads, min_prop),])
}

# -----------------------------------------------------------------------------
# Methods: Transcript-level TPM estimates using sleuth
# Last Modified: 28/03/17
# -----------------------------------------------------------------------------
tx2gene <- function(t2g) {
   colnames(t2g) <- c("target_id", "ens_gene", "ext_gene")
   #rownames(t2g) <- t2g$target_id
   
   return(t2g)
}

tx2Ens <- function(ensGene) {
   return(tx2gene(ensGene[,c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name")]))
}

## Gene-level TPMs with patches/scaffold sequences (*_PATCH)   ## Please refer to line 64 (in guide-to-the/hg19.R)
## https://www.ncbi.nlm.nih.gov/grc/help/patches
list2Matrix <- function(list, kallisto.table) {
   genes <- unique(kallisto.table$target_id)
   samples <- unique(kallisto.table$sample)  
 
   list.matrix <- data.frame(matrix(list, nrow=length(genes), byrow=T))
   rownames(list.matrix) <- genes
   colnames(list.matrix) <- samples
   
   return(list.matrix)
}

## Remove patches (*_PATCH)
getTPMGene <- function(tpm.gene.patch) {
   overlaps <- intersect(rownames(tpm.gene.patch), rownames(ensGene))
   tpm.gene <- tpm.gene.patch[overlaps,]
   
   return(tpm.gene)
}

getLog2andMedian <- function(tpm.gene) {
   tpm.gene.log2 <- log2(tpm.gene + 0.01)
   tpm.gene.log2$MEDIAN <- mapply(x = 1:nrow(tpm.gene.log2), function(x) median(as.numeric(tpm.gene.log2[x,])))
 
   return(tpm.gene.log2)
}

# =============================================================================
# Inner Class  : Extentions of guide-to-the/hg19.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 10/04/18
# =============================================================================
# -----------------------------------------------------------------------------
# Method: Get non-redundant gene list (ensGene)
# Last Modified: 10/04/18
# -----------------------------------------------------------------------------
getEnsGeneFiltered <- function(tpm.gene, ensGene, autosomeOnly, proteinCodingOnly, proteinCodingNonRedundantOnly) {   ## ADD 11/04/18; ADD 01/02/18
   if (autosomeOnly)
      tpm.gene <- tpm.gene[intersect(rownames(tpm.gene), rownames(subset(ensGene, chromosome_name %in% paste0("chr", 1:22)))),]
   if (proteinCodingOnly)
      tpm.gene <- tpm.gene[intersect(rownames(tpm.gene), rownames(subset(ensGene, gene_biotype == "protein_coding"))),]
   if (proteinCodingNonRedundantOnly) {
      tpm.gene <- tpm.gene[intersect(rownames(tpm.gene), rownames(subset(ensGene, protein_coding_non_redundant == T))),]
   }
 
   return(tpm.gene)
}

# -----------------------------------------------------------------------------
# Methods: Add ensGene$protein_coding_non_redundant (in guide-to-the/hg19.R)
# Last Modified: 24/04/17
# -----------------------------------------------------------------------------
isNonRedundant <- function(ensembl_gene_id, ensGene) {
   gene <- ensGene[ensembl_gene_id,]
 
   ensGene.chr <- subset(ensGene, chromosome_name == gene$chromosome_name)
   ensGene.chr.start <- subset(ensGene.chr, end_position >= gene$start_position)
   ensGene.chr.start.end <- subset(ensGene.chr.start, start_position <= gene$end_position)
 
   if (nrow(ensGene.chr.start.end) == 1)
      return(T)
   return(F)
}
#ensGene.pcg <- subset(ensGene, gene_biotype == "protein_coding")
#ensGene.pcg$non_redundant <- F
#ensGene.pcg$non_redundant <- mapply(x = 1:nrow(ensGene.pcg), function(x) isNonRedundant(ensGene.pcg$ensembl_gene_id[x], ensGene.pcg))
## > nrow(ensGene.pcg)
## [1] 20327
## > nrow(subset(ensGene.pcg, non_redundant == T))
## [1] 13487
## > nrow(subset(ensGene.pcg, non_redundant == F))
## [1] 6840

#ensGene$protein_coding_non_redundant <- NA
#ensGene[rownames(ensGene.pcg),]$protein_coding_non_redundant <- ensGene.pcg$non_redundant

# =============================================================================
# Inner Class: Collections of test/obsolete/deprecated methods
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 01/02/18
# =============================================================================
# -----------------------------------------------------------------------------
# Methods: GTEx threshold to determine a detected/expressed gene
# Link: https://www.ncbi.nlm.nih.gov/pubmed/25954002
# Last Modified: 03/02/17
# -----------------------------------------------------------------------------
# getDetected <- function(expr, value) {   ## Genes being expressed at TPM > 0.1 in at least one sample
#    detected <- mapply(x = 1:nrow(expr), function(x) length(which(as.numeric(expr[x,]) > value)))
#  
#    return(which(detected >= 1))
# }

# getExpressed <- function(expr) {   ## Not expressed (TPM = 0) genes in any of the samples 
#    return(mapply(x = 1:nrow(expr), function(x) !any(as.numeric(expr[x,]) == 0)))
# }

# -----------------------------------------------------------------------------
# Method: Fisher's combined probability test
# Last Modified: 06/11/17
# -----------------------------------------------------------------------------
# fishers <- function(x, y) {
#    return(pchisq(-2*log(x)-2*log(y), 4, low=F))
# }

# -----------------------------------------------------------------------------
# Method: Residuals of expression
# Last Modified: 06/02/17
# -----------------------------------------------------------------------------
# residualsOf <- function(expr, pheno, covariate) {
#    pheno.expr <- getFinalPhenotype(expr, pheno, covariate)
#    expr.pheno <- getFinalExpression(expr, pheno.expr)
#  
#    trait <- as.factor(pheno.expr[,covariate])
#  
#    expr.res <- expr.pheno
#    for (x in 1:nrow(expr.pheno))
#       expr.res[x,] <- as.vector(as.numeric(resid(lm(as.numeric(expr.pheno[x,]) ~ trait))))
#  
#    return(expr.res)
# }
