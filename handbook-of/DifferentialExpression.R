# =============================================================================
# Library: Handbook of Differential Gene Expression
# Name: handbook-of/DifferentialExpression.R
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 24/04/17
# =============================================================================

# -----------------------------------------------------------------------------
# Methods: GTEx threshold on determining a detected/expressed gene
# Link: https://www.ncbi.nlm.nih.gov/pubmed/25954002
# Last Modified: 03/02/17
# -----------------------------------------------------------------------------
getDetectedGTEx <- function(expr, value) {   ## Genes being expressed at FPKM > 0.1 in at least one sample
   detected <- mapply(x = 1:nrow(expr), function(x) length(which(as.numeric(expr[x,]) > value)))
   
   return(which(detected >= 1))
}

getExpressedGTEx <- function(expr) {   ## Not expressed (FPKM = 0) genes in any of the samples 
   return(mapply(x = 1:nrow(expr), function(x) !any(as.numeric(expr[x,]) == 0)))
}

getNotExpressed <- function(expr) {   ## Not expressed (TPM = 0) genes in any of the samples 
   return(mapply(x = 1:nrow(expr), function(x) any(as.numeric(expr[x,]) == 0)))
}

getNumberOfNotExpressed <- function(expr) {   ## Not expressed (TPM = 0) genes in half of the samples 
   return(mapply(x = 1:nrow(expr), function(x) length(which(expr[x,] == 0))))
}

getNotExpressedGTEx <- function(expr) {   ## Not expressed (TPM = 0) genes in half of the samples 
   numbers <- mapply(x = 1:nrow(expr), function(x) length(which(expr[x,] == 0)))
   
   return(which(numbers > ncol(expr) - 10))
}

getNotExpressedInHalf <- function(expr) {   ## Not expressed (TPM = 0) genes in half of the samples 
   numbers <- mapply(x = 1:nrow(expr), function(x) length(which(expr[x,] == 0)))
 
   return(which(numbers > ncol(expr)/2))
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

## http://sape.inf.usi.ch/quick-reference/ggplot2/colour
plotPCA <- function(x, y, pca, trait, wd.pca, variable, file.main, legend.x, legend.y, cols) {
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
   
   pdf(paste0(wd.pca, "pca_", variable, "_", names(scores)[x], "-", names(scores)[y], ".pdf"))
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
# Methods: Differential expression
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

## source("https://bioconductor.org/biocLite.R")
## biocLite("qvalue")
testFDR <- function(P, test.fdr) {
   if (test.fdr == "Q") {
      library(qvalue)
      return(qvalue(P)$qvalue)
   
   } else if (test.fdr == "BH")
      return(p.adjust(P, "BH"))
}

sortP <- function(de) {
   return(de[order(de$P),])
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
   
   ## Effect size
   de[,3] <- median00(expr.pheno, samples.expr.wt)
   de[,4] <- median00(expr.pheno, samples.expr.mut)
   de$Effect <- de[,4] - de[,3]
 
   ## NOTE: Must sort AFTER effect size and BEFORE annotation!!
   de <- sortP(de)
   #return(cbind(refGene[rownames(de), c("name2", "chrom", "strand", "txStart", "txEnd")], de))   ## CHANGE 08/03/17: For MAGIC pipeline
   return(de)
}

onlyVariable <- function(de, effect) {
   de.effect <- rbind(subset(de, Effect >= effect), subset(de, Effect <= -effect))

   return(sortP(de.effect))
}

significantAndVariable <- function(de, effect, fdr) {
   return(subset(onlyVariable(de, effect), FDR <= fdr))
}

# -----------------------------------------------------------------------------
# Method: Residuals of expression
# Last Modified: 06/02/17
# -----------------------------------------------------------------------------
residualsOf <- function(expr, pheno, covariate) {
   pheno.expr <- getFinalPhenotype(expr, pheno, covariate)
   expr.pheno <- getFinalExpression(expr, pheno.expr)
   
   trait <- as.factor(pheno.expr[,covariate])
   
   expr.res <- expr.pheno
   for (x in 1:nrow(expr.pheno))
      expr.res[x,] <- resid(lm(as.numeric(expr.pheno[x,]) ~ trait))
 
   return(expr.res)
}

# -----------------------------------------------------------------------------
# Method: Volcano plots
# Last Modified: 13/02/17
# -----------------------------------------------------------------------------
plotVolcano <- function(de, fdr, effect, file.de, file.main, legend.x) {
   de$log10P <- -log10(de$P)
   xmax <- max(de$Effect)
   ymax <- max(de$log10P)
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

# -----------------------------------------------------------------------------
# Method: Heatmap
# Last Modified: 17/02/17
# -----------------------------------------------------------------------------
## install.packages("gplots")
plotHeatmap <- function(P, test.fdr) {}

# -----------------------------------------------------------------------------
# Methods: Miscellaneous/shortcuts
# Last Modified: 28/04/17
# -----------------------------------------------------------------------------
getEnsGene <- function(name) {
   return(subset(ensGene, external_gene_name == name))
}

getEnsGeneID <- function(name) {
   return(subset(ensGene.gene, external_gene_name == name)$ensembl_gene_id)
}


# =============================================================================
# Inner Class: kallisto/Sleuth FileReader
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 24/04/17
# =============================================================================

# -----------------------------------------------------------------------------
# Methods: Transcript TPM estimates/aboundants from kallisto
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

tx2Ref <- function(refGene) {
   return(tx2gene(refGene[,c("name", "name2", "name2")]))
}

orderBySamples <- function(tpm) {
   return(tpm[,order(colnames(tpm), decreasing=F)])
}

# -----------------------------------------------------------------------------
# Method: Gene-level TPM estimates/aboundants from kallisto
# Last Modified: 28/03/17
# -----------------------------------------------------------------------------
getGeneTPM <- function(t2g, tpm) {
   t2g.tpm <- t2g[rownames(tpm),]      ## Keep only filtered transcripts
   genes <- unique(t2g.tpm$ens_gene)
 
   tpm.gene <- toTable(0, ncol(tpm), length(genes), colnames(tpm))
   rownames(tpm.gene) <- genes
   for (g in 1:length(genes)) {
      transcripts <- rownames(subset(t2g.tpm, ens_gene == genes[g]))   ## Find corresponding transcripts
      tpm.transcripts <- tpm[transcripts,]                             ## and sum up TPM estimates of each gene
  
      if (length(transcripts) == 1)
         tpm.gene[g,] <- tpm.transcripts
      else
         tpm.gene[g,] <- mapply(s = 1:ncol(tpm.gene), function(s) sum(tpm.transcripts[,s]))
   }
 
   return(tpm.gene)
}


# =============================================================================
# Inner Class: Collections of test/obsolete/deprecated methods
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 24/04/17
# =============================================================================

# -----------------------------------------------------------------------------
# Methods: GTEx threshold to determine a gene to be detected/expressed
# Link: https://www.ncbi.nlm.nih.gov/pubmed/25954002
# Last Modified: 03/02/17
# -----------------------------------------------------------------------------
reportZeroTPM <- function(tpm) {
   zeros <- toTable(0, 2, nrow(tpm), c("ID", "Zero"))
   rownames(zeros) <- rownames(tpm)
   zeros$Zero <- mapply(x = 1:nrow(tpm), function(x) length(which(as.numeric(tpm[x,]) == 0)))
 
   zeros <- zeros[order(zeros$Zero, decreasing=T),]
   return(zeros)
}

reportGTExThreshold <- function(genes) {
   reports <- toTable(0, 5, length(genes), c("RefGene", "Raw", "Final", "Not_Detected", "Not_Expressed"))
   rownames(reports) <- genes
   notDetected <- setdiff(rownames(expr), rownames(expr.d))
   notExpressed <- setdiff(rownames(expr.d), rownames(expr.d.e))
 
   for (g in 1:length(genes)) {
      gene <- genes[g]
      transcripts <- rownames(getRefGene("name2", gene))
      reports$RefGene[g] <- length(transcripts)
  
      transcripts.raw <- intersect(transcripts, rownames(expr))
      reports$Raw[g] <- length(transcripts.raw)
  
      reports$Final[g] <- length(intersect(transcripts.raw, rownames(expr.d.e)))
      reports$Not_Detected[g] <- length(intersect(transcripts.raw, notDetected))
      reports$Not_Expressed[g] <- length(intersect(transcripts.raw, notExpressed))
   }
 
   return(reports)
}

# -----------------------------------------------------------------------------
# Method: Modification of sleuth_to_matrix
# Tricks: To be able to output obs_norm_filt
# Last Modified: 28/03/17
# -----------------------------------------------------------------------------
sleuth_to_matrix0 <- function(obj, which_df, which_units) {
   data <- as.data.frame(obj[[which_df]])
 
   res <- list()
 
   s_data <- data %>%
   select_("target_id", "sample", which_units) %>% tidyr::spread_("sample", which_units)
   rownames(s_data) <- s_data$target_id
   s_data$target_id <- NULL
   s_data <- as.matrix(s_data)
   s_data <- s_data[, sample(1:ncol(s_data))]
   res[["data"]] <- s_data
 
   condition_order <- match(colnames(s_data), as.character(obj$sample_to_condition$sample))
   res[["condition"]] <- obj$sample_to_condition$condition[condition_order]
 
   res
}
