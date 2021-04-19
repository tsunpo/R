# =============================================================================
# Library      : Transcription
# Name         : handbook-of/Transcription.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 07/01/19
# =============================================================================

# -----------------------------------------------------------------------------
# Methods: Principal component analysis
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
plotPCA <- function(x, y, pca, trait, wd.de.data, file.name, size, file.main, legend, cols, samples, flip.x, flip.y, legend.title=NA) {
   scores <- pcaScores(pca)
   trait[is.na(trait)] <- "NA"
   #trait.v <- sort(unique(trait), decreasing=F)
   trait.v <- c("B", "X", "N")
   
   if (isNA(cols))
      cols <- c("red", "deepskyblue", "forestgreen", "purple3", "blue", "gold", "lightsalmon", "turquoise1", "limegreen")   #, "salmon", "tomato", "steelblue2", "cyan")
   if (length(trait.v) > length(cols))
      cols <- rainbow(length(trait.v))
   else
      cols <- cols[1:length(trait.v)]
   trait.col <- mapply(x = 1:length(trait), function(x) return(cols[which(trait[x] == trait.v)]))   ## Assign colours to each subtypes
   
   cols[which(trait.v == "NA")] <- "lightgray"
   trait.col[which(trait == "NA")] <- "lightgray"
   xlab.txt <- paste0("PC", x, " (", pcaProportionofVariance(pca, x), "%)")
   ylab.txt <- paste0("PC", y, " (", pcaProportionofVariance(pca, y), "%)")
   
   pdf(file.path(wd.de.data, paste0(file.name, "_", names(scores)[x], "-", names(scores)[y], ".pdf")), height=size, width=size)
   #plot(scores[, x]*flip.x, scores[, y]*flip.y, col=trait.col, pch=16, cex=1.5, main=file.main[1], xlab=xlab.txt, ylab=ylab.txt)
   plot(scores[, x]*flip.x, scores[, y]*flip.y, col=NA, pch=19, cex=1.5, main=file.main[1], xlab=xlab.txt, ylab="", cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
   idx <- which(trait == "B")
   points(scores[idx, x]*flip.x, scores[idx, y]*flip.y, col=trait.col[idx], pch=19, cex=1.7)
   idx <- which(trait == "X")
   points(scores[idx, x]*flip.x, scores[idx, y]*flip.y, col=trait.col[idx], pch=19, cex=1.7)
   idx <- which(trait == "N")
   points(scores[idx, x]*flip.x, scores[idx, y]*flip.y, col=trait.col[idx], pch=19, cex=1.7)

   #if (!is.null(samples))
   #   for (s in 1:length(samples)) {
   #      sample <- samples[s]
   #      if (sample == "S125650" || sample == "S125652")
   #        text(scores[sample, x]*flip.x, scores[sample, y]*flip.y, sample, col="black", adj=c(0, -0.5), cex=1)
   #   }
   
   mtext(file.main[2], cex=1.3, line=0.3)
   for (l in 1:length(trait.v))
      trait.v[l] <- trait.v[l]
      #trait.v[l] <- paste0(trait.v[l], " (n=", length(which(trait == trait.v[l])), ")")
   
   #trait.v[5] <- "Others"   ## For PCA ALL
   #if (BASE != "") {
   #   trait.v <- trait.v[1:4]
   #   trait.v <- paste0(BASE, " ", trait.v)
   #   cols <- cols[1:4]
   #}
   if (is.na(legend.title))
      legend(legend, trait.v, col=cols, pch=19, pt.cex=2.5, cex=1.8)   ##bty="n")
   else
      legend(legend, title=legend.title, trait.v, col=cols, pch=19, pt.cex=2.5, cex=1.8)
   mtext(ylab.txt, side=2, line=2.75, cex=1.8)
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

testT <- function(q3, q4) {
   return(t.test(as.numeric(q3), as.numeric(q4))$p.value)
}

testU <- function(q3, q4) {
   trait <- rep(0, length(q3))
   trait <- c(trait, rep(1, length(q4)))
   trait <- as.factor(trait)
 
   expr <- as.numeric(c(q3, q4))
   return(wilcox.test(expr ~ trait, exact=F)$p.value)
}

## Student's t-test
testStudents <- function(expr, samples.mut, samples.wt) {
   expr.mut <- expr[,samples.mut]
   expr.wt  <- expr[,samples.wt]
 
   return(mapply(x = 1:nrow(expr), function(x) t.test(expr.mut[x,], expr.wt[x,])$p.value))
}

## Wilcoxon rank-sum test (Mann-Whitney U test) with continuity correction (logical indicator exact=F)
## http://gdac.broadinstitute.org/runs/analyses__latest/reports/cancer/LUSC-TP/Correlate_Clinical_vs_mRNAseq/nozzle.html
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

getExpressed <- function(expr) {   ## Not expressed (TPM = 0) genes in any of the samples 
   return(mapply(x = 1:nrow(expr), function(x) !any(as.numeric(expr[x,]) == 0)))
}

## Differential expression analysis methods
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
   if (test == "Wilcoxon" || test == "Mann–Whitney" || test == "U" || test == "wilcox.test") {
      de$P <- testWilcoxon(expr.pheno, pheno.expr, predictor)
   } else if (test == "Student's" || test == "t.test") {
      de$P <- testStudents(expr, samples.expr.mut, samples.expr.wt)
   }
   
   ## FDR
   de$FDR <- testFDR(de$P, test.fdr)
   
   ## Log2 fold change
   de[,3] <- median00(expr.pheno, samples.expr.wt)
   de[,4] <- median00(expr.pheno, samples.expr.mut)
   de$LOG2_FC <- de[,4] - de[,3]
 
   ## NOTE: Sort AFTER fold change and BEFORE annotation!!
   de <- de[order(de$P),]
   return(de)
}

## https://www.biostars.org/p/186368/
## https://www.biostars.org/p/157240/
## https://www.biostars.org/p/143458/#157303                            ## estimates vs. count-based methods; transcript level quantification tools
## https://liorpachter.wordpress.com/2015/08/17/a-sleuth-for-rna-seq/   ## estimates vs. count-based methods
## https://liorpachter.wordpress.com/2013/08/26/magnitude-of-effect-vs-statistical-significance/  ## Read counts accumulated across a gene cannot be used directly to estimate fold change 

# -----------------------------------------------------------------------------
# Method: Volcano plot
# Last Modified: 14/04/18
# -----------------------------------------------------------------------------
fdrToP <- function(fdr, de) {
   de.sig <- subset(de, FDR <= fdr)
 
   return(max(de.sig$P))
}

pvalueToFDR <- function(pvalue, de) {
   de.sig <- subset(de, P <= pvalue)
 
   return(round0(max(de.sig$FDR)*100, digits=0))
}

rhoToP <- function(rho, de) {
   de.rho <- rbind(subset(de, RHO >= rho), subset(de, RHO <= rho*-1))
 
   return(max(de.rho$P))
}

getVolcanoGenes <- function(file.tab, de) {
   genes <- readTable(file.tab, header=T, rownames=F, sep="\t")
   rownames(genes) <- genes$GENE
   
   return(genes[intersect(genes$GENE, de$external_gene_name),])
}

# -----------------------------------------------------------------------------
# Method: Fisher's combined probability test
# Last Modified: 06/11/17
# -----------------------------------------------------------------------------
fishers <- function(x, y) {
   return(pchisq(-2*log(x)-2*log(y), 4, low=F))
}

# -----------------------------------------------------------------------------
# Method: Gene Set Enrichment Analysis (GSEA)
# Link(s): https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/FAQ#Where_are_the_GSEA_statistics_.28ES.2C_NES.2C_FDR.2C_FWER.2C_nominal_p_value.29_described.3F
# Last Modified: 16/08/18
# -----------------------------------------------------------------------------
## GSEAPreranked tool
## https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/FAQ#Can_I_use_GSEA_to_analyze_SNP.2C_SAGE.2C_ChIP-Seq_or_RNA-Seq_data.3F
## https://gsea-msigdb.github.io/gsea-gpmodule/v19/index.html

## http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#RNK:_Ranked_list_file_format_.28.2A.rnk.29
writeRNKformat <- function(de.tpm.gene, wd.de.data, file.name) {
   #de.sorted <- de.tpm.gene[order(de.tpm.gene$LOG2_FC, decreasing=T),]
   de.tpm.gene$WEIGHT <- de.tpm.gene$LOG2_FC*(-log10(de.tpm.gene$P))
   de.sorted <- de.tpm.gene[order(de.tpm.gene$WEIGHT, decreasing=T),]
   
   file <- file.path(wd.de.data, paste0(file.name, ".rnk"))
   writeTable(de.sorted[,c("ensembl_gene_id", "WEIGHT")], file, colnames=F, rownames=F, sep="\t")
}

## http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GRP
writeGRPformat <- function(list, wd.de.data, file.name) {
   writeTable(list, file.path(wd.de.data, paste0(file.name, ".grp")), colnames=F, rownames=F, sep="\t")
}

# =============================================================================
# Inner Class: Collections of test/obsolete/deprecated methods
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 01/02/18
# =============================================================================

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
