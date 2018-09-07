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
plotPCA <- function(x, y, pca, trait, wd.de.data, file.name, size, file.main, legend, cols, samples, flip.x, flip.y) {
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
   
   pdf(file.path(wd.de.data, paste0(file.name, "_", names(scores)[x], "-", names(scores)[y], ".pdf")), height=size, width=size)
   plot(scores[,x]*flip.x, scores[,y]*flip.y, col=trait.col, pch=16, cex=1.5, main=file.main, xlab=xlab, ylab=ylab)
   
   if (!is.null(samples)) {
      for (s in 1:length(samples))
         text(scores[s, x]*flip.x, scores[s, y]*flip.y, samples[s], col="black", adj=c(0, -0.75), cex=0.75)
   }

   legend(legend, trait.v, col=cols, pch=16, cex=1)   ##bty="n")
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
   de$LOG2_FC <- de[,4] - de[,3]
 
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

## https://www.biostars.org/p/186368/

## https://www.biostars.org/p/157240/
## https://www.biostars.org/p/143458/#157303                            ## estimates vs. count-based methods; transcript level quantification tools
## https://liorpachter.wordpress.com/2015/08/17/a-sleuth-for-rna-seq/   ## estimates vs. count-based methods
## https://liorpachter.wordpress.com/2013/08/26/magnitude-of-effect-vs-statistical-significance/  ## Read counts accumulated across a gene cannot be used directly to estimate fold change 
 
# -----------------------------------------------------------------------------
# Method: Volcano plots
# Last Modified: 13/02/17
# -----------------------------------------------------------------------------
fdrToP <- function(fdr, de) {
   de.sig <- subset(de, FDR <= fdr)
   de.sig$log10P <- -log10(de.sig$P)
 
   return(max(de.sig$P))
}

plotVolcano <- function(de, fdr, genes, file.de, file.main) {
   de.sig <- subset(de, FDR <= fdr)
   de.sig$log10P <- -log10(de.sig$P)
 
   de$log10P <- -log10(de$P)
   xmax <- max(de$LOG2_FC)
   ymax <- max(de$log10P)
   p <- max(de.sig$P)
 
   ##
   pdf(file.de, height=7, width=7)
   plot(de$LOG2_FC, de$log10P, pch=16, xlim=c(-xmax, xmax), ylim=c(0, ymax), xlab="Median fold change (log2FC RB1/WT)", ylab="Significance (-log10 P-value)", col="darkgray", main=file.main)
 
   abline(h=c(-log10(fdrToP(fdr, de))), lty=5)
   text(xmax*-1 + 2*xmax/80, -log10(fdrToP(fdr, de)) + ymax/42, "Q<0.05", cex=0.85)
   abline(h=c(-log10(fdrToP(0.1, de))), lty=5, col="darkgray")
   text(xmax*-1 + 2*xmax/200, -log10(fdrToP(0.1, de)) + ymax/42, "Q<0.1", col="darkgray", cex=0.85)
   
   de.up   <- subset(de.sig, LOG2_FC > 0)
   points(de.up$LOG2_FC, de.up$log10P, pch=16, col="red")
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
 
   legend("topleft", legend=c("Up-regulation", "Down-regulation"), col=c("red", "dodgerblue"), pch=19)
   dev.off()
}

# -----------------------------------------------------------------------------
# Method: Gene Set Enrichment Analysis (GSEA)
# Last Modified: 16/08/18
# -----------------------------------------------------------------------------
## GSEAPreranked tool
## https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/FAQ#Can_I_use_GSEA_to_analyze_SNP.2C_SAGE.2C_ChIP-Seq_or_RNA-Seq_data.3F
## https://gsea-msigdb.github.io/gsea-gpmodule/v19/index.html

## http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#RNK:_Ranked_list_file_format_.28.2A.rnk.29
writeRNKformat <- function(de.tpm.gene, wd.de.data, file.name) {
   de.sorted <- de.tpm.gene[order(de.tpm.gene$LOG2_FC, decreasing=T),]
   file <- file.path(wd.de.data, paste0(file.name, ".rnk"))
   writeTable(de.sorted[,c("ensembl_gene_id", "LOG2_FC")], file, colnames=F, rownames=F, sep="\t")
}

## http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GRP
writeGRPformat <- function(list, wd.de.data, file.name) {
   writeTable(list, file.path(wd.de.data, paste0(file.name, ".grp")), colnames=F, rownames=F, sep="\t")
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

getExpressed <- function(expr) {   ## Not expressed (TPM = 0) genes in any of the samples 
   return(mapply(x = 1:nrow(expr), function(x) !any(as.numeric(expr[x,]) == 0)))
}

# =============================================================================
# Inner Class: Collections of test/obsolete/deprecated methods
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 01/02/18
# =============================================================================

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
