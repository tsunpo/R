# =============================================================================
# Library      : Differential Gene Expression
# Name         : handbook-of/DifferentialExpression.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 07/01/19
# =============================================================================

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

## Transcript-level estimates with patches/scaffold sequences (*_PATCH)   ## See line 30 in guide-to-the/hg19.R
## https://www.ncbi.nlm.nih.gov/grc/help/patches
list2Matrix <- function(list, kallisto.table) {
   genes <- unique(kallisto.table$target_id)
   samples <- unique(kallisto.table$sample)  
 
   list.matrix <- data.frame(matrix(list, nrow=length(genes), byrow=T))
   rownames(list.matrix) <- genes
   colnames(list.matrix) <- samples
 
   return(list.matrix)
}

# -----------------------------------------------------------------------------
# Method: tpm.gene file manipulations
# Last Modified: 10/04/18
# -----------------------------------------------------------------------------
## Remove patches (*_PATCH)
getGeneTPM <- function(tpm.gene.patch, ensGene) {
   overlaps <- intersect(rownames(tpm.gene.patch), rownames(ensGene))
   tpm.gene <- tpm.gene.patch[overlaps,]
 
   return(tpm.gene)
}

## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5143225/
## https://www.nature.com/articles/nmeth.3580/figures/7
getLog2andMedian <- function(tpm.gene, pseudocount) {
   tpm.gene.log2 <- log2(tpm.gene + pseudocount)
   tpm.gene.log2$MEDIAN <- mapply(x = 1:nrow(tpm.gene.log2), function(x) median(as.numeric(tpm.gene.log2[x,])))
 
   return(tpm.gene.log2)
}

getEnsGeneFiltered <- function(tpm.gene, ensGene, autosomeOnly, proteinCodingOnly) {   ## ADD 11/04/18; ADD 01/02/18
   if (autosomeOnly)
      tpm.gene <- tpm.gene[intersect(rownames(tpm.gene), rownames(subset(ensGene, chromosome_name %in% paste0("chr", 1:22)))),]
   if (proteinCodingOnly)
      tpm.gene <- tpm.gene[intersect(rownames(tpm.gene), rownames(subset(ensGene, gene_biotype == "protein_coding"))),]
 
   return(tpm.gene)
}

# =============================================================================
# Methods: Density plot and histogram
# Last Modified: 25/11/18
# =============================================================================
## http://www.sthda.com/english/wiki/abline-r-function-an-easy-way-to-add-straight-lines-to-a-plot-using-r-software
## http://www.sthda.com/english/wiki/line-types-in-r-lty
plotDensity <- function(medians, BASE, file.name, detected, pseudocount, ymax) {
   xlab.text <- paste0("log2(TPM+", pseudocount, ")")
   ylab.text <- "Density"
   d <- density(medians)
   #d$y <- d$n/sum(d$y) * d$y   ## Convert to counts
   q <- as.numeric(quantile(medians))
   if (is.null(ymax))
      ymax <- max(d$y)
   numbers <- formatC(length(medians), format="f", big.mark=",", digits=0)
   conditioned <- "Expressed"
   if (detected)
      conditioned <- "Detected"
   
   pdf(file.name, height=6, width=6)
   plot(d, xlab=xlab.text, ylab=ylab.text, main=paste0(conditioned, " genes in ", BASE), ylim=c(0, ymax), lwd=1.5)
   abline(v=q, col=c("red", "blue", "blue", "blue", "blue"), lty=c(1, 5, 1, 5, 1), lwd=c(0.85, 0.85, 0.85, 0.85, 0.85))
   for (x in 2:5)
      text((q[x] + q[x-1])/2, (ymax + min(d$y))/2, paste0("Q", (x-1)), cex=0.85, col="blue")
   text(q[1], ymax, "TPM=0", cex=0.85, col="red") 
   text(q[3], ymax, "Median", cex=0.85, col="blue") 
   text(q[5], ymax, "Maximum", cex=0.85, col="blue") 
   
   mtext(paste0("(n=", numbers, ")"), cex=1.2, line=0.3)
   rug(jitter(medians))
   dev.off()
}

## https://homepage.divms.uiowa.edu/~luke/classes/STAT4580/histdens.html
## https://www.r-graph-gallery.com/190-mirrored-histogram/
## https://www.statmethods.net/graphs/density.html
plotHistogram <- function(medians, BASE, file.name, detected, pseudocount, ymax, breaks=15) {
   h <- hist(medians, breaks=breaks) 
   if (is.null(ymax))
      ymax <- max(h$counts)
   numbers <- formatC(length(medians), format="f", big.mark=",", digits=0)
   conditioned <- "Expressed"
   if (detected)
      conditioned <- "Detected"
   
   pdf(file.name, height=6, width=6)
   hist(medians, ylab="Frequency", xlab=paste0("log2(TPM+", pseudocount, ")"), main=paste0(conditioned, " genes in ", BASE), breaks=breaks, ylim=c(0, ymax)) 
   mtext(paste0("(n=", numbers, ")"), cex=1.2, line=0.3)
   dev.off()
}

plotTxDensityHistogram <- function(gene, BASE, wd.rt.plots, tpm.gene.log2, ensGene.rt.tx, tx.q4.fix.all, pseudocount) {
   ensembl_gene_id <- rownames(subset(ensGene.rt.tx, external_gene_name == gene))
   ensGene.rt.tx.gene <- ensGene.rt.tx[ensembl_gene_id,]
   main.text <- paste0(gene, " (", ensembl_gene_id, ") in ", BASE, " (n=", ncol(tpm.gene.log2) - 1, ")")
   xlab.text <- paste0("log2(TPM+", pseudocount, ")")
   for (q in 1:4)
      if (length(intersect(ensembl_gene_id, tx.q4.fix.all[[q]])) == 1)
         q.text <- paste0("Q", q)
   
   if (ensGene.rt.tx.gene$CONSIST == 1) {
      if (ensGene.rt.tx.gene$CD == 1) {
         mtext <- paste0("Co-directional (", q.text, "); RT(")
         if (ensGene.rt.tx.gene$RT == 1)
            mtext <- paste0(mtext, "R), Tx(+)")
         else
            mtext <- paste0(mtext, "L), Tx(-)")
      } else {
         mtext <- paste0("Head-on (", q.text, "); RT(")  
         if (ensGene.rt.tx.gene$RT == 1)
            mtext <- paste0(mtext, "R), Tx(-)")
         else
            mtext <- paste0(mtext, "L), Tx(+)")
      }
   } else {
      mtext <- paste0("Inconsistent RT (", q.text, "); ")
      if (ensGene.rt.tx.gene$strand == 1) {
         mtext <- paste0(mtext, "Tx(+), TSS(")
         if (ensGene.rt.tx.gene$RT_1 == 1)
            mtext <- paste0(mtext, "R), TES(L)")
         else
            mtext <- paste0(mtext, "L), TES(R)")
      } else {
         if (ensGene.rt.tx.gene$RT_1 == 1)
            mtext <- paste0(mtext, "TSS(R), TES(L)")
         else
            mtext <- paste0(mtext, "TSS(L), TES(R)")
         mtext <- paste0(mtext, ", Tx(-)")
      }
   }
   
   medians <- as.numeric(tpm.gene.log2[ensembl_gene_id, -ncol(tpm.gene.log2)])
   d <- density(medians)
   #d$y <- d$n/sum(d$y) * d$y   ## Convert to counts
 
   ## Density plots
   file.name  <- file.path(wd.rt.plots, paste0("density_", base, "_tx_", gene, ".pdf"))
   pdf(file.name, height=6, width=6)
   plot(d, xlab=xlab.text, ylab="Density", main=main.text, lwd=1.5)
   mtext(mtext, cex=1.2, line=0.3)
   rug(jitter(medians))
   dev.off()
 
   ## Histogram
   file.name  <- file.path(wd.rt.plots, paste0("hist_", base, "_tx_", gene, ".pdf"))
   pdf(file.name, height=6, width=6)
   hist(medians, xlab=xlab.text, ylab="Frequency", main=main.text) 
   mtext(mtext, cex=1.2, line=0.3)
   dev.off()
}

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
   trait.v <- sort(unique(trait), decreasing=T)
   
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
   plot(scores[, x]*flip.x, scores[, y]*flip.y, col="lightgray", pch=16, cex=1.5, main=file.main[1], xlab=xlab.txt, ylab="", cex.axis=1.2, cex.lab=1.3, cex.main=1.5)
   idx <- which(trait != "NA")
   points(scores[idx, x]*flip.x, scores[idx, y]*flip.y, col=trait.col[idx], pch=16, cex=1.5, main=file.main[1], xlab=xlab.txt, ylab=ylab.txt)

   if (!is.null(samples))
      for (s in 1:length(samples)) {
         sample <- samples[s]
         
         if (sample == "NGP" || sample == "SKNAS" || sample == "CLBGA" || sample == "LS")
            text(scores[sample, x]*flip.x, scores[sample, y]*flip.y, sample, col="black", adj=c(1, -0.75), cex=1.3)
         else if (sample == "TR14")
            text(scores[sample, x]*flip.x, scores[sample, y]*flip.y, sample, col="black", adj=c(-0.25, 1.25), cex=1.3)
         else
            text(scores[sample, x]*flip.x, scores[sample, y]*flip.y, sample, col="black", adj=c(0, -0.75), cex=1.3)
      }
   
   mtext(file.main[2], cex=1.2, line=0.3)
   for (l in 1:length(trait.v))
      trait.v[l] <- trait.v[l]
      #trait.v[l] <- paste0(trait.v[l], " (n=", length(which(trait == trait.v[l])), ")")
   
   if (is.na(legend.title))
      legend(legend, trait.v, col=cols, pch=16, cex=1.25)   ##bty="n")
   else
      legend(legend, title=legend.title, trait.v, col=cols, pch=16, cex=1.25)
   mtext(ylab.txt, side=2, line=2.85, cex=1.3)
   dev.off()
}
plotPCA(1, 2, pca.de, trait, wd.rt.plots, "PCA_NBL-CL_chrs", size=6, file.main, "topright", c("red", "lightcoral", "lightskyblue3", "blue"), samples1, flip.x=-1, flip.y=1, legend.title=NA)


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
 
   return(round0(max(de.sig$FDR)*100, digits=1))
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
