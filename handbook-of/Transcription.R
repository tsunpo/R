# =============================================================================
# Library      : Transcription
# Name         : handbook-of/Transcription.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 07/01/19
# =============================================================================
library(qvalue)

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

removeMedian0 <- function(tpm.gene, tpm=0) {
   medians <- mapply(x = 1:nrow(tpm.gene), function(x) median(as.numeric(tpm.gene[x,])))
   tpm.gene <- tpm.gene[which(medians > tpm),]
 
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
plotDensity <- function(medians, BASE, file.name, title, pseudocount=1, ymax) {
   xlab.text <- paste0("log2(TPM+", pseudocount, ")")
   ylab.text <- "Density"
   d <- density(medians)
   #d$y <- d$n/sum(d$y) * d$y   ## Convert to counts
   q <- as.numeric(quantile(medians))
   if (is.null(ymax))
      ymax <- max(d$y)
   numbers <- formatC(length(medians), format="f", big.mark=",", digits=0)

   pdf(file.name, height=6, width=6)
   plot(d, xlab=xlab.text, ylab=ylab.text, main=paste0(title, " genes in ", BASE), ylim=c(0, ymax), lwd=1.5)
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

plotDensity0 <- function(medians, BASE, file.name, title, pseudocount=1, ymax, tpm=NA) {
   xlab.text <- expression("log" * ""[2] * "(TPM + 1)")
   ylab.text <- "Density"
   d <- density(medians)
   #d$y <- d$n/sum(d$y) * d$y   ## Convert to counts
   q <- as.numeric(quantile(medians))
   if (is.null(ymax))
      ymax <- max(d$y)
   numbers <- formatC(length(medians), format="f", big.mark=",", digits=0)
 
   pdf(file.name, height=6, width=6)
   par(mar=c(5.1, 4.7, 4.1, 1.4))
   plot(d, xlab=xlab.text, ylab=ylab.text, main=paste0(numbers, " ", title, " genes in ", BASE), ylim=c(0, ymax), cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
   abline(v=q, col=c("black", "black", "black", "black", "black"), lty=c(1, 5, 1, 5, 1), lwd=c(1, 1, 1, 1, 1))
   #for (x in 2:5)
   #   text((q[x] + q[x-1])/2, (ymax + min(d$y))/2, paste0("Q", (x-1)), cex=0.85, col="blue")
   if (!is.na(tpm)) {
      if (tpm == "r5p47")
   	     text(q[4] + 0.5, ymax, ">5 reads in >47% samples", cex=1.8, col="black")
      else
      	  text(q[3], ymax, paste0("Median TPM > ", tpm), cex=1.8, col="black")
   }
   #text(q[3], ymax, "Median", cex=1, col="black") 
   #text(q[5], ymax, "Maximum", cex=0.85, col="blue") 
 
   #mtext(paste0("n=", numbers), cex=1.4, line=0.25)
   rug(jitter(medians))
   dev.off()
}

## https://homepage.divms.uiowa.edu/~luke/classes/STAT4580/histdens.html
## https://www.r-graph-gallery.com/190-mirrored-histogram/
## https://www.statmethods.net/graphs/density.html
plotHistogram <- function(medians, BASE, file.name, title, pseudocount, ymax, breaks=15) {
   h <- hist(medians, breaks=breaks) 
   if (is.null(ymax))
      ymax <- max(h$counts)
   numbers <- formatC(length(medians), format="f", big.mark=",", digits=0)
   
   pdf(file.name, height=6, width=6)
   hist(medians, ylab="Frequency", xlab=paste0("log2(TPM+", pseudocount, ")"), main=paste0(title, " genes in ", BASE), breaks=breaks, ylim=c(0, ymax)) 
   mtext(paste0("(n=", numbers, ")"), cex=1.2, line=0.3)
   dev.off()
}

plotDensityHistogram <- function(tpm.gene, file.main, title, tpm=0) {
   pcs <- c(1, 0.1, 0.01)
   for (c in 1:1) {
      pc <- pcs[c]
  
      tpm.gene.log2 <- getLog2andMedian(tpm.gene, pseudocount=pc)
      plotDensity0(tpm.gene.log2$MEDIAN, BASE, paste0(file.main, "_density_pc", pc, ".pdf"), title, pseudocount=pc, NULL, tpm=tpm)
      #plotHistogram(tpm.gene.log2$MEDIAN, BASE, paste0(file.main, "_hist_pc", pc, ".pdf"), title, pseudocount=pc, NULL)
   }
}

##
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
plotPCA <- function(x, y, pca, trait, wd.de.data, file.name, size, file.main, legend, legend.title, cols, flip.x, flip.y) {
   scores <- pcaScores(pca)
   trait[is.na(trait)] <- "NA"
   trait.v <- legend.title
   
   #if (isNA(cols))
   #   cols <- c("red", "deepskyblue", "forestgreen", "purple3", "blue", "gold", "lightsalmon", "turquoise1", "limegreen")   #, "salmon", "tomato", "steelblue2", "cyan")
   #if (length(trait.v) > length(cols))
   #   cols <- rainbow(length(trait.v))
   #else
   #   cols <- cols[1:length(trait.v)]
   trait.col <- mapply(x = 1:length(trait), function(x) return(cols[which(trait[x] == trait.v)]))   ## Assign colours to each subtypes
   
   cols[which(trait.v == "NA")] <- "lightgray"
   trait.col[which(trait == "NA")] <- "lightgray"
   xlab.txt <- paste0("PC", x, " (", pcaProportionofVariance(pca, x), "%)")
   ylab.txt <- paste0("PC", y, " (", pcaProportionofVariance(pca, y), "%)")
   
   png(file.path(wd.de.data, paste0(file.name, "_", names(scores)[x], "-", names(scores)[y], ".png")), units="in", res=300, height=size, width=size)
   par(mar=c(5.1, 4.7, 4.1, 1.4))
   plot(scores[, x]*flip.x, scores[, y]*flip.y, col=trait.col, pch=19, cex=2, main=file.main[1], xlab=xlab.txt, ylab=ylab.txt, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
   #plot(scores[, x]*flip.x, scores[, y]*flip.y, col=NA, pch=19, cex=1.5, main=file.main[1], xlab=xlab.txt, ylab="", cex.axis=1.8, cex.lab=1.9, cex.main=2.1)

   #idx <- which(trait == "Non-FS")
   #points(scores[idx, x]*flip.x, scores[idx, y]*flip.y, col=trait.col[idx], pch=19, cex=2)
   #idx <- which(trait == "PCAWG_FS")
   #points(scores[idx, x]*flip.x, scores[idx, y]*flip.y, col=trait.col[idx], pch=19, cex=2)
   #idx <- which(trait == "GEL_FS")
   #points(scores[idx, x]*flip.x, scores[idx, y]*flip.y, col=trait.col[idx], pch=19, cex=2)
   
   #idx <- which(trait == "Non-FS")
   #points(scores[idx, x]*flip.x, scores[idx, y]*flip.y, col=trait.col[idx], pch=19, cex=2)
   #idx <- which(trait == "PCAWG_FS_Late")
   #points(scores[idx, x]*flip.x, scores[idx, y]*flip.y, col=trait.col[idx], pch=19, cex=2)
   #idx <- which(trait == "PCAWG_FS_Very_Late")
   #points(scores[idx, x]*flip.x, scores[idx, y]*flip.y, col=trait.col[idx], pch=19, cex=2)

   #idx <- which(trait == "Non-FS")
   #points(scores[idx, x]*flip.x, scores[idx, y]*flip.y, col=trait.col[idx], pch=19, cex=2)
   #idx <- which(trait == "FS")
   #points(scores[idx, x]*flip.x, scores[idx, y]*flip.y, col=trait.col[idx], pch=19, cex=2)
   
   idx <- which(trait == "Early")
   points(scores[idx, x]*flip.x, scores[idx, y]*flip.y, col=trait.col[idx], pch=19, cex=2)
   idx <- which(trait == "Late")
   points(scores[idx, x]*flip.x, scores[idx, y]*flip.y, col=trait.col[idx], pch=19, cex=2)
   
   #mtext(file.main[2], cex=1.2, line=0.3)
   for (l in 1:length(trait.v))
   	  #trait.v[l] <- paste0(trait.v[l], " (n=", length(which(trait == trait.v[l])), ")")
      trait.v[l] <- trait.v[l]
   
   #trait.v[5] <- "Others"   ## For PCA ALL
   #if (BASE != "") {
   #   trait.v <- trait.v[1:4]
   #   trait.v <- paste0(BASE, " ", trait.v)
   #   cols <- cols[1:4]
   #}
   #if (is.na(legend.title))
   #   legend(legend, trait.v, col=cols, pch=19, pt.cex=2.5, cex=1.8)   ##bty="n")
   #else
      legend(legend, trait.v, col=cols, pch=19, pt.cex=2.5, cex=1.6, bg="transparent")
   #mtext(ylab.txt, side=2, line=2.75, cex=1.9)
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
## BiocManager::install("qvalue")
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

removeMissing <- function(expr) {
   return(mapply(x = 1:nrow(expr), function(x) !any(is.na(expr[x,]))))
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

adjustcolor.red   <- adjustcolor(red, alpha.f=0.15)
adjustcolor.blue  <- adjustcolor(blue, alpha.f=0.15)
adjustcolor.green <- adjustcolor(green, alpha.f=0.25)
adjustcolor.yellow <- adjustcolor(yellow, alpha.f=0.25)
adjustcolor.skyblue <- adjustcolor("skyblue", alpha.f=0.25)
adjustcolor.black  <- adjustcolor("black", alpha.f=0.15)
adjustcolor.dimgray  <- adjustcolor("dimgray", alpha.f=0.25)
adjustcolor.darkgray  <- adjustcolor("darkgray", alpha.f=0.25)
adjustcolor.lightgray  <- adjustcolor("lightgray", alpha.f=0.25)
adjustcolor.white   <- adjustcolor("white", alpha.f=0)

## Bottom-right	-0.12	1.25
##  Level-right	-0.05	
plotVolcano <- function(de, pvalue, genes, file.de, file.main, xlab.text, ylab.text, legend, legends, cols, cols2, fold=1, ymax=0) {
   de.sig <- subset(de, P <= pvalue)
   de.sig$log10P <- -log10(de.sig$P)
 
   de$log10P <- -log10(de$P)
   xmax <- max(de$LOG2_FC)
   xmin <- min(de$LOG2_FC)
   if (ymax == 0) ymax <- max(de$log10P)
   #ymax <- 7
 
   pdf(file.de, height=6, width=6)
   par(mar=c(5.1, 4.7, 4.1, 1.4), xpd=F)
   plot(de$LOG2_FC, de$log10P, pch=16, xlim=c(xmin, xmax), ylim=c(0, ymax), xlab=xlab.text, ylab=ylab.text, col="lightgray", main=file.main[1], cex=1.7, cex.axis=1.6, cex.lab=1.7, cex.main=1.8)
   abline(h=c(-log10(pvalue)), lty=5, lwd=2)
 
   de.up   <- subset(de.sig, LOG2_FC > fold)
   points(de.up$LOG2_FC, de.up$log10P, pch=16, col=cols[1], cex=1.7)
   de.down <- subset(de.sig, LOG2_FC < -fold)
   points(de.down$LOG2_FC, de.down$log10P, pch=16, col=cols[2], cex=1.7)
   abline(v=fold, lty=5, col=red, lwd=2)
   abline(v=-fold, lty=5, col=blue, lwd=2)
 
   par(xpd=T)
   if (nrow(genes) != 0) {
      for (g in 1:nrow(genes)) {
         gene <- subset(de, external_gene_name == genes[g,]$GENE)
         gene <- cbind(gene, genes[g,])
   
         if (nrow(gene) > 0) {
            points(gene$LOG2_FC, gene$log10P, pch=1, col="black", cex=1.7)
    
            #if (!is.na(gene$ADJ_1))
               #if (is.na(gene$ADJ_2))
                  #text(gene$LOG2_FC, gene$log10P, paste(genes[g,]$GENE, genes[g,]$ADJ_3), col="black", adj=gene$ADJ_1, cex=1.7)
               #else
                  #text(gene$LOG2_FC, gene$log10P, paste(genes[g,]$GENE, genes[g,]$ADJ_3), col="black", adj=c(gene$ADJ_1, gene$ADJ_2), cex=1.7)
            #else
               #if (gene$LOG2_FC > 0)
                  #text(gene$LOG2_FC, gene$log10P, paste(genes[g,]$GENE, genes[g,]$ADJ_3), col="black", adj=c(0, -0.5), cex=1.7)
               #else
                  #text(gene$LOG2_FC, gene$log10P, paste(genes[g,]$GENE, genes[g,]$ADJ_3), col="black", adj=c(1, -0.3), cex=1.7)
         } else
            print(genes[g,])
      }
   }
 
   #legends[1] <- paste0(nrow(de.up), " ", legends[1])
   #legends[2] <- paste0(nrow(de.down), " ", legends[2])
   #axis(side=1, at=seq(-8, 8, by=1), labels=c(-3, -2, -1, 0, 1, 2, 3), cex.axis=1.2)
   #legend(legend, legend=legends, col=cols2, pch=19, pt.cex=2, cex=1.6)
   dev.off()
}

# -----------------------------------------------------------------------------
# Method: Cannoli plot
# Last Modified: 01/07/22
# -----------------------------------------------------------------------------
getSRC0 <- function(wd.de.data, BASE, tpm.gene.log2, samples.tpm, COR, n) {
   colnames <- c("RHO", "P", "Q")
   src <- toTable(0, length(colnames), nrow(tpm.gene.log2), colnames)
   rownames(src) <- rownames(tpm.gene.log2)
 
   ## SRC
   src$RHO <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples.tpm[, COR], method="spearman", exact=F)[[4]])
   src$P   <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples.tpm[, COR], method="spearman", exact=F)[[3]])
   src <- src[!is.na(src$P),]
 
   ## FDR
   library(qvalue)
   src$Q <- qvalue(src$P)$qvalue
   src <- src[order(src$P),]
 
   ## Ensembl gene annotations
   annot <- ensGene[,c("ensembl_gene_id", "external_gene_name")]
   src.tpm.gene <- cbind(annot[rownames(src),], src)   ## BE EXTRA CAREFUL!!
 
   #writeTable(src.tpm.gene, file.path(wd.de.data, paste0("SRC_", BASE, "_tpm-gene_TPM-vs-", COR, "_n", n, ".txt")), colnames=T, rownames=F, sep="\t")
   save(src.tpm.gene, file=file.path(wd.de.data, paste0("SRC_", BASE, "_tpm-gene_TPM-vs-", COR, "_n", n, ".RData")))
}

getSRC <- function(wd.de.data, BASE, tpm.gene.log2, samples.tpm, COR, n, Q, q) {
   colnames <- c("RHO", "P", "Q")
   src <- toTable(0, length(colnames), nrow(tpm.gene.log2), colnames)
   rownames(src) <- rownames(tpm.gene.log2)
 
   ## SRC
   src$RHO <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples.tpm[, COR], method="spearman", exact=F)[[4]])
   src$P   <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples.tpm[, COR], method="spearman", exact=F)[[3]])
   src <- src[!is.na(src$P),]
   
   ## FDR
   library(qvalue)
   src$Q <- qvalue(src$P)$qvalue
   src <- src[order(src$P),]
   
   ## Ensembl gene annotations
   annot <- ensGene[,c("ensembl_gene_id", "external_gene_name")]
   src.tpm.gene <- cbind(annot[rownames(src),], src)   ## BE EXTRA CAREFUL!!
   
   #writeTable(src.tpm.gene, file.path(wd.de.data, paste0("SRC_", BASE, "_tpm-gene_TPM-vs-", COR, "_n", n, "_", Q, q, ".txt")), colnames=T, rownames=F, sep="\t")
   save(src.tpm.gene, file=file.path(wd.de.data, paste0("SRC_", BASE, "_tpm-gene_TPM-vs-", COR, "_n", n, "_", Q, q, ".RData")))
}

getCannoli <- function(de.data, BASE, n, gene.list=NULL, TEST="COR", TEST2="Purity", M2="_M2") {
   load(file=file.path(de.data, paste0("SRC_", BASE, "_tpm-gene_TPM-vs-", TEST, "_n", n, M2, ".RData")))
   #load(file=file.path(de.data, paste0("SRC_", BASE, "_tpm-gene_", TEST, "-vs-TPM_q_n", n, M2, ".RData")))
   de1 <- src.tpm.gene
   load(file=file.path(de.data, paste0("SRC_", BASE, "_tpm-gene_TPM-vs-", TEST2, "_n", n, M2, ".RData")))
   #load(file=file.path(de.data, paste0("SRC_", BASE, "_tpm-gene_", TEST2, "-vs-TPM_q_n", n, M2, ".RData")))
   de2 <- src.tpm.gene
 
   overlaps <- intersect(rownames(de1), rownames(de2))
   de <- cbind(de1[overlaps, c("P", "Q", "RHO")], de2[overlaps, c("P", "Q", "RHO")])
   rownames(de) <- overlaps
   de <- cbind(ensGene[overlaps, c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")], de)
   colnames(de) <- c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype", "P1", "Q1", "Effect1", "P2", "Q2", "Effect2")
 
   if (!is.null(gene.list))
      de <- de[intersect(rownames(de), gene.list),]
   
   #de.positive <- subset(subset(de, Effect1 > 0), Effect2 > 0)
   #de.positive.sig <- subset(subset(de.positive, P1 <= 0.001), P2 <= 0.001)
   #de.negative <- subset(subset(de, Effect1 < 0), Effect2 < 0)
   #de.negative.sig <- subset(subset(de.negative, P1 <= 0.001), P2 <= 0.001)
   
   return(de)
}

getCannoliGenes <- function(de, pvalue, genes0) {
   de.pos.pos <- subset(subset(de, Effect2 > 0), Effect1 > 0)
   de.pos.pos.sig <- subset(subset(de.pos.pos, P1 <= pvalue), P2 <= pvalue)
   de.pos.neg <- subset(subset(de, Effect2 > 0), Effect1 < 0)
   de.pos.neg.sig <- subset(subset(de.pos.neg, P1 <= pvalue), P2 <= pvalue)
 
   if (nrow(de.pos.pos.sig) > 0 && nrow(de.pos.pos.sig) < 5)
      genes0 <- c(genes0, as.vector(de.pos.pos.sig$external_gene_name))
   if (nrow(de.pos.neg.sig) > 0 && nrow(de.pos.neg.sig) < 5)
      genes0 <- c(genes0, as.vector(de.pos.neg.sig$external_gene_name))
 
   colnames <- c("GENE", "ADJ_1", "ADJ_2")
   genes <- toTable(NA, length(colnames), length(genes0), colnames)
   genes$GENE <- genes0
 
   return(genes)
}

plotCannoli<- function(de, pvalue, genes, file.de, file.main, xlab.text, ylab.text, legend, legends, cols, cols2, fold=1, ymax=0, col="black", pos="bottomright") {
   plot.de <- as.vector(unlist(strsplit0(file.de, ".pdf"))[1])
 
   pdf(file.de, height=6, width=6)
   par(mar=c(5.1, 4.7, 4.1, 1.4), xpd=F)
   plot(de$Effect1, de$Effect2, pch=16, xlab=xlab.text, ylab=ylab.text, xaxt="n", yaxt="n", col="white", main=file.main[1], cex=1.7, cex.axis=1.6, cex.lab=1.7, cex.main=1.8)
   #abline(h=c(-log10(pvalue)), lty=5, lwd=2)
   
   de.neg <- subset(de, Effect2 < 0)
   points(de.neg$Effect1, de.neg$Effect2, pch=16, col="lightgray", cex=1.7)

   de.pos.pos <- subset(subset(de, Effect2 > 0), Effect1 > 0)
   de.pos.pos.sig <- subset(subset(de.pos.pos, P1 <= pvalue), P2 <= pvalue)
   points(de.pos.pos$Effect1,     de.pos.pos$Effect2,     pch=16, col=cols[1], cex=1.7)
   points(de.pos.pos.sig$Effect1, de.pos.pos.sig$Effect2, pch=16, col=cols2[1], cex=1.7)
   if (nrow(de.pos.pos.sig) > 0)
      writeTable(de.pos.pos.sig, file.path(paste0(plot.de, "_pos.pos.txt")), colnames=T, rownames=F, sep="\t")
   
   de.pos.neg <- subset(subset(de, Effect2 > 0), Effect1 < 0)
   de.pos.neg.sig <- subset(subset(de.pos.neg, P1 <= pvalue), P2 <= pvalue)
   points(de.pos.neg$Effect1,     de.pos.neg$Effect2,     pch=16, col=cols[2], cex=1.7)
   points(de.pos.neg.sig$Effect1, de.pos.neg.sig$Effect2, pch=16, col=cols2[2], cex=1.7)
   if (nrow(de.pos.neg.sig) > 0)
      writeTable(de.pos.neg.sig, file.path(paste0(plot.de, "_pos.neg.txt")), colnames=T, rownames=F, sep="\t")
   
   #text(max(de.pos.pos$Effect1)/2, max(de.pos.pos$Effect2)/2, paste0(round0(nrow(de.pos.pos)/nrow(de)*100, 1), " %"), col="white", font=2, cex=1.8)
   #text(min(de.pos.neg$Effect1)/2, max(de.pos.pos$Effect2)/2, paste0(round0(nrow(de.pos.neg)/nrow(de)*100, 1), " %"), col="white", font=2, cex=1.8)
   
   abline(v=0, lty=5, col="black", lwd=2)
   abline(h=0, lty=5, col="black", lwd=2)
 
   lm.fit <- lm(de$Effect1 ~ de$Effect2)
   abline(lm.fit, col=col, lwd=7)
   
   cor <- cor.test(de$Effect1, de$Effect2, method="spearman", exact=F)
   #legend(pos, c(paste0("N = ", separator(nrow(de))), paste0("rho = ", round0(cor[[4]], digits=2)), paste0("P = 1.00E-00")), text.col=c(col, col, adjustcolor.white), text.font=2, bty="n", cex=1.8)
   legend(pos, c("", paste0("rho = ", round0(cor[[4]], digits=2)), paste0("P = 1.00E-00")), text.col=c(col, col, adjustcolor.white), text.font=2, bty="n", cex=1.8)
   legend(pos, expression(bolditalic('P')~"                   "), text.col=col, text.font=2, bty="n", cex=1.8)
   legend(pos, paste0("   = ", scientific(cor[[3]])), text.col=col, text.font=2, bty="n", cex=1.8)
   
   par(xpd=T)
   if (nrow(genes) != 0) {
      for (g in 1:nrow(genes)) {
         gene <- subset(de, external_gene_name == genes[g,]$GENE)
         gene <- cbind(gene, genes[g,])
     
         if (nrow(gene) > 0) {
            points(gene$Effect1, gene$Effect2, pch=1, col="black", cex=1.7)
      
            if (!is.na(gene$ADJ_1))
               if (is.na(gene$ADJ_2))
                  text(gene$Effect1, gene$Effect2, paste(genes[g,]$GENE, genes[g,]$ADJ_3), col="black", adj=gene$ADJ_1, cex=1.7)
               else
                  text(gene$Effect1, gene$Effect2, paste(genes[g,]$GENE, genes[g,]$ADJ_3), col="black", adj=c(gene$ADJ_1, gene$ADJ_2), cex=1.7)
            else
               if (gene$Effect1 > 0)
                  text(gene$Effect1, gene$Effect2, paste(genes[g,]$GENE, genes[g,]$ADJ_3), col="black", adj=c(0, -0.5), cex=1.7)
               else
                  text(gene$Effect1, gene$Effect2, paste(genes[g,]$GENE, genes[g,]$ADJ_3), col="black", adj=c(1, -0.3), cex=1.7)
         } else
            print(genes[g,])
      }
   }
   
   legends[1] <- paste0(nrow(de.pos.pos.sig), " ", legends[1])
   legends[2] <- paste0(nrow(de.pos.neg.sig), " ", legends[2])
   axis(side=1, at=seq(-0.5, 0.5, by=0.5), labels=c(-0.5, 0, 0.5), cex.axis=1.6)
   axis(side=2, at=seq(-0.5, 0.5, by=0.5), labels=c(-0.5, 0, 0.5), cex.axis=1.6)
   legend(legend, legend=legends, col=cols2, pch=19, pt.cex=2, cex=1.6)
   dev.off()
}

plotCannoliGenes<- function(de, genes, file.de, file.main, xlab.text, ylab.text, col="black", pos="bottomright") {
   plot.de <- as.vector(unlist(strsplit0(file.de, ".pdf"))[1])
 
   pdf(file.de, height=6, width=6)
   par(mar=c(5.1, 4.7, 4.1, 1.4), xpd=F)
   plot(de$Effect1, de$Effect2, pch=16, xlab=xlab.text, ylab=ylab.text, xaxt="n", yaxt="n", col="lightgray", main=file.main[1], cex=1.7, cex.axis=1.6, cex.lab=1.7, cex.main=1.8)
   #abline(h=c(-log10(pvalue)), lty=5, lwd=2)

   overlaps <- intersect(rownames(de), genes)
   de.genes <- de[overlaps, ]
   
   de.pos.pos <- subset(subset(de.genes, Effect2 > 0), Effect1 > 0)
   de.pos.neg <- subset(subset(de.genes, Effect2 > 0), Effect1 < 0)
   text( 0.5, 0.5, paste0(round0(nrow(de.pos.pos)/nrow(de.genes)*100, 1), " %"), col=col, font=2, cex=1.8)
   text(-0.5, 0.5, paste0(round0(nrow(de.pos.neg)/nrow(de.genes)*100, 1), " %"), col=col, font=2, cex=1.8)

   de.neg.pos <- subset(subset(de.genes, Effect2 < 0), Effect1 > 0)
   de.neg.neg <- subset(subset(de.genes, Effect2 < 0), Effect1 < 0)
   text( 0.5, -0.5, paste0(round0(nrow(de.neg.pos)/nrow(de.genes)*100, 1), " %"), col=col, font=2, cex=1.8)
   text(-0.5, -0.5, paste0(round0(nrow(de.neg.neg)/nrow(de.genes)*100, 1), " %"), col=col, font=2, cex=1.8)
     
   abline(v=0, lty=5, col="black", lwd=2)
   abline(h=0, lty=5, col="black", lwd=2)
 
   lm.fit <- lm(de$Effect1 ~ de$Effect2)
   abline(lm.fit, col="black", lwd=7)
 
   cor <- cor.test(de$Effect1, de$Effect2, method="spearman", exact=F)
   #legend(pos, c(paste0("N = ", separator(nrow(de))), paste0("rho = ", round0(cor[[4]], digits=2)), paste0("P = 1.00E-00")), text.col=c(col, col, adjustcolor.white), text.font=2, bty="n", cex=1.8)
   legend(pos, c("", paste0("rho = ", round0(cor[[4]], digits=2)), paste0("P = 1.00E-00")), text.col=c("black", "black", adjustcolor.white), text.font=2, bty="n", cex=1.8)
   legend(pos, expression(bolditalic('P')~"                   "), text.col="black", text.font=2, bty="n", cex=1.8)
   legend(pos, paste0("   = ", scientific(cor[[3]])), text.col="black", text.font=2, bty="n", cex=1.8)
 
   if (length(genes) != 0) {
      for (g in 1:length(genes)) {
         gene <- subset(de, ensembl_gene_id == genes[g])
         
         if (nrow(gene) > 0) {
            points(gene$Effect1, gene$Effect2, pch=1, col=col, cex=1.7)
         }
      }
   }

   axis(side=1, at=seq(-0.5, 0.5, by=0.5), labels=c(-0.5, 0, 0.5), cex.axis=1.6)
   axis(side=2, at=seq(-0.5, 0.5, by=0.5), labels=c(-0.5, 0, 0.5), cex.axis=1.6)
   dev.off()
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

## RNK format
## http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#RNK:_Ranked_list_file_format_.28.2A.rnk.29
writeRNKformat <- function(de.tpm.gene, wd.de.data, file.name) {
   #de.sorted <- de.tpm.gene[order(de.tpm.gene$LOG2_FC, decreasing=T),]
   de.tpm.gene$WEIGHT <- de.tpm.gene$LOG2_FC * (-log10(de.tpm.gene$P))
   #de.tpm.gene$WEIGHT <- de.tpm.gene$RHO
   de.sorted <- de.tpm.gene[order(de.tpm.gene$WEIGHT, decreasing=T),]
   
   file <- file.path(wd.de.data, paste0(file.name, ".rnk"))
   #writeTable(de.sorted[,c("ensembl_gene_id", "WEIGHT")], file, colnames=F, rownames=F, sep="\t")
   writeTable(de.sorted[,c("external_gene_name", "WEIGHT")], file, colnames=F, rownames=F, sep="\t")
}

writeRNKformatCNA <- function(de, wd.de.data, file.name, pos) {
   de$WEIGHT <- de$Effect2 * de$Effect1
  
   de.sorted <- de[order(de$WEIGHT, decreasing=T),]
 
   file <- file.path(wd.de.data, paste0(file.name, ".rnk"))
   writeTable(de.sorted[,c("external_gene_name", "WEIGHT")], file, colnames=F, rownames=F, sep="\t")
}

writeRNKformatNOCNA <- function(de, wd.de.data, file.name, pos) {
   de$WEIGHT <- de$Effect1
 
   de.sorted <- de[order(de$WEIGHT, decreasing=T),]
 
   file <- file.path(wd.de.data, paste0(file.name, ".rnk"))
   writeTable(de.sorted[,c("external_gene_name", "WEIGHT")], file, colnames=F, rownames=F, sep="\t")
}

## http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GRP
writeGRPformat <- function(list, wd.de.data, file.name) {
   writeTable(list, file.path(wd.de.data, paste0(file.name, ".grp")), colnames=F, rownames=F, sep="\t")
}

## Gene sets in GMT format
## https://www.gsea-msigdb.org/gsea/downloads.jsp (MSigDB 7.5.1)
## https://www.gsea-msigdb.org/gsea/msigdb/mouse_geneset_resources.jsp

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
