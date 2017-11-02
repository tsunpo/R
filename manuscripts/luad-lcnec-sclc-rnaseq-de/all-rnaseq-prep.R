# =============================================================================
# Name: lcnec_rnaseq_de_kallisto_ensembl.R
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 26/03/17
# =============================================================================

# -----------------------------------------------------------------------------
# Sleuth installation (v0.29.0)
# -----------------------------------------------------------------------------
source("http://bioconductor.org/biocLite.R")
biocLite("rhdf5", "tximport")
install.packages("devtools", "stringr", "cluster", "survival", "tidyr", "dplyr")   ## ADD

library("devtools")
devtools::install_github("pachterlab/sleuth")

# -----------------------------------------------------------------------------
# Transcript TPM estimates/aboundants from kallisto
# -----------------------------------------------------------------------------
wd.source <- "/Users/tpyang/Work/local/R"
sources <- c("DailyMeal.R", "DifferentialAnalysis.R")
invisible(sapply(sources, function(s) source(file.path(wd.source, "handbookof", s))))
load(file.path(wd.source, "guidetothe", "hg19.transcript.RData"))

###
##
library("sleuth")
wd <- "/Users/tpyang/Work/uni-koeln/tyang2/LCNEC/"
wd.de <- paste0(wd, "analysis/expression/kallisto/lcnec-rnaseq-de-kallisto-ensembl/")
wd.kallisto <- paste0(wd, "ngs/RNA/kallisto_hg19.ensembl_quant-b100--bias/")
setwd(wd)

#load("/Users/tpyang/Work/uni-koeln/tyang2/LCNEC/analysis/expression/kallisto/kallisto_refseq_quant-b100--bias_n69.RData")
load("/Users/tpyang/Work/uni-koeln/tyang2/LCNEC/analysis/expression/FPKM/lcnec_expr.RData")
samples.expr <- readTable(paste0(wd, "ngs/RNA/lcnec_n69.list"), header=F, rownames=F, sep="\t")
pheno.expr <- pheno[samples.expr,]

tsvs <- sapply(samples.expr, function(s) file.path(wd.kallisto, s))
s2c <- data.frame(path=tsvs, sample=samples.expr, stringsAsFactors=F)
t2g <- tx2Ens(ensGene)

so <- sleuth_prep(s2c, ~1, target_mapping=t2g)
# reading in kallisto results
# .....................................................................
# normalizing est_counts
# 88335 targets passed the filter
# normalizing tpm
# merging in metadata
# normalizing bootstrap samples
# summarizing bootstraps

sleuth_to_matrix0 <- function(obj, which_df, which_units) {
   data <- as.data.frame(obj[[which_df]])
 
   res <- list()
 
   s_data <- data %>%
   select_("target_id", "sample", which_units) %>%
   tidyr::spread_("sample", which_units)
   rownames(s_data) <- s_data$target_id
   s_data$target_id <- NULL
   s_data <- as.matrix(s_data)
   s_data <- s_data[, sample(1:ncol(s_data))]
   res[["data"]] <- s_data
 
   condition_order <- match(colnames(s_data), as.character(obj$sample_to_condition$sample))
   res[["condition"]] <- obj$sample_to_condition$condition[condition_order]
 
   res
}

orderBySamples <- function(tpm) {
   return(tpm[,order(colnames(tpm), decreasing=F)])
}

##
tpm.norm.filt <- sleuth_to_matrix0(so, "obs_norm_filt", "tpm")$data
tpm.norm.filt <- orderBySamples(tpm.norm.filt)
tpm.norm <- sleuth_to_matrix(so, "obs_norm", "tpm")$data
tpm.norm <- orderBySamples(tpm.norm)
save(tpm.norm.filt, tpm.norm, file=paste0(wd.de, "data/lcnec_kallisto_tpm.RData"))

tpm2.norm.filt <- kallisto_table(so, use_filtered=T, normalized=T)$tpm
tpm2.raw.filt <- kallisto_table(so, use_filtered=T, normalized=F)$tpm
tpm2.norm <- kallisto_table(so, use_filtered=F, normalized=T)$tpm
tpm2.raw <- kallisto_table(so, use_filtered=F, normalized=F)$tpm
save(tpm2.norm, tpm2.norm.filt, tpm2.raw, tpm2.raw.filt, file=paste0(wd.de, "data/lcnec_kallisto_tpm2.RData"))

# -----------------------------------------------------------------------------
# Gene-level TPM estimates/aboundants from kallisto
# -----------------------------------------------------------------------------
sumGeneTPM <- function(t2g, tpm) {
   t2g.tpm <- t2g[rownames(tpm),]      ## Keep only filtered transcripts
   genes <- unique(t2g.tpm$ens_gene)

   tpm.gene <- toTable(0, ncol(tpm), length(genes), colnames(tpm))
   rownames(tpm.gene) <- genes
   for (g in 1:length(genes)) {
      transcripts <- rownames(subset(t2g.tpm, ens_gene == genes[g]))   ## Find corresponding transcripts
      tpm.transcripts <- tpm[transcripts,]                                    ## and sum up TPM estimates of each gene

      if (length(transcripts) == 1)
         tpm.gene[g,] <- tpm.transcripts
      else
         tpm.gene[g,] <- mapply(s = 1:ncol(tpm.gene), function(s) sum(tpm.transcripts[,s]))
   }
 
   return(tpm.gene)
}

##
overlaps <- intersect(rownames(t2g), rownames(tpm.norm.filt))
tpm <- tpm.norm.filt[overlaps,]

overlaps <- intersect(rownames(subset(ensGene, transcript_biotype == "protein_coding")), rownames(tpm.norm.filt))
tpm.pct <- tpm.norm.filt[overlaps,]
#tpm.exp <- tpm[getExpressedGTEx(tpm),]
# > nrow(tpm.exp)
# [1] 34011

tpm.gene <- sumGeneTPM(t2g, tpm)
overlaps <- intersect(rownames(subset(ensGene.gene, gene_biotype == "protein_coding")), rownames(tpm.gene))
tpm.gene.pcg <- tpm.gene[overlaps,]
# > nrow(tpm.norm.filt)
# [1] 88335
# > nrow(tpm)
# [1] 83255
# > nrow(tpm.pct)
# [1] 46730

# > nrow(tpm.gene)
# [1] 18917
# > nrow(tpm.gene.pcg)
# [1] 16921

## Filter for expressed genes in ALL (n=69) samples
zeros <- reportZeroTPM(tpm.gene)
table <- as.data.frame(table(zeros))

tpm.gene.exp <- tpm.gene[getExpressedGTEx(tpm.gene),]
overlaps <- intersect(rownames(subset(ensGene.gene, gene_biotype == "protein_coding")), rownames(tpm.gene.exp))
tpm.gene.exp.pcg <- tpm.gene.exp[overlaps,]
# > nrow(tpm.gene.exp)
# [1] 15670
# > nrow(tpm.gene.exp.pcg)
# [1] 14749

## Filter for expressed genes in FINAL ANALYSIS (n=54) samples 
# > nrow(de.tpm.gene.exp2)
# [1] 16271
# > nrow(de.tpm.gene.exp2.pcg)
# [1] 15132

save(tpm.norm.filt, tpm, tpm.gene, file=paste0(wd.de, "data/lcnec_kallisto_tpm.norm.filt.RData"))











# -----------------------------------------------------------------------------
# Determine detected/expressed genes (TPM; n=54)
# Last Modified: 28/03/17
# -----------------------------------------------------------------------------
## For QC
tpm.norm.log2 <- log2(tpm.norm + 0.01)
tpm.norm.filt.log2 <- log2(tpm.norm.filt + 0.01)   ## Filtered transcript-level TPMs

## For downstream D.E.
tpm.log2 <- log2(tpm + 0.01)                     ## Transcript-level TPMs (matching with ensGene; with 0 values)
tpm.pct.log2 <- log2(tpm.pct + 0.01)             ## Transcript-level protein-coding-transcript TPMs 
tpm.gene.log2 <- log2(tpm.gene + 0.01)           ## Gene-level TPMs (with 0 values)
tpm.gene.pcg.log2 <- log2(tpm.gene.pcg + 0.01)   ## Gene-level, protein-coding-gene TPMs (with 0 values)

#tpm.gene.exp.log2 <- log2(tpm.gene.exp + 0.01)           ## Gene-level, expressed TPMs (with 0 values in up to 2 samples)
#tpm.gene.exp.pcg.log2 <- log2(tpm.gene.exp.pcg + 0.01)   ## Gene-level, expressed, protein-coding-gene TPMs (with 0 values in up to 2 samples)

###
## Filtered/un-filtered transcript-level TPMs
ymax <- 20000   #length(which(medians.magic < 6))

pdf(paste0(wd.de, "data/hist_kallisto_tpm.norm.log2.pdf"))
hist(median0(tpm.norm.log2), main=c("Median Transcript Expression", "Ensembl cDNAs (n=180,253)"), xlab="log2(TPM + 0.01)", ylim=c(0, ymax))   #, xlim=c(7, 24), ylim=c(0, ymax))
dev.off()

pdf(paste0(wd.de, "data/hist_kallisto_tpm.norm.filt.log2.pdf"))
hist(median0(tpm.norm.filt.log2), main=c("Median Transcript Expression", "Ensembl cDNAs (Filtered; n=88,335)"), xlab="log2(TPM + 0.01)", ylim=c(0, ymax))   #, xlim=c(7, 24), ylim=c(0, ymax))
dev.off()

## Transcript-level protein-coding-transcript TPMs
ymax <- 17000   #length(which(medians.magic < 6))

pdf(paste0(wd.de, "data/hist_kallisto_tpm.log2.pdf"))
hist(median0(tpm.log2), main=c("Median Transcript Expression", "Ensembl cDNAs (Filtered; n=83,255)"), xlab="log2(TPM + 0.01)", ylim=c(0, ymax))   #, xlim=c(7, 24), ylim=c(0, ymax))
dev.off()

pdf(paste0(wd.de, "data/hist_kallisto_tpm.coding.log2.pdf"))
hist(median0(tpm.pct.log2), main=c("Median Transcript Expression", "Ensembl cDNAs (Protein Coding; n=46,730)"), xlab="log2(TPM + 0.01)", ylim=c(0, ymax))   #, xlim=c(7, 24), ylim=c(0, ymax))
dev.off()

## Gene-level, protein-coding-gene TPMs (with 0 values)
ymax <- 3200   #length(which(medians.magic < 6))

pdf(paste0(wd.de, "data/hist_kallisto_tpm.gene.log2.pdf"))
hist(median0(tpm.gene.log2), main=c("Median Gene Expression", "Ensembl cDNAs (Filtered; n=18,917)"), xlab="log2(TPM + 0.01)", ylim=c(0, ymax))   #, xlim=c(7, 24), ylim=c(0, ymax))
dev.off()

pdf(paste0(wd.de, "data/hist_kallisto_tpm.gene.pcg.log2.pdf"))
hist(median0(tpm.gene.pcg.log2), main=c("Median Gene Expression", "Ensembl cDNAs (Protein Coding; n=16,921)"), xlab="log2(TPM + 0.01)", ylim=c(0, ymax))   #, xlim=c(7, 24), ylim=c(0, ymax))
dev.off()

## Gene-level, expressed, protein-coding-gene TPMs (without 0 values in ALL)
ymax <- 3200   #length(which(medians.magic < 6))

pdf(paste0(wd.de, "data/hist_kallisto_tpm.gene.exp.log2.pdf"))
hist(median0(tpm.gene.exp.log2), main=c("Median Gene Expression", "Ensembl cDNAs (Expressed; n=15,670)"), xlab="log2(TPM + 0.01)", ylim=c(0, ymax))   #, xlim=c(7, 24), ylim=c(0, ymax))
dev.off()

pdf(paste0(wd.de, "data/hist_kallisto_tpm.gene.exp.pcg.log2.pdf"))
hist(median0(tpm.gene.exp.pcg.log2), main=c("Median Gene Expression", "Ensembl cDNAs (Expressed Coding; n=14,749)"), xlab="log2(TPM + 0.01)", ylim=c(0, ymax))   #, xlim=c(7, 24), ylim=c(0, ymax))
dev.off()

## Gene-level, expressed, protein-coding-gene TPMs (without 0 values in FINAL ANALYSIS samples)
ymax <- 3200   #length(which(medians.magic < 6))

pdf(paste0(wd.de, "data/hist_kallisto_tpm.gene.exp2.log2.pdf"))
hist(median0(tpm.gene.exp2.log2), main=c("Median Gene Expression", "Ensembl cDNAs (Expressed; n=16,271)"), xlab="log2(TPM + 0.01)", ylim=c(0, ymax))   #, xlim=c(7, 24), ylim=c(0, ymax))
dev.off()

pdf(paste0(wd.de, "data/hist_kallisto_tpm.gene.exp2.pcg.log2.pdf"))
hist(median0(tpm.gene.exp2.pcg.log2), main=c("Median Gene Expression", "Ensembl cDNAs (Expressed Coding; n=15,132)"), xlab="log2(TPM + 0.01)", ylim=c(0, ymax))   #, xlim=c(7, 24), ylim=c(0, ymax))
dev.off()

# -----------------------------------------------------------------------------
# Differential expression (RB1NEW vs WT; n=69-15NA)
# Last Modified: 28/03/17
# -----------------------------------------------------------------------------
pipeDE <- function(expr, pheno.expr, dao, file.de, file.main, annot) {
   if (dao$exp.only) {
      expr <- getFinalExpression(expr, pheno.expr)
      expr <- expr[getExpressedGTEx(expr),]
   }
   if (!dao$is.log2)
      expr <- log2(expr + 0.01) 
 
   ## Differential expression
   de <- differentialAnalysis(expr, pheno.expr, dao$predictor, dao$predictor.wt, dao$test, dao$test.fdr)
   de <- cbind(annot[rownames(de),], de)
   
   writeTable(de, paste0(file.de, ".txt"), colnames=T, rownames=T, sep="\t")
   writeTable(onlyVariable(de, dao$effect), paste0(file.de, "_Effect>", dao$effect, ".txt"), colnames=T, rownames=T, sep="\t")   ## Only varible genes
   
   ## Volcano plot
   plotVolcano(de, dao$fdr, dao$effect, file.de, file.main, legend.x=NA)
   
   return(de)
}

filenameDE <- function(wd.de, title, dao) {
   return(paste0(wd.de, title, "_", dao$predictor, "_", dao$test, "_", dao$test.fdr))
}

## Data access object for test parameters
## Test: Wilcoxon/Wilcox/U/Students/ttest
## FDR: Q/BH
## DE:  RB1NEW vs WT(0) as factor
dao <- data.frame(is.log2=T, exp.only=F, test="Wilcox", test.fdr="Q", fdr=0.05, effect=1, predictor="RB1NEW", predictor.wt=0, stringsAsFactors=F)

###
## Transcript-level TPMs (matching with ensGene; with 0 values)
annot <- cbind(t2g, ensGene[rownames(t2g), "transcript_biotype"])

file.de <- filenameDE(wd.de, "LCNEC_DE_kallisto_ens_tpm", dao)
file.main <- c("Differential Transcript-Level Expression in LCNEC", "Mutant RB1 (n=20) vs WT (n=34)")
#expr <- tpm.log2
#de.tpm <- pipeDE(expr, pheno.expr, dao, file.de, file.main, annot)
plotVolcano(de.tpm, dao$fdr, dao$effect, file.de, file.main, legend.x=NA)

###
## Gene-level TPMs (with 0 values)
annot.gene <- ensGene.gene[,c("ensembl_gene_id", "external_gene_name", "gene_biotype")]

file.de <- filenameDE(wd.de, "LCNEC_DE_kallisto_ens_tpm-gene", dao)
file.main <- c("Differential Gene-Level Expression in LCNEC", "Mutant RB1 (n=20) vs WT (n=34)")
#expr <- tpm.gene.log2
#de.tpm.gene <- pipeDE(expr, pheno.expr, dao, file.de, file.main, annot.gene)
plotVolcano(de.tpm.gene, dao$fdr, dao$effect, file.de, file.main, legend.x=NA)

## Gene-level, protein-coding-gene TPMs (with 0 values)
file.de <- filenameDE(wd.de, "LCNEC_DE_kallisto_ens_tpm-gene-pcg", dao)
file.main <- c("Differential Gene Expression in LCNEC", "Mutant RB1 (n=20) vs WT (n=34)")
#expr <- tpm.gene.pcg.log2
#de.tpm.gene.pcg <- pipeDE(expr, pheno.expr, dao, file.de, file.main, annot.gene)
plotVolcano(de.tpm.gene.pcg, dao$fdr, dao$effect, file.de, file.main, legend.x=NA)

###
## Gene-level, expressed TPMs (without 0 values in ALL (n=69) samples)
file.de <- filenameDE(wd.de, "LCNEC_DE_kallisto_ens_tpm-gene-exp", dao)
file.main <- c("Differential Gene-Level Expression in LCNEC", "Mutant RB1 (n=20) vs WT (n=34)")
#expr <- tpm.gene.exp.log2
#de.tpm.gene.exp <- pipeDE(expr, pheno.expr, dao, file.de, file.main, annot.gene)
plotVolcano(de.tpm.gene.exp, dao$fdr, dao$effect, file.de, file.main, legend.x=NA)

## Gene-level, expressed, protein-coding-gene TPMs (without 0 values)
file.de <- filenameDE(wd.de, "LCNEC_DE_kallisto_ens_tpm-gene-exp-pcg", dao)
file.main <- c("Differential Gene Expression in LCNEC", "Mutant RB1 (n=20) vs WT (n=34)")
#expr <- tpm.gene.exp.pcg.log2
#de.tpm.gene.exp.pcg <- pipeDE(expr, pheno.expr, dao, file.de, file.main, annot.gene)
plotVolcano(de.tpm.gene.exp.pcg, dao$fdr, dao$effect, file.de, file.main, legend.x=NA)

###
## Gene-level, expressed TPMs (without 0 values in FINAL ANALYSIS (n=54) samples)
dao$is.log2=F
dao$exp.only=T

file.de <- filenameDE(wd.de, "LCNEC_DE_kallisto_ens_tpm-gene-exp2", dao)
file.main <- c("Differential Gene-Level Expression in LCNEC", "Mutant RB1 (n=20) vs WT (n=34)")
expr <- tpm.gene
de.tpm.gene.exp2 <- pipeDE(expr, pheno.expr, dao, file.de, file.main, annot.gene)
#plotVolcano(de.tpm.gene.exp, dao$fdr, dao$effect, file.de, file.main, legend.x=NA)

## Gene-level, expressed, protein-coding-gene TPMs (without 0 values)
file.de <- filenameDE(wd.de, "LCNEC_DE_kallisto_ens_tpm-gene-exp2-pcg", dao)
file.main <- c("Differential Gene Expression in LCNEC", "Mutant RB1 (n=20) vs WT (n=34)")
expr <- tpm.gene.pcg
de.tpm.gene.exp2.pcg <- pipeDE(expr, pheno.expr, dao, file.de, file.main, annot.gene)
#plotVolcano(de.tpm.gene.exp.pcg, dao$fdr, dao$effect, file.de, file.main, legend.x=NA)

save(de.tpm, de.tpm.gene, de.tpm.gene.pcg, file=paste0(wd.de, "lcnec_kallisto_de.tpm.gene_RB1NEW.RData"))
save(de.tpm.gene.exp, de.tpm.gene.exp.pcg, file=paste0(wd.de, "lcnec_kallisto_de.tpm.gene.exp_RB1NEW.RData"))
save(de.tpm.gene.exp2, de.tpm.gene.exp2.pcg, file=paste0(wd.de, "lcnec_kallisto_de.tpm.gene.exp2_RB1NEW.RData"))

# -----------------------------------------------------------------------------
# Differential expression on residuals of expression (RB1NEW vs WT; n=69-15NA)
# Last Modified: 29/03/17
# -----------------------------------------------------------------------------
copyEffectSize <- function(de.res, de.tpm) {
   de.res.temp <- de.tpm[rownames(de.res),]       ## Use the original annotation and TPM effect size
   de.res.temp[,c("P", "FDR")] <- de.res[,c("P", "FDR")]   ## Only replace the D.E. P and FDR
   
   return(de.res.temp)
}

###
## Gene-level TPMs (without 0 values in FIANL ANALYSIS)
tpm.gene.exp2.log2 <- tpm.gene.log2[rownames(de.tpm.gene.exp2),]

covariate <- "Classification2"   ## Two-Subtype Residual
file.de <- filenameDE(wd.de, paste0("LCNEC_DE_kallisto_ens_tpm-gene-exp2-res2-", covariate), dao)
file.main <- c("Mutant RB1 Differential Expression in LCNEC", "Two-Subtype Residual")
expr.res <- residualsOf(tpm.gene.exp2.log2, pheno.expr, covariate)

de.tpm.gene.exp2.res2 <- differentialAnalysis(expr.res, pheno.expr, dao$predictor, dao$predictor.wt, dao$test, dao$test.fdr)
de.tpm.gene.exp2.res2 <- copyEffectSize(de.tpm.gene.exp2.res2, de.tpm.gene)
writeTable(de.tpm.gene.exp2.res2, paste0(file.de, ".txt"), colnames=T, rownames=T, sep="\t")
writeTable(onlyVariable(de.tpm.gene.exp2.res2, dao$effect), paste0(file.de, "_Effect>", dao$effect, ".txt"), colnames=T, rownames=T, sep="\t")   ## Only varible genes
plotVolcano(de.tpm.gene.exp2.res2, dao$fdr, dao$effect, file.de, file.main, legend.x=NA)

##
covariate <- "Classification"   ## Four-Subtype Residual
file.de <- filenameDE(wd.de, paste0("LCNEC_DE_kallisto_ens_tpm-gene-exp2-res-", covariate), dao)
file.main <- c("Mutant RB1 Differential Expression in LCNEC", "Four-Subtype Residual")
expr.res <- residualsOf(tpm.gene.exp2.log2, pheno.expr, covariate)

de.tpm.gene.exp2.res <- differentialAnalysis(expr.res, pheno.expr, dao$predictor, dao$predictor.wt, dao$test, dao$test.fdr)
de.tpm.gene.exp2.res <- copyEffectSize(de.tpm.gene.exp2.res, de.tpm.gene)
writeTable(de.tpm.gene.exp2.res, paste0(file.de, ".txt"), colnames=T, rownames=T, sep="\t")
writeTable(onlyVariable(de.tpm.gene.exp2.res, dao$effect), paste0(file.de, "_Effect>", dao$effect, ".txt"), colnames=T, rownames=T, sep="\t")   ## Only varible genes
plotVolcano(de.tpm.gene.exp2.res, dao$fdr, dao$effect, file.de, file.main, legend.x=NA)

save(de.tpm.gene.exp2.res, de.tpm.gene.exp2.res2, file=paste0(wd.de, "lcnec_kallisto_de.tpm.gene.exp2.res2_RB1NEW_Classification.RData"))

###
## Gene-level, protein-coding-gene TPMs (without 0 values in FIANL ANALYSIS)
tpm.gene.exp2.pcg.log2 <- tpm.gene.pcg.log2[rownames(de.tpm.gene.exp2.pcg),]

covariate <- "Classification2"   ## Two-Subtype Residual
file.de <- filenameDE(wd.de, paste0("LCNEC_DE_kallisto_ens_tpm-gene-exp2-pcg-res2-", covariate), dao)
file.main <- c("Mutant RB1 Differential Expression in LCNEC", "Two-Subtype Residual")
expr.res <- residualsOf(tpm.gene.exp2.pcg.log2, pheno.expr, covariate)

de.tpm.gene.exp2.pcg.res2 <- differentialAnalysis(expr.res, pheno.expr, dao$predictor, dao$predictor.wt, dao$test, dao$test.fdr)
de.tpm.gene.exp2.pcg.res2 <- copyEffectSize(de.tpm.gene.exp2.pcg.res2, de.tpm.gene.pcg)
writeTable(de.tpm.gene.exp2.pcg.res2, paste0(file.de, ".txt"), colnames=T, rownames=T, sep="\t")
writeTable(onlyVariable(de.tpm.gene.exp2.pcg.res2, dao$effect), paste0(file.de, "_Effect>", dao$effect, ".txt"), colnames=T, rownames=T, sep="\t")   ## Only varible genes
plotVolcano(de.tpm.gene.exp2.pcg.res2, dao$fdr, dao$effect, file.de, file.main, legend.x=NA)

##
covariate <- "Classification"   ## Four-Subtype Residual
file.de <- filenameDE(wd.de, paste0("LCNEC_DE_kallisto_ens_tpm-gene-exp2-pcg-res-", covariate), dao)
file.main <- c("Mutant RB1 Differential Expression in LCNEC", "Four-Subtype Residual")
expr.res <- residualsOf(tpm.gene.exp2.pcg.log2, pheno.expr, covariate)

de.tpm.gene.exp2.pcg.res <- differentialAnalysis(expr.res, pheno.expr, dao$predictor, dao$predictor.wt, dao$test, dao$test.fdr)
de.tpm.gene.exp2.pcg.res <- copyEffectSize(de.tpm.gene.exp2.pcg.res, de.tpm.gene.pcg)
writeTable(de.tpm.gene.exp2.pcg.res, paste0(file.de, ".txt"), colnames=T, rownames=T, sep="\t")
writeTable(onlyVariable(de.tpm.gene.exp2.pcg.res, dao$effect), paste0(file.de, "_Effect>", dao$effect, ".txt"), colnames=T, rownames=T, sep="\t")   ## Only varible genes
plotVolcano(de.tpm.gene.exp2.pcg.res, dao$fdr, dao$effect, file.de, file.main, legend.x=NA)

save(de.tpm.gene.exp2.pcg.res, de.tpm.gene.exp2.pcg.res2, file=paste0(wd.de, "lcnec_kallisto_de.tpm.gene.exp2.cpg.res2_RB1NEW_Classification.RData"))




###
## Gene-level, protein-coding-gene TPMs (without 0 values in ALL)
covariate <- "Classification2"   ## Two-Subtype Residual
file.de <- filenameDE(wd.de, paste0("LCNEC_DE_kallisto_ens_tpm-gene-exp-pcg-res-", covariate), dao)
file.main <- c("Mutant RB1 Differential Expression in LCNEC", "Two-Subtype Residual")
expr.res <- residualsOf(tpm.gene.pcg.log2, pheno.expr, covariate)

de.res2 <- differentialAnalysis(expr.res, pheno.expr, dao$predictor, dao$predictor.wt, dao$test, dao$test.fdr)
de.res2 <- copyEffectSize(de.res2, de.tpm.gene.pcg)
writeTable(de.res2, paste0(file.de, ".txt"), colnames=T, rownames=T, sep="\t")
writeTable(onlyVariable(de.res2, dao$effect), paste0(file.de, "_Effect>", dao$effect, ".txt"), colnames=T, rownames=T, sep="\t")   ## Only varible genes
plotVolcano(de.res2, dao$fdr, dao$effect, file.de, file.main, legend.x=NA)

##
covariate <- "Classification"   ## Four-Subtype Residual
file.de <- filenameDE(wd.de, paste0("LCNEC_DE_kallisto_ens_tpm-gene-exp-pcg-res-", covariate), dao)
file.main <- c("Mutant RB1 Differential Expression in LCNEC", "Four-Subtype Residual")
expr.res <- residualsOf(tpm.gene.pcg.log2, pheno.expr, covariate)

de.res <- differentialAnalysis(expr.res, pheno.expr, dao$predictor, dao$predictor.wt, dao$test, dao$test.fdr)
de.res <- copyEffectSize(de.res, de.tpm.gene.pcg)
writeTable(de.res, paste0(file.de, ".txt"), colnames=T, rownames=T, sep="\t")
writeTable(onlyVariable(de.res, dao$effect), paste0(file.de, "_Effect>", dao$effect, ".txt"), colnames=T, rownames=T, sep="\t")   ## Only varible genes
plotVolcano(de.res, dao$fdr, dao$effect, file.de, file.main, legend.x=NA)

save(de.res, de.res2, file=paste0(wd.de, "lcnec_kallisto_de.tpm.gene.exp.cpg.res_RB1NEW_Classification.RData"))

###
## Gene-level, EXPRESSED, protein-coding-gene TPMs (without 0 values)
covariate <- "Classification2"   ## Two-Subtype Residual
file.de <- filenameDE(wd.de, paste0("LCNEC_DE_kallisto_ens_tpm-gene-exp-pcg-res-", covariate), dao)
file.main <- c("Mutant RB1 Differential Expression in LCNEC", "Two-Subtype Residual")
expr.res <- residualsOf(tpm.gene.exp.pcg.log2, pheno.expr, covariate)

de.res2 <- differentialAnalysis(expr.res, pheno.expr, dao$predictor, dao$predictor.wt, dao$test, dao$test.fdr)
de.res2 <- copyEffectSize(de.res2, de.tpm.gene.exp.pcg)
writeTable(de.res2, paste0(file.de, ".txt"), colnames=T, rownames=T, sep="\t")
writeTable(onlyVariable(de.res2, dao$effect), paste0(file.de, "_Effect>", dao$effect, ".txt"), colnames=T, rownames=T, sep="\t")   ## Only varible genes
plotVolcano(de.res2, dao$fdr, dao$effect, file.de, file.main, legend.x=NA)

##
covariate <- "Classification"   ## Four-Subtype Residual
file.de <- filenameDE(wd.de, paste0("LCNEC_DE_kallisto_ens_tpm-gene-exp-pcg-res-", covariate), dao)
file.main <- c("Mutant RB1 Differential Expression in LCNEC", "Four-Subtype Residual")
expr.res <- residualsOf(tpm.gene.exp.pcg.log2, pheno.expr, covariate)

de.res <- differentialAnalysis(expr.res, pheno.expr, dao$predictor, dao$predictor.wt, dao$test, dao$test.fdr)
de.res <- copyEffectSize(de.res, de.tpm.gene.exp.pcg)
writeTable(de.res, paste0(file.de, ".txt"), colnames=T, rownames=T, sep="\t")
writeTable(onlyVariable(de.res, dao$effect), paste0(file.de, "_Effect>", dao$effect, ".txt"), colnames=T, rownames=T, sep="\t")   ## Only varible genes
plotVolcano(de.res, dao$fdr, dao$effect, file.de, file.main, legend.x=NA)

save(de.res, de.res2, file=paste0(wd.de, "lcnec_kallisto_de.tpm.gene.exp.cpg.res_RB1NEW_Classification.RData"))

# -----------------------------------------------------------------------------
# Replace Ensembl Gene IDs to gene name in Reactome results
# -----------------------------------------------------------------------------
wd.reactome <- "/Users/tpyang/Work/uni-koeln/tyang2/LCNEC/analysis/expression/kallisto/lcnec-rnaseq-de-kallisto-ensembl/de.tpm.gene.exp2.pcg/pathway/"
#wd.reactome <- "/Users/tpyang/Work/uni-koeln/tyang2/LCNEC/analysis/expression/kallisto/lcnec-rnaseq-de-kallisto-ensembl/de.tpm.gene.exp/pathway/"
list <- readTable(paste0(wd.reactome, "untitled.txt"), header=T, rownames=T, sep="\t")
colnames(list) <- c("ensembl_gene_id",	"external_gene_name")

reactome <- read.csv(paste0(wd.reactome, "result.csv"))
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
writeTable(reactome, paste0(wd.reactome, "result.tsv"), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Principal component analysis (PCA)
# http://www.bioconductor.org/help/workflows/rnaseqGene/#eda
# https://onlinecourses.science.psu.edu/stat505/node/55
# Last Modified: 29/03/17
# -----------------------------------------------------------------------------
pca <- getPCA(t(expr.magic))
#summary(pca)

##
wd.pca <- paste0(wd, "analysis/expression/data/")
#samples.expr <- names(expr.d.e)
#pheno.expr <- pheno[samples.expr,]

traits <- pheno.expr[,c("Sex", "Smoking_Hx", "Stage_UICC")]
colnames(traits) <- c("Sex", "Smoking History", "UICC Stage")
plotPCAs(1, 2, pca, traits, wd.pca, NA, 75, -50)   ## TO-DO

traits <- pheno.expr[,c("Classification", "Survival_censor")]
colnames(traits) <- c("Histological Subtype", "Survival")
plotPCAs(1, 2, pca, traits, wd.pca, NA, 75, -50, c("deepskyblue","purple3","red","blue"))   ## TO-DO

traits <- pheno.expr[,c("Age", "Classification2")]
colnames(traits) <- c("Age", "Histological Subtype")
plotPCAs(1, 2, pca, traits, wd.pca, NA, 60, -50)   ## TO-DO





# -----------------------------------------------------------------------------
# Molecular classification / PCA
# Last Modified: 10/02/17
# -----------------------------------------------------------------------------
pca.de <- getPCA(t(tpm.gene.exp2.pcg.log2[rownames(significantAndVariable(de.tpm.gene.exp2.pcg, effect, fdr)),]))   ## BUG FIX 13/02/17: Perform PCA using normalised data
pca.de.res2 <- getPCA(t(tpm.gene.exp2.pcg.log2[rownames(significantAndVariable(de.tpm.gene.exp2.pcg.res2, effect, fdr)),]))
pca.de.res <- getPCA(t(tpm.gene.exp2.pcg.log2[rownames(significantAndVariable(de.tpm.gene.exp2.pcg.res, effect, fdr)),]))

###
## RB1 status on D.E genes
pheno.expr <- pheno[samples.expr,]

trait <- pheno.expr[,"RB1NEW"]
trait[which(trait == 0)] <- "WT"
trait[which(trait == 1)] <- "RB1"

file.main <- c("LCNEC RB1 Status on 54 D.E. Genes", "TPM Expression")
plotPCA(1, 2, pca.de, trait, wd.de, "RB1_54DE", file.main, NA, NA, c("purple", "red", "deepskyblue", "blue"))

file.main <- c("LCNEC RB1 Status on 20 D.E. Genes", "Two-Subtype Residual")
plotPCA(1, 2, pca.de.res2, trait, wd.de, "RB1_20DE", file.main, NA, NA, c("purple", "red", "deepskyblue", "blue"))

file.main <- c("LCNEC RB1 Status on 7 D.E. Genes", "Four-Subtype Residual")
plotPCA(1, 2, pca.de.res, trait, wd.de, "RB1_7DE", file.main, NA, NA, c("purple", "red", "deepskyblue", "blue"))




## Find outliers
outliers <- subset(pca.de[,c("PC1", "PC2")], PC1 > 2.5)
subset(cbind(outliers, pheno[rownames(outliers),]), RB1NEW == 0)

outliers <- subset(subset(pca.de[,c("PC1", "PC2")], PC2 > 2), PC1 > 0.5)
subset(cbind(outliers, pheno[rownames(outliers),]), RB1NEW != "NA")

outliers <- subset(pca.de[,c("PC1", "PC2")], PC1 < 0.5)
subset(cbind(outliers, pheno[rownames(outliers),]), RB1NEW == 1)

outliers <- pca.de[c("S00714", "S02229"), c("PC1", "PC2")]
subset(cbind(outliers, pheno[rownames(outliers),]))

## Find outliers
outliers <- subset(pca.de.resid[,c("PC1", "PC2")], PC1 < 1)
subset(cbind(outliers, pheno[rownames(outliers),]), RB1NEW == 1)

outliers <- subset(subset(pca.de.resid[,c("PC1", "PC2")], PC1 > 0), PC1 < 1)
subset(cbind(outliers, pheno[rownames(outliers),]), RB1NEW == 0)

outliers <- pca.de.resid[c("S00714", "S02229"), c("PC1", "PC2")]
subset(cbind(outliers, pheno[rownames(outliers),]))

###
## LCNEC subtypes on D.E genes
trait <- pheno.expr[,"Classification"]

file.main <- c("LCNEC subtypes on 51 D.E. Genes", "FPKM Expression")
plotPCA(1, 2, pca.de, trait, wd.de, "Classification_51DE", file.main, 5, -2.2, c("deepskyblue","gold","blue","red"))

file.main <- c("LCNEC subtypes on 18 D.E. Genes", "2-Subtype Residuals")
plotPCA(1, 2, pca.de.resid2, trait, wd.de, "Classification_18DE", file.main, 2.97, -2.87, c("deepskyblue","gold","blue","red"))

file.main <- c("LCNEC subtypes on 4 D.E. Genes", "4-Subtype Residuals")
plotPCA(1, 2, pca.de.resid, trait, wd.de, "Classification_4DE", file.main, 1.427, -1.365, c("deepskyblue","gold","blue","red"))

## Find outliers
outliers <- subset(subset(pca.de[,c("PC1", "PC2")], PC2 > 2.5), PC1 > 0)
subset(cbind(outliers, pheno[rownames(outliers),]), Classification != "LCNEC")

# -----------------------------------------------------------------------------
# Fisher's exact (RB1NEW vs WT)
# Last Modified: 08/02/17
# -----------------------------------------------------------------------------
trait1 <- "RB1NEW"
trait2 <- "Classification2"
table <- getFishersTable(pheno.expr, trait1, trait2) 
table
#                           RB1NEW_0 RB1NEW_1
# Classification2_LCNEC           30       14
# Classification2_LCNEC+Mix        4        6
testFishers(table)
# [1] 0.1470132

nrow(subset(subset(pheno.expr, Classification == "LCNEC"), RB1NEW != 1))
# [1] 30
nrow(subset(subset(pheno.expr, Classification == "LCNEC"), RB1NEW == 1))
# [1] 14
nrow(subset(subset(pheno.expr, Classification != "LCNEC"), RB1NEW != 1))
# [1] 4
nrow(subset(subset(pheno.expr, Classification != "LCNEC"), RB1NEW == 1))
# [1] 6
fisher.test(rbind(c(30, 14), c(4, 6)))$p.value
# [1] 0.1470132

# -----------------------------------------------------------------------------
# Differential expression (LCNEC vs LCNEC-Mix)
# Last Modified: 28/02/17
# -----------------------------------------------------------------------------
samples <- sort(rownames(pheno))
samples.expr <- intersect(colnames(expr.d.e.log2), samples)
pheno.expr <- pheno[samples.expr,]

wd.de <- paste0(wd, "analysis/expression/de/")
test <- "Wilcox"   ## Wilcoxon/Wilcox/U/Students/ttest
test.fdr <- "Q"    ## Q/BH
effect <- 1
fdr <- 0.05

## DE (LCNEC vs Mix)
predictor <- "Classification2"
predictor.wt <- "LCNEC"
file.de.sub <- paste0(wd.de, "LCNEC_DE_SUB_RPKM_", predictor, "_", test, "_", test.fdr)

results.de.sub <- differentialAnalysis(expr.d.e.log2, pheno.expr, predictor, predictor.wt, test, test.fdr)
writeTable(results.de.sub, paste0(file.de.sub, ".txt"), colnames=T, rownames=T, sep="\t")

## Only varible genes
writeTable(onlyVariable(results.de.sub, effect), paste0(file.de, "_Effect>", effect, ".txt"), colnames=T, rownames=T, sep="\t")





# -----------------------------------------------------------------------------
# Differential expression on the residuals of expression (LCNEC vs Mix; n=69-15NA)
# Last Modified: 01/03/17
# -----------------------------------------------------------------------------
covariate <- "RB1NEW"

## DE on the residuals (LCNEC vs Mix)
expr.resid <- residualsOf(expr.d.e.log2, pheno.expr, covariate)
results.de.resid <- differentialAnalysis(expr.resid, pheno.expr, predictor, predictor.wt, test, test.fdr)

## Use the original FPKM effect size
results.de.resid[,8:10] <- results.de.sub[rownames(results.de.resid), 8:10]
writeTable(results.de.resid, paste0(file.de, "_Resid.txt2"), colnames=T, rownames=T, sep="\t")

## Only varible genes
writeTable(onlyVariable(results.de.resid, effect), paste0(file.de, "_Resid_Effect>", effect, ".txt2"), colnames=T, rownames=T, sep="\t")





# -----------------------------------------------------------------------------
# Sandbox: Plot un-detected, un-expressed genes
# Last Modified: 07/03/17
# -----------------------------------------------------------------------------
expr.log2     <- log2(1 + expr)
expr.d.log2   <- log2(1 + expr.d)
expr.d.e.log2

undetected <- setdiff(rownames(expr.log2), rownames(expr.d.log2))
unexpressed <- setdiff(rownames(expr.d.log2), rownames(expr.d.e.log2))

expr.unexpressed.log2 <- expr.d.log2[unexpressed,]
medians.unexpressed.log2 <- median0(expr.unexpressed.log2)
unexpressed.1 <- unexpressed[which(medians.unexpressed.log2 > 1)]
unexpressed.2 <- unexpressed[which(medians.unexpressed.log2 < 1)]

expr.removed.log2 <- expr.log2[c(undetected, unexpressed.1, unexpressed.2),]










# -----------------------------------------------------------------------------
# Plot differential analysis (RB1NEW vs WT)
# Last Modified: 09/02/17
# -----------------------------------------------------------------------------
plotDE <- function(expr, samples.expr.mut, name, name2, legend, legend.y) {
   red  <- rgb(1, 0, 0, 0.5)
   blue <- rgb(0, 0, 1, 0.5)
   samples.expr.wt <- setdiff(colnames(expr), samples.expr.mut)
   xmin <- min(expr[name,])
   xmax <- max(as.numeric(expr[name,]))
   
   pdf(paste0("fpkm-d-e-log2_", name2, "_", name, ".pdf"))
   hist(as.numeric(expr[name, samples.expr.mut]), col=red, xlim=c(xmin, xmax+0.1), ylim=c(0, legend.y), main=paste0(name2, " (", name, ")"), xlab="log2(1 + FPKM)")
   hist(as.numeric(expr[name, samples.expr.wt]), col=blue, add=T)
 
   box()
   legend(x=xmin, y=legend.y, legend=legend, fill=c(red, blue), border=c(red, blue))
   dev.off()
}

legend <- c("RB1 (n=20)", "WT (n=49)")
plotDE(expr.d.e.log2, samples.expr.mut, "NM_004526", "MCM2", legend, 12)
plotDE(expr.d.e.log2, samples.expr.mut, "NM_053056", "CCND1", legend, 14)

plotDE(expr.d.e.log2, samples.expr.mut, "NM_000465", "BARD1", legend, 11)
plotDE(expr.d.e.log2, samples.expr.mut, "NM_203394", "E2F7", legend, 13)
plotDE(expr.d.e.log2, samples.expr.mut, "NM_001204192", "TP73", legend, 30)
plotDE(expr.d.e.log2, samples.expr.mut, "NM_058195", "CDKN2A", legend, 11)

# -----------------------------------------------------------------------------
# Manual validation
# Last Modified: 09/02/17
# -----------------------------------------------------------------------------
samples.expr.mut
#  [1] "S00550" "S00567" "S00587" "S00609" "S00710" "S00718" "S00732" "S00739" "S00794" "S00795" "S01407" "S01496" "S01499" "S01505" "S01508" "S01574" "S01581" "S02236"
# [19] "S02259" "S02266"

getRefGene("name2", "MCM2")
# name chrom strand   txStart     txEnd  cdsStart    cdsEnd exonCount
# NR_073375 NR_073375  chr3      + 127317199 127341278 127341278 127341278        16
# NM_004526 NM_004526  chr3      + 127317199 127341278 127317309 127340616        16

wilcox.test(as.numeric(expr.d.e.log2["NM_004526",]) ~ pheno.expr[,"RB1NEW"], exact=F)$p.value
# [1] 2.583733e-07
median(as.numeric(expr.d.e.log2["NM_004526", samples.expr.wt]))
# [1] 4.515555
median(as.numeric(expr.d.e.log2["NM_004526", samples.expr.mut]))
# [1] 5.656553

# > wilcox.test(as.numeric(expr["NM_004526",]) ~ pheno.expr[,trait], exact=F)
# Wilcoxon rank sum test with continuity correction
# data:  as.numeric(expr["NM_004526", ]) by pheno.expr[, trait]
# W = 100, p-value = 2.584e-07
# alternative hypothesis: true location shift is not equal to 0
# 
# > wilcox.test(as.numeric(expr["NM_004526",]) ~ pheno.expr[,trait])
# Wilcoxon rank sum test
# data:  as.numeric(expr["NM_004526", ]) by pheno.expr[, trait]
# W = 100, p-value = 1.593e-08
# alternative hypothesis: true location shift is not equal to 0

# -----------------------------------------------------------------------------
# Differential expression (LCNEC vs LCNEC-Mix; n=69)
# Last Modified: 06/02/17
# -----------------------------------------------------------------------------
wd.de <- paste0(wd, "analysis/expression/de/")
test <- "Wilcox"   ## Wilcoxon/Wilcox/U/Students/ttest
test.fdr <- "Q"    ## Q/BH
effect <- 1
fdr <- 0.05

###
## DE (LCNEC vs LCNEC-Mix)
pheno.expr <- pheno[samples.expr,]   ## NOTE
predictor <- "Classification2"
predictor.wt <- "LCNEC"
file.de <- paste0(wd.de, "LCNEC_DE_RPKM_", predictor, "_", test, "_", test.fdr)

results.de <- differentialAnalysis(expr.d.e.log2, pheno.expr, predictor, predictor.wt, test, test.fdr)
writeTable(results.de, paste0(file.de, ".txt"), colnames=T, rownames=T, sep="\t")

## Only varible genes
pdf(paste0(wd.de, "hist_expr-d-e_effect.pdf"))
hist(abs(results.de$Effect))
dev.off()

writeTable(onlyVariable(results.de, effect), paste0(file.de, "_Effect>", effect, ".txt"), colnames=T, rownames=T, sep="\t")

# -----------------------------------------------------------------------------
# Volcano plots (LCNEC vs LCNEC-Mix)
# Last Modified: 07/02/17
# -----------------------------------------------------------------------------
file.vol <- paste0(file.de, "_Effect>", effect, ".pdf")
file.main <- "55 LCNEC vs 14 LCNEC-Mix"
plotVolcano(results.de, fdr, effect, file.vol, file.main, 1.8)
