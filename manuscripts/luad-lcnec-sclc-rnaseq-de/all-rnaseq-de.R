# =============================================================================
# Manuscript: Journey to the centre of the cancer cell
# Chapter I:
# Name:   all-rnaseq-de.R
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 31/10/17 (Reformation Day)
# =============================================================================
#wd.source <- "/projects/cangen/tyang2/dev/R"              ## tyang2@cheops
#wd.source <- "/ngs/cangen/tyang2/dev/R"                   ## tyang2@gauss
wd.source <- "/Users/tpyang/Work/dev/R"                    ## tpyang@localhost

handbooks <- c("Common.R", "DifferentialExpression.R")     ## Required handbooks/libraries for the manuscript
invisible(sapply(handbooks, function(b) source(file.path(wd.source, "handbook-of", b))))

load(file.path(wd.source, "guide-to-the", "hg19.RData"))   ## The bioinformatician's guide to the human genome

# -----------------------------------------------------------------------------
# Transcript TPM estimates/aboundants from kallisto (v0.43.0)
# -----------------------------------------------------------------------------
getTSVS <- function(wd, samples, BASE) {
   wd.kallisto <- paste0(wd, "ngs/RNA/kallisto_hg19.ensembl_quant-b100--bias/")
   wd.de <- paste0(wd, "analysis/expression/kallisto/", tolower(BASE), "-rnaseq-de-kallisto-ensembl/")
   
   return(sapply(samples, function(s) file.path(wd.kallisto, s)))
}

wd <- paste0("/ngs/cangen/tyang2/LUAD/")
samples.luad <- readTable(paste0(wd, "ngs/RNA/luad_rna_n48.list"), header=F, rownames=F, sep="\t")
tsvs.luad <- getTSVS(wd, samples.luad, "LUAD")

wd <- paste0("/ngs/cangen/tyang2/LCNEC/")
samples.lcnec <- readTable(paste0(wd, "ngs/RNA/lcnec_rna_n69.list"), header=F, rownames=F, sep="\t")
tsvs.lcnec <- getTSVS(wd, samples.lcnec, "LCNEC")

wd <- paste0("/ngs/cangen/tyang2/SCLC/")
samples.sclc <- readTable(paste0(wd, "ngs/RNA/sclc_rna_n81.list"), header=F, rownames=F, sep="\t")[,1]
tsvs.sclc <- getTSVS(wd, samples.sclc, "SCLC")

###
##
library("sleuth")
wd.all <- "/Users/tpyang/Work/uni-koeln/tyang2/ALL/"
wd.all.de <- paste0(wd.all, "analysis/expression/kallisto/all-rnaseq-de/")
setwd(wd.all)

tsvs <- c(tsvs.luad, tsvs.lcnec, tsvs.sclc)
samples <- c(samples.luad, samples.lcnec, samples.sclc)
s2c <- data.frame(path=tsvs, sample=samples, stringsAsFactors=F)
t2g <- tx2Ens(ensGene)

so <- sleuth_prep(s2c, ~1, target_mapping=t2g)
# reading in kallisto results
# dropping unused factor levels
# .....................................................................
# normalizing est_counts
# 88145 targets passed the filter
# normalizing tpm
# merging in metadata
# normalizing bootstrap samples
# summarizing bootstraps

tpm.norm.filt <- sleuth_to_matrix0(so, "obs_norm_filt", "tpm")$data
tpm.norm.filt <- orderBySamples(tpm.norm.filt)
save(tpm.norm.filt, file=paste0(wd.all.de, "data/all_kallisto_tpm.norm.filt.RData"))

#kallisto.table <- kallisto_table(so, use_filtered=T, normalized=T)
#save(kallisto.table, file=paste0(wd.de, "data/sclc_kallisto_kallisto.table.RData"))

# -----------------------------------------------------------------------------
# Gene-level TPM estimates/aboundants from kallisto
# Last Modified: 30/06/17
# -----------------------------------------------------------------------------
transcripts <- intersect(rownames(t2g), rownames(tpm.norm.filt))
tpm <- tpm.norm.filt[transcripts,]   ## Only keep transcripts in Ensembl BioMart annotation (see guide-to-the/hg19.R)

tpm.gene <- getGeneTPM(t2g, tpm)
codings <- intersect(rownames(subset(ensGene.gene, gene_biotype == "protein_coding")), rownames(tpm.gene))
tpm.gene.pcg <- tpm.gene[codings,]
# > nrow(tpm.norm.filt)
# [1] 88145
# > nrow(tpm)
# [1] 83059
# > nrow(tpm.gene)
# [1] 18898
# > nrow(tpm.gene.pcg)       ## Gene-level, protein-coding-gene TPMs (with 0 values)
# [1] 16897
# > nrow(tpm.gene.pcg.exp3)
# [1] 15368
# > nrow(tpm.gene.pcg.exp)   ## Gene-level, protein-coding-gene TPMs (without 0 values)
# [1] 13759
save(tpm, tpm.gene, file=paste0(wd.de, "data/all_kallisto_tpm.gene.RData"))

##
rownames <- rownames(tpm.gene.pcg)
unexpressed.0 <- rownames[getNotExpressed(tpm.gene.pcg[, rownames(subset(pheno.all, Cancer_Type == 0))])]
unexpressed.1 <- rownames[getNotExpressed(tpm.gene.pcg[, rownames(subset(pheno.all, Cancer_Type == 1))])]
unexpressed.2 <- rownames[getNotExpressed(tpm.gene.pcg[, rownames(subset(pheno.all, Cancer_Type == 2))])]

#unexpressed <- intersect(unexpressed.0, intersect(unexpressed.1, unexpressed.2))
#tpm.gene.pcg.exp3 <- tpm.gene.pcg[setdiff(rownames, unexpressed),]
#tpm.gene.pcg.exp3.non <- tpm.gene.pcg[unexpressed,]

#test <- mapply(x = 1:nrow(tpm.gene.pcg.exp3.non), function(x) length(which(tpm.gene.pcg.exp3.non[x,] == 0)))

##
rownames <- rownames(tpm.gene.pcg)
unexpressed.half.0 <- rownames[getNotExpressedInHalf(tpm.gene.pcg[, rownames(subset(pheno.all, Cancer_Type == 0))])]
unexpressed.half.1 <- rownames[getNotExpressedInHalf(tpm.gene.pcg[, rownames(subset(pheno.all, Cancer_Type == 1))])]
unexpressed.half.2 <- rownames[getNotExpressedInHalf(tpm.gene.pcg[, rownames(subset(pheno.all, Cancer_Type == 2))])]

unexpressed.half <- intersect(unexpressed.half.0, intersect(unexpressed.half.1, unexpressed.half.2))

##
pdf("/Users/tpyang/Work/uni-koeln/tyang2/ALL/analysis/expression/kallisto/all-rnaseq-de-kallisto-ensembl/data/hist_tpm.gene.pcg.exp_n198.pdf")
plot(hist(median0(tpm.gene.pcg.exp.log2)))
dev.off()

##
unexpressed <- mapply(x = 1:nrow(tpm.gene.pcg.log2), function(x) length(which(tpm.gene.pcg.log2[x,] == log2(0.01))))
pdf("/Users/tpyang/Work/uni-koeln/tyang2/ALL/analysis/expression/kallisto/all-rnaseq-de-kallisto-ensembl/data/hist_tpm.gene.pcg_unexpressed_n198.pdf")
plot(hist(unexpressed))
dev.off()

pdf("/Users/tpyang/Work/uni-koeln/tyang2/ALL/analysis/expression/kallisto/all-rnaseq-de-kallisto-ensembl/data/hist_tpm.gene.pcg_unexpressed-20_n198.pdf")
plot(hist(unexpressed[which(unexpressed <= 20)]))
dev.off()

# -----------------------------------------------------------------------------
# Phenotypes
# Last Modified: 17/08/17
# -----------------------------------------------------------------------------
#load("/ngs/cangen/tyang2/ALL/analysis/expression/kallisto/all-rnaseq-de-kallisto-ensembl/data/all_kallisto_tpm.gene.RData")
#pheno.all <- readTable("/ngs/cangen/tyang2/ALL/analysis/expression/kallisto/all-rnaseq-de-kallisto-ensembl/data/pheno.all_n198.txt", header=T, rownames=T, sep="\t")
#codings <- intersect(rownames(subset(ensGene.gene, gene_biotype == "protein_coding")), rownames(tpm.gene))
#tpm.gene.pcg <- tpm.gene[codings,]

load("/Users/tpyang/Work/uni-koeln/tyang2/ALL/analysis/expression/kallisto/all-rnaseq-de-kallisto-ensembl/data/all_kallisto_tpm.gene.pcg.exp.log2_n198.RData")#writeTable(pheno.all, file="/Users/tpyang/Work/uni-koeln/tyang2/ALL/analysis/expression/kallisto/all-rnaseq-de-kallisto-ensembl/data/pheno.all.txt", colnames=T, rownames=F, sep="\t")
pheno.all <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/ALL/analysis/expression/kallisto/all-rnaseq-de-kallisto-ensembl/data/pheno.all_n198.txt", header=T, rownames=T, sep="\t")

## SCLC
sclc.rb1 <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/SCLC/ngs/WGS/sclc_wgs_n101_rb92.list", header=F, rownames=F, sep="")
sclc.rna <- rownames(subset(pheno.all, Cancer_Type == 2))
sclc.rna.rb1 <- intersect(sclc.rna, sclc.rb1)
sclc.rna.wt  <- setdiff(sclc.rna, sclc.rna.rb1)
pheno.all[sclc.rna.rb1,]$RB1 <- "RB1"
pheno.all[sclc.rna.wt,]$RB1 <- "WT"

## LUAD
luad <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/LUAD/ngs/WGS/luad_ngs_n230.list", header=F, rownames=F, sep="")
luad.rb1 <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/LUAD/ngs/WGS/luad_ngs_n230_rb10.list", header=F, rownames=F, sep="") 

writeTable(pheno.all, "/Users/tpyang/Work/uni-koeln/tyang2/ALL/analysis/expression/kallisto/all-rnaseq-de-kallisto-ensembl/data/pheno.all_n198.txt", colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Find D.E genes between SCLC, LUAD and LCNEC (with "additive" effects)
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
#tpm.gene.pcg.exp <- tpm.gene.pcg[getExpressedGTEx(tpm.gene.pcg),]   ## Gene-level, protein-coding-gene TPMs (without 0 values)
#tpm.gene.pcg.exp.log2 <- log2(tpm.gene.pcg.exp + 0.01)

###
## !!!
expr.pheno.log2 <- tpm.gene.pcg.exp3.log2[,rownames(pheno.all)]   ## VERY VERY VERY IMPORTANT!!!
#expr.pheno.log2 <- tpm.gene.pcg.log2[,rownames(pheno.all)]   ## VERY VERY VERY IMPORTANT!!!

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
annot.gene <- ensGene.gene[,c("ensembl_gene_id", "external_gene_name")]
#overlaps <- intersect(rownames(de), rownames(annot.gene))         ## BE EXTRA CAREFUL!! NOT intersect(rownames(annot.gene), rownames(de)) 
de.all.tpm.gene.pcg.exp3 <- cbind(annot.gene[rownames(de),], de[,-1])   ## BE EXTRA CAREFUL!!

save(de.all.tpm.gene.pcg.exp3, pheno.all, file="/Users/tpyang/Work/uni-koeln/tyang2/ALL/analysis/expression/kallisto/all-rnaseq-de-kallisto-ensembl/de-all-tpm-gene-pcg-exp3/de_all_kallisto_tpm.gene.pcg.exp3_anova_src_n198.RData")
writeTable(de.all.tpm.gene.pcg.exp3, file="/Users/tpyang/Work/uni-koeln/tyang2/ALL/analysis/expression/kallisto/all-rnaseq-de-kallisto-ensembl/de-all-tpm-gene-pcg-exp3/de_all_kallisto_tpm.gene.pcg.exp3_anova_src_n198.txt", colnames=T, rownames=F, sep="\t")

##
annot.gene <- ensGene.gene[,c("ensembl_gene_id", "external_gene_name")]
#overlaps <- intersect(rownames(de), rownames(annot.gene))         ## BE EXTRA CAREFUL!! NOT intersect(rownames(annot.gene), rownames(de)) 
de.tpm.gene.pcg.exp <- cbind(annot.gene[rownames(de),], de[,-1])   ## BE EXTRA CAREFUL!!

save(de.tpm.gene.pcg.exp, pheno.all, file="/Users/tpyang/Work/uni-koeln/tyang2/ALL/analysis/expression/kallisto/all-rnaseq-de-kallisto-ensembl/de-tpm-gene-pcg-exp/de_all_kallisto_tpm.gene.pcg.exp_anova_src_n198.RData")
writeTable(de.tpm.gene.pcg.exp, file="/Users/tpyang/Work/uni-koeln/tyang2/ALL/analysis/expression/kallisto/all-rnaseq-de-kallisto-ensembl/de-tpm-gene-pcg-exp/de_all_kallisto_tpm.gene.pcg.exp_anova_src_n198.txt", colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Not expressed genes
# Last Modified: 17/09/17
# -----------------------------------------------------------------------------
notexp <- setdiff(rownames(tpm.gene.pcg.log2), rownames(tpm.gene.pcg.exp.log2))

sig <- rownames(rbind(subset(de.tpm.gene.pcg, SRC_RHO >= 0.6), subset(de.tpm.gene.pcg, SRC_RHO <= -0.6)))
# > length(sig)
# [1] 1008
sig.exp <- rownames(rbind(subset(de.tpm.gene.pcg.exp, SRC_RHO >= 0.6), subset(de.tpm.gene.pcg.exp, SRC_RHO <= -0.6)))
# > length(sig.exp)
# [1] 783
notsig <- rownames(rbind(subset(de.tpm.gene.pcg, SRC_RHO < 0.6), subset(de.tpm.gene.pcg, SRC_RHO > -0.6)))
notsig.exp <- rownames(rbind(subset(de.tpm.gene.pcg.exp, SRC_RHO < 0.6), subset(de.tpm.gene.pcg.exp, SRC_RHO > -0.6)))

length(intersect(notexp, sig))
length(intersect(notexp, notsig))

length(intersect(notexp, sig.exp))
length(intersect(notexp, notsig.exp))

# -----------------------------------------------------------------------------
# Not expressed genes in all cancer/tissue type
# Last Modified: 28/10/17
# -----------------------------------------------------------------------------
notexp <- setdiff(rownames(tpm.gene.pcg.log2), rownames(tpm.gene.pcg.exp3.log2))

sig <- rownames(rbind(subset(de.tpm.gene.pcg, SRC_RHO >= 0.6), subset(de.tpm.gene.pcg, SRC_RHO <= -0.6)))
# > length(sig)
# [1] 1008
sig.exp <- rownames(rbind(subset(de.all.tpm.gene.pcg.exp3, SRC_RHO >= 0.6), subset(de.all.tpm.gene.pcg.exp3, SRC_RHO <= -0.6)))
# > length(sig.exp)
# [1] 946
notsig <- rownames(rbind(subset(de.tpm.gene.pcg, SRC_RHO < 0.6), subset(de.tpm.gene.pcg, SRC_RHO > -0.6)))
notsig.exp <- rownames(rbind(subset(de.all.tpm.gene.pcg.exp3, SRC_RHO < 0.6), subset(de.all.tpm.gene.pcg.exp3, SRC_RHO > -0.6)))

length(intersect(notexp, sig))
length(intersect(notexp, notsig))

length(intersect(notexp, sig.exp))
length(intersect(notexp, notsig.exp))


# -----------------------------------------------------------------------------
# Comparisions between SCLC, LUAD and LCNEC
# Last Modified: 18/08/17
# -----------------------------------------------------------------------------
#wd <- "/ngs/cangen/tyang2/ALL/"
wd <- "/Users/tpyang/Work/uni-koeln/tyang2/ALL/"
wd.de <- paste0(wd, "analysis/expression/kallisto/all-rnaseq-de-kallisto-ensembl/de-all-tpm-gene-pcg-exp3/")

#install.packages('beeswarm')
library(beeswarm)

plotBeeswarm <- function(gene, wd.de, expr.pheno.log2, pheno.all, position) {
   ensembl_gene_id <- subset(ensGene.gene, external_gene_name == gene)$ensembl_gene_id
   gene.tpms <- cbind(t(expr.pheno.log2)[rownames(pheno.all), ensembl_gene_id], pheno.all)
   colnames(gene.tpms)[1] <- "LOG2_TPM"
 
   pdf(paste0(wd.de, "plots/beeswarm/beeswarm_kallisto_tpm.gene.pcg.exp.log2_", gene, ".pdf"), height=6, width=4)
   ymin <- min(gene.tpms$LOG2_TPM)
   ymax <- max(gene.tpms$LOG2_TPM)
   boxplot(LOG2_TPM ~ Cancer_Type, data=gene.tpms, outline=F, names=c("LUAD", "LCNEC", "SCLC"), ylim=c(ymin, ymax), ylab="log2(TPM + 0.01)", main=paste0(gene, " (", ensembl_gene_id, ")"))
 
   beeswarm(LOG2_TPM ~ Cancer_Type, data=subset(gene.tpms, RB1=="NA"), col="lightgray", pch=16, add=T)
   beeswarm(LOG2_TPM ~ Cancer_Type, data=subset(gene.tpms, RB1=="WT"), col="blue", pch=16, add=T)
   beeswarm(LOG2_TPM ~ Cancer_Type, data=subset(gene.tpms, RB1=="RB1"), col="red", pch=16, add=T)
 
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
   ensembl_gene_id <- subset(ensGene.gene, external_gene_name == gene)$ensembl_gene_id[1]   ## To avoid ENSG00000269846 (one of RBL1)
   gene.tpms <- cbind(t(expr.pheno.log2)[rownames(pheno.all), ensembl_gene_id], pheno.all)
   colnames(gene.tpms)[1] <- "LOG2_TPM"
 
   pdf(paste0(wd.de, "plots/boxplot/boxplot_kallisto_tpm.gene.pcg.exp.log2_", gene, ".pdf"), height=6, width=4)
   ymin <- min(gene.tpms$LOG2_TPM)
   ymax <- max(gene.tpms$LOG2_TPM)
   boxplot(LOG2_TPM ~ Cancer_Type, data=gene.tpms, outline=T, names=c("LUAD", "LCNEC", "SCLC"), ylim=c(ymin, ymax), ylab="log2(TPM + 0.01)", main=paste0(gene, " (", ensembl_gene_id, ")"))
 
   dev.off()
}
genes <- c("CDKN2A", "PSIP1")
for (g in 1:length(genes)) {
 #plotBox(genes[g], wd.de, tpm.gene.pcg.log2[setdiff(rownames(tpm.gene.pcg.log2), "ENSG00000269846"),], pheno.all)
 plotBeeswarm(genes[g], wd.de, tpm.gene.pcg.exp.log2, pheno.all)
}

genes <- c("TP53BP1", "XRCC1", "XRCC4", "FEN1", "EXO1", "WRN")
genes <- c("UCHL1", "REST", "HES1")
genes <- c("RBL1", "RBL2", "XRCC3")
genes <- c("RAD51B", "RAD51C", "RAD51D")
genes <- c("AHRR", "AHR", "NFKB1", "NFKB2", "TLR2")
genes <- c("TRADD", "TNFRSF1A", "CHD2")
genes <- c("RNF138", "RNF168")
genes <- c("FP15737", "LDOC1", "TMEM179", "DPYSL5", "CRMP1", "DPYSL3")
genes <- c("NOTCH2", "NOTCH3", "NOTCH4")
genes <- c("KRAS", "XRCC2", "DLL3")
genes <- c("MAX", "LIG4", "NHEJ1", "MUM1", "MYD88", "PARP3")
genes <- c("CD274", "PDCD1LG2", "IL4R", "TGFBR2")
genes <- c("XPC", "ERCC6L2", "ERCC6L", "ERCC2", "ERCC3", "ERCC6", "ERCC5", "ERCC8", "ERCC1", "ERCC4", "ATRIP")
genes <- c("MYC", "CDK6", "RB1", "TRADD", "MYD88")
for (g in 1:length(genes)) {
   #plotBox(genes[g], wd.de, tpm.gene.pcg.log2[setdiff(rownames(tpm.gene.pcg.log2), "ENSG00000269846"),], pheno.all)
   plotBeeswarm(genes[g], wd.de, tpm.gene.pcg.exp.log2, pheno.all)
}

genes <- c("CIITA", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-DRA", "HLA-DMA", "HLA-DMB", "HLA-DPA1", "HLA-DPB1", "HLA-DOA", "HLA-DOB")
genes <- c("RBL1", "RBL2")
genes <- c("TP53", "TP53BP1", "TP53BP2")
for (g in 1:length(genes)) {
   #plotBox(genes[g], wd.de, tpm.gene.pcg.log2[setdiff(rownames(tpm.gene.pcg.log2), "ENSG00000269846"),], pheno.all)
   plotBeeswarm(genes[g], wd.de, tpm.gene.pcg.log2[setdiff(rownames(tpm.gene.pcg.log2), "ENSG00000269846"),], pheno.all, "bottomleft")
}

genes <- c("CD40", "IL6", "IL10", "POLQ", "ASF1A", "ASF1B", "SMARCAL1", "TIMELESS")
genes <- c("RAG1")
genes <- c("E2F4")
genes <- c("LIG4")
genes <- c("E2F1", "E2F2", "E2F3", "E2F5", "E2F6", "E2F7", "E2F8")
for (g in 1:length(genes)) {
   plotBox(genes[g], wd.de, tpm.gene.pcg.log2[setdiff(rownames(tpm.gene.pcg.log2), "ENSG00000269846"),], pheno.all)
}

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

plotVolcanoALL <- function(de, rho, genes, file.de, file.main) {
   p <- rhoToP(rho, de)
   de.sig <- subset(de, SRC_P <= p)
   de.sig$log10P <- -log10(de.sig$SRC_P)
   #rho <- pToRHO(p, de)
 
   de$log10P <- -log10(de$SRC_P)
   xmax <- max(de$SCLC_LUAD)
   ymax <- max(de$log10P)
   #p <- max(significantAndVariable(de, effect, fdr)$P)
 
   pdf(file.de, height=7, width=7)
   plot(de$SCLC_LUAD, de$log10P, pch=16, xlim=c(-xmax, xmax), ylim=c(0, ymax), xlab=c("Effect size (log2FC SCLC/LUAD)", "", ""), ylab="Significance (-log10 P-value)", col="darkgray", main=file.main)
   abline(h=c(-log10(rhoToP(0.4, de))), lty=5, col="darkgray")
   text(xmax*-1 + 2*xmax/29, -log10(rhoToP(0.4, de)) + ymax/42, paste0("rho=±", 0.4), col="darkgray")
   
   de.up <- subset(de.sig, SCLC_LUAD > 0)
   points(de.up$SCLC_LUAD, de.up$log10P, pch=16, col="red")

   de.down <- subset(de.sig, SCLC_LUAD < 0)
   points(de.down$SCLC_LUAD, de.down$log10P, pch=16, col="dodgerblue")
   
   for (g in 1:length(genes)) {
      gene <- subset(de, external_gene_name == genes[g])
      if (nrow(gene) > 0) {
         points(gene$SCLC_LUAD, gene$log10P, pch=1, col="black")
       
         if (genes[g] == "TP53BP1" || genes[g] == "DNMT1")
            text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=c(-0.07, 1.2), cex=0.75)
         else if (genes[g] == "ATM" || genes[g] == "DLL3")
            text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=c(-0.05, 1.2), cex=0.75)
         else if (genes[g] == "")
            text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=c(-0.02, 1.2), cex=0.75)
         else if (genes[g] == "SMARCA4")
            text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=c(-0.06, 0.9), cex=0.75)
         else if (genes[g] == "WRN" || genes[g] == "KAT5")
            text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=c(-0.12, 1.2), cex=0.75)
         else if (genes[g] == "LIG4" || genes[g] == "NHEJ1")
            text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=1.15, cex=0.75)
         else if (genes[g] == "TLR6" || genes[g] == "TLR3" || genes[g] == "CD28")
            text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=1.12, cex=0.75)
         else if (genes[g] == "CDC7" || genes[g] == "XRCC3" || genes[g] == "HMGB1" || genes[g] == "CTLA4" || genes[g] == "CD86")
            text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=1.1, cex=0.75)
         else if (genes[g] == "XRCC1" || genes[g] == "RUNX1")
            text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=1.07, cex=0.75)  
         else if (genes[g] == "SMARCB1")
            text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=c(1.05, 0), cex=0.75)
         else if (genes[g] == "CD274")
            text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=c(1.11, 0.4), cex=0.75)
         else if (genes[g] == "TLR7")
            text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=c(1.14, 0.8), cex=0.75)
         else if (genes[g] == "HES1" || genes[g] == "REST" || genes[g] == "NFKB1" || genes[g] == "NFKB2" || genes[g] == "RUNX2"  || genes[g] == "TLR4")
            text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=c(0, -0.5), cex=0.75)
         else if (genes[g] == "BRIP1")
            text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=c(-0.05, -0.45), cex=0.75)
         else if (genes[g] == "BRCA2")
            text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=c(-0.1, 0), cex=0.75)
         else if (genes[g] == "TOP2A" || genes[g] == "ERCC6L")
            text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=-0.09, cex=0.75)
         else if (genes[g] == "BRCA1" || genes[g] == "WRN" || genes[g] == "TP53" || genes[g] == "EXO1" || genes[g] == "XRCC2" || genes[g] == "AHRR" || genes[g] == "MYD88" || genes[g] == "DNMT3A")
            text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=-0.1, cex=0.75)
         else if (genes[g] == "POLQ")
            text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=-0.11, cex=0.75)
         else if (genes[g] == "BLM")
            text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=-0.15, cex=0.75)
         else if (genes[g] == "XRCC4" || genes[g] == "NOTCH1" || genes[g] == "HMGB2" || genes[g] == "CDK6" || genes[g] == "ATR" || genes[g] == "RAD51B" || genes[g] == "SETD2" || genes[g] == "PARP3" || genes[g] == "TLR8" || genes[g] == "TGFBR3")
            text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=c(1.05, 1.3), cex=0.75)
         else if (genes[g] == "DNA2" || genes[g] == "RAD51D" || genes[g] == "EZH2" || genes[g] == "MAD2L2" || genes[g] == "ASCL1" || genes[g] == "DPYSL5" || genes[g] == "FBXO5" || genes[g] == "RNF138" || genes[g] == "RNASEH2A" || genes[g] == "PARP1" || genes[g] == "CDK2" || genes[g] == "PCNA" || genes[g] == "PSIP1" || genes[g] == "TMPO" || genes[g] == "MSH2" || genes[g] == "SIRT1")
            text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=c(1, -0.5), cex=0.75)
         else if (genes[g] == "RBBP8" || genes[g] == "RBL1" || genes[g] == "ERCC6L2")
            text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=c(1.05, -0.4), cex=0.75)
         else if (genes[g] == "PARP2")
            text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=c(1.06, -0.43), cex=0.75)
         else
            if (gene$SCLC_LUAD > 0)
               text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=c(0, -0.5), cex=0.75)
            else
               text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=c(1, -0.5), cex=0.75)
      } else
         print(genes[g])
   }
   
   abline(h=c(-log10(p)), lty=5)
   text(xmax*-1 + 2*xmax/29, -log10(p) + ymax/42, paste0("rho=±", 0.6))
   legend("topleft", legend=c("Upregulation", "Downregulation"), col=c("red", "dodgerblue"), pch=19)
   dev.off()
}

###
## Volcano plot
colnames <- c("Ensembl_Gene_ID", "external_gene_name", "SRC_RHO", "SRC_P", "SRC_FDR", "ANOVA_P", "ANOVA_FDR", "LUAD", "LCNEC", "SCLC", "LCNEC_LUAD", "SCLC_LUAD")
colnames(de.all.tpm.gene.pcg.exp3) <- colnames
de.all <- de.all.tpm.gene.pcg.exp3
#de.all[setdiff(rownames(de.all), "ENSG00000269846"),]

wd <- "/Users/tpyang/Work/uni-koeln/tyang2/ALL/"
wd.de <- paste0(wd, "analysis/expression/kallisto/all-rnaseq-de-kallisto-ensembl/de-all-tpm-gene-pcg-exp3/")
file.de <- paste0(wd.de, "volcanoplot_all_tpm-gene-pcg-exp3_rho6")

## Cycle   ##"ERCC6L2", "ERCC6L"
genes <- c("CDK6", "CDK4", "MCM2", "RNASEH2A", "CDKN2A", "CDKN2C", "PARP1", "PARP2", "BARD1", "BRIP1", "BRCA1", "BRCA2", "CDC7", "CCNE2", "CCND1", "E2F2", "E2F7", "E2F1", "ATM", "ATR", "CHEK1", "CHEK2", "CDK2", "POLD4", "TOPBP1", "CLSPN", "TOP2A", "TOP1", "RB1", "RBL1", "RBL2", "TP73", "TP53", "TP53BP1", "PCNA")
genes <- c("STMN1", "FBXO5", "TUBA1A", "TMPO", "MAD2L2", genes)
file.main <- c("RB1mut-associated differential expression in LUAD, LCNEC and SCLC", "", "Cell cycle pathway", "")
plot.de <- paste0(file.de, "_Cycle.pdf")
plotVolcanoALL(de.all, 0.6, genes, plot.de, file.main)

## 53BP1   ##"UBE2N", "RNF8", "DNMT3A", "DNMT3B", "DNMT3A", "KAT5", "SIRT1", "HELLS", "MAX"
genes <- c("MSH2", "MSH6", "DNMT1", "E2F2", "E2F7", "E2F1", "SMARCA4", "SMARCB1", "ARID1A", "NHEJ1", "LIG1", "CHEK2", "LIG4", "POLQ", "DNA2", "SETD2", "RNF138", "RNF168", "PSIP1", "HMGB1", "HMGB2", "PARP1", "PARP2", "PARP3", "BARD1", "BRIP1", "BRCA1", "BRCA2", "RB1", "RBL1", "RBL2", "TP53", "TP53BP1", "XRCC1", "EZH2", "XRCC2", "XRCC3", "XRCC4", "RAD51AP1", "RBBP8", "FEN1", "EXO1", "BLM", "WRN")
file.main <- c("RB1mut-associated gene expression profiling in SCLC", "", "DNA repair pathway", "")
plot.de <- paste0(file.de, "_53BP1.pdf")
plotVolcanoALL(de.all, 0.6, genes, plot.de, file.main)

## XRCC   ##"UBE2N", "RNF8" 
genes <- c("SMARCA4", "SMARCB1", "NHEJ1", "LIG4", "POLQ", "SETD2", "RNF138", "RNF168", "PSIP1", "HMGB1", "HMGB2", "PARP1", "PARP2", "PARP3", "BARD1", "BRIP1", "BRCA1", "BRCA2", "RB1", "RBL1", "RBL2", "TP53", "TP53BP1", "XRCC1", "EZH2", "XRCC2", "XRCC3", "XRCC4", "RAD51AP1", "RBBP8", "FEN1", "EXO1", "BLM", "WRN")
file.main <- c("RB1mut-associated gene expression profiling in SCLC", "", "DNA repair pathway", "")
plot.de <- paste0(file.de, "_XRCC.pdf")
plotVolcanoALL(de.all, 0.6, genes, plot.de, file.main)

## NOTCH
genes <- c("CDK9", "KRAS", "MYC", "MAX", "MYCL", "MYCN", "NEUROD1", "ASCL1", "DLL3", "PSIP1", "HMGB2", "BRD4", "EZH2", "UCHL1", "HES1", "REST", "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "EGFR", "ERBB2", "ERBB3", "DPYSL5", "CRMP1", "DPYSL3")
file.main <- c("RB1mut-associated gene expression profiling in SCLC", "", "Neuroendocrine marker", "")
plot.de <- paste0(file.de, "_NOTCH.pdf")
plotVolcanoALL(de.all, 0.6, genes, plot.de, file.main)

#file.main <- c("RB1mut-associated gene expression profiling in SCLC", "(with genes not-expressed in LUAD)", "Neuroendocrine marker", "")
#file.de <- paste0(wd.de, "volcanoplot_all_tpm-gene-pcg_rho6_NOTCH.pdf")
#plotVolcanoALL(de.tpm.gene.pcg[setdiff(rownames(de.tpm.gene.pcg), "ENSG00000269846"),], 0.6, genes, file.de, file.main)

## Immune response   ##"IL6R", "AHRR", "AHR", "IDH1"
#genes <- c("CD28", "CD86", "CD274", "PDCD1LG2", "CTLA4", "PARP3", "MYD88", "TNFSF13", "CYP1A2", "CYP1B1", "TRADD", "TNFRSF1A", "TNFRSF10B", "NFKB1", "NFKB2", "TLR2", "TLR5", "TLR3", "TLR6", "TLR7", "TLR8", "TLR4", "BCL2", "BCL2L11", "BCL3", "RUNX1", "RUNX2")
#file.main <- c("RB1mut-associated gene expression profiling in SCLC", "", "Immune response", "")
#file.de <- paste0(wd.de, "_Immune_TLRX.pdf")
#plotVolcanoALL(de.all, 0.6, genes, file.de, file.main)

genes <- c("CD274", "PDCD1LG2", "MYD88", "NFKB1", "NFKB2", "TLR2", "BCL3", "IL6", "IL6ST", "IL6R")
file.main <- c("RB1mut-associated gene expression profiling in SCLC", "", "Immune response", "")
plot.de <- paste0(file.de, "_Immune.pdf")
plotVolcanoALL(de.all, 0.6, genes, plot.de, file.main)

# -----------------------------------------------------------------------------
# Volcano plots
# Last Modified: 22/10/17
# -----------------------------------------------------------------------------
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
   xmax <- max(de$SCLC_LUAD)
   ymax <- max(de$log10P)
   #p <- max(significantAndVariable(de, effect, fdr)$P)
 
   pdf(file.de, height=7, width=7)
   plot(de$SCLC_LUAD, de$log10P, pch=16, xlim=c(-xmax, xmax), ylim=c(0, ymax), xlab=c("Effect size (log2FC SCLC/LUAD)", "", ""), ylab="Significance (-log10 P-value)", col="darkgray", main=file.main)
   abline(h=c(-log10(rhoToP(0.4, de))), lty=5, col="darkgray")
   text(xmax*-1 + 2*xmax/29, -log10(rhoToP(0.4, de)) + ymax/42, paste0("rho=±", 0.4), col="darkgray")
 
   de.up <- subset(de.sig, SCLC_LUAD > 0)
   points(de.up$SCLC_LUAD, de.up$log10P, pch=16, col="red")
 
   de.down <- subset(de.sig, SCLC_LUAD < 0)
   points(de.down$SCLC_LUAD, de.down$log10P, pch=16, col="dodgerblue")
 
   for (g in 1:length(genes)) {
      gene <- subset(de, external_gene_name == genes[g])
      if (nrow(gene) > 0) {
         points(gene$SCLC_LUAD, gene$log10P, pch=1, col="black")
   
         if (genes[g] == "TP53BP1" || genes[g] == "DNMT1")
            text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=c(-0.07, 1.2), cex=0.75)
         else if (genes[g] == "NHEJ1")
            text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=1.08, cex=0.75)
         else if (genes[g] == "LIG4" || genes[g] == "RBL2")
            text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=1.15, cex=0.75)
         else if (genes[g] == "PARP1")
            text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=c(1, -0.5), cex=0.75)
         else if (genes[g] == "BRCA2")
            text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=c(-0.1, 0), cex=0.75)
         else if (genes[g] == "TP53BP2")
            text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=-0.08, cex=0.75)
         else if (genes[g] == "XRCC4")
            text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=c(1.04, 1.4), cex=0.75)
         else
            if (gene$SCLC_LUAD > 0)
               text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=c(0, -0.5), cex=0.75)
            else
               text(gene$SCLC_LUAD, gene$log10P, genes[g], col="black", adj=c(1, -0.5), cex=0.75)
      } else
         print(genes[g])
   }
 
   abline(h=c(-log10(rhoToP(0.6, de))), lty=5)
   text(xmax*-1 + 2*xmax/29, -log10(rhoToP(0.6, de)) + ymax/42, paste0("rho=±", 0.6))
   legend("topleft", legend=c("Upregulation", "Downregulation"), col=c("red", "dodgerblue"), pch=19)
   dev.off()
}

###
## Volcano plot
de.all <- de.all.tpm.gene.pcg.exp3
#de.all[setdiff(rownames(de.all), "ENSG00000269846"),]

wd <- "/Users/tpyang/Work/uni-koeln/tyang2/ALL/"
wd.de <- paste0(wd, "analysis/expression/kallisto/all-rnaseq-de-kallisto-ensembl/de-all-tpm-gene-pcg-exp3/")
file.de <- paste0(wd.de, "volcanoplot_all_tpm-gene-pcg-exp3_rho6")

## NHEJ and alternative pathways to DSB repair
genes <- c("XRCC1", "XRCC2", "XRCC3", "XRCC4", "NHEJ1", "LIG1", "LIG3", "LIG4", "RB1", "RBL1", "RBL2", "E2F4", "TP53", "TP53BP1", "TP53BP2", "MAX")  #"BRCA1", "BRCA2", "PARP1", "PARP2", "CHEK1", "CHEK2", "MSH2", "FEN1", "BLM", "EXO1", "WRN")
file.main <- c("RB1mut-associated differential expression in LUAD, LCNEC and SCLC", "", "NHEJ and alternative pathways to DSB repair", "")
plot.de <- paste0(file.de, "_53BP1_Lite5.pdf")
plotVolcanoALLLite(de.all, 0.6, genes, plot.de, file.main)

## Proliferation / Apoptosis inhibition
genes <- c("TRADD", "TNFRSF1A", "PARP3", "BCL2", "BCL2L11", "BBC3", "PMAIP1")
file.main <- c("RB1mut-associated differential expression in LUAD, LCNEC and SCLC", "", "Apoptosis inhibition (Proliferation)", "")
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
# Differential expression (LUAD vs SCLC; n=99)
# Last Modified: 22/10/17
# -----------------------------------------------------------------------------
wd <- "/Users/tpyang/Work/uni-koeln/data/ALL/"
wd.de <- paste0(wd, "analysis/expression/")
setwd(wd)

filenameDE <- function(wd.de, title, dao) {
   return(paste0(wd.de, title, "_", dao$predictor, "_", dao$test, "_", dao$test.fdr))
}

pipeDE <- function(expr, pheno.expr, dao, file.de, annot) {
   ## Differential expression
   de <- differentialAnalysis(expr, pheno.expr, dao$predictor, dao$predictor.wt, dao$test, dao$test.fdr)
   de <- cbind(annot[rownames(de),], de)
 
   writeTable(de[-c(1)], paste0(file.de, ".txt"), colnames=T, rownames=T, sep="\t")
 
   return(de)
}

## Data access object for test parameters
## Test: Wilcoxon/Wilcox/U/Students/ttest
## FDR:  Q/BH
## DE:   small_cell_carcinoma vs adenocarcinoma as factor
dao <- data.frame(test="Wilcox", test.fdr="Q", fdr=0.05, effect=0, predictor="Cancer_Type", predictor.wt=0, stringsAsFactors=F)

pheno.all2 <- subset(pheno.all, Cancer_Type != 1)
all2.tpm.gene.pcg.log2 <- tpm.gene.pcg.log2[, rownames(pheno.all2)]
all2.tpm.gene.pcg.exp.log2 <- tpm.gene.pcg.exp.log2[, rownames(pheno.all2)]

file.de <- filenameDE(wd.de, "de.all2.tpm.pcg.", dao)
de.all2.tpm.pcg <- pipeDE(all2.tpm.gene.pcg.log2, pheno.all2, dao, file.de, annot.gene)
# Samples with MUT Cancer_Type: 81
# Samples with  WT Cancer_Type: 48
writeTable(de.all2.tpm.pcg[-c(1,3)], paste0(file.de, ".txt"), colnames=T, rownames=T, sep="\t")
save(de.all2.tpm.pcg, file=paste0(wd.de, "de.all2.tpm.pcg_Cancer_Type.RData"))

file.de <- filenameDE(wd.de, "de.all2.tpm.pcg.exp.", dao)
de.all2.tpm.pcg.exp <- pipeDE(all2.tpm.gene.pcg.exp.log2, pheno.all2, dao, file.de, annot.gene)
# Samples with MUT Cancer_Type: 81
# Samples with  WT Cancer_Type: 48
writeTable(de.all2.tpm.pcg.exp[-c(1,3)], paste0(file.de, ".txt"), colnames=T, rownames=T, sep="\t")
save(de.all2.tpm.pcg.exp, file=paste0(wd.de, "de.all2.tpm.pcg.exp_Cancer_Type.RData"))

plotVolcanoALL2 <- function(de, p, genes, file.de, file.main) {
   de.sig <- subset(de, P <= p)
   de.sig$log10P <- -log10(de.sig$P)

   de$log10P <- -log10(de$P)
   xmax <- max(de$Effect)
   ymax <- max(de$log10P)
   #p <- max(significantAndVariable(de, effect, fdr)$P)
 
   pdf(file.de, height=7, width=7)
   plot(de$Effect, de$log10P, pch=16, xlim=c(-xmax, xmax), ylim=c(0, ymax), xlab=c("Effect size (log2FC SCLC/LUAD)", "", ""), ylab="Significance (-log10 P-value)", col="darkgray", main=file.main)

   de.up <- subset(de.sig, Effect > 0)
   points(de.up$Effect, de.up$log10P, pch=16, col="red")
 
   de.down <- subset(de.sig, Effect < 0)
   points(de.down$Effect, de.down$log10P, pch=16, col="dodgerblue")
 
   for (g in 1:length(genes)) {
      gene <- subset(de, external_gene_name == genes[g])
      if (nrow(gene) > 0) {
         points(gene$Effect, gene$log10P, pch=1, col="black")
   
         if (genes[g] == "TP53BP1" || genes[g] == "DNMT1")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(-0.07, 1.2), cex=0.75)
         else if (genes[g] == "NHEJ1")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=1.07, cex=0.75)  
         else
            if (gene$Effect > 0)
               text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(0, -0.5), cex=0.75)
            else
               text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(1, -0.5), cex=0.75)
      } else
         print(genes[g])
   }
 
   abline(h=c(-log10(p)), lty=5)
   text(xmax*-1 + 2*xmax/20, -log10(p) + ymax/42, "P=1.00E-6")
   legend("topleft", legend=c("Upregulation", "Downregulation"), col=c("red", "dodgerblue"), pch=19)
   dev.off()
}
## Volcano plot
## MHC I
genes <- c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F")
file.main <- c("Differential expression between LUAD and SCLC", "", "MHC class I", "")
file.de <- paste0(wd.de, "volcanoplot_all2_tpm-gene-pcg-exp_p6_MHCI.pdf")
plotVolcanoALL2(de.all2.tpm.pcg.exp[setdiff(rownames(de.all2.tpm.pcg.exp), "ENSG00000269846"),], 1E-06, genes, file.de, file.main)

## MHC II
genes <- c("CIITA", "HLA-DRA", "HLA-DMA", "HLA-DMB", "HLA-DPA1", "HLA-DPB1", "HLA-DOA", "HLA-DOB", "HLA-DRB1", "HLA-DRB5", "HLA-DQA1", "HLA-DQB1")
file.main <- c("Differential expression between LUAD and SCLC", "", "MHC class II", "")
file.de <- paste0(wd.de, "volcanoplot_all2_tpm-gene-pcg-exp_p6_MHCII.pdf")
plotVolcanoALL2(de.all2.tpm.pcg.exp[setdiff(rownames(de.all2.tpm.pcg.exp), "ENSG00000269846"),], 1E-06, genes, file.de, file.main)





file.de <- filenameDE(wd.de, "de.ccle.rpkm.lung.pcg", dao)
de.ccle.rpkm.lung.pcg <- pipeDE(ccle.rpkm.lung.pcg.log2, ccle.pheno.lung, dao, file.de, annot.gene)
de.ccle.rpkm.lung.pcg <- de.ccle.rpkm.lung.pcg[!is.na(de.ccle.rpkm.lung.pcg$P),]
de.ccle.rpkm.lung.pcg$FDR <- testFDR(de.ccle.rpkm.lung.pcg$P, dao$test.fdr)
# > dim(de.ccle.rpkm.lung.pcg[!is.na(de.ccle.rpkm.lung.pcg$P),])
# [1] 19436     8
#writeTable(de.ccle.rpkm.lung.pcg[-c(1,3)], paste0(file.de, ".txt"), colnames=T, rownames=T, sep="\t")
save(de.ccle.rpkm.lung.pcg, file=paste0(wd.de, "de.ccle.rpkm.lung.pcg_Hist.Subtype1.RData"))








# -----------------------------------------------------------------------------
# Permutations on expression
# Last Modified: 09/09/17
# -----------------------------------------------------------------------------
load("/re/home/tyang2/ALL/analysis/expression/kallisto/all-rnaseq-de-kallisto-ensembl/data/all_kallisto_tpm.gene.pcg.log2_n198.RData")
pheno.all <- readTable("/re/home/tyang2/ALL/analysis/expression/kallisto/all-rnaseq-de-kallisto-ensembl/data/pheno.all_n198.txt", header=T, rownames=T, sep="\t")

expr.pheno.log2 <- tpm.gene.pcg.log2[,rownames(pheno.all)]   ## VERY VERY VERY IMPORTANT!!!

###
##
colnames <- c(1:1000)
permutations <- toTable(0, length(colnames), nrow(expr.pheno.log2), colnames)
rownames(permutations) <- rownames(expr.pheno.log2)

for (p in 1:1000) {
 samples <- colnames(expr.pheno.log2)
 idx <- sample(1:ncol(expr.pheno.log2), ncol(expr.pheno.log2), replace=F)
 expr.pheno.log2 <- expr.pheno.log2[,idx]
 colnames(expr.pheno.log2) <- samples
 
 #permutations[,p] <- mapply(x = 1:nrow(expr.pheno.log2), function(x) cor.test(as.numeric(expr.pheno.log2[x,]), pheno.all$RB1_RATE, method="spearman", exact=F)[[3]])
 permutations[,p] <- mapply(x = 1:nrow(expr.pheno.log2), function(x) testANOVA(x, expr.pheno.log2, pheno.all))
}

save(permutations, file="/re/home/tyang2/ALL/analysis/expression/kallisto/all-rnaseq-de-kallisto-ensembl/pheno.all.perm.RData")
#writeTable(permutations, file="/re/home/tyang2/ALL/analysis/expression/kallisto/all-rnaseq-de-kallisto-ensembl/pheno.all.perm.txt", colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Permutation results
# Last Modified: 08/09/17
# -----------------------------------------------------------------------------
load("/Users/tpyang/Work/uni-koeln/tyang2/ALL/analysis/expression/kallisto/all-rnaseq-de-kallisto-ensembl/data/all_kallisto_tpm.gene.pcg.log2_n198.RData")#writeTable(pheno.all, file="/Users/tpyang/Work/uni-koeln/tyang2/ALL/analysis/expression/kallisto/all-rnaseq-de-kallisto-ensembl/data/pheno.all.txt", colnames=T, rownames=F, sep="\t")
pheno.all <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/ALL/analysis/expression/kallisto/all-rnaseq-de-kallisto-ensembl/data/pheno.all_n198.txt", header=T, rownames=T, sep="\t")
load("/Users/tpyang/Work/uni-koeln/tyang2/ALL/analysis/expression/kallisto/all-rnaseq-de-kallisto-ensembl/de-tpm-gene-pcg/de_all_kallisto_tpm.gene.pcg_anova_src_n198.RData")

load("/Users/tpyang/Work/uni-koeln/tyang2/ALL/analysis/expression/kallisto/all-rnaseq-de-kallisto-ensembl/de-tpm-gene-pcg/permutations/all.expr.perm_1000.RData")

empiricalP <- function(p, permutations) {
 sig <- 0
 for (i in 1:length(permutations))
  if (permutations[i] <= p)
   sig++
  
  if (sig == 0)
   return(0);
 return(sig/length(permutations)) 
}

permutations <- permutations[rownames(de.tpm.gene.pcg),]
de.tpm.gene.pcg$SRC_P1000 <- NA
de.tpm.gene.pcg$SRC_P1000 <- mapply(x = 1:nrow(de.tpm.gene.pcg), function(x) empiricalP(de.tpm.gene.pcg$SRC_P[x], permutations[x,]))

##
de.tpm.gene.pcg <- de.tpm.gene.pcg[order(de.tpm.gene.pcg$SRC_P),]
length(which(de.tpm.gene.pcg$SRC_P1000 == 0))

##
expr.pheno.log2 <- expr.pheno.log2[rownames(de.tpm.gene.pcg),]

# -----------------------------------------------------------------------------
# TPM~iCN plot
# Last Modified: 19/08/17
# -----------------------------------------------------------------------------
gene <- "MYC"
cn <- readTable(paste0("/Users/tpyang/Work/uni-koeln/tyang2/TMP/SCLC/tmp/backup/analysis/8p12/data/expr.log2.", gene, ".txt"), header=T, rownames=T, sep="\t")
ensembl <- subset(ensGene.gene, external_gene_name == gene)$ensembl_gene_id
cn$TPM <- as.numeric(tpm.gene.pcg.log2[ensembl, rownames(cn)])

rho <- cor.test(as.numeric(cn$TPM), log10(cn$iCN_median), method="spearman", exact=F)[[4]]
p   <- cor.test(as.numeric(cn$TPM), log10(cn$iCN_median), method="spearman", exact=F)[[3]]

pdf(paste0("/Users/tpyang/Work/uni-koeln/tyang2/SCLC/analysis/expression/kallisto/sclc-rnaseq-de-kallisto-ensembl/plots/plot_tpm~icn_", gene, "_n67.pdf"), height=5, width=4.5)
plot(as.numeric(cn$TPM) ~ cn$iCN_median, xlab="iCN", ylab="log2(TPM + 0.01)", main=paste0(gene, " (", ensembl, ")"))
text(26.8, 0, paste0("rho=", round(rho, 2)))
text(28, -0.5, paste0("p=", formatC(p, format="e", digits=2)))
dev.off()

##
gene <- "MYCL"
cn <- readTable(paste0("/Users/tpyang/Work/uni-koeln/tyang2/TMP/SCLC/tmp/backup/analysis/8p12/data/expr.log2.", gene, ".txt"), header=T, rownames=T, sep="\t")
ensembl <- subset(ensGene.gene, external_gene_name == gene)$ensembl_gene_id
cn$TPM <- as.numeric(tpm.gene.pcg.log2[ensembl, rownames(cn)])

rho <- cor.test(as.numeric(cn$TPM), log10(cn$iCN_median), method="spearman", exact=F)[[4]]
p   <- cor.test(as.numeric(cn$TPM), log10(cn$iCN_median), method="spearman", exact=F)[[3]]

pdf(paste0("/Users/tpyang/Work/uni-koeln/tyang2/SCLC/analysis/expression/kallisto/sclc-rnaseq-de-kallisto-ensembl/plots/plot_tpm~icn_", gene, "_n67.pdf"), height=5, width=4.5)
plot(as.numeric(cn$TPM) ~ cn$iCN_median, xlab="iCN", ylab="log2(TPM + 0.01)", main=paste0(gene, " (", ensembl, ")"))
text(41.35, 1.7, paste0("rho=", round(rho, 2)))
text(43.18, 1.18, paste0("p=", formatC(p, format="e", digits=2)))
dev.off()

##
gene <- "MYCN"
cn <- readTable(paste0("/Users/tpyang/Work/uni-koeln/tyang2/TMP/SCLC/tmp/backup/analysis/8p12/data/expr.log2.", gene, ".txt"), header=T, rownames=T, sep="\t")
ensembl <- subset(ensGene.gene, external_gene_name == gene)$ensembl_gene_id
cn$TPM <- as.numeric(tpm.gene.pcg.log2[ensembl, rownames(cn)])

rho <- cor.test(as.numeric(cn$TPM), log10(cn$iCN_median), method="spearman", exact=F)[[4]]
p   <- cor.test(as.numeric(cn$TPM), log10(cn$iCN_median), method="spearman", exact=F)[[3]]

pdf(paste0("/Users/tpyang/Work/uni-koeln/tyang2/SCLC/analysis/expression/kallisto/sclc-rnaseq-de-kallisto-ensembl/plots/plot_tpm~icn_", gene, "_n67.pdf"), height=5, width=4.5)
plot(as.numeric(cn$TPM) ~ cn$iCN_median, xlab="iCN", ylab="log2(TPM + 0.01)", main=paste0(gene, " (", ensembl, ")"))
text(41.35, 1.7, paste0("rho=", round(rho, 2)))
text(43.18, 1.18, paste0("p=", formatC(p, format="e", digits=2)))
dev.off()










# -----------------------------------------------------------------------------
# PCA
# Last Modified: 10/02/17
# -----------------------------------------------------------------------------
pca.de <- getPCA(t(tpm.gene.pcg.exp.log2[rownames(de.tpm.gene.pcg.exp.p6.cd),]))   ## BUG FIX 13/02/17: Perform PCA using normalised data

###
## ALL
trait <- pheno.all[,"Cancer_Type"]

file.main <- "LCNEC RB1 Status on 54 D.E. Genes"
plotPCA(1, 2, pca.de, trait, wd.de, "RB1_54DE", file.main, NA, NA, c("purple", "red", "deepskyblue", "blue"))

# -----------------------------------------------------------------------------
# Replace Ensembl Gene IDs to gene name in Reactome results
# -----------------------------------------------------------------------------
wd.reactome <- "/Users/tpyang/Work/uni-koeln/tyang2/ALL/analysis/expression/kallisto/all-rnaseq-de-kallisto-ensembl/de-tpm-gene-pcg-exp/pathway/296DE_P1E-6/"
list <- de.tpm.gene.pcg.exp.p6.cd[,1:2]
wd.reactome <- "/Users/tpyang/Work/uni-koeln/tyang2/LCNEC/analysis/expression/kallisto/lcnec-rnaseq-de-kallisto-ensembl/de-tpm-gene-pcg-exp/pathway/129DE_FDR5%/"
list <- de.tpm.gene.pcg.exp2[,1:2]
#list <- readTable(paste0(wd.reactome, "untitled.txt"), header=T, rownames=T, sep="\t")
#colnames(list) <- c("ensembl_gene_id",	"external_gene_name")

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
# Comparisions between SCLC, LUAD and LCNEC
# Last Modified: 30/06/17
# -----------------------------------------------------------------------------
#wd <- "/ngs/cangen/tyang2/ALL/"
wd <- "/Users/tpyang/Work/uni-koeln/tyang2/ALL/"
wd.de <- paste0(wd, "analysis/expression/kallisto/all-rnaseq-de-kallisto-ensembl/")

#install.packages('beeswarm')
library(beeswarm)

plotBeeswarm <- function(gene, wd.de, expr.pheno.log2, pheno.all) {
   ensembl_gene_id <- subset(ensGene.gene, external_gene_name == gene)$ensembl_gene_id
   gene.tpms <- cbind(t(expr.pheno.log2)[rownames(pheno.all), ensembl_gene_id], pheno.all)
   colnames(gene.tpms)[1] <- "LOG2_TPM"
 
   pdf(paste0(wd.de, "plots/beeswarm/beeswarm_kallisto_tpm.gene.pcg.log2_", gene, ".pdf"), height=6, width=4)
   ymin <- min(gene.tpms$LOG2_TPM)
   ymax <- max(gene.tpms$LOG2_TPM)
   boxplot(LOG2_TPM ~ Cancer_Type, data=gene.tpms, outline=F, names=c("LUAD", "LCNEC", "SCLC"), ylim=c(ymin, ymax), ylab="log2(TPM + 0.01)", main=paste0(gene, " (", ensembl_gene_id, ")"))
 
   #beeswarm(LOG2_TPM ~ Cancer_Type, data=subset(gene.tpms, RB1=="NA"), col="lightgray", pch=16, add=T)
   #beeswarm(LOG2_TPM ~ Cancer_Type, data=subset(gene.tpms, RB1=="WT"), col="blue", pch=16, add=T)
   #beeswarm(LOG2_TPM ~ Cancer_Type, data=subset(gene.tpms, RB1=="RB1"), col="red", pch=16, add=T)
 
   #gene.tpms[intersect(rownames(gene.tpms), unique(muts.gene.sil$PAT_ID)),]$RB1 <- "RB1 (silent)"
   #gene.tpms["S00396",]$RB1 <- "RB1 (silent)"
   #gene.tpms["S00396",]$LOG2_TPM <- -7
   #gene.tpms["S00006",]$RB1 <- "RB1 (silent)"
   #gene.tpms["S00006",]$LOG2_TPM <- -7
   #beeswarm(LOG2_TPM ~ Cancer_Type, data=subset(gene.tpms, RB1=="RB1 (silent)"), col="yellow", pch=16, add=T)
 
   legend("bottomright", legend = c("RB1", "RB1 (silent)", "WT"), pch=16, col=c("red", "yellow", "blue"))
   dev.off()
}

plotBox <- function(gene, wd.de, expr.pheno.log2, pheno.all) {
   ensembl_gene_id <- subset(ensGene.gene, external_gene_name == gene)$ensembl_gene_id
   gene.tpms <- cbind(t(expr.pheno.log2)[rownames(pheno.all), ensembl_gene_id], pheno.all)
   colnames(gene.tpms)[1] <- "LOG2_TPM"
 
   pdf(paste0(wd.de, "plots/boxplot/boxplot_kallisto_tpm.gene.pcg.exp.log2_", gene, ".pdf"), height=6, width=4)
   ymin <- min(gene.tpms$LOG2_TPM)
   ymax <- max(gene.tpms$LOG2_TPM)
   boxplot(LOG2_TPM ~ Cancer_Type, data=gene.tpms, outline=T, names=c("LUAD", "LCNEC", "SCLC"), ylim=c(ymin, ymax), ylab="log2(TPM + 0.01)", main=paste0(gene, " (", ensembl_gene_id, ")"))
 
   dev.off()
}

#genes <- c("PSIP1", "PRIM1", "CD151", "PVRL2", "CDC7", "CCNE2", "MCM2", "CDKN2A", "BARD1", "MYCL", "MYC", "MYCN", "SYP")   ## Not "expressed": CHGA (in LUAD); NCAM1 (in LCNEC and LUAD)
#genes <- c("TMSB15A", "TNFRSF21", "HELLS", "ZNF519", )
#genes <- c("MURC", "CHGA", "NCAM1")   ## Not "expressed"
#genes <- c("CHEK1", "IDH1", "NOTCH2NL")
#genes <- c("TOP1", "TOPBP1", "POLA1", "POLE2", "POLD3", POLD4", "NMNAT1", "RB1", "CDK1", "CDK2", "CDK4", "CDK6", "EXO1")
genes <- c("CHEK1", "CHEK2", "RAD51", "RAD51AP1", "BLM", "MSH2", "EXO1", "TOP3A", "TOP3B")
genes <- c("HELLS", "SMARCA4", "ARID1A")
genes <- c("BCL2", "BCL2L11", "BCL3")
for (g in 1:length(genes)) {
   plotBox(genes[g], wd.de, tpm.gene.pcg.log2, pheno.all)
   #plotBeeswarm(genes[g], wd.de, tpm.gene.pcg.exp.log2, pheno.all)
}




x <- "ENSG00000100604"
t.test(as.numeric(expr.pheno.log2[x, rownames(subset(pheno.all, Cancer_Type == 0))]), as.numeric(expr.pheno.log2[x, rownames(subset(pheno.all, Cancer_Type == 1))]))$p.value


 
gene <- "MCM2"
ensembl_gene_id <- subset(ensGene.gene, external_gene_name == gene)$ensembl_gene_id
gene.tpms <- cbind(t(expr.pheno.log2)[,ensembl_gene_id], pheno.all)
colnames(gene.tpms)[1] <- "LOG2_TPM"

pdf(paste0(wd.de, "plots/boxplot_kallisto_tpm.log2_", gene, ".pdf"), height=6, width=4)
ymin <- min(gene.tpms$LOG2_TPM)
ymax <- max(gene.tpms$LOG2_TPM)
boxplot(LOG2_TPM ~ Cancer_Type, data=gene.tpms, outline=T, names=c("LUAD", "LCNEC", "SCLC"), ylim=c(ymin, ymax), ylab="log2(TPM + 0.01)", main=paste0(gene, " (", ensembl_gene_id, ")"))

#beeswarm(LOG2_TPM ~ Cancer_Type, data=subset(gene.tpms, RB1=="NA"), col="lightgray", pch=16, add=T)
#beeswarm(LOG2_TPM ~ Cancer_Type, data=subset(gene.tpms, RB1=="WT"), col="blue", pch=16, add=T)
#beeswarm(LOG2_TPM ~ Cancer_Type, data=subset(gene.tpms, RB1=="RB1"), col="red", pch=16, add=T)

#gene.tpms[intersect(rownames(gene.tpms), unique(muts.gene.sil$PAT_ID)),]$RB1 <- "RB1 (silent)"
#gene.tpms["S00396",]$RB1 <- "RB1 (silent)"
#gene.tpms["S00396",]$LOG2_TPM <- -7
#gene.tpms["S00006",]$RB1 <- "RB1 (silent)"
#gene.tpms["S00006",]$LOG2_TPM <- -7
#beeswarm(LOG2_TPM ~ Cancer_Type, data=subset(gene.tpms, RB1=="RB1 (silent)"), col="yellow", pch=16, add=T)

legend("bottomright", legend = c("RB1", "RB1 (silent)", "WT"), pch=16, col=c("red", "yellow", "blue"))
dev.off()

dev.off()





# -----------------------------------------------------------------------------
# Permutations on RB1 status (WRONG!!)
# Last Modified: 22/08/17
# -----------------------------------------------------------------------------
colnames <- c(1:1000)
permutations <- toTable(0, length(colnames), nrow(expr.pheno.log2), colnames)
rownames(permutations) <- rownames(expr.pheno.log2)

for (p in 1:1000) {
 idx <- sample(1:nrow(pheno.all), nrow(pheno.all), replace=F)
 pheno.all.perm <- pheno.all[idx,]
 
 permutations[,p] <- mapply(x = 1:nrow(expr.pheno.log2), function(x) cor.test(as.numeric(expr.pheno.log2[x,]), pheno.all.perm$RB1_RATE, method="spearman", exact=F)[[3]])
}