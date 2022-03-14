# =============================================================================
# Manuscript   :
# Chapter      :
# Name         : manuscripts/expression/nbl-tpm-de.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 29/05/20
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"              ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"             ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Commons.R", "Transcription.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))

# -----------------------------------------------------------------------------
# https://apple.stackexchange.com/questions/254380/why-am-i-getting-an-invalid-active-developer-path-when-attempting-to-use-git-a
# https://github.com/kharchenkolab/conos/wiki/Installing-Conos-for-Mac-OS
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
#wd <- "/projects/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE <- "ICGC"
base <- tolower(BASE)

wd.rna <- file.path(wd, BASE, "ngs/RNA")
#wd.rna.raw <- file.path(wd.rna, "kallisto_hg19.ensembl_quant-b100--bias")

wd.anlys <- file.path(wd, BASE, BASE, "analysis")
wd.de    <- file.path(wd.anlys, "expression", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

t2g <- tx2Ens(ensGene.transcript)

# =============================================================================
# ICGC TPM calculated using Kallisto (version 0.42.1) and Gencode v19 protein coding transcripts
# Link(s)      : https://dcc.icgc.org/releases/PCAWG/transcriptome/transcript_expression
# Last Modified: 12/03/22
# =============================================================================
tpm <- readTable(file.path(wd.rna, "pcawg.rnaseq.transcript.expr.tpm.tsv.gz"), header=T, rownames=T, sep="")[,-1]
ids <- colnames(tpm)
ids <- gsub("\\.", "-", ids)
ids <- gsub("X", "", ids)
dim(tpm)
# [1] 95309  1359

aliquots <- readTable(file.path(wd.rna, "rnaseq.extended.metadata.aliquot_id.V4.txt"), header=T, rownames=T, sep="\t")
overlaps <- intersect(ids, rownames(aliquots))
aliquots <- aliquots[overlaps,]
colnames(tpm) <- aliquots$icgc_specimen_id

## Gene-level TPMs (All)
genes <- unique(t2g$ens_gene)
tpm.gene <- toTable(NA, ncol(tpm), 0, colnames(tpm))
for (g in 1:length(genes)) {
   transcripts <- subset(t2g, ens_gene== genes[g])
   
   tpm.transcript <- tpm[grep(paste(transcripts$target_id, collapse="|"), rownames(tpm)),]
   if (nrow(tpm.transcript) > 0) {
      tpm.g <- toTable(NA, ncol(tpm), 1, colnames(tpm))
      rownames(tpm.g) <- genes[g]
      
      tpm.g[1,] <- mapply(x = 1:ncol(tpm.transcript), function(x) sum(tpm.transcript[, x]))
      tpm.gene <- rbind(tpm.gene, tpm.g)
   }
}
save(tpm.gene, file=file.path(wd.de.data, paste0(base, "_kallisto_0.42.1_tpm.gene.RData")))
#writeTable(tpm.gene, gzfile(file.path(wd.de.data, paste0(base, "_kallisto_0.42.1_tpm.gene.txt.gz"))), colnames=T, rownames=T, sep="\t")
dim(tpm.gene)
# [1] 20720  1359

## Gene-level TPMs (Expressed)
#load(file=file.path(wd.de.data, paste0(base, "_kallisto_0.42.1_tpm.gene.RData")))
tpm.gene <- removeMedian0(tpm.gene)
save(tpm.gene, file=file.path(wd.de.data, paste0(base, "_kallisto_0.42.1_tpm.gene.median0.RData")))
#writeTable(tpm.gene, gzfile(file.path(wd.de.data, paste0(base, "_kallisto_0.42.1_tpm.gene.median0.txt.gz"))), colnames=T, rownames=T, sep="\t")
nrow(tpm.gene)
# [1] 18502

# =============================================================================
# Density plot and histograms (See DifferentialExpression.R)
# Figure(s)    : Figure S2 (A and B)
# Last Modified: 06/09/20; 29/05/20
# =============================================================================
## All genes
file.main <- file.path(wd.de.data, paste0(base, "_kallisto_0.42.1_tpm.gene"))
load(paste0(file.main, ".RData"))
plotDensityHistogram(tpm.gene, file.main, "Total Ensembl")

## Expressed genes
file.main <- file.path(wd.de.data, paste0(base, "_kallisto_0.42.1_tpm.gene.median0"))
load(paste0(file.main, ".RData"))
plotDensityHistogram(tpm.gene, file.main, "Expressed")

# =============================================================================
# ICGC FPKM aligned with TopHat2 and STAR aligners
# Link(s)      : https://dcc.icgc.org/releases/PCAWG/transcriptome/transcript_expression
# Last Modified: 12/03/22
# =============================================================================
fpkm <- readTable(file.path(wd.rna, "tophat_star_fpkm.v2_aliquot_gl.tsv.gz"), header=T, rownames=T, sep="")[,-1]
ids <- colnames(fpkm)
ids <- gsub("\\.", "-", ids)
ids <- gsub("X", "", ids)
colnames(fpkm) <- ids
dim(fpkm)
# [1] 57820  1521

aliquots <- readTable(file.path(wd.rna, "rnaseq.extended.metadata.aliquot_id.V4.txt"), header=T, rownames=T, sep="\t")
overlaps <- intersect(ids, rownames(aliquots))
fpkm.gene <- fpkm[, overlaps]
aliquots <- aliquots[overlaps,]
colnames(fpkm.gene) <- aliquots$icgc_specimen_id
test <- mapply(x = 1:nrow(fpkm.gene), function(x) unlist(strsplit(rownames(fpkm.gene)[x], "\\."))[1])
rownames(fpkm.gene) <- test
dim(fpkm.gene)
# [1] 57820  1359

###
##
overlaps <- intersect(rownames(tpm.gene), rownames(fpkm.gene))
tpm.gene.x <- tpm.gene[overlaps,]
x.median <- mapply(x = 1:nrow(tpm.gene.x), function(x) median(as.numeric(tpm.gene.x[x,])))
fpkm.gene.y <- fpkm.gene[overlaps,]
y.median <- mapply(x = 1:nrow(fpkm.gene.y), function(x) median(as.numeric(fpkm.gene.y[x,])))

file.name <- file.path(wd.de.data, paste0(base, "_FPKM-vs-TPM"))
plotCorrelation(file.name, "ICGC gene expression", "TPM", "FPKM", x.median, y.median, "topright", line=2.4)
file.name <- file.path(wd.de.data, paste0(base, "_FPKM-vs-TPM_log2+1"))
plotCorrelation(file.name, "ICGC gene expression", "log2(TPM + 1)", "log2(FPKM + 1)", log2(x.median+1), log2(y.median+1), "topright", line=2.4)
