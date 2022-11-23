# =============================================================================
# Manuscript   :
# Name         : manuscripts/expression/esad-tpm.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 22/08/20
# =============================================================================
wd.src <- "/ngs/cangen/tyang2/dev/R"              ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"             ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Commons.R", "Transcription.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Human Genome
load(file.path(wd.src.ref, "hg19.RData"))

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
BASE <- "ESAD2"
base <- tolower(BASE)

#wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
wd.rna   <- file.path(wd, BASE, "ngs/3RNA")
wd.rna.raw <- file.path(wd.rna, "kallisto_hg19.ensembl_quant-b100--single--bias")
wd.meta  <- file.path(wd, BASE, "metadata")

wd.anlys <- file.path(wd, BASE, "analysis")
wd.de    <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

samples <- readTable(file.path(wd.rna, "Sample_IDs_batch2_aim2_2021_final_RNA.txt"), header=T, rownames=3, sep="\t")
f1 <- rownames(subset(samples, Cheops_Folder_1 != ""))
f2 <- rownames(subset(samples, Cheops_Folder_2 != ""))
samples$Cheops_Folder <- ""
samples[f1,]$Cheops_Folder <- samples[f1,]$Cheops_Folder_1
samples[f2,]$Cheops_Folder <- samples[f2,]$Cheops_Folder_2
samples <- subset(samples, Cheops_Folder != "")
writeTable(samples, file.path(wd.rna, "esad2_3rna_n97.txt"), colnames=F, rownames=F, sep="\t")

#fastqs <- readTable(file.path(wd.rna, "esad2_3rna_n101.list"), header=T, rownames=T, sep="")
fastqs <- readTable(file.path(wd.rna, "esad2_3rna_n101.txt"), header=T, rownames=F, sep="\t")
fastqs$ID <- mapply(x = 1:nrow(fastqs), function(x) paste(unlist(strsplit(fastqs[x,]$FILE, "_"))[2:3], collapse="_"))
rownames(fastqs) <- fastqs$ID
overlaps <- intersect(rownames(fastqs), rownames(samples))
writeTable(fastqs[overlaps,], file.path(wd.rna, "esad2_3rna_n97.list"), colnames=F, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# From transcript-level estimates to gene-level TPMs using sleuth (v0.29.0)
# Based on https://pachterlab.github.io/sleuth_walkthroughs/boj/analysis.html
#
# With sleuth's default quality-control filters: minimum 5 reads in at least 47% of the samples
# https://pachterlab.github.io/sleuth/docs/basic_filter.html
# https://groups.google.com/forum/#!topic/kallisto-sleuth-users/QrKxxEEFnE0
# -----------------------------------------------------------------------------
library("sleuth")   ## R version 3.3.2 (on MacBook)

tsv <- file.path(wd.rna.raw, overlaps)
s2c <- data.frame(path=tsv, sample=, stringsAsFactors=F)
t2g <- tx2Ens(ensGene.transcript)   ## Use full Ensembl transcripts with patches/scaffold sequences (*_PATCH) to avoid warning messages
                                    ## Please refer to line 49 (in guide-to-the/hg19.R)
so <- sleuth_prep(s2c, target_mapping=t2g, gene_mode=T, aggregation_column="ens_gene", extra_bootstrap_summary=T, min_reads=5, min_prop=0.47)   ## Default quality-control filters
# reading in kallisto results
# dropping unused factor levels
# .....................................................................
# normalizing est_counts
# 23899 targets passed the filter
# normalizing tpm
# merging in metadata
# aggregating by column: ens_gene
# 13714 genes passed the filter
# summarizing bootstraps

## Transcript-level estimates with patches/scaffold sequences (*_PATCH)   ## See line 30 in guide-to-the/hg19.R
## https://www.ncbi.nlm.nih.gov/grc/help/patches
tpm.norm      <- kallisto_table(so, use_filtered=F, normalized=T, include_covariates=F)
tpm.norm.filt <- kallisto_table(so, use_filtered=T, normalized=T, include_covariates=F)
save(tpm.norm,      file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.norm.RData")))
save(tpm.norm.filt, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.norm.filt_r5p47.RData")))

## Gene-level TPMs (All)
tpm.gene <- getGeneTPM(list2Matrix(tpm.norm$tpm, tpm.norm), ensGene)             ## Gene-level TPMs (without filtering)
save(tpm.gene, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.RData")))
writeTable(tpm.gene, gzfile(file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.txt.gz"))), colnames=T, rownames=T, sep="\t")
nrow(tpm.gene)
# [1] 34908

## Gene-level TPMs (Expressed)
load(file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.RData")))
tpm.gene <- removeMedian0(tpm.gene)
save(tpm.gene, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
nrow(tpm.gene)
# [1] 15381

## Gene-level TPMs with default filters
tpm.gene <- getGeneTPM(list2Matrix(tpm.norm.filt$tpm, tpm.norm.filt), ensGene)   ## Gene-level TPMs (with default filters)
save(tpm.gene, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.r5p47.RData")))
nrow(tpm.gene)
# [1] 13004

# =============================================================================
# Density plot and histograms (See DifferentialExpression.R)
# Figure(s)    : Figure S2 (A and B)
# Last Modified: 13/09/20; 29/05/20
# =============================================================================
## All genes
file.main <- file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene"))
load(paste0(file.main, ".RData"))
plotDensityHistogram(tpm.gene, file.main, "Total Ensembl")

## Expressed genes
file.main <- file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.median0"))
load(paste0(file.main, ".RData"))
plotDensityHistogram(tpm.gene, file.main, "Expressed Ensembl")

## Consistently expressed genes
file.main <- file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.r5p47"))
load(paste0(file.main, ".RData"))
plotDensityHistogram(tpm.gene, file.main, "Filtered")

# =============================================================================
# RefSeq IDs
# Last Modified: 04/05/22
# =============================================================================
samples.wes <- readTable(file.path(wd.meta, "Patienten_Follow_Up_joined_final_3.0_tyang2.txt"), header=T, rownames=T, sep="\t")
samples.wes <- subset(samples.wes, IS_WES_QC == F)
samples.wes <- subset(samples.wes, IS_RNA_SEQ == T)
rownames(samples.wes) <- paste0(samples.wes$Patient_ID, "_B")
nrow(samples.wes)
# [1] 64
# [1] 49
# [1] 48

samples.aim1 <- samples.aim1[colnames(tpm.gene.1),]
colnames(tpm.gene.1) <- samples.aim1$PATIENT_ID2

overlaps.1 <- intersect(rownames(samples.wes), colnames(tpm.gene.1))
overlaps.2 <- intersect(rownames(samples.wes), colnames(tpm.gene.2))

overlaps <- intersect(rownames(tpm.gene.1), rownames(tpm.gene.2))
tpm.gene <- cbind(tpm.gene.1[overlaps, overlaps.1], tpm.gene.2[overlaps, overlaps.2])
save(tpm.gene, file=file.path(wd.de.data, paste0("esad1+2", "_kallisto_0.43.1_tpm.gene.RData")))
file.main <- file.path(wd.de.data, paste0("esad1+2", "_kallisto_0.43.1_tpm.gene_34908x64"))
writeTable(tpm.gene, gzfile(paste0(file.main, ".txt.gz")), colnames=T, rownames=T, sep="\t")

##
tpm.gene.median0 <- removeMedian0(tpm.gene, 0)
#save(tpm.gene.median0, file=file.path(wd.de.data, paste0("esad1+2", "_kallisto_0.43.1_tpm.gene.median0.RData")))
writeTable(tpm.gene.median0, gzfile(file.path(wd.de.data, paste0("esad1+2", "_kallisto_0.43.1_tpm.gene.median0_15911x48.txt.gz"))), colnames=T, rownames=T, sep="\t")
nrow(tpm.gene.median0)
# [1] 15911

tpm.gene.median1 <- removeMedian0(tpm.gene, 1)
nrow(tpm.gene.median1)
# [1] 13850

# =============================================================================
# Density plot and histograms (See DifferentialExpression.R)
# Figure(s)    : Figure S2 (A and B)
# Last Modified: 04/05/22; 06/09/20; 29/05/20
# =============================================================================
BASE <- "EAC 1+2"

## All genes
file.main <- file.path(wd.de.data, paste0("esad1+2", "_kallisto_0.43.1_tpm.gene"))
#load(paste0(file.main, ".RData"))
plotDensityHistogram(tpm.gene, file.main, "Total Ensembl")

## Expressed genes
file.main <- file.path(wd.de.data, paste0("esad1+2", "_kallisto_0.43.1_tpm.gene.median0"))
#load(paste0(file.main, ".RData"))
plotDensityHistogram(tpm.gene.median0, file.main, "Expressed Ensembl")

# =============================================================================
# Reference    : HGNC (hg19)
# Last Modified: 20/06/22; 04/05/22
# =============================================================================
## All genes
overlaps <- intersect(rownames(tpm.gene), rownames(hgnc.ref.ens))
length(overlaps)
# [1] 17924

file.main <- file.path(wd.de.data, paste0("esad1+2", "_kallisto_0.43.1_tpm.gene_17924x48"))
tpm.gene.refgene <- tpm.gene[overlaps,]
rownames <- rownames(tpm.gene.refgene)
rownames(tpm.gene.refgene) <- hgnc.ref.ens[rownames,]$symbol
writeTable(tpm.gene.refgene, gzfile(paste0(file.main, ".txt.gz")), colnames=T, rownames=T, sep="\t")

file.main <- file.path(wd.de.data, paste0("esad1+2", "_kallisto_0.43.1_tpm.gene_17924x48"))
#load(paste0(file.main, ".RData"))
plotDensityHistogram(tpm.gene.refgene, file.main, "Total Ensembl / RefGene")

## Expressed genes
overlaps <- intersect(rownames(tpm.gene.median0), rownames(hgnc.ref.ens))
length(overlaps)
# [1] 13608

file.main <- file.path(wd.de.data, paste0("esad1+2", "_kallisto_0.43.1_tpm.gene.median0_13608x48"))
tpm.gene.median0.refgene <- tpm.gene.median0[overlaps,]
rownames <- rownames(tpm.gene.median0.refgene)
rownames(tpm.gene.median0.refgene) <- hgnc.ref.ens[rownames,]$symbol
writeTable(tpm.gene.median0.refgene, gzfile(paste0(file.main, ".txt.gz")), colnames=T, rownames=T, sep="\t")

file.main <- file.path(wd.de.data, paste0("esad1+2", "_kallisto_0.43.1_tpm.gene.median0_13608x48"))
#load(paste0(file.main, ".RData"))
plotDensityHistogram(tpm.gene.median0.refgene, file.main, "Expressed Ensembl / RefGene")

###
## Expressed genes (median TPM > 1)
overlaps <- intersect(rownames(tpm.gene.median1), rownames(hgnc.ref.ens))
length(overlaps)
# [1] 11992

file.main <- file.path(wd.de.data, paste0("esad1+2", "_kallisto_0.43.1_tpm.gene.median1_11992x48"))
tpm.gene.median1.refgene <- tpm.gene.median1[overlaps,]
rownames <- rownames(tpm.gene.median1.refgene)
rownames(tpm.gene.median1.refgene) <- hgnc.ref.ens[rownames,]$symbol
writeTable(tpm.gene.median1.refgene, gzfile(paste0(file.main, ".txt.gz")), colnames=T, rownames=T, sep="\t")










overlaps <- intersect(rownames(tpm.gene.median1), rownames(hgnc.ref.ens.vcf))
length(overlaps)

file.main <- file.path(wd.de.data, paste0("esad1+2", "_kallisto_0.43.1_tpm.gene.median1_mut=6209"))
tpm.gene.median1.vcf <- tpm.gene.median1[overlaps,]
rownames <- rownames(tpm.gene.median1.vcf)
rownames(tpm.gene.median1.vcf) <- hgnc.ref.ens.vcf[rownames,]$symbol
writeTable(tpm.gene.median1.vcf, gzfile(paste0(file.main, ".txt.gz")), colnames=T, rownames=T, sep="\t")

##
overlaps <- intersect(rownames(tpm.gene.median0), rownames(hgnc.ref.ens.vcf))
length(overlaps)

file.main <- file.path(wd.de.data, paste0("esad1+2", "_kallisto_0.43.1_tpm.gene.median0_mut=7185"))
tpm.gene.median0.vcf <- tpm.gene.median0[overlaps,]
rownames <- rownames(tpm.gene.median0.vcf)
rownames(tpm.gene.median0.vcf) <- hgnc.ref.ens.vcf[rownames,]$symbol
writeTable(tpm.gene.median0.vcf, gzfile(paste0(file.main, ".txt.gz")), colnames=T, rownames=T, sep="\t")

###
##
overlaps <- intersect(rownames(tpm.gene.median0), rownames(hgnc.ref.ens))
length(overlaps)
# [1] 13608

file.main <- file.path(wd.de.data, paste0("esad1+2", "_kallisto_0.43.1_tpm.gene.median0_refgene=13608"))
tpm.gene.median0.refgene <- tpm.gene.median0[overlaps,]
rownames <- rownames(tpm.gene.median0.refgene)
rownames(tpm.gene.median0.refgene) <- hgnc.ref.ens[rownames,]$symbol
writeTable(tpm.gene.median0.refgene, gzfile(paste0(file.main, ".txt.gz")), colnames=T, rownames=T, sep="\t")

###
## 20/06/22





## Expressed genes
file.main <- file.path(wd.de.data, paste0("esad1+2", "_kallisto_0.43.1_tpm.gene.median0"))
#load(paste0(file.main, ".RData"))
plotDensityHistogram(tpm.gene.refgene, file.main, "Total Ensembl & RefGene")





expressed.genes <- rownames(tpm.gene.median1)
ensGene.expressed <- ensGene[expressed.genes,]

hgnc   <- readTable(file.path(wd.reference, "hgnc/hgnc_complete_set_less.txt"), header=T, rownames=F, sep="\t")
refseq <- readTable(file.path(wd.reference, "hgnc/hg19u_refSeq_genes.txt"), header=F, rownames=F, sep="\t")
refseq.genes <- unique(refseq$V10)
refseq.genes.overlaps <- intersect(ensGene.expressed$external_gene_name, refseq.genes)
ensGene.expressed.overlaps <- subset(ensGene.expressed, external_gene_name %in% refseq.genes.overlaps)

hgnc <- subset(hgnc, ensembl_gene_id %in% ensGene$ensembl_gene_id)
