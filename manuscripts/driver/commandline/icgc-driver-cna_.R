# =============================================================================
# Manuscript   : 
# Chapter      : 
# Name         : manuscripts/driver/icgc-driver-cna.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 22/03/22
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "Asymmetry.R", "Mutation.R", "Survival.R", "Transcription.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))

# -----------------------------------------------------------------------------
# Set up working directory
# Last Modified: 22/03/22
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
#wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE <- "ICGC"
base <- tolower(BASE)

wd.icgc     <- file.path(wd, BASE, "consensus")
wd.icgc.cna <- file.path(wd.icgc, "copy_number_alterations")

wd.meta     <- file.path(wd, BASE, "metadata")

wd.anlys  <- file.path(wd, BASE, "analysis")
wd.driver <- file.path(wd.anlys, "driver", paste0(base, "-driver"))
wd.driver.data  <- file.path(wd.driver, "data")
wd.driver.plots <- file.path(wd.driver, "plots")

# -----------------------------------------------------------------------------
# Load data
# Last Modified: 22/03/22
# -----------------------------------------------------------------------------
samples <- readTable(file.path(wd.meta, paste0("samples.txt")), header=T, rownames=T, sep="\t")[, -1]
cor.rho <- readTable(file.path(wd.meta, paste0("cor.txt")), header=F, rownames=F, sep="")

release <- readTable(file.path(wd.meta, "data_release", "release_may2016.v1.4.tsv"), header=T, rownames=F, sep="")
rownames(release) <- release$tumor_wgs_aliquot_id
nrow(release)
# [1] 2834

list <- strsplit0(readTable(file.path(wd.icgc.cna, "copy_number_alterations.list"), header=F, rownames=F, sep=""), "_mutcall_filtered.vcf", 1)
length(list)
# [1] 2777

overlaps <- intersect(rownames(release), list)
length(overlaps)
# [1] 2595
release <- release[overlaps,]
rownames(release) <- release$tumor_wgs_icgc_specimen_id

overlaps <- intersect(rownames(samples), rownames(release))
length(overlaps)
# [1] 2443
release     <- release[overlaps,]
samples.cna <- samples[overlaps,]
samples.cna$tumor_wgs_aliquot_id <- release$tumor_wgs_aliquot_id
samples.cna <- setProliferation(samples.cna, cor.rho)
nrow(samples.cna)
# [1] 2443
#writeTable(samples.cna[, c("icgc_specimen_id", "tumor_wgs_aliquot_id")], file.path(wd.meta, paste0("copy_number_alterations.txt")), colnames=F, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Assigning copy numbers for each Ensembl gene in each sample
# Last Modified: 22/03/22
# -----------------------------------------------------------------------------
## See cmd_icgc-driver-cna.R

# -----------------------------------------------------------------------------
# Differential copy number alterations (CNA)
# Last Modified: 22/03/22
# -----------------------------------------------------------------------------
colnames <- rownames(samples.cna)
rownames <- rownames(ensGene)
cna.gene <- toTable(NA, length(colnames), length(rownames), colnames)
rownames(cna.gene) <- rownames

for (s in 1:nrow(samples.cna)) {
   sample.cna <- samples.cna[s,]
 
   cna.dupl <- readTable(file.path(wd.driver.data, "copy_number_alterations", paste0(sample.cna$icgc_specimen_id, "_ens.iCN.seg.gz")), header=T, rownames=F, sep="")
   cna <- as.data.frame(sort(table(cna.dupl$ensembl_gene_id), decreasing=T))
   rownames(cna) <- cna$Var1
 
   cna.1 <- subset(cna, Freq == 1)
   cna.1$CN <- cna.dupl[which(cna.dupl$ensembl_gene_id %in% rownames(cna.1)),]$CN
 
   cna.2 <- subset(cna, Freq > 1)
   cna.2$CN <- NA
   for (c in 1:nrow(cna.2)) {
      ens <- rownames(cna.2[c,])
      cna.2$CN[c] <- mean(subset(cna.dupl, ensembl_gene_id == ens)$CN)
   }
   cna <- rbind(cna.1[,c(1,3)], cna.2[,c(1,3)])
   colnames(cna) <- c("ensembl_gene_id", "CN")
 
   overlaps <- intersect(rownames, cna$ensembl_gene_id)
   cna.gene[overlaps, s] <- cna[overlaps,]$CN
}

cna.gene.nona <- cna.gene[removeMissing(cna.gene),]
save(samples.cna, cna.gene, cna.gene.nona, file=file.path(wd.driver.data, "icgc-driver-cna_iCN.seg.RData"))

# -----------------------------------------------------------------------------
# Wilcoxon rank sum test
# Last Modified: 22/03/22
# -----------------------------------------------------------------------------
## Test: Wilcoxon/Mann–Whitney/U/wilcox.test
## FDR : Q
argv      <- data.frame(predictor="SG1", predictor.wt="G1", test="Wilcoxon", test.fdr="Q", stringsAsFactors=F)
file.name <- paste0("de_", base, "_cna-gene_SG1_wilcox_q_n2443_iCN.seg")

de <- differentialAnalysis(cna.gene.nona, samples.cna, argv$predictor, argv$predictor.wt, argv$test, argv$test.fdr)
# Samples with MUT SG1: 1536
# Samples with  WT SG1: 907

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.cna.gene <- cbind(annot[rownames(de),], de)

save(de.cna.gene,  file=file.path(wd.driver.data, paste0(file.name, ".RData")))
writeTable(de.cna.gene, file.path(wd.driver.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")
nrow(de.cna.gene)
# [1] 18523
