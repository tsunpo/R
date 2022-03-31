# =============================================================================
# Manuscript   : 
# Chapter      : Chromosome replication timing of the human genome
# Name         : manuscripts/asymmetries/cll-asym-rt.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 20/06/18
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "Asymmetry.R", "Mutation.R", "Survival.R")
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
wd.icgc.vcf <- file.path(wd.icgc, "consensus_snv_indel")

wd.meta     <- file.path(wd, BASE, "metadata", "data_release")

wd.anlys  <- file.path(wd, BASE, "analysis")
wd.driver <- file.path(wd.anlys, "driver", paste0(base, "-driver"))
wd.driver.data  <- file.path(wd.driver, "data")
wd.driver.plots <- file.path(wd.driver, "plots")

# -----------------------------------------------------------------------------
# Differential point mutations (SNV)
# Last Modified: 23/03/22
# -----------------------------------------------------------------------------
vcf.gene <- readTable(file.path(wd.icgc.vcf, "final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz"), header=T, rownames=F, sep="\t")

colnames <- unique(vcf.gene$t_alt_count)
rownames <- unique(vcf.gene$Hugo_Symbol)
mut.gene <- toTable(0, length(colnames), length(rownames), colnames)
rownames(mut.gene) <- rownames

for (s in 1:nrow(samples.mut)) {
   sample.mut <- samples.mut[s,]
 
   mut.dupl <- readTable(file.path(wd.driver.data, "point_mutations", paste0(sample.mut$icgc_specimen_id, "_mutcall_filtered_ens.vcf.gz")), header=T, rownames=F, sep="")
   mut <- as.data.frame(sort(table(mut.dupl$ensembl_gene_id), decreasing=T))
   rownames(mut) <- mut$Var1
 
   overlaps <- intersect(rownames, rownames(mut))
   mut.gene[overlaps, s] <- mut[overlaps,]$Freq
}

save(mut.gene, file=file.path(wd.driver.data, "icgc-consensus-mut.RData"))
