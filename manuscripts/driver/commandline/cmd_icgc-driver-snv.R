#!/usr/bin/env Rscript
args  <- commandArgs(TRUE)
START <- as.numeric(args[1])
END   <- as.numeric(args[2])

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
# Set working directory
# Last Modified: 22/03/22
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
#wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE <- "ICGC"
base <- tolower(BASE)

wd.icgc     <- file.path(wd, BASE, "consensus")
wd.icgc.vcf <- file.path(wd.icgc, "point_mutations", "Converted_Data")

wd.meta     <- file.path(wd, BASE, "metadata")

wd.anlys  <- file.path(wd, BASE, "analysis")
wd.driver <- file.path(wd.anlys, "driver", paste0(base, "-driver"))
wd.driver.data  <- file.path(wd.driver, "data")
wd.driver.plots <- file.path(wd.driver, "plots")

wd.rt <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

load(file=file.path(wd.rt.data, paste0("icgc_wgs.RData")))
load(file=file.path(wd.rt.data, paste0("icgc_wgs_samples_n2612.RData")))

# -----------------------------------------------------------------------------
# Load data
# Last Modified: 22/03/22
# -----------------------------------------------------------------------------
list <- strsplit0(readTable(file.path(wd.icgc.vcf, "../point_mutations.list"), header=F, rownames=F, sep=""), "_mutcall_filtered.vcf", 1)
length(list)
# [1] 2703

#overlaps <- intersect(rownames(release), list)
#length(overlaps)
# [1] 2521
overlaps <- intersect(totals$wgs_id, list)
length(overlaps)
# [1] 2674

#nrow(samples.mut)
# [1] 2373
overlaps <- intersect(rownames(samples), totals[overlaps,]$specimen_id)
length(overlaps)
# [1] 2542
samples.vcf <- samples[overlaps,]
rownames(totals) <- totals$specimen_id
samples.vcf$tumor_wgs_aliquot_id <- totals[overlaps,]$wgs_id

# -----------------------------------------------------------------------------
# Step 1: Finding mutations locate within Ensembl genes
# Last Modified: 23/01/18
# -----------------------------------------------------------------------------
for (s in START:END) {
   sample.vcf <- samples.vcf[s,]
 
   vcf <- read.peiflyne.mutcall.filtered.vcf(file.path(wd.icgc.vcf, paste0(sample.vcf$tumor_wgs_aliquot_id, "_mutcall_filtered.vcf.gz")), pass=T, rs=F)
   vcf.gene <- getSNVinEnsGene(vcf, ensGene)
 
   writeTable(vcf.gene, gzfile(file.path(wd.driver.data, "point_mutations", paste0(sample.vcf$icgc_specimen_id, "_mutcall_filtered_ens.vcf.gz"))), colnames=T, rownames=F, sep="\t")
}
