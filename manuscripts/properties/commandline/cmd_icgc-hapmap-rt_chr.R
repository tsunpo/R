#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
CHR  <- as.numeric(args[1])

# =============================================================================
# Manuscript   : 
# Chapter      : 
# Name         : manuscripts/driver/icgc-driver-cna.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 22/03/22
# =============================================================================
wd.src <- "/nfs/users/nfs_t/ty2/dev/R"            ## @nfs
#wd.src <- "/Users/ty2/Work/dev/R"                ## @localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "GenomicProperty.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.bed.gc.1kb.RData"))

# -----------------------------------------------------------------------------
# Set up working directory
# Last Modified: 22/03/22
# -----------------------------------------------------------------------------
wd <- "/lustre/scratch127/casm/team294rr/ty2"   ## @lustre
#wd <- "/Users/ty2/Work/sanger/ty2"             ## @localhost
BASE <- "ICGC"
base <- tolower(BASE)

wd.meta <- file.path(wd, BASE, "metadata/genetic_map_HapMapII_GRCh37")

wd.anlys  <- file.path(wd, BASE, "analysis")
wd.rt <- file.path(wd.anlys, "properties", paste0(base, "-hapmap-rt"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

# -----------------------------------------------------------------------------
# 
# Last Modified: 03/07/23
# -----------------------------------------------------------------------------
chr <- chrs[CHR]
bed.gc.chr <- subset(bed.gc, CHR == chr)
bed.gc.chr$BED <- rownames(bed.gc.chr)

map.chr <- readTable(file.path(wd.meta, paste0("genetic_map_GRCh37_", chr, ".txt")), header=T, rownames=F, sep="")
colnames(map.chr) <- c("CHR", "START", "RATE", "MAP")
map.chr$BED <- mapply(x = 1:nrow(map.chr), function(x) getGenomicProperty(map.chr$CHR[x], map.chr$START[x], bed.gc.chr))

save(map.chr, file=file.path(wd.rt.data, paste0("genetic_map_GRCh37_", chr, ".RData")))
