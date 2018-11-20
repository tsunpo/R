# =============================================================================
# Library      : Phase-specific gene sets of cell cycle from multiple studies
# Name         : guide-to-the/cycle.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 07/09/18
# =============================================================================
source("/Users/tpyang/Work/dev/R/handbook-of/Common.R")

wd.src.ref <- "/Users/tpyang/Work/dev/R/guide-to-the"
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"

# -----------------------------------------------------------------------------
# Tirosh et al., Nature 2016
# -----------------------------------------------------------------------------
wd.meta  <- file.path(wd, "HeLa", "metadata/Tirosh 2016")
core.G1S <- readTable(file.path(wd.meta, "Tirosh_G1-S.list"), header=F, rownames=F, sep="")
core.G2M <- readTable(file.path(wd.meta, "Tirosh_G2-M.list"), header=F, rownames=F, sep="")
core.Stemness <- readTable(file.path(wd.meta, "Tirosh_Stemness.list"), header=F, rownames=F, sep="")

## Core G1-S genes
core.G1S <- gsub("MLF1IP", "CENPU", core.G1S)
core.G1S <- subset(ensGene, external_gene_name %in% core.G1S)$ensembl_gene_id             ## 43/43

## Core G2-M genes
core.G2M <- unique(core.G2M)   ## HJURP is duplicates as it has two ensembl gene IDs
core.G2M <- subset(ensGene, external_gene_name %in% core.G2M)$ensembl_gene_id             ## 54/55

## Stemness genes
core.Stemness <- subset(ensGene, external_gene_name %in% core.Stemness)$ensembl_gene_id   ## 58/63

# -----------------------------------------------------------------------------
# Dominguez et al., Cell Res. 2016
# -----------------------------------------------------------------------------
wd.meta  <- file.path(wd, "HeLa", "metadata/Dominguez 2016")
periodic.G1S <- readTable(file.path(wd.meta, "Dominguez_G1-S.list"), header=F, rownames=F, sep="")   ## 304
periodic.G2M <- readTable(file.path(wd.meta, "Dominguez_G2-M.list"), header=F, rownames=F, sep="")   ## 876

save(core.G1S, core.G2M, core.Stemness, periodic.G1S, periodic.G2M, file=file.path(wd.src.ref, "cycle.RData"))
