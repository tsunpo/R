# =============================================================================
# Library      : The bioinformatician's guide to the human genome
# Name         : guide-to-the/hg38.R
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 02/11/23
# =============================================================================
source("/Users/ty2/Work/dev/R/handbook-of/Commons.R")

wd.reference <- "/Users/ty2/Work/local/reference/hg38"
wd.src.ref   <- "/Users/ty2/Work/dev/R/guide-to-the"

# =============================================================================
# Reference    : UCSC Genome Browser (Dec 2013/GRCh38/hg38)
# Link(s)      : http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/
# Last Modified: 02/11/23
# =============================================================================

# -----------------------------------------------------------------------------
# File: cytoBand.txt.gz / Download Version: 2022-10-28
# -----------------------------------------------------------------------------
cytoBand <- readTable(file.path(wd.reference, "ucsc/cytoBand.txt.gz"), header=F, rownames=F, sep="")
colnames(cytoBand) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")
cytoBand$name2 <- paste0(gsub("chr", "", cytoBand$chrom), cytoBand$name)
cytoBand$size <- cytoBand$chromEnd - cytoBand$chromStart
# > nrow(cytoBand)
# [1] 1549

# -----------------------------------------------------------------------------
# File: chromInfo.txt.gz / Download Version: 2022-10-26
# -----------------------------------------------------------------------------
chrs <- paste0("chr", c(1:22, "X", "Y", "M"))

chromInfo <- readTable(file.path(wd.reference, "ucsc/chromInfo.txt.gz"), header=F, rownames=F, sep="")
colnames(chromInfo) <- c("chrom", "size")
rownames(chromInfo) <- chromInfo$chrom
chromInfo <- subset(chromInfo, chrom %in% chrs)[, 1:2]
chromInfo <- chromInfo[chrs,]   ## ADD 02/11/23

## The Bioinformatician's Guide to the Human Genome
save(cytoBand, chromInfo, chrs, file=file.path(wd.src.ref, "hg38.RData"))
