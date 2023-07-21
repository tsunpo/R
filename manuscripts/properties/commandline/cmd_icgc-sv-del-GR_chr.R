#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
CHR  <- as.numeric(args[1])
FILE <- args[2]
VALUE <- as.numeric(args[3])

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
handbooks  <- c("Commons.R", "GenomicProperty.R", "TranscriptionReplicationInteraction.R")
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

wd.icgc     <- file.path(wd, BASE, "consensus")
wd.icgc.sv  <- file.path(wd.icgc, "sv")
wd.icgc.sv.plots <- file.path(wd.icgc.sv, "del", "plots")

wd.nr3  <- file.path(wd, "dev/nr3/data/genome_properties")
wd.meta <- file.path(wd, BASE, "metadata")

wd.anlys  <- file.path(wd, BASE, "analysis")
wd.rt <- file.path(wd.anlys, "properties", paste0(base, "-sv-del"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

# -----------------------------------------------------------------------------
# 
# Last Modified: 03/07/23
# -----------------------------------------------------------------------------
sv <- readTable(file.path(wd.nr3, "results/2017_06_positions/svpos_with_hg19_props.csv"), header=T, rownames=F, sep=",")
colnames(sv)[c(1,2)] <- c("CHR", "START")
sv.del <- subset(sv, class1 == "Del")

fs <- data.frame(get(load(file.path(wd.nr3, "results/GRanges", paste0(FILE, "_GR.RData")))))
colnames(fs)[c(1:3,VALUE)] <- c("CHR", "START", "END", "VALUE")
fs$BED <- mapply(x = 1:nrow(fs), function(x) paste0("P", x))
rownames(fs) <- fs$BED

# -----------------------------------------------------------------------------
# 
# Last Modified: 05/07/23
# -----------------------------------------------------------------------------
chr <- chrs[CHR]
sv.del.chr <- subset(sv.del, CHR == chr)
sv.del.chr$SIZE <- NA
fs.chr <- subset(fs, CHR == chr)

removed <- c()
ids <- unique(sv.del.chr$breakpoint_id)
sv.del.chr.1 <- sv.del.chr[0,]
sv.del.chr.2 <- sv.del.chr[0,]
for (d in 1:length(ids)) {
	  bps <- subset(sv.del.chr, breakpoint_id == ids[d])
	  
	  if (nrow(bps) == 2) {
	  	  bps.1 <- bps[1,]
	  	  bps.2 <- bps[2,]
	  	  if (bps.1$START > bps.2$START) {
	  	  	  bps.1 <- bps[2,]
	  	  	  bps.2 <- bps[1,]
	  	  }
	  	  
	  	  bps.1$SIZE <- bps.2$START - bps.1$START
	  	  bps.2$SIZE <- bps.2$START - bps.1$START
	  	  
	  	  sv.del.chr.1 <- rbind(sv.del.chr.1, bps.1)
	  	  sv.del.chr.2 <- rbind(sv.del.chr.2, bps.2)
	  } else {
	  	  removed <- c(removed, ids[d])
	  }
}

## Breakpoint 1
sv.del.chr.1$BED <- mapply(x = 1:nrow(sv.del.chr.1), function(x) getGenomicProperty(sv.del.chr.1$CHR[x], sv.del.chr.1$START[x], fs))
sv.del.chr.1$VALUE <- fs[sv.del.chr.1$BED,]$VALUE
sv.del.chr.1 <- sv.del.chr.1[, -49]

## Breakpoint 2
sv.del.chr.2$BED <- mapply(x = 1:nrow(sv.del.chr.2), function(x) getGenomicProperty(sv.del.chr.2$CHR[x], sv.del.chr.2$START[x], fs))
sv.del.chr.2$VALUE <- fs[sv.del.chr.2$BED,]$VALUE
sv.del.chr.2 <- sv.del.chr.2[, -49]

save(removed, sv.del.chr.1, sv.del.chr.2, file=file.path(wd.rt.data, paste0("icgc_wgs_sv_del_", FILE, "_chr", CHR, ".RData")))
