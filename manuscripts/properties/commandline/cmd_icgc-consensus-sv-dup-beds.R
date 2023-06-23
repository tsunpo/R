#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
h  <- as.numeric(args[1])
rt <- as.numeric(args[2])

# =============================================================================
# Name         : manuscripts/cmd_icgc-consensus-sv-del.R
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 25/04/23
# =============================================================================
wd.src <- "/nfs/users/nfs_t/ty2/dev/R"           ## @nfs
#wd.src <- "/Users/ty2/Work/dev/R"                 ## @localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R","GenomicProperty.R", "ReplicationForkDirectionality.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.bed.gc.1kb.RData"))

# -----------------------------------------------------------------------------
# Set up working directory
# Last Modified: 22/03/22
# -----------------------------------------------------------------------------
wd <- "/lustre/scratch127/casm/team294rr/ty2/"   ## @lustre
#wd <- "/Users/ty2/Work/sanger/ty2"               ## @localhost
BASE <- "ICGC"
base <- tolower(BASE)

wd.icgc     <- file.path(wd, BASE, "consensus")
wd.icgc.sv  <- file.path(wd.icgc, "sv")
wd.icgc.sv.plots <- file.path(wd.icgc.sv, "dupl", "plots")

wd.meta     <- file.path(wd, BASE, "metadata")

wd.anlys  <- file.path(wd, BASE, "analysis")
wd.rt <- file.path(wd.anlys, "replication", paste0(base, "-sv-dup"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

load(file=file.path(wd.icgc.sv, "dup", "data", paste0("final.sv.dup.RData")))

kb <- 20
if (rt == 0) {
	  load(file=file.path(wd.rt.data, "bstrps", paste0("sclc-nl_rpkm.gc.cn.d.rt.log2s.nrfd.", kb, "kb_", "m2-m1", ".RData")))
} else if (rt == 1) {
	  load(file=file.path(wd.rt.data, "bstrps", paste0("sclc_rpkm.gc.cn.d.rt.log2s.nrfd.", kb, "kb_", "m2-m1", ".RData")))
} else if (rt == 2) {
	  load(file=file.path(wd.rt.data, "bstrps", paste0("nbl_rpkm.gc.cn.d.rt.log2s.nrfd.", kb, "kb_", "m2-m1", ".RData")))
} else if (rt == 3) {
	  load(file=file.path(wd.rt.data, "bstrps", paste0("cll_rpkm.gc.cn.d.rt.log2s.nrfd.", kb, "kb_", "m2-m1", ".RData")))
}

# -----------------------------------------------------------------------------
# RT vs. Del size
# Last Modified: 21/04/23
# -----------------------------------------------------------------------------
#getDELRTNRFD <- function(del, nrds.RT.NRFD, bed.gc) {
#	  chr <- del$chrom1
#	  bed.gc.chr <- subset(bed.gc, CHR == paste0("chr", chr))
#	  
#	  bed.gc.chr.start <- subset(bed.gc.chr, START <= del$end2)
#	  bed.gc.chr.start.end <- subset(bed.gc.chr.start, END >= del$start1)
#	  
#	  return(intersect(rownames(bed.gc.chr.start.end), rownames(nrds.RT.NRFD)))
#}

#for (h in 1:nrow(table.sv.del.20)) {
	  hist <- as.vector(table.sv.dup$Histology[h])
	  totals.sv.hist <- subset(totals.sv, histology_abbreviation == hist)
	
	  dups <- c()
	  for (s in 1:nrow(totals.sv.hist)) {
		    id2 <- totals.sv.hist$specimen_id[s]
		    dup <- readTable(file.path(wd.icgc.sv, "dup", "data", hist, paste0(id2, ".somatic.sv.dup.bedpe.gz")), header=T, rownames=F, sep="\t")
      
		    if (is.null(dups))
		    	  dups <- dup
		    else
		    	  dups <- rbind(dups, dup)
	  }
	  dups$size <- dups$end2 - dups$start1

	  dups$BED <- mapply(x = 1:nrow(dups), function(x) getRandomBreakpointBED(dups[x,], nrds.RT.NRFD, bed.gc))
	  dups <- dups[!is.na(dups$BED), ]
	  #dups <- cbind(dups, nrds.RT.NRFD[dups$BED, c("SPLINE", "RFD", "NRFD")])
	  
	  save(dups, file=file.path(wd.rt.data, "beds", paste0(hist, ".sv.dup.rt", rt, ".beds.RData")), version=2)
#}