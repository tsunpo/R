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
wd.icgc.sv.plots <- file.path(wd.icgc.sv, "del", "plots")

wd.meta     <- file.path(wd, BASE, "metadata")

wd.anlys  <- file.path(wd, BASE, "analysis")
wd.rt <- file.path(wd.anlys, "replication", paste0(base, "-sv-del"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

load(file=file.path(wd.icgc.sv, "del", "data", paste0("final.sv.del.RData")))

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
	  hist <- as.vector(table.sv.del.20$Histology[h])
	  totals.sv.hist <- subset(totals.sv, histology_abbreviation == hist)
	
	  dels <- c()
	  for (s in 1:nrow(totals.sv.hist)) {
		    id2 <- totals.sv.hist$specimen_id[s]
		    del <- readTable(file.path(wd.icgc.sv, "del", "data", hist, paste0(id2, ".somatic.sv.del.bedpe.gz")), header=T, rownames=F, sep="\t")
      
		    if (is.null(dels))
		    	  dels <- del
		    else
		       dels <- rbind(dels, del)
	  }
	  dels$size <- dels$end2 - dels$start1

	  dels$BED <- mapply(x = 1:nrow(dels), function(x) getRandomBreakpointBED(dels[x,], nrds.RT.NRFD, bed.gc))
	  dels <- dels[!is.na(dels$BED), ]
	  dels <- cbind(dels, nrds.RT.NRFD[dels$BED, c("SPLINE", "RFD", "NRFD")]
	  
	  save(pixels, file=file.path(wd.icgc.sv, "del", "data", "beds", paste0(hist, ".sv.del.rt", rt, ".beds.RData")), version=2)
#}