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
sv.del <- subset(sv, class1 == "Del")

load(file=file.path(wd.rt.data, paste0("icgc_wgs_sv.RData")))

# -----------------------------------------------------------------------------
# 
# Last Modified: 05/07/23
# -----------------------------------------------------------------------------
chr <- chrs[CHR]
sv.del <- subset(sv.del, seqnames == chr)

removed.chr <- c()
samples.chr <- c()
icgc.dels.chr <- NULL
for (d in 1:nrow(sv.del)) {
	  id <- rownames(subset(release.sv, tumor_wgs_icgc_sample_id == sv.del$icgc_sample_id[d]))
	  
	  if (length(id) != 0) {
	  	  icgc <- readTable(file.path(wd.icgc.sv, "final", paste0(id, ".pcawg_consensus_1.6.161116.somatic.sv.bedpe.gz")), header=T, rownames=F, sep="\t")
	  	  icgc.del <- subset(icgc, svclass == "DEL")
	  	  
	  	  if (nrow(icgc.del) > 0) {
	  		    samples.chr <- c(id, samples.chr)
	  		
	  		    sv_id <- unlist(strsplit(sv.del$breakpoint_id[d], ":"))[4]
	  		    icgc.del.id <- subset(icgc.del, sv_id == sv_id)
	  		    if (nrow(icgc.del.id) > 0)
	  		    	  if (is.null(icgc.dels.chr))
	  		    	  	  icgc.dels.chr <- icgc[0,]
	  		    	  else
	  			         icgc.dels.chr <- rbind(icgc.del.id, icgc.dels.chr)
	  	  }
	  } else {
	  	  removed.chr <- c(d, removed.chr)
	  }
}

save(removed.chr, samples.chr, icgc.dels.chr, file=file.path(wd.rt.data, paste0("icgc_wgs_sv_chr", CHR, ".RData")))