#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
CHR   <- as.numeric(args[1])

# =============================================================================
# Name: cmd-rt_2a_nrd.gc.cn.d_chr.R (commandline mode)
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 23/04/19; 22/10/18
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"         ## tyang2@cheops
#wd.src <- "/re/home/tyang2/dev/R"                ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"             ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "ReplicationTiming.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.bed.gc.1kb.RData"))

# -----------------------------------------------------------------------------
# Add BED information for the ensGenes
# Last Modified: 20/10/19
# -----------------------------------------------------------------------------
getPosTSS <- function(ensGene.g) {
   if (ensGene.g$strand > 0)
      return(ensGene.g$start_position)
   else
      return(ensGene.g$end_position)
}

getPosTTS <- function(ensGene.g) {
   if (ensGene.g$strand > 0)
      return(ensGene.g$end_position)
   else
      return(ensGene.g$start_position)
}

getEnsGeneBED <- function(pos, bed.gc.chr) {
   bed.gc.chr.start <- subset(bed.gc.chr, pos >= START)
   bed.gc.chr.start.end <- subset(bed.gc.chr.start, pos <= END)    ## ADD "=" 19/11/18
 
   if (nrow(bed.gc.chr.start.end) != 0)
      return(rownames(bed.gc.chr.start.end))
   else
      return("")
}

chr <- paste0("chr", CHR)
ensGene.chr <- subset(ensGene, chromosome_name == chr)
bed.gc.chr <- subset(bed.gc, CHR == chr)

ensGene.bed.chr <- toTable(NA, 4, nrow(ensGene.chr), c("TS", "TT", "TSS", "TTS"))
rownames(ensGene.bed.chr) <- rownames(ensGene.chr)
ensGene.bed.chr$TS <- mapply(x = 1:nrow(ensGene.chr), function(x) getPosTSS(ensGene.chr[x,]))
ensGene.bed.chr$TT <- mapply(x = 1:nrow(ensGene.chr), function(x) getPosTTS(ensGene.chr[x,]))
ensGene.bed.chr$TSS <- mapply(x = 1:nrow(ensGene.chr), function(x) getEnsGeneBED(ensGene.bed.chr$TS[x], bed.gc.chr))
ensGene.bed.chr$TTS <- mapply(x = 1:nrow(ensGene.chr), function(x) getEnsGeneBED(ensGene.bed.chr$TT[x], bed.gc.chr))

save(ensGene.bed.chr, file=file.path(wd.src.ref, paste0("hg19.ensGene.bed.1kb.", chr, ".RData")))
