#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
START <- as.numeric(args[1])
END   <- as.numeric(args[2])

# =============================================================================
# Name: 1b_cmd-rt_bam.rpkm.corr.gc.R (commandline mode)
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 11/06/17
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"         ## tyang2@cheops
#wd.src <- "/re/home/tyang2/dev/R"                ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"             ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "ReplicationTiming.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))

# -----------------------------------------------------------------------------
# Calculate absolute RPKM
# Last Modified: 14/05/17
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"   ## tyang2@cheops
BASE <- "ICGC"
base <- tolower(BASE)

wd.anlys  <- file.path(wd, BASE, "analysis")
wd.driver <- file.path(wd.anlys, "driver", paste0(base, "-driver"))
wd.driver.data  <- file.path(wd.driver, "data", "copy_number_alterations")
wd.driver.plots <- file.path(wd.driver, "plots")

wd.meta     <- file.path(wd, BASE, "metadata")
samples.cna <- readTable(file.path(wd.meta, paste0("copy_number_alterations.txt")), header=F, rownames=F, sep="\t")

for (c in START:END) {
   sample <- samples.cna[c,]

   colnames <- c("ensembl_gene_id", "CN")
   ensGene.cna <- toTable("", length(colnames), 0, colnames)
   segs <- read.peiflyne.cn.seg(file.path("/projects/cangen/PCAWG-Repository/PCAWG.raw.data/copy_number/sclust_final_copy_number_analysis_files", sample$V2, paste0(sample$V2, "_iCN.seg")))   ## See ReplicationTiming.R / Inner Class / PeifLyne File Reader
   for (s in 1:nrow(segs)) {
      seg <- segs[s,]
      ensGene.chr.end.start <- getEnsGenesFromSegment(seg)
      
      if (nrow(ensGene.chr.end.start) > 0) {
         ens <- toTable("", length(colnames), nrow(ensGene.chr.end.start), colnames)
         ens$ensembl_gene_id <- ensGene.chr.end.start$ensembl_gene_id
         ens$CN              <- seg$CN
         
         ensGene.cna <- rbind(ensGene.cna, ens)
      }
   }
   
   rm(segs)
   rm(ensGene.chr.end.start)
   writeTable(ensGene.cna, gzfile(file.path(wd.driver.data, paste0(sample$V1, "_ens.cn.seg.gz"))), colnames=T, rownames=F, sep="\t")
}
