# =============================================================================
# Manuscript   : 
# Chapter I    : RB1-loss differential gene expression in neuroendocrine tumours
# Name         : manuscripts/expression/lcnec-tpm-de-path.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 08/08/18
# =============================================================================
#wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Common.R", "DifferentialExpression.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Human Genome
load(file.path(wd.src.ref, "hg19.RData"))

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost

BASE <- "LCNEC"
base <- tolower(BASE)
wd.anlys <- file.path(wd, BASE, "analysis")

wd.de       <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")

# -----------------------------------------------------------------------------
# Replace Ensembl Gene IDs to gene name in Reactome results (Up-regulation)
# -----------------------------------------------------------------------------
wd.de.data.reactome <- file.path(wd.de.data, "pathway_q0.05_up")

list <- readTable(file.path(wd.de.data.reactome, "genes_rb1_q0.05_n117_up.txt"), header=F, rownames=F, sep="")
list <- ensGene[list, c("ensembl_gene_id",	"external_gene_name")]

reactome <- read.csv(file.path(wd.de.data.reactome, "result.csv"))
colnames(reactome) <- gsub("X.", "", colnames(reactome))
reactome$Submitted.entities.found <- as.vector(reactome$Submitted.entities.found)
for (r in 1:nrow(reactome)) {
   ids <- as.vector(reactome$Submitted.entities.found[r])
   ids <- unlist(strsplit(ids, ";"))
 
   for (i in 1:length(ids))
      if (nrow(list[ids[i],]) != 0)
         ids[i] <- list[ids[i],]$external_gene_name
 
   reactome$Submitted.entities.found[r] <- paste(ids, collapse=";")
}
writeTable(reactome, file.path(wd.de.data.reactome, "result.tsv"), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Replace Ensembl Gene IDs to gene name in Reactome results (Down-regulation)
# -----------------------------------------------------------------------------
wd.de.data.reactome <- file.path(wd.de.data, "pathway_q0.05_down")

list <- readTable(file.path(wd.de.data.reactome, "genes_rb1_q0.05_n28_down.txt"), header=F, rownames=F, sep="")
list <- ensGene[list, c("ensembl_gene_id",	"external_gene_name")]

reactome <- read.csv(file.path(wd.de.data.reactome, "result.csv"))
colnames(reactome) <- gsub("X.", "", colnames(reactome))
reactome$Submitted.entities.found <- as.vector(reactome$Submitted.entities.found)
for (r in 1:nrow(reactome)) {
   ids <- as.vector(reactome$Submitted.entities.found[r])
   ids <- unlist(strsplit(ids, ";"))
 
   for (i in 1:length(ids))
      if (nrow(list[ids[i],]) != 0)
      ids[i] <- list[ids[i],]$external_gene_name
 
   reactome$Submitted.entities.found[r] <- paste(ids, collapse=";")
}
writeTable(reactome, file.path(wd.de.data.reactome, "result.tsv"), colnames=T, rownames=F, sep="\t")
