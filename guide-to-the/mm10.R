# =============================================================================
# Library      : The Bioinformatician's Guide to the Mouse Genome
# Name         : guide-to-the/mm10.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 03/01/22
# =============================================================================
source("/Users/tpyang/Work/dev/R/handbook-of/Commons.R")

wd.reference <- "/Users/tpyang/Work/local/reference/mm10"
wd.src.ref   <- "/Users/tpyang/Work/dev/R/guide-to-the"

# =============================================================================
# Reference    : Ensembl BioMart
# Link(s)      : https://nov2020.archive.ensembl.org/biomart/martview/
# Last Modified: 03/01/22
# =============================================================================

# -----------------------------------------------------------------------------
# Database: Ensembl Gene
# Table: ensGene / Count: 55,184 (out of 56,305; Unique results only) / Download Version: 2022-01-03 (GRCm38.p6)
# -----------------------------------------------------------------------------
ensGene <- readTable(file.path(wd.reference, "ensembl/BioMart.GRCm38.p6.EnsemblGene.ensembl_gene_id_20220103.txt.gz"), header=T, rownames=F, sep="\t")
colnames(ensGene) <- c("ensembl_gene_id", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype", "external_gene_name")
rownames(ensGene) <- ensGene$ensembl_gene_id
# > nrow(ensGene)
# [1] 56305
# > nrow(subset(ensGene, gene_biotype == "protein_coding"))
# [1] 22287

## Remove patches (*_PATCH)
## https://www.ncbi.nlm.nih.gov/grc/help/patches
ensGene <- subset(ensGene, chromosome_name %in% c(1:19, "X", "Y", "MT"))   ## ADD 23/02/17
ensGene$chromosome_name <- paste0("chr", ensGene$chromosome_name)          ## ADD 27/04/17
# > nrow(ensGene)
# [1] 55401
# > nrow(subset(ensGene, gene_biotype == "protein_coding"))
# [1] 21859

freq <- as.data.frame(table(ensGene$external_gene_name))   ## ADD 20/08/17
freq <- freq[order(freq$Freq, decreasing=T),]
# > nrow(subset(freq, Freq == 1))
# [1] 55184

# -----------------------------------------------------------------------------
# Database: Ensembl Gene / Transcripts (with patches/scaffold sequences)
# Table: ensGene / Count: 215,170 / Download Version: 2017-03-25
# -----------------------------------------------------------------------------
ensGene.transcript <- readTable(file.path(wd.reference, "ensembl/BioMart.GRCm38.p6.EnsemblGene.ensembl_transcript_id_20220103.txt.gz"), header=T, rownames=F, sep="\t")
colnames(ensGene.transcript) <- c("ensembl_transcript_id", "ensembl_gene_id", "chromosome_name", "strand", "transcript_start", "transcript_end", "start_position", "end_position", "transcript_biotype", "gene_biotype", "external_gene_name")   ##"hugo_gene", "refseq_mrna", "refseq_mrna"
rownames(ensGene.transcript) <- ensGene.transcript$ensembl_transcript_id
#ensGene.transcript$chromosome_name <- paste0("chr", ensGene.transcript$chromosome_name)   ## REMOVED 15/03/18
ensGene.transcript <- ensGene.transcript[, -c(7,8,10)]   ## CHANGED 14/06/19 Now with strand information (column 4); ADD 14/03/18
# > nrow(ensGene.transcript)
# [1] 144778
# > length(unique(ensGene.transcript$ensembl_gene_id))
# [1] 56305
# > length(unique(ensGene.transcript$external_gene_name))
# [1] 55494

## REMOVED 02/11/17: Keep all transcripts (including *_PATCH) to run "sleuth_prep" to avoid warning messages
#ensGene.transcript <- subset(ensGene.transcript, chromosome_name %in% c(1:19, "X", "Y", "MT"))   ## ADD 23/02/17

# =============================================================================
# Reference    : UCSC Genome Browser (Dec 2011/GRCm38/mm10)
# Link(s)      : http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/
# Last Modified: 03/01/22
# =============================================================================

# -----------------------------------------------------------------------------
# Database: Ensembl Gene / Transcripts / Exons
# Table: ensGene / Download Version: 06-Apr-2014
# -----------------------------------------------------------------------------
chrs <- paste0("chr", c(1:19, "X", "Y", "M"))

## The Bioinformatician's Guide to the Mouse Genome
save(ensGene, ensGene.transcript, chrs, file=file.path(wd.src.ref, "mm10.RData"))
