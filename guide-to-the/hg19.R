# =============================================================================
# Library      : The bioinformatician's guide to the reference genome
# Name         : guide-to-the/hg19.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 15/03/18
# =============================================================================
wd.reference <- "/Users/tpyang/Work/local/reference/hg19/"
wd.src.ref   <- "/Users/tpyang/Work/dev/R/guide-to-the/"

source("/Users/tpyang/Work/dev/R/handbook-of/Common.R")

# =============================================================================
# Reference    : Ensembl BioMart
# Link(s)      : http://grch37.ensembl.org/biomart/martview/
# Last Modified: 02/11/17
# =============================================================================
# -----------------------------------------------------------------------------
# Database: Ensembl Gene
# Table: ensGene / Count: 57,773 (out of 63,677; Unique results only) / Download Version: 2017-03-27
# -----------------------------------------------------------------------------
ensGene <- readTable(paste0(wd.reference, "ensembl/BioMart.GRCh37.p13.EnsemblGene.ensembl_gene_id_20170327.txt.gz"), header=T, rownames=F, sep="\t")
colnames(ensGene) <- c("ensembl_gene_id", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype", "external_gene_name")
rownames(ensGene) <- ensGene$ensembl_gene_id
# > nrow(ensGene)
# [1] 63677
# > nrow(subset(ensGene, gene_biotype == "protein_coding"))
# [1] 22810

## Remove patches (*_PATCH)
## https://www.ncbi.nlm.nih.gov/grc/help/patches
ensGene <- subset(ensGene, chromosome_name %in% c(1:22, "X", "Y", "MT"))   ## ADD 23/02/17
ensGene$chromosome_name <- paste0("chr", ensGene$chromosome_name)          ## ADD 27/04/17
# > nrow(ensGene)
# [1] 57773
# > nrow(subset(ensGene, gene_biotype == "protein_coding"))
# [1] 20327

## Add non-redundant protein coding gene list in handbook-of/DifferentialExpression.R from line 345

freq <- as.data.frame(table(ensGene$external_gene_name))   ## ADD 20/08/17
freq <- freq[order(freq$Freq, decreasing=T),]
# > nrow(subset(freq, Freq == 1))
# [1] 55539
# > subset(freq, Var1 == "RBL1")   ## NOTE 20/08/17
#       Var1 Freq
# 29322 RBL1    2

# -----------------------------------------------------------------------------
# Database: Ensembl Gene / Transcripts (with patches/scaffold sequences)
# Table: ensGene / Count: 215,170 / Download Version: 2017-03-25
# -----------------------------------------------------------------------------
ensGene.transcript <- readTable(paste0(wd.reference, "ensembl/BioMart.GRCh37.p13.EnsemblGene.ensembl_transcript_id_20170325.txt.gz"), header=T, rownames=F, sep="\t")
colnames(ensGene.transcript) <- c("ensembl_transcript_id", "ensembl_gene_id", "chromosome_name", "strand", "transcript_start", "transcript_end", "start_position", "end_position", "transcript_biotype", "gene_biotype", "external_gene_name")   ##"hugo_gene", "refseq_mrna", "refseq_mrna"
rownames(ensGene.transcript) <- ensGene.transcript$ensembl_transcript_id
#ensGene.transcript$chromosome_name <- paste0("chr", ensGene.transcript$chromosome_name)   ## REMOVED 15/03/18
ensGene.transcript <- ensGene.transcript[, -c(4,7,8,10)]   ## ADD 14/03/18
# > nrow(ensGene.transcript)
# [1] 215170
# > length(unique(ensGene.transcript$ensembl_gene_id))
# [1] 63677
# > length(unique(ensGene.transcript$external_gene_name))
# [1] 56638

## REMOVED 02/11/17: Keep all transcripts (including *_PATCH) to run "sleuth_prep" to avoid warning messages
#ensGene.transcript <- subset(ensGene.transcript, chromosome_name %in% c(1:22, "X", "Y", "MT"))   ## ADD 23/02/17

# =============================================================================
# Reference    : UCSC Genome Browser (Feb 2009/GRCh37/hg19)
# Link(s)      : http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/
# Last Modified: 15/03/18
# =============================================================================
# -----------------------------------------------------------------------------
# Database: Ensembl Gene / Transcripts / Exons
# Table: ensGene / Download Version: 06-Apr-2014
# -----------------------------------------------------------------------------
chrs <- paste0("chr", c(1:22, "X", "Y", "M"))

ensGene.tmp <- readTable(paste0(wd.reference, "ucsc/ensGene.txt.gz"), header=T, rownames=F, sep="")[, -1]   ## Get rid of first column "bin"
rownames(ensGene.tmp) <- ensGene.tmp$name
ensGene.tmp <- subset(ensGene.tmp, chrom %in% chrs)
# > nrow(ensGene.tmp)
# [1] 196354
# > length(unique(ensGene.tmp$name2))
# [1] 57773

overlaps <- intersect(ensGene.tmp$name, ensGene.transcript$ensembl_transcript_id)
length(overlaps)
# > length(overlaps)
# [1] 204940

colnames <- c("name", "name2", "chrom", "strand", "exonStart", "exonEnd", "cdsStartStat", "cdsEndStat", "exonFrame")
ensGene.transcript.exon <- toTable(0, length(colnames), 0, colnames)
for (t in 1:length(overlaps)) {
   transcript <- ensGene.tmp[overlaps[t],]
   exonStarts <- as.numeric(unlist(strsplit(transcript$exonStarts, ",")))
   exonEnds   <- as.numeric(unlist(strsplit(transcript$exonEnds, ",")))
   exonFrames <- as.numeric(unlist(strsplit(transcript$exonFrames, ",")))
   
   exons <- toTable(0, length(colnames), length(exonStarts), colnames)
   exons$name   <- transcript$name
   exons$name2  <- transcript$name2
   exons$chrom  <- transcript$chrom
   exons$strand <- transcript$strand
   exons$cdsStartStat <- transcript$cdsStartStat
   exons$cdsEndStat   <- transcript$cdsEndStat
   for (e in 1:length(exonStarts))
      exons[e, c(5,6,9)] <- c(exonStarts[e], exonEnds[e], exonFrames[e])
   
   ensGene.transcript.exon <- rbind(ensGene.transcript.exon, exons)
}
ensGene.transcript.exon$exonStart <- ensGene.transcript.exon$exonStart + 1   ## Convert 0-based start coordinates (in 0-start 1-end UCSC format) to 1-based start coordinates (Ensembl)
                                                                             ## https://www.biostars.org/p/6131/
# > nrow(ensGene.transcript.exon)
# [1] 1195408
# > nrow(subset(ensGene.transcript.exon, exonFrame != -1))   ## Note that an exonFrames value of -1 means that the exon is entirely UTR
# [1] 724012                                                 ## https://groups.google.com/a/soe.ucsc.edu/forum/#!topic/genome/U-w4b_ZS2j0

#ensGene.tmp <- subset(subset(ensGene.tmp, cdsStartStat == "cmpl"), cdsEndStat == "cmpl")   ## Only keep transcripts with "complete" CDS start and end information (like in RefGene)
                                                                                            ## https://groups.google.com/a/soe.ucsc.edu/forum/#!topic/genome/Uz5ozC9vkCQ
# -----------------------------------------------------------------------------
# File: cytoBand.txt.gz / Download Version: 14-Jun-2009
# -----------------------------------------------------------------------------
cytoBand <- readTable(paste0(wd.reference, "ucsc/cytoBand.txt.gz"), header=F, rownames=F, sep="")
colnames(cytoBand) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")
# > nrow(cytoBand)
# [1] 862

cytoBand$name2 <- paste0(gsub("chr", "", cytoBand$chrom), cytoBand$name)
cytoBand$size <- cytoBand$chromEnd - cytoBand$chromStart

# -----------------------------------------------------------------------------
# File: chromInfo.txt.gz / Download Version: 27-Apr-2009
# -----------------------------------------------------------------------------
chromInfo <- readTable(paste0(wd.reference, "ucsc/chromInfo.txt.gz"), header=F, rownames=F, sep="")
colnames(chromInfo) <- c("chrom", "size")
rownames(chromInfo) <- chromInfo$chrom
chromInfo <- subset(chromInfo, chrom %in% chrs)[, 1:2]

###
## The Bioinformatician's Guide to the Human Genome
save(ensGene, ensGene.transcript, ensGene.transcript.exon, cytoBand, chromInfo, chrs, file=paste0(wd.src.ref, "hg19.RData"))

# =============================================================================
# Reference    : Miscellaneous
# Last Modified: 13/05/17
# =============================================================================
# -----------------------------------------------------------------------------
# File: human-genome.1kb-grid.bed
# Link(s): https://github.com/andrej-fischer/cloneHD
#          ftp://ftp.sanger.ac.uk/pub/teams/153/bed/
# -----------------------------------------------------------------------------
#bed.raw <- readTable(paste0(wd.reference, "collections/human-genome.1kb-grid.bed.gz"), header=F, rownames=F, sep="")
#colnames(bed.raw) <- c("CHR", "START", "END")
#bed.raw$BED <- paste0("P", rownames(bed.raw))
#rownames(bed.raw) <- bed.raw$BED
bed <- readTable(paste0(wd.ngs, "coverage/S00022/S00022_ANALYSIS/S00022_cn.txt"), header=F, rownames=T, sep="")[,c(1:4,10)]
colnames(bed) <- c("BED", "CHR", "START", "END", "GC")
rownames(bed) <- bed$BED
bed <- bed[,-1]

save(bed, file=paste0(wd.src.ref, "hg19.1kb.gc.RData"))
