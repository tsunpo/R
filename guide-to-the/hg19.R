# =============================================================================
# Library      : The bioinformatician's guide to the reference genome
# Name         : guide-to-the/hg19.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 02/11/17
# =============================================================================
wd.reference <- "/Users/tpyang/Work/local/reference/hg19/"
wd.guidetothe <- "/Users/tpyang/Work/dev/R/guide-to-the/"

source("/Users/tpyang/Work/dev/R/handbook-of/Common.R")

# =============================================================================
# Reference    : Ensembl BioMart
# Link(s)      : http://grch37.ensembl.org/biomart/martview/
#                ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.75.cdna.all.fa.gz
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
# [1] 20327

## Remove patches (*_PATCH)
## https://www.ncbi.nlm.nih.gov/grc/help/patches
ensGene <- subset(ensGene, chromosome_name %in% c(1:22, "X", "Y", "MT"))   ## ADD 23/02/17
ensGene$chromosome_name <- paste0("chr", ensGene$chromosome_name)          ## ADD 27/04/17
# > nrow(ensGene)
# [1] 57773

freq <- as.data.frame(table(ensGene$external_gene_name))   ## ADD 20/08/17
freq <- freq[order(freq$Freq, decreasing=T),]
# > nrow(subset(freq, Freq == 1))
# [1] 55539
# > subset(freq, Var1 == "RBL1")   ## NOTE 20/08/17
#       Var1 Freq
# 29322 RBL1    2

# -----------------------------------------------------------------------------
# Database: Ensembl Gene / Transcripts
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

## REMOVED 02/11/17: Keep all transcripts (including *_PATCH) to run "sleuth_prep" without warning message
#ensGene.transcript <- subset(ensGene.transcript, chromosome_name %in% paste0("chr", c(1:22, "X", "Y", "MT")))   ## ADD 23/02/17
# > nrow(ensGene.transcript)                             ## Same as line 92
# [1] 196354
# > length(unique(ensGene.transcript$ensembl_gene_id))   ## Same as line 95
# [1] 57773
# > length(unique(ensGene.transcript$external_gene_name))
# [1] 55773

# > length(unique(subset(ensGene.transcript, transcript_biotype == "protein_coding")$ensembl_gene_id))
# [1] 22642
# > ensGene.transcript <- subset(ensGene.transcript, chromosome_name %in% paste0("chr", c(1:22, "X", "Y", "MT")))
# > length(unique(subset(ensGene.transcript, transcript_biotype == "protein_coding")$ensembl_gene_id))
# [1] 20167
# > length(unique(ensGene.transcript.exon$name2))   ## After line 122
# [1] 20167

# =============================================================================
# Reference    : UCSC Genome Browser (Feb 2009/GRCh37/hg19)
# Link(s)      : http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/
# Last Modified: 15/03/18
# =============================================================================

# -----------------------------------------------------------------------------
# Database: Ensembl Gene / Transcripts / Exons
# Table: ensGene / Download Version: 06-Apr-2014
# -----------------------------------------------------------------------------
ensGene.tmp <- readTable(paste0(wd.reference, "ucsc/ensGene.txt.gz"), header=T, rownames=F, sep="")[, -1]   ## Get rid of first column "bin"
rownames(ensGene.tmp) <- ensGene.tmp$name
ensGene.tmp <- subset(ensGene.tmp, chrom %in% paste0("chr", c(1:22, "X", "Y", "M")))
# > nrow(ensGene.tmp)
# [1] 196354
# > length(unique(ensGene.tmp$name2))
# [1] 57773

ensGene.transcript.protein_coding <- subset(ensGene.transcript, transcript_biotype == "protein_coding")
codings <- intersect(ensGene.tmp$name, ensGene.transcript.protein_coding$ensembl_transcript_id)
# > length(codings)
# [1] 81745

colnames <- c("name", "name2", "chrom", "strand", "exonStart", "exonEnd", "cdsStartStat", "cdsEndStat", "exonFrame")
ensGene.transcript.exon <- toTable(0, length(colnames), 0, colnames)
for (t in 1:length(codings)) {
   transcript <- ensGene.tmp[codings[t],]
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
# [1] 603154
# > nrow(subset(ensGene.transcript.exon, exonFrame != -1))                   ## Note that an exonFrames value of -1 means that the exon is entirely UTR
# [1] 555597                                                                 ## https://groups.google.com/a/soe.ucsc.edu/forum/#!topic/genome/U-w4b_ZS2j0

#ensGene.tmp <- subset(subset(ensGene.tmp, cdsStartStat == "cmpl"), cdsEndStat == "cmpl")  ## Only keep transcripts with "complete" CDS start and end information (like in RefGene)
                                                                                           ## https://groups.google.com/a/soe.ucsc.edu/forum/#!topic/genome/Uz5ozC9vkCQ
# -----------------------------------------------------------------------------
# Database: RefSeq Gene
# Table: refGene / Download Version: 04-Mar-2018
# -----------------------------------------------------------------------------
refGene <- readTable(paste0(wd.reference, "ucsc/refGene.txt.gz"), header=T, rownames=F, sep="")[, -1]   ## Get rid of first column "bin"
# > nrow(refGene)
# [1] 70019

refGene <- subset(refGene, chrom %in% paste0("chr", c(1:22, "X", "Y", "M")))
# > nrow(refGene)
# [1] 66872
# > nrow(refGene[grep("NM_", refGene$name),])
# [1] 50799
# > nrow(refGene[grep("NR_", refGene$name),])
# [1] 16073

#table <- as.data.frame(table(refGene$name))
#refGene <- subset(refGene, name %in% subset(table, Freq == 1)$Var1)
#rownames(refGene) <- refGene$name
# Should NOT remove duplicate names first, otherwise e.g. NM_001191005 will be excluded in the final refGene
# > subset(refGene, name == "NM_001191005")
#               name                chrom strand  txStart    txEnd cdsStart   cdsEnd exonCount
# 2325  NM_001191005                 chr1      - 24290836 24306953 24297631 24306745         6
# 62667 NM_001191005 chr1_gl000191_random      -    34172    50281    40974    50073         6

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
chromInfo <- chromInfo[1:24, 1:2]

###
## The Bioinformatician's Guide to the Human Genome
save(ensGene, ensGene.transcript, ensGene.transcript.exon, cytoBand, chromInfo, file=paste0(wd.guidetothe, "hg19.RData"))
save(refGene, file=paste0(wd.guidetothe, "hg19.refGene.RData"))

# -----------------------------------------------------------------------------
# File: rmsk.txt.gz / Download Version: 27-Apr-2009
# -----------------------------------------------------------------------------
rmsk <- readTable(paste0(wd.reference, "ucsc/rmsk.txt.gz"), header=F, rownames=F, sep="")[,-1]
colnames(rmsk) <- c("bin", "milliDiv", "milliDel", "milliIns", "genoName", "genoStart", "genoEnd", "genoLeft", "strand", "repName", "repClass", "repFamily", "repStart", "repEnd", "repLeft", "id")
rmsk <- rmsk[,c("genoName", "genoStart", "genoEnd", "strand", "repName", "repClass", "repFamily")]

save(rmsk, file=paste0(wd.guidetothe, "hg19.rmsk.RData"))

# =============================================================================
# Reference: Miscellaneous
# Last Modified: 13/05/17
# =============================================================================

# -----------------------------------------------------------------------------
# File: human-genome.1kb-grid.bed
# Link(s): https://github.com/andrej-fischer/cloneHD
#          ftp://ftp.sanger.ac.uk/pub/teams/153/bed/
# -----------------------------------------------------------------------------
chrs <- paste0("chr", c(1:22, "X", "Y"))

#bed.raw <- readTable(paste0(wd.reference, "collections/human-genome.1kb-grid.bed.gz"), header=F, rownames=F, sep="")
#colnames(bed.raw) <- c("CHR", "START", "END")
#bed.raw$BED <- paste0("P", rownames(bed.raw))
#rownames(bed.raw) <- bed.raw$BED
bed <- readTable(paste0(wd.ngs, "coverage/S00022/S00022_ANALYSIS/S00022_cn.txt"), header=F, rownames=T, sep="")[,c(1:4,10)]
colnames(bed) <- c("BED", "CHR", "START", "END", "GC")
rownames(bed) <- bed$BED
bed <- bed[,-1]

save(bed, chrs, file=paste0(wd.guidetothe, "hg19.1kb.gc.RData"))

# =============================================================================
# Reference    : UCSC Table Browser (Feb 2009/GRCh37/hg19)
# Link(s)      : https://genome.ucsc.edu/cgi-bin/hgTables
#                http://genome.ucsc.edu/FAQ/FAQdownloads#download10
# Last Modified: 14/03/18
# =============================================================================
