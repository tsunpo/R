# =============================================================================
# Reference: The Bioinformatician's Guide to the Human Genome
# Name: guide-to-the/hg19.R
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 26/04/17
# =============================================================================
wd.reference <- "/Users/tpyang/Work/local/reference/hg19/"
wd.guidetothe <- "/Users/tpyang/Work/local/R/guide-to-the/"

source("/Users/tpyang/Work/local/R/handbook-of/Analytics.R")

# =============================================================================
# Source: Ensembl BioMart
# Link: http://grch37.ensembl.org/biomart/martview/
#       ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.75.cdna.all.fa.gz
# Last Modified: 27/04/17
# =============================================================================

# -----------------------------------------------------------------------------
# Database: Ensembl Gene / Transcripts
# Table: ensGene / Count: 196,354 (out of 215,171; Unique results only) / Date: 2017-03-25
# -----------------------------------------------------------------------------
ensGene <- readTable(paste0(wd.reference, "ensembl/BioMart.GRCh37.p13.EnsemblGene_20170325.txt.gz"), header=T, rownames=F, sep="")
colnames(ensGene) <- c("ensembl_transcript_id", "ensembl_gene_id", "chromosome_name", "strand", "transcript_start", "transcript_end", "start_position", "end_position", "transcript_biotype", "gene_biotype", "external_gene_name")   ##"hugo_gene", "refseq_mrna", "refseq_mrna"
rownames(ensGene) <- ensGene$ensembl_transcript_id
# > nrow(ensGene)
# [1] 215170

ensGene <- subset(ensGene, chromosome_name %in% c(1:22, "X", "Y", "MT"))   ## ADD 23/02/17
ensGene$chromosome_name <- paste0("chr", ensGene$chromosome_name)          ## ADD 27/04/17
# > nrow(ensGene)
# [1] 196354
# > length(unique(ensGene$ensembl_gene_id))
# [1] 57773
# > length(unique(ensGene$external_gene_name))
# [1] 55773

# -----------------------------------------------------------------------------
# Database: Ensembl Gene / Genes
# Table: ensGene / Count: 57,773 (out of 63,677; Unique results only) / Date: 2017-03-27
# -----------------------------------------------------------------------------
ensGene.gene <- readTable(paste0(wd.reference, "ensembl/BioMart.GRCh37.p13.EnsemblGene.ensembl_gene_id_20170327.txt.gz"), header=T, rownames=F, sep="")
colnames(ensGene.gene) <- c("ensembl_gene_id", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype", "external_gene_name")
rownames(ensGene.gene) <- ensGene.gene$ensembl_gene_id
# > nrow(ensGene.gene)
# [1] 63677

ensGene.gene <- subset(ensGene.gene, chromosome_name %in% c(1:22, "X", "Y", "MT"))   ## ADD 23/02/17
ensGene.gene$chromosome_name <- paste0("chr", ensGene.gene$chromosome_name)          ## ADD 27/04/17
# > nrow(ensGene.gene)
# [1] 57773
# > subset(freq, Var1 == "RBL1")   ## NOTE 20/08/17
# Var1 Freq
# 30073 RBL1    2

freq <- as.data.frame(table(ensGene$external_gene_name))   ## ADD 20/08/17
freq <- freq[order(freq$Freq, decreasing=T),]
dim(subset(freq, Freq == 1))


# =============================================================================
# Source: UCSC Table Browser
# Link: https://genome.ucsc.edu/cgi-bin/hgTables
#       http://genome.ucsc.edu/FAQ/FAQdownloads#download10
# Last Modified: 28/01/17
# =============================================================================

# -----------------------------------------------------------------------------
# Database: RefSeq Gene
# Table: refGene / Row Count: 57,214 (out of 62,815) / Date: 2017-01-25
# -----------------------------------------------------------------------------
refGene <- readTable(paste0(wd.reference, "ucsc/ucsc.hg19.refGene_20170125.txt.gz"), header=T, rownames=F, sep="")[, -1]   ## Rid first column "bin"
# > nrow(refGene)
# [1] 62815

table <- as.data.frame(table(refGene$name))
refGene <- subset(refGene, name %in% subset(table, Freq == 1)$Var1)
rownames(refGene) <- refGene$name
# Remove duplicate names (not random chrs) first, otherwise e.g. NM_001191005 will be included in the final refGene
# > subset(refGene.raw, name == "NM_001191005")
#               name                chrom strand  txStart    txEnd cdsStart   cdsEnd exonCount
# 2325  NM_001191005                 chr1      - 24290836 24306953 24297631 24306745         6
# 62667 NM_001191005 chr1_gl000191_random      -    34172    50281    40974    50073         6
# > nrow(refGene)
# [1] 57262

refGene <- subset(refGene, chrom %in% paste0("chr", c(1:22, "X", "Y", "M")))
# > nrow(refGene)
# [1] 57214
# > nrow(refGene[grep("NM_", refGene$name),])
# [1] 44353
# > nrow(refGene[grep("NR_", refGene$name),])
# [1] 12861

###
## The Bioinformatician's Guide to the Transcription
#save(refGene, ensGene, ensGene.gene, file=paste0(wd.guidetothe, "hg19.transcript.RData"))


# =============================================================================
# Source: UCSC Genome Browser (Feb 2009/GRCh37/hg19)
# Link: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/
# Last Modified: 26/04/17
# =============================================================================

# -----------------------------------------------------------------------------
# File: cytoBand.txt.gz / Date: 14-Jun-2009
# -----------------------------------------------------------------------------
cytoBand <- readTable(paste0(wd.reference, "ucsc/cytoBand.txt.gz"), header=F, rownames=F, sep="")
colnames(cytoBand) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")
# > nrow(cytoBand)
# [1] 862

cytoBand$name2 <- paste0(gsub("chr", "", cytoBand$chrom), cytoBand$name)
cytoBand$size <- cytoBand$chromEnd - cytoBand$chromStart

# -----------------------------------------------------------------------------
# File: chromInfo.txt.gz / Date: 27-Apr-2009
# -----------------------------------------------------------------------------
chromInfo <- readTable(paste0(wd.reference, "ucsc/chromInfo.txt.gz"), header=F, rownames=F, sep="")
colnames(chromInfo) <- c("chrom", "size")
chromInfo <- chromInfo[1:24, 1:2]

###
## The Bioinformatician's Guide to the Genome
#save(chromInfo, cytoBand, file=paste0(wd.guidetothe, "hg19.genome.RData"))
save(ensGene, ensGene.gene, refGene, chromInfo, cytoBand, file=paste0(wd.guidetothe, "hg19.RData"))

# -----------------------------------------------------------------------------
# File: rmsk.txt.gz / Date: 27-Apr-2009
# -----------------------------------------------------------------------------
rmsk <- readTable(paste0(wd.reference, "ucsc/rmsk.txt.gz"), header=F, rownames=F, sep="")[,-1]
colnames(rmsk) <- c("bin", "milliDiv", "milliDel", "milliIns", "genoName", "genoStart", "genoEnd", "genoLeft", "strand", "repName", "repClass", "repFamily", "repStart", "repEnd", "repLeft", "id")
rmsk <- rmsk[,c("genoName", "genoStart", "genoEnd", "strand", "repName", "repClass", "repFamily")]

save(rmsk, file=paste0(wd.guidetothe, "hg19.rmsk.RData"))


# =============================================================================
# Source: Miscellaneous
# Last Modified: 13/05/17
# =============================================================================

# -----------------------------------------------------------------------------
# File: human-genome.1kb-grid.bed
# Link: https://github.com/andrej-fischer/cloneHD
#       ftp://ftp.sanger.ac.uk/pub/teams/153/bed/
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
