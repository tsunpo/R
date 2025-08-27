# =============================================================================
# Manuscript   :
# Chapter      :
# Name         : 
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 23/07/25; 15/04/25
# =============================================================================

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
wd.de.data  <- file.path("/lustre/scratch125/casm/staging/team294/ty2/SSC/ngs/pacbio")

# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
samplesheet <- read.table(
   file = file.path(wd.de.data, paste0("samplesheet.csv")),
   sep = ",",
   header = T
)
pd_ids <- sort(unique(samplesheet$donor_id))

# Prepare a data.frame
samplesheet_new <- data.frame(
   id = character(),
   bam = character(),
   cell_barcodes = character(),
   celltypes = character(),
   mutations = character(),
   stringsAsFactors = FALSE
)

# Loop through pd_ids and construct the paths
for (pd_id in pd_ids) {
   new_row <- data.frame(
      id = paste0(pd_id),
      bam = file.path(wd.de.data, "mn7", pd_id, paste0(pd_id, ".bam")),
      cell_barcodes = file.path(wd.de.data, "mn7", pd_id, paste0(pd_id, "_trencadis_cell_barcodes.txt")),
      celltypes = file.path(wd.de.data, "mn7", pd_id, paste0(pd_id, "_trencadis_celltypes.csv")),
      mutations = file.path(wd.de.data, "mn7", pd_id, paste0(pd_id, "_trencadis_mutations.tsv")),
      stringsAsFactors = FALSE
   )
   
   # Add to main samplesheet
   samplesheet_new <- rbind(samplesheet_new, new_row)
}

# Now, save it to CSV
write.table(
   samplesheet_new,
   file = file.path(wd.de.data, "nf-trencadis-seq", "samplesheet_trencadis.csv"),
   sep = ",",
   row.names = FALSE,
   quote = FALSE
)

# -----------------------------------------------------------------------------
# *_trencadis_mutations.tsv (from mn7)
# -----------------------------------------------------------------------------
for (pd_id in pd_ids) {
   vcf <- read.delim(file.path(wd.de.data, "mn7", "mn7.hg38.bed.40.coding.unique.vcf.bed"), header=F)
 
   mutations <- vcf[, c(1, 3, 4, 5)]
   colnames(mutations) <- c("chr", "pos", "ref", "alt")
 
   write.table(
      mutations,
      file = file.path(wd.de.data, "mn7", pd_id, paste0(pd_id, "_trencadis_mutations.tsv")),
      sep = "\t",
      col.names = TRUE,
      row.names = FALSE,
      quote = FALSE
   )
}

# -----------------------------------------------------------------------------
# *_trencadis_mutations.tsv (from Scomatic)
# -----------------------------------------------------------------------------
for (pd_id in pd_ids) {
   lines <- readLines(file.path(wd.de.data, "out", pd_id, paste0(pd_id, ".calling.step2.pass.tsv")))
   header_line <- grep("^#CHROM", lines)
   header <- unlist(strsplit(sub("^#", "", lines[header_line]), "\t"))
   vcf <- read.table(
      file = file.path(wd.de.data, "out", pd_id, paste0(pd_id, ".calling.step2.pass.tsv")),
      skip = header_line,
      sep = "\t",
      col.names = header,
      stringsAsFactors = FALSE,
      check.names = FALSE
   )
   
   mutations <- vcf[, c(1, 2, 4, 5)]
   colnames(mutations) <- c("chr", "pos", "ref", "alt")
   
   write.table(
      mutations,
      file = file.path(wd.de.data, "out", pd_id, paste0(pd_id, "_trencadis_mutations.tsv")),
      sep = "\t",
      col.names = TRUE,
      row.names = FALSE,
      quote = FALSE
   )
}

# -----------------------------------------------------------------------------
# *_trencadis_celltypes.csv and *_trencadis_barcodes.txt
# -----------------------------------------------------------------------------
for (pd_id in pd_ids) {
   celltypes <- read.table(
      file = file.path(wd.de.data, paste0("celltypes_progenitor.tsv")),
      sep = "\t",
      header = TRUE
   )
   
   colnames(celltypes) <- c("sample_id", "barcode", "celltype")
   celltypes_new <- subset(celltypes, sample_id %in% subset(samplesheet, donor_id == pd_id)$sample_id)
 
   write.table(
      celltypes_new[, -1],
      file = file.path(wd.de.data, "mn7", pd_id, paste0(pd_id, "_trencadis_celltypes.csv")),
      sep = ",",
      col.names = TRUE,
      row.names = FALSE,
      quote = FALSE
   )
   
   write.table(
      celltypes_new[, 2],
      file = file.path(wd.de.data, "mn7", pd_id, paste0(pd_id, "_trencadis_cell_barcodes.txt")),
      col.names = FALSE,
      row.names = FALSE,
      quote = FALSE
   )
}

# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
mn7 <- c("KDM5B", "PTPN11", "NF1", "SMAD6", "CUL3", "MIB1", "RASA2", "PRRC2A", "PTEN", "RIT1", "ROBO1", "DDX3X", "CSNK2B", "KRAS", "FGFR3", "PPM1D", "ARID1A", "BRAF", "HRAS", "KMT2E", "EP300", "SCAF4", "BMPR2", "TCF12", "CCAR2", "DHX9", "NSD1", "LZTR1", "FGFR2", "ARHGAP35", "CBL", "SSX1", "RBM12", "FAM222B", "SMAD4", "AR", "KDM5C", "KMT2D", "CTNNB1", "RAF1")

write.table(
   mn7,
   file = file.path(wd.de.data, paste0("genes_mn7.txt")),
   col.names = FALSE,
   row.names = FALSE,
   quote = FALSE
)

