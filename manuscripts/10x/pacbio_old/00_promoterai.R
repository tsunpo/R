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
wd.de.data  <- file.path("/lustre/scratch125/casm/staging/team294/ty2/SSC/ngs/pacbio/")

vcf <- read.delim(file.path(wd.de.data, "mn7", "mn7.hg38.bed.unique.vcf"), header=F)

# -----------------------------------------------------------------------------
# Ensembl's consequence severity ranking
# -----------------------------------------------------------------------------
# Ensembl canonical severity order
vep_severity <- c(
   "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
   "stop_gained", "frameshift_variant", "stop_lost", "start_lost", "transcript_amplification",
   "inframe_insertion", "inframe_deletion", "missense_variant", "protein_altering_variant",
   "splice_region_variant", "incomplete_terminal_codon_variant", "start_retained_variant",
   "stop_retained_variant", "synonymous_variant", "coding_sequence_variant", "mature_miRNA_variant",
   "5_prime_UTR_variant", "3_prime_UTR_variant", "non_coding_transcript_exon_variant",
   "intron_variant", "NMD_transcript_variant", "non_coding_transcript_variant",
   "upstream_gene_variant", "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification",
   "TF_binding_site_variant", "regulatory_region_ablation", "regulatory_region_amplification",
   "feature_elongation", "regulatory_region_variant", "feature_truncation", "intergenic_variant"
)

# Extract main consequence term
vcf$Main_Consequence <- sapply(
   strsplit(vcf$V7, ","), 
   function(x) x[1]
)

# Assign severity rank
vcf$Severity_Rank <- match(vcf$Main_Consequence, vep_severity)

# Plot histogram
library(ggplot2)
library(ggtext)

# Ensure proper factor ordering by severity
vcf$Main_Consequence <- factor(
   sapply(strsplit(vcf$V7, ","), function(x) x[1]),
   levels = rev(vep_severity)
)

bold_these <- c(
   "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
   "stop_gained", "frameshift_variant", "stop_lost", "start_lost", "transcript_amplification",
   "inframe_insertion", "inframe_deletion", "missense_variant",
   "protein_altering_variant", "splice_region_variant", "incomplete_terminal_codon_variant",
   "start_retained_variant", "stop_retained_variant", "synonymous_variant", "coding_sequence_variant", "mature_miRNA_variant"
)

vcf$Main_Consequence_label <- sapply(
   sapply(strsplit(vcf$V7, ","), function(x) x[1]),
   function(term) {
      if (term %in% bold_these) paste0("**", term, "**") else term
   }
)

vcf$Main_Consequence_label <- factor(
   vcf$Main_Consequence_label,
   levels = rev(sapply(vep_severity, function(term) {
      if (term %in% bold_these) paste0("**", term, "**") else term
   }))
)

# Plot with transcript_ablation on top
plot <- ggplot(
   subset(vcf, V6 != ""),
   aes(x = Main_Consequence_label)
) +
   geom_bar(fill = "#00BFC4") +
   scale_x_discrete(drop = FALSE) +
   coord_flip() +
   theme_minimal(base_size = 12) +
   theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 17),
      axis.text.y = ggtext::element_markdown(size = 11.5),
      axis.text.x = element_text(size = 12),
      axis.title = element_text(size = 15)
 ) +
 ylab("Count") +
 xlab("Consequence") +
 ggtitle("61,142 SNVs")

# Save the plot as a PNG file
png(file.path(wd.de.data, "mn7/promoterai", "mn7.unique.png"), width = 600, height = 600)
print(plot)  # necessary to render the saved plot
dev.off()

# -----------------------------------------------------------------------------
# Convert to PromoterAI's format
# -----------------------------------------------------------------------------
vcf.promoter <- subset(vcf, Main_Consequence %in% c("upstream_gene_variant", "5_prime_UTR_variant", "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant", "regulatory_region_variant", "regulatory_region_ablation", "regulatory_region_amplification", "feature_truncation", "feature_elongation"))

vcf.promoter$chrom <- paste0("chr", mapply(x = 1:nrow(vcf.promoter), function(x) unlist(strsplit(vcf.promoter$V2[x], "_"))[1]))
vcf.promoter$pos <- mapply(x = 1:nrow(vcf.promoter), function(x) unlist(strsplit(vcf.promoter$V2[x], ":"))[2])
vcf.promoter$ref <- mapply(x = 1:nrow(vcf.promoter), function(x) unlist(strsplit(vcf.promoter$V1[x], "_"))[3])
vcf.promoter$alt <- mapply(x = 1:nrow(vcf.promoter), function(x) unlist(strsplit(vcf.promoter$V1[x], "_"))[4])

# Now, save it to CSV
write.table(
   vcf.promoter[, c("chrom", "pos", "ref", "alt")],
      file = file.path(wd.de.data, "mn7", "promoterai", "mn7.promoter.tsv"),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
)

# -----------------------------------------------------------------------------
# Convert to PrimateAI's format
# -----------------------------------------------------------------------------
vcf.missense <- subset(vcf, Main_Consequence == "missense_variant")

vcf.missense$CHROM <- paste0("chr", mapply(x = 1:nrow(vcf.missense), function(x) unlist(strsplit(vcf.missense$V1[x], "_"))[1]))
vcf.missense$POS <- mapply(x = 1:nrow(vcf.missense), function(x) unlist(strsplit(vcf.missense$V2[x], ":"))[2])
vcf.missense$REF <- mapply(x = 1:nrow(vcf.missense), function(x) unlist(strsplit(vcf.missense$V1[x], "_"))[3])
vcf.missense$ALT <- mapply(x = 1:nrow(vcf.missense), function(x) unlist(strsplit(vcf.missense$V1[x], "_"))[4])
vcf.missense$ID <- "."
vcf.missense$QUAL <- "."
vcf.missense$FILTER <- "."
vcf.missense$INFO <- "."

# Now, save it to CSV
write.table(
 vcf.missense[, c("CHROM", "POS", "REF", "ALT", "ID", "QUAL", "FILTER", "INFO")],
   file = file.path(wd.de.data, "mn7", "promoterai", "mn7.missense.vcf"),
   sep = "\t",
   row.names = FALSE,
   quote = FALSE
)




mutations <- vcf[, c(1, 3, 4, 5)]
colnames(mutations) <- c("chr", "pos", "ref", "alt")




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

