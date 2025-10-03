# =============================================================================
# Name         : 
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 12/09/25; 20/02/25; 25/07/24; 14/03/24
# =============================================================================
wd.src <- "/nfs/users/nfs_t/ty2/dev/R"            ## ty2@farm
#wd.src <- "/Users/ty2/Work/dev/R"                ## ty2@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Commons.R", "SingleCellTranscriptomics.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
wd <- "/lustre/scratch125/casm/staging/team294/ty2"   ## ty2@farm22
#wd <- "/Users/ty2/Work/sanger/ty2"                   ## ty2@localhost
BASE <- "SSC"
base <- tolower(BASE)

wd.rna <- file.path(wd, BASE)
wd.rna.raw <- file.path(wd.rna, "ngs", "10x")

wd.de    <- file.path(wd.rna, "analysis", paste0(base, ""))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

load(file.path(wd.de.data, paste0("ssc_filtered_normalised_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_ordered_annotated.RData")))

so.integrated$orig.ident <- ifelse(
   so.integrated$orig.ident == "SeuratProject",
   "VL00297",
   so.integrated$orig.ident
)
unique(so.integrated$orig.ident)
so.integrated$patient_id <- str_extract(so.integrated$orig.ident, "^[^_]+")

so.integrated$cell_type <- Idents(so.integrated)

so.integrated@meta.data <- so.integrated@meta.data %>%
   mutate(age_group = cut(
      age,
      breaks = c(20, 30, 40, 50, 60, 70, 80),
      labels = c("20-30", "30-40", "40-50", "50-60", "60-70", "70-80"),
      right = FALSE
   ))

# -----------------------------------------------------------------------------
# Prepare for PacBio
# -----------------------------------------------------------------------------
library(Biostrings)

# Step 1: Remove the "-X_Y" suffix from cell barcodes
so.integrated@meta.data$clean_barcodes <- gsub("-\\d+_\\d+$", "", rownames(so.integrated@meta.data))  # Cleaned barcode without suffix
# Step 2: Compute the reverse complement of the cleaned barcodes
so.integrated@meta.data$clean_barcodes_rc <- as.character(reverseComplement(DNAStringSet(so.integrated@meta.data$clean_barcodes)))

# Step 2: Create a table of barcode frequencies
barcode_counts <- table(so.integrated@meta.data$clean_barcodes)

# Step 3: Identify duplicated barcodes (frequency > 1)
duplicated_barcodes <- as.data.frame(barcode_counts[barcode_counts > 1])
colnames(duplicated_barcodes) <- c("barcodes", "Frequency")
duplicated_barcodes[order(-duplicated_barcodes$Frequency),]

##
Clean_Barcode = gsub("-\\d+_\\d+$", "", rownames(so.integrated@meta.data))  # Cleaned barcode without suffix

# Step 3: Extract relevant information from metadata
barcode_table <- data.frame(
	  Seurat_Barcode = rownames(so.integrated@meta.data),  # Original barcode with "-X_Y" suffix
	  Clean_Barcode = Clean_Barcode,
	  Reverse_Complement = as.character(reverseComplement(DNAStringSet(Clean_Barcode))),
	  Sample_ID = so.integrated@meta.data$sample.id,  # Sample ID
	  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# Outputs for PacBio
# Last Modified: 20/02/25
# -----------------------------------------------------------------------------
samples.1 <- c("PD53624b_M", "PD53623b_2N", "PD53623b_4N", "PD53623b_ST2")  # PacBio_B1
samples.2 <- c("PD53623b_M", "PD53621b_M", "PD53626b_ST1")                  # PacBio_B2
samples <- c(samples.1, samples.2)

# Extract metadata from the Seurat object
metadata <- so.integrated@meta.data
# Filter metadata for the specified samples
filtered_metadata <- metadata[metadata$sample.id %in% samples, ]
# Extract cell barcodes (rownames are typically cell barcodes in Seurat metadata)
filtered_metadata$cell_barcode <- rownames(filtered_metadata)

# Select the desired columns
output_table <- filtered_metadata[, c("sample.id", "cell_barcode", "clean_barcodes", "clean_barcodes_rc", "cell_type", "age", "age_group")]
# Save the table to a CSV file
writeTable(output_table, file=file.path(wd, "SSC/ngs/pacbio",  "10x_pacbio_cell_barcode_rc.txt"), rownames=F, colnames=T, sep="\t")

# Step 2: Create a table of barcode frequencies
barcode_counts <- table(output_table$clean_barcodes_rc)
# Step 3: Identify duplicated barcodes (frequency > 1)
duplicated_barcodes <- as.data.frame(barcode_counts[barcode_counts > 1])
colnames(duplicated_barcodes) <- c("barcodes", "Frequency")
duplicated_barcodes[order(-duplicated_barcodes$Frequency),]

# Ensure the 'Sample ID' column respects the specified order
output_table$`sample.id` <- factor(output_table$`sample.id`, levels = samples)
# Summarise the number of cells per sample ID
cell_counts <- output_table %>%
	  dplyr::group_by(`sample.id`) %>%
	  dplyr::summarise(cell_count = n())

# Save the table to a CSV file
writeTable(cell_counts, file=file.path(wd, "SSC/ngs/pacbio", "10x_pacbio_cell_counts.txt"), rownames=F, colnames=T, sep="\t")

# -----------------------------------------------------------------------------
# Output an .h5ad file from a Seurat object 
# Last Modified: 28/03/25
# -----------------------------------------------------------------------------
load(file.path(wd.de.data, paste0("ssc_filtered_normalised_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_ordered_annotated_monocle3+phase.RData")))

library(SeuratDisk)

h5Seurat <- file.path(wd, "SSC/ngs/10x", "ssc_annotated_monocle3+phase_clean.h5Seurat")
SaveH5Seurat(so.integrated, filename = h5Seurat)
Convert(h5Seurat, dest = "h5ad")
