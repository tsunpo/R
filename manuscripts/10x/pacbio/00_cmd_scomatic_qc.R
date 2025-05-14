args <- commandArgs(trailingOnly = TRUE)
pd_id <- as.character(args[1])

# =============================================================================
# Manuscript   :
# Chapter      :
# Name         : 
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 15/04/25
# =============================================================================
library(ggplot2)
library(dplyr)

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
wd.de.data  <- file.path("/lustre/scratch127/casm/team294rr/ty2/SSC/ngs/pacbio/out_max13", pd_id)
wd.de.plots <- file.path(wd.de.data, "00_QC")
dir.create(wd.de.plots, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Pe-processing
# -----------------------------------------------------------------------------
# VCF file
lines <- readLines(file.path(wd.de.data, paste0(pd_id, ".calling.step2.pass.tsv")))
header_line <- grep("^#CHROM", lines)
header <- unlist(strsplit(sub("^#", "", lines[header_line]), "\t"))
vcf <- read.table(
   file = file.path(wd.de.data, paste0(pd_id, ".calling.step2.pass.tsv")),
   skip = header_line,
   sep = "\t",
   col.names = header,
   stringsAsFactors = FALSE,
   check.names = FALSE
)

# MAF file
maf <- read.table(file.path(wd.de.data, paste0(pd_id, ".calling.step2.pass.nochr.maf")), sep = "\t", header = T, comment.char = "#", quote = "", fill = TRUE)

vcf$key <- paste0(gsub("^chr", "", vcf$CHROM), ":", vcf$Start)
maf$key <- paste0(maf$Chromosome, ":", maf$Start_Position)
table(vcf$key %in% maf$key)

merged <- merge(vcf, maf,
                by.x = "key",
                by.y = "key",
                all.x = TRUE)  # use all = TRUE if you want outer/full join

mn7 <- c("KDM5B", "PTPN11", "NF1", "SMAD6", "CUL3", "MIB1", "RASA2", "PRRC2A", "PTEN", "RIT1", "ROBO1", "DDX3X", "CSNK2B", "KRAS", "FGFR3", "PPM1D", "ARID1A", "BRAF", "HRAS", "KMT2E", "EP300", "SCAF4", "BMPR2", "TCF12", "CCAR2", "DHX9", "NSD1", "LZTR1", "FGFR2", "ARHGAP35", "CBL", "SSX1", "RBM12", "FAM222B", "SMAD4", "AR", "KDM5C", "KMT2D", "CTNNB1", "RAF1")
merged$Gene_Class <- ifelse(merged$Hugo_Symbol %in% mn7, "In List", "Not in List")
# Set factor levels to control plotting order
merged$Gene_Class <- factor(merged$Gene_Class, levels = c("Not in List", "In List"))

# Filter
merged <- merged[merged$dbSNP_RS == "" | merged$dbSNP_RS == "novel", ]
write.table(merged, file = file.path(wd.de.data, paste0(pd_id, ".calling.step2.pass.nochr.vep.maf.vcf")), sep = "\t", quote = FALSE, row.names = FALSE)
prefix <- "novel+COSV_"
#prefix <- ""

inList <- length(intersect(unique(merged$SYMBOL), mn7))
notInList <- length(unique(merged$SYMBOL)) - inList
inList.txt <- paste0("In (n=", inList, ")")
notInList.txt <- paste0("Not in (n=", notInList, ")")

##
table <- as.data.frame(sort(table(merged$Cell_types), decreasing=T))
write.table(table, file = file.path(wd.de.data, paste0(pd_id, "_table.txt")), sep = "\t", quote = FALSE, row.names = FALSE)


# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
# Step 1: Full counts for all Cell_types
celltype_counts <- as.data.frame(table(merged$Cell_types))
colnames(celltype_counts) <- c("Cell_type", "Total")

# Step 2: Subset counts for "In List" genes
inlist_subset <- subset(merged, Gene_Class == "In List")

inlist_counts <- as.data.frame(table(inlist_subset$Cell_types))
colnames(inlist_counts) <- c("Cell_type", "In_List")

inlist_symbols <- inlist_subset %>%
   dplyr::group_by(Cell_types) %>%
   dplyr::summarise(Hugo_Symbols = paste(unique(Hugo_Symbol), collapse = ", ")) %>%
   dplyr::rename(Cell_type = Cell_types)

# Merge Hugo symbols
inlist_counts <- merge(inlist_counts, inlist_symbols, by = "Cell_type", all.x = TRUE)

# Step 3: Merge the two
celltype_counts <- merge(celltype_counts, inlist_counts, by = "Cell_type", all.x = TRUE)
celltype_counts$In_List[is.na(celltype_counts$In_List)] <- 0

# Step 5 (optional): Order by Total
celltype_counts <- celltype_counts[order(celltype_counts$In_List, decreasing = TRUE), ]
write.table(celltype_counts, file = file.path(wd.de.plots, paste0(pd_id, "_celltype_mutation_counts.txt")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(inlist_subset, file = file.path(wd.de.plots, paste0(pd_id, "_VCF.txt")), sep = "\t", quote = FALSE, row.names = FALSE)

# -----------------------------------------------------------------------------
# Depth of Coverage (DP) Distribution
# -----------------------------------------------------------------------------
plot <- ggplot(merged, aes(x = Dp)) +
   geom_histogram(bins = 50) +
   labs(title = "All Genes", x = "Read Depth (DP)", y = "Count") +
   theme(plot.title = element_text(hjust = 0.5))
# Print the volcano plot
#pdf(file = file.path(wd.de.plots, paste0(prefix, "hist_DP", "", ".pdf")), width = 5, height = 5)
#print(plot)
#dev.off()

# Define your 40 genes
plot <- ggplot(merged, aes(x = Dp, fill = Gene_Class)) +
   geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
   scale_fill_manual(
      values = c("In List" = "blue", "Not in List" = "red"),
      labels = c("Not in List" = notInList.txt, "In List" = inList.txt)
    ) +
   labs(title = "In 40 Gene List", x = "Read Depth (DP)", y = "Count", fill = "Hugo Symbol") +
   theme(plot.title = element_text(hjust = 0.5))
# Print
pdf(file = file.path(wd.de.plots, paste0(prefix, "hist_DP_mn7", "", ".pdf")), width = 5, height = 5)
print(plot)
dev.off()

# -----------------------------------------------------------------------------
# Depth of Coverage (DP) Distribution
# -----------------------------------------------------------------------------
wilcox_result <- wilcox.test(Dp ~ Gene_Class, data = merged)

# Boxplot
boxplot <- ggplot(merged, aes(x = Gene_Class, y = Dp, fill = Gene_Class)) +
   geom_boxplot(alpha = 0.6) +
   scale_fill_manual(values = c("In List" = "blue", "Not in List" = "red")) +
   labs(
      title = "Read Depth by Gene Category",
      x = "Gene Category",
      y = "Read Depth (DP)",
      fill = "Hugo Symbol"
   ) +
   theme_minimal() +
   theme(plot.title = element_text(hjust = 0.5)) +
   annotate("text", x = 1.5, y = max(merged$Dp, na.rm = TRUE) * 0.95,
            label = paste0("Wilcoxon p = ", signif(wilcox_result$p.value, 3)))
# Print the volcano plot
pdf(file = file.path(wd.de.plots, paste0(prefix, "boxplot_DP_mn7", "", ".pdf")), width = 5, height = 5)
print(boxplot)
dev.off()

# -----------------------------------------------------------------------------
# Variant Allele Frequency (VAF) Distribution
# -----------------------------------------------------------------------------
plot <- ggplot(merged, aes(x = VAF)) +
   geom_histogram(bins = 50) +
   labs(title = "All Genes", x = "Variant Allele Frequency (VAF)", y = "Count") +
   theme(plot.title = element_text(hjust = 0.5))
# Print the volcano plot
#pdf(file = file.path(wd.de.plots, paste0(prefix, "hist_VAF", "", ".pdf")), width = 5, height = 5)
#print(plot)
#dev.off()

# Define your 40 genes
plot <- ggplot(merged, aes(x = VAF, fill = Gene_Class)) +
   geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
   scale_fill_manual(
      values = c("In List" = "blue", "Not in List" = "red"),
      labels = c("Not in List" = notInList.txt, "In List" = inList.txt)
   ) +
   labs(title = "In 40 Gene List", x = "Variant Allele Frequency (VAF)", y = "Count", fill = "Hugo Symbol") +
   theme(plot.title = element_text(hjust = 0.5))
# Print the volcano plot
pdf(file = file.path(wd.de.plots, paste0(prefix, "hist_VAF_mn7", "", ".pdf")), width = 5, height = 5)
print(plot)
dev.off()

# -----------------------------------------------------------------------------
# DP vs. (VAF)
# -----------------------------------------------------------------------------
plot <- ggplot(merged, aes(x = Dp, y = VAF, color = Gene_Class)) +
   geom_point(alpha = 0.7, size = 2) +
   scale_color_manual(values = c("In List" = "blue", "Not in List" = "red")) +
   labs(
      title = "DP vs. VAF",
      x = "Read Depth (DP)",
      y = "Variant Allele Frequency (VAF)",
      color = "Hugo Symbol"
   ) +
 theme_minimal() +
 theme(plot.title = element_text(hjust = 0.5))

# Print the volcano plot
pdf(file = file.path(wd.de.plots, paste0(prefix, "hist_DP-vs-VAF_mn7", "", ".pdf")), width = 5, height = 5)
print(plot)
dev.off()

# -----------------------------------------------------------------------------
# Allele Balance: ALT / (ALT + REF)
# Shows how often ALT alleles are observed relative to REF.
# In single-cell, this can reveal allelic dropout or technical noise.
# -----------------------------------------------------------------------------
merged$allele_balance <- merged$Cc / (merged$Cc + merged$Nc)

plot <- ggplot(merged, aes(x = allele_balance)) +
   geom_histogram(bins = 50) +
   labs(title = "All Genes", x = "Allele Balance (ALT / (ALT + REF))", y = "Count") +
   theme(plot.title = element_text(hjust = 0.5))
# Print the volcano plot
#pdf(file = file.path(wd.de.plots, paste0(prefix, "hist_Allele Balance", "", ".pdf")), width = 5, height = 5)
#print(plot)
#dev.off()

# Define your 40 genes
plot <- ggplot(merged, aes(x = allele_balance, fill = Gene_Class)) +
   geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
   scale_fill_manual(
      values = c("In List" = "blue", "Not in List" = "red"),
      labels = c("Not in List" = notInList.txt, "In List" = inList.txt)
   ) +
   labs(title = "In 40 Gene List", x = "Allele Balance (ALT / (ALT + REF))", y = "Count", fill = "Hugo Symbol") +
   theme(plot.title = element_text(hjust = 0.5))
# Print the volcano plot
pdf(file = file.path(wd.de.plots, paste0(prefix, "hist_Allele Balance_mn7", "", ".pdf")), width = 5, height = 5)
print(plot)
dev.off()
