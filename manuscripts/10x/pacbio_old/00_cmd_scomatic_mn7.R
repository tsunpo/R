# =============================================================================
# Manuscript   :
# Chapter      :
# Name         : 
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 26/06/25
# =============================================================================

all <- readTable(file.path("/Users/ty2/Work/sanger/ty2/SSC/ngs/pacbio", "unique_genes.txt"), sep = "", header = F, rownames=F)

# -----------------------------------------------------------------------------
# Pe-processing
# -----------------------------------------------------------------------------
# VCF file
nano <- read.table(file.path("/Users/ty2/Work/sanger/ty2/SSC/ngs/pacbio", "out_filteredTargVars.vcf.txt"), sep = "\t", header = T, comment.char = "#", quote = "", fill = TRUE)

nano.mn7 <- subset(nano, SYMBOL %in% mn7)
nano.mn7.soma <- subset(nano.mn7, !grepl("rs[0-9]+", Existing_variation))

nano.mn7.soma.coding <- nano.mn7.soma[sapply(
   strsplit(nano.mn7.soma$Consequence, ","),
   function(x) any(x %in% bold_these)
), ]

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
#wd.de.data  <- file.path("/lustre/scratch127/casm/team294rr/ty2/SSC/ngs/pacbio/out_max13", pd_id)
wd.de.data  <- file.path("/Users/ty2/Work/sanger/ty2/SSC/ngs/pacbio/out_max13", pd_id)
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
# 
# Last Modified: 16/06/25
# -----------------------------------------------------------------------------
on_target <- merged[merged$Gene_Class == "In List", ]
off_target <- merged[merged$Gene_Class == "Not in List", ]

extract_main_consequence <- function(x) {
   # pick the most impactful from multi-annotation strings
   primary <- unlist(strsplit(x, ","))[1]
   return(primary)
}

on_target$Main_Consequence <- sapply(on_target$Consequence, extract_main_consequence)
off_target$Main_Consequence <- sapply(off_target$Consequence, extract_main_consequence)

# Compare distributions
table(on_target$Main_Consequence)
table(off_target$Main_Consequence)

# Or normalized
prop.table(table(on_target$Main_Consequence))
prop.table(table(off_target$Main_Consequence))

# Visualize
library(ggplot2)
combined <- rbind(
   data.frame(Consequence = on_target$Main_Consequence, Group = "On-target"),
   data.frame(Consequence = off_target$Main_Consequence, Group = "Off-target")
)

ggplot(combined, aes(x = Consequence, fill = Group)) +
   geom_bar(position = "fill") + coord_flip() +
   ylab("Proportion") + xlab("Consequence") +
   ggtitle("Relative Consequence Distribution: On-target vs Off-target")

# -----------------------------------------------------------------------------
# 
# Last Modified: 16/06/25
# -----------------------------------------------------------------------------
df <- merged

df$Effect_Class <- case_when(
   df$Variant_Classification %in% c("Silent", "Synonymous") ~ "syn",
   df$Variant_Classification %in% c("Missense_Mutation", "Nonsense_Mutation", "Splice_Site", "Frame_Shift", "Nonstop_Mutation") ~ "non-syn",
   TRUE ~ "non-coding"
)
table(df$Effect_Class, df$Gene_Class)

table(df$Consequence, df$Effect_Class)
table(df$Consequence, df$Gene_Class)
table(df$Consequence, df$Effect_Class, df$Gene_Class)  # for 3-way view

# -----------------------------------------------------------------------------
# Check exon coverage
# Last Modified: 16/06/25
# -----------------------------------------------------------------------------
# Count how many cell types each mutation was seen in
library(dplyr)

merged$n_celltypes <- mapply(function(dp) {
   # Handle NA safely
   if (is.na(dp) || dp == "") {
      return(0)
   } else {
      return(length(strsplit(dp, ",")[[1]]))
   }
}, merged$Cell_types)

dp_single_celltype <- subset(merged, n_celltypes == 1)

on_target <- dp_single_celltype %>% filter(Gene_Class == "In List")
off_target <- dp_single_celltype %>% filter(Gene_Class == "Not in List")

on_dp <- as.numeric(on_target$Dp)
off_dp <- as.numeric(off_target$Dp)

on_dp <- on_dp[!is.na(on_dp)]
off_dp <- off_dp[!is.na(off_dp)]

# Wilcoxon test
wilcox.test(on_dp, off_dp)

# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
# Boxplot
dp_data <- rbind(
   data.frame(Dp = as.numeric(on_target$Dp), Group = "On-target"),
   data.frame(Dp = as.numeric(off_target$Dp), Group = "Off-target")
)

# Compute Wilcoxon p-value
wilcox_res <- wilcox.test(Dp ~ Group, data = dp_data)
pval <- wilcox_res$p.value

# Create the plot with scientific notation for p-value
ggplot(dp_data, aes(x = Group, y = Dp, fill = Group)) +
   geom_boxplot(outlier.shape = NA, alpha = 0.6) +
   geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
   theme_minimal() +
   theme(legend.position = "none") +
   ylab("Read Depth") +
   xlab("") +
   annotate("text",
          x = 1.5,
          y = max(dp_data$Dp, na.rm = TRUE) * 1.05,
          label = paste0("Wilcoxon p = ", format(pval, scientific = TRUE, digits = 3)),
          size = 5,
          fontface = "italic")

# Boxplot
dp_data <- rbind(
   data.frame(Nc = as.numeric(on_target$Nc), Group = "On-target"),
   data.frame(Nc = as.numeric(off_target$Nc), Group = "Off-target")
)

# Compute Wilcoxon p-value
wilcox_res <- wilcox.test(Nc ~ Group, data = dp_data)
pval <- wilcox_res$p.value

# Create the plot with scientific notation for p-value
ggplot(dp_data, aes(x = Group, y = Nc, fill = Group)) +
   geom_boxplot(outlier.shape = NA, alpha = 0.6) +
   geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
   theme_minimal() +
   theme(legend.position = "none") +
   ylab("Number of Cells") +
   xlab("") +
   annotate("text",
          x = 1.5,
          y = max(dp_data$Nc, na.rm = TRUE) * 1.05,
          label = paste0("Wilcoxon p = ", format(pval, scientific = TRUE, digits = 3)),
          size = 5,
          fontface = "italic")

# Boxplot
dp_data <- rbind(
   data.frame(Nc = as.numeric(on_target$VAF), Group = "On-target"),
   data.frame(Nc = as.numeric(off_target$VAF), Group = "Off-target")
)

# Compute Wilcoxon p-value
wilcox_res <- wilcox.test(Nc ~ Group, data = dp_data)
pval <- wilcox_res$p.value

# Create the plot with scientific notation for p-value
ggplot(dp_data, aes(x = Group, y = Nc, fill = Group)) +
   geom_boxplot(outlier.shape = NA, alpha = 0.6) +
   geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
   theme_minimal() +
   theme(legend.position = "none") +
   ylab("VAF") +
   xlab("") +
   annotate("text",
          x = 1.5,
          y = max(dp_data$Nc, na.rm = TRUE) * 1.05,
          label = paste0("Wilcoxon p = ", format(pval, scientific = TRUE, digits = 3)),
          size = 5,
          fontface = "italic")

# Boxplot
dp_data <- rbind(
   data.frame(Nc = as.numeric(on_target$Bc), Group = "On-target"),
   data.frame(Nc = as.numeric(off_target$Bc), Group = "Off-target")
)

# Compute Wilcoxon p-value
wilcox_res <- wilcox.test(Nc ~ Group, data = dp_data)
pval <- wilcox_res$p.value

# Create the plot with scientific notation for p-value
ggplot(dp_data, aes(x = Group, y = Nc, fill = Group)) +
   geom_boxplot(outlier.shape = NA, alpha = 0.6) +
   geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
   theme_minimal() +
   theme(legend.position = "none") +
   ylab("Number of ALT reads") +
   xlab("") +
   annotate("text",
          x = 1.5,
          y = max(dp_data$Nc, na.rm = TRUE) * 1.05,
          label = paste0("Wilcoxon p = ", format(pval, scientific = TRUE, digits = 3)),
          size = 5,
          fontface = "italic")

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
dp_single_celltype$Main_Consequence <- sapply(
   strsplit(dp_single_celltype$Consequence, ","), 
   function(x) x[1]
)

# Assign severity rank
dp_single_celltype$Severity_Rank <- match(dp_single_celltype$Main_Consequence, vep_severity)

# Label on/off-target group
dp_single_celltype$Target_Group <- ifelse(
   dp_single_celltype$Gene_Class == "In List", 
   "On-target", 
   ifelse(dp_single_celltype$Gene_Class == "Not in List", "Off-target", NA)
)

# Plot histogram
library(ggplot2)
library(ggtext)

# Ensure proper factor ordering by severity
dp_single_celltype$Main_Consequence <- factor(
   sapply(strsplit(dp_single_celltype$Consequence, ","), function(x) x[1]),
   levels = rev(vep_severity)
)

bold_these <- c(
   "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
   "stop_gained", "frameshift_variant", "stop_lost", "start_lost", "transcript_amplification",
   "inframe_insertion", "inframe_deletion", "missense_variant",
   "protein_altering_variant", "splice_region_variant", "incomplete_terminal_codon_variant",
   "start_retained_variant", "stop_retained_variant", "synonymous_variant", "coding_sequence_variant", "mature_miRNA_variant"
)

dp_single_celltype$Main_Consequence_label <- sapply(
   sapply(strsplit(dp_single_celltype$Consequence, ","), function(x) x[1]),
   function(term) {
      if (term %in% bold_these) paste0("**", term, "**") else term
   }
)

dp_single_celltype$Main_Consequence_label <- factor(
   dp_single_celltype$Main_Consequence_label,
   levels = rev(sapply(vep_severity, function(term) {
      if (term %in% bold_these) paste0("**", term, "**") else term
   }))
)

# Plot with transcript_ablation on top
ggplot(
   subset(dp_single_celltype, Target_Group == "On-target"),
   aes(x = Main_Consequence_label)
) +
   geom_bar(fill = "#00BFC4") +
   scale_x_discrete(drop = FALSE) +
   coord_flip() +
   theme_minimal() +
   theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = ggtext::element_markdown()  # Enables bold text on y-axis
   ) +
   ylab("Count") +
   xlab("Consequence") +
   ggtitle("On-target")

# Plot with transcript_ablation on top
ggplot(
   subset(dp_single_celltype, Target_Group == "Off-target"),
   aes(x = Main_Consequence_label)
) +
   geom_bar(fill = "#F8766D") +
   scale_x_discrete(drop = FALSE) +
   coord_flip() +
   theme_minimal() +
   theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = ggtext::element_markdown()  # Enables bold text on y-axis
   ) +
   ylab("Count") +
   xlab("Consequence") +
   ggtitle("Off-target")

# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
# Filter on-target mutations in coding consequences
on_target_coding <- subset(dp_single_celltype,
                           Target_Group == "On-target" &
                            sapply(strsplit(Consequence, ","), `[`, 1) %in% bold_these
)
nrow(on_target_coding)  # Expect this to be zero

# -----------------------------------------------------------------------------
# dN/dS on off-target mutations
# -----------------------------------------------------------------------------
install.packages("remotes")
remotes::install_github("im3sanger/dndscv")

library(dndscv)

# Step 1: prepare mutation table
# Use only mutations from off-target genes and only coding consequences
df <- subset(dp_single_celltype, Target_Group == "Off-target" & Consequence %in% bold_these)
genes <- unique(df$SYMBOL)

df <- subset(dp_single_celltype, Target_Group == "Off-target" & SYMBOL %in% genes)

# Example input structure (you’ll need to build this from your current table):
mutations <- data.frame(
   sampleID = pd_id,
   chr = df$CHROM,
   pos = df$Start,
   ref = df$REF,
   mut = df$ALT
)
mutations$chr <- gsub("^chr", "", mutations$chr)

dndsout <- dndscv(mutations, refdb = "hg38")
dndsout$globaldnds

# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
library(ggplot2)

# Manually enter your dNdScv results
dnds_df <- data.frame(
   class = c("Missense", "Nonsense", "Splice", "Truncating", "All"),
   mle = c(0.9218095, 3.8325060, 1.5182295, 2.2202898, 1.0344325),
   cilow = c(0.3769065, 0.5937707, 0.2608016, 0.5510450, 0.4338408),
   cihigh = c(2.254492, 24.736994, 8.838216, 8.946068, 2.466459)
)

# Calculate error bars
dnds_df$lower_error <- dnds_df$mle - dnds_df$cilow
dnds_df$upper_error <- dnds_df$cihigh - dnds_df$mle

# Set severity-based factor levels
dnds_df$class <- factor(dnds_df$class, levels = rev(c("Truncating", "Nonsense", "Splice", "Missense", "All")))

# Plot with centered title
ggplot(dnds_df, aes(x = class, y = mle)) +
   geom_bar(stat = "identity", fill = "skyblue", color = "black") +
   geom_errorbar(aes(ymin = mle - lower_error, ymax = mle + upper_error), width = 0.3) +
   geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
   coord_flip() +
   labs(y = "dN/dS (MLE)", x = NULL, title = "Global dN/dS Estimates") +
   theme_minimal() +
   theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
celltypes <- c("Stage_0", "Stage_0A", "Stage_0B", "Stage_1", "Stage_2", "Stage_3", "19", "Leptotene", "Zygotene", "Pachytene", "Diplotene", "6",
               "Meiotic_division", "Early_spermatid", "Late_spermatid", "23", "Sertoli", "Leydig", "PMC", "Fibrotic_PMC", "Macrophage", "Endothelial")
info <- c("DP", "NC", "BC", "CC")










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
