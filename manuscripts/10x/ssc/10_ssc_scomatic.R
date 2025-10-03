# =============================================================================
# Manuscript   :
# Chapter      :
# Name         : 
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 15/04/25
# =============================================================================

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
wd.de.data  <- file.path("/Users/ty2/Work/sanger/ty2/SSC/ngs/pacbio/nf-scomatic/PD53626b")
#wd.de.plots <- file.path(wd.de.data, "10_mn7")
#dir.create(wd.de.plots, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Pe-processing
# -----------------------------------------------------------------------------
infile  <- file.path(wd.de.data, "mn7.hg38.bed.40.coding.unique.out")
outfile <- file.path(wd.de.data, "mn7.hg38.bed.40.coding.unique.out.common")

# Read as plain TSV; avoid factors and keep header as-is
mn7 <- read.table(infile,
						sep = "\t", header = TRUE, comment.char = "#",
						quote = "", fill = TRUE, check.names = FALSE,
						stringsAsFactors = FALSE)

# Coerce VAF safely to numeric (handles character columns)
VAF_num <- suppressWarnings(as.numeric(mn7$VAF))

# Keep rows with 0.5 <= VAF <= 1 and finite
keep <- is.finite(VAF_num) & VAF_num >= 0.5 & VAF_num <= 1
mn7_common <- mn7[keep, , drop = FALSE]

# Write result (TSV) with the same columns; NA as "."
write.table(mn7_common, outfile,
				sep = "\t", quote = FALSE, row.names = FALSE, na = ".")

cat(sprintf("[ok] Read %d rows, kept %d (VAF in [0.5, 1]).\nOutput: %s\n",
				nrow(mn7), nrow(mn7_common), outfile))

# -----------------------------------------------------------------------------
# Pe-processing
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Add VAF_SOMA & CCF_SOMA from 'soma' column and write filtered .common file
# -----------------------------------------------------------------------------

if (!exists("wd.de.data")) wd.de.data <- getwd()
infile  <- file.path(wd.de.data, "mn7.hg38.bed.40.coding.unique.out")
outfile <- file.path(wd.de.data, "mn7.hg38.bed.40.coding.unique.out.common")

# Read
mn7 <- read.table(infile, sep = "\t", header = TRUE, comment.char = "#",
						quote = "", fill = TRUE, check.names = FALSE,
						stringsAsFactors = FALSE)

stopifnot("soma" %in% names(mn7))

# Split 'soma' like "DP|NC|CC|BC|BQ|BCf|BCr" and extract needed fields
parts <- strsplit(mn7$soma, "\\|")

get_num <- function(lst, idx) {
	vapply(lst, function(x) {
		if (length(x) >= idx) suppressWarnings(as.numeric(x[idx])) else NA_real_
	}, numeric(1))
}

soma_DP <- get_num(parts, 1L)  # DP
soma_NC <- get_num(parts, 2L)  # NC
soma_CC <- get_num(parts, 3L)  # CC
soma_BC <- get_num(parts, 4L)  # BC

mn7$VAF_SOMA <- ifelse(is.finite(soma_DP) & soma_DP > 0, soma_BC / soma_DP, NA_real_)
mn7$CCF_SOMA <- ifelse(is.finite(soma_NC) & soma_NC > 0, soma_CC / soma_NC, NA_real_)

# Filter by overall VAF ∈ [0.5, 1] and write out (includes new columns)
VAF_num <- suppressWarnings(as.numeric(mn7$VAF))
keep <- is.finite(VAF_num) & VAF_num >= 0.5 & VAF_num <= 1

mn7_common <- mn7[keep, , drop = FALSE]
write.table(mn7_common, outfile, sep = "\t", quote = FALSE, row.names = FALSE, na = ".")

cat(sprintf("[ok] Augmented %d rows; wrote %d rows (VAF in [0.5,1]) to: %s\n",
				nrow(mn7), nrow(mn7_common), outfile))

# -----------------------------------------------------------------------------
# Write VAF==0 to *.zero, and VAF==0|NA to *.zero_or_na
# -----------------------------------------------------------------------------
VAF_num <- suppressWarnings(as.numeric(mn7$VAF))

# Exactly zero  (vectorized '&')
keep_zero <- is.finite(VAF_num) & VAF_num == 0
mn7_zero <- mn7[keep_zero, , drop = FALSE]
write.table(mn7_zero, out_zero, sep = "\t", quote = FALSE, row.names = FALSE, na = ".")

# Zero OR NA
keep_zero_or_na <- is.na(VAF_num) | (is.finite(VAF_num) & VAF_num == 0)
mn7_zero_or_na <- mn7[keep_zero_or_na, , drop = FALSE]
write.table(mn7_zero_or_na, out_zero_or_na, sep = "\t", quote = FALSE, row.names = FALSE, na = ".")

cat(sprintf("[ok] total=%d  VAF==0=%d  VAF==0|NA=%d\n",
				nrow(mn7), nrow(mn7_zero), nrow(mn7_zero_or_na)))
cat(sprintf("Wrote:\n  %s\n  %s\n", out_zero, out_zero_or_na))

# -----------------------------------------------------------------------------
# Build mn7_rare = mn7 - mn7_common - mn7_zero, and write to file
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({ library(dplyr) })

# Paths
if (!exists("wd.de.data")) wd.de.data <- getwd()
infile_all   <- file.path(wd.de.data, "mn7.hg38.bed.40.coding.unique.out")
infile_common <- file.path(wd.de.data, "mn7.hg38.bed.40.coding.unique.out.common")
infile_zero   <- file.path(wd.de.data, "mn7.hg38.bed.40.coding.unique.out.zero")
outfile_rare  <- file.path(wd.de.data, "mn7.hg38.bed.40.coding.unique.out.rare")

# Load mn7_common if not present
if (!exists("mn7_common")) {
	if (file.exists(infile_common)) {
		mn7_common <- read.table(infile_common, sep = "\t", header = TRUE, comment.char = "#",
										 quote = "", fill = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
	} else {
		mn7_common <- mn7[0, ]  # empty (no rows to subtract)
	}
}

# Load mn7_zero if not present
if (!exists("mn7_zero")) {
	if (file.exists(infile_zero)) {
		mn7_zero <- read.table(infile_zero, sep = "\t", header = TRUE, comment.char = "#",
									  quote = "", fill = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
	} else {
		mn7_zero <- mn7[0, ]
	}
}

# Use shared columns for reliable matching
by_cols <- Reduce(intersect, list(names(mn7), names(mn7_common), names(mn7_zero)))

# Subtractions: mn7 - mn7_common - mn7_zero
mn7_rare <- mn7 %>%
	dplyr::anti_join(mn7_common[, by_cols, drop = FALSE], by = by_cols) %>%
	dplyr::anti_join(mn7_zero[,   by_cols, drop = FALSE], by = by_cols)

# Save
write.table(mn7_rare, outfile_rare, sep = "\t", quote = FALSE, row.names = FALSE, na = ".")

cat(sprintf("[ok] mn7 rows: %d | common: %d | zero: %d | rare: %d\nOutput: %s\n",
				nrow(mn7), nrow(mn7_common), nrow(mn7_zero), nrow(mn7_rare), outfile_rare))

# -----------------------------------------------------------------------------
# Annotation
# -----------------------------------------------------------------------------
# MAF file
vcf <- subset(subset(subset(mn7_rare, Bc >= 5), Cc >=2), VAF_SOMA == 0)
vcf <- vcf[vcf$Bc != vcf$Cc,]

maf <- read.table(file.path("/Users/ty2/Work/sanger/ty2/SSC/ngs/pacbio/mn7/mn7.hg38.bed.40.coding.unique.vcf.maf"), sep = "\t", header = T, comment.char = "#", quote = "", fill = TRUE)

vcf$key <- paste0(gsub("^chr", "", vcf$CHROM), ":", vcf$Start)
maf$key <- maf$KEY
table(vcf$key %in% maf$key)

merged <- merge(vcf, maf,
					 by.x = "key",
					 by.y = "key",
					 all.x = TRUE)  # use all = TRUE if you want outer/full join

mn7_rare_snv_23b <- merged

intersect(mn7_rare_snv_23b$key, mn7_rare_snv_26b$key)
