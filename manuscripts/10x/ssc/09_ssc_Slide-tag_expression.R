# =============================================================================
# Name         : 
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 15/09/25
# =============================================================================

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
wd <- "/lustre/scratch125/casm/staging/team294/ty2"   ## ty2@farm22
#wd <- "/Users/ty2/Work/sanger/ty2"                   ## ty2@localhost
BASE <- "SSC"
base <- tolower(BASE)

wd.rna <- file.path(wd, BASE)

wd.de    <- file.path(wd.rna, "analysis", paste0(base, ""))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

load(file.path(wd.de.data, paste0("ssc_filtered_normalised_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_ordered_annotated.RData")))

so.integrated$cell_type <- Idents(so.integrated)

# -----------------------------------------------------------------------------
# After running SComatic genotyping
# # Last Modified: 17/09/25
# -----------------------------------------------------------------------------
wd.nf.data  <- file.path(wd.rna, "ngs/pacbio/SComatic/results/Step4_VariantCalling/new_beta")
barcodes <- read.table(
   file = file.path(wd.nf.data, "scomatics_PD53626b.tsv"),
   sep = "\t",
   header = T,
   stringsAsFactors = FALSE
)
barcodes <- subset(barcodes, Cell_type_observed == "germ")

meta <- read.table(
   file = file.path(wd.nf.data, "celltypes.tsv"),
   sep = "\t",
   header = T,
   stringsAsFactors = FALSE
)

barcodes.st <- subset(subset(barcodes, CHROM == "chr11"), Base_observed == "A")
meta.st <- subset(meta, Index %in% barcodes.st$CB)

# -----------------------------------------------------------------------------
# Standard Seurat re-processing workflow
# 01_QC
# https://satijalab.org/seurat/articles/pbmc3k_tutorial
# -----------------------------------------------------------------------------
library(Seurat)

so.integrated.st <- subset(so.integrated, sample.id == "PD53626b_ST1")

# -----------------------------------------------------------------------------
# MAST
# Last Modified: 19/09/25
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
	library(Seurat)
	library(dplyr)
	library(ggplot2)
	library(readr)
})

# -----------------------------
# Inputs
# -----------------------------
obj  <- so.integrated.st
gene <- "CBL"

# Your variant calls (exact cell IDs; no suffix changes here)
barcodes.ref <- subset(subset(barcodes, CHROM == "chr11"), Base_observed == "T")
barcodes.alt <- subset(subset(barcodes, CHROM == "chr11"), Base_observed == "A")

# -----------------------------------------------------------------------------
# To quantify per-cell allelic imbalance (AI) at chr6:31631400 (per cell type)
# Last Modified: 20/09/25
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(readr)
	library(Seurat)
	library(grid) # not strictly needed here, but harmless if already attached elsewhere
})

# -----------------------------
# CONFIG
# -----------------------------
locus_chr <- "chr11"
locus_pos <- 119284966          # if your table has Start column
gene      <- "CBL"

# choose the Seurat object you’re using
obj <- if (exists("so.integrated.st")) so.integrated.st else so.st
DefaultAssay(obj) <- "RNA"

# Helpers
get_cb <- function(df){
	if ("CB" %in% names(df))    return(as.character(df$CB))
	if ("Index" %in% names(df)) return(as.character(df$Index))
	if ("cell" %in% names(df))  return(as.character(df$cell))
	if (!is.null(rownames(df)))  return(rownames(df))
	stop("Couldn't find a barcode column (tried CB/Index/cell).")
}
has_col <- function(df, nm) nm %in% colnames(df)

# -----------------------------
# 1) ALT cells (from your prepared table)
# -----------------------------
alt_cells_use <- unique(get_cb(barcodes.alt))

# -----------------------------
# 2) Collapse reads per ALT cell at the locus using RAW 'barcodes'
# -----------------------------
# ensure we can address the barcode column as 'CB'
bar_raw <- barcodes
if (!"CB" %in% names(bar_raw) && "Index" %in% names(bar_raw)) {
	bar_raw$CB <- bar_raw$Index
}

# filter to locus (CHR + Start if present) and only the ALT cells of interest
bar_loc <- bar_raw %>%
	dplyr::filter(.data$CHROM == !!locus_chr) %>%
	{ if (has_col(., "Start")) dplyr::filter(., .data$Start == !!locus_pos) else . } %>%
	dplyr::filter(.data$CB %in% alt_cells_use) %>%
	dplyr::mutate(Base_observed = toupper(.data$Base_observed))

# collapse A/G read counts per cell
alt_ai <- bar_loc %>%
	dplyr::group_by(.data$CB) %>%
	dplyr::summarise(
		REF_reads   = sum(Num_reads[Base_observed == "T"], na.rm = TRUE),
		ALT_reads   = sum(Num_reads[Base_observed == "A"], na.rm = TRUE),
		total_reads = REF_reads + ALT_reads,
		.groups = "drop"
	) %>%
	dplyr::mutate(
		AF_alt = dplyr::if_else(total_reads > 0, ALT_reads / total_reads, NA_real_),
		AI     = abs(AF_alt - 0.5),
		log2_ALT_REF = dplyr::case_when(
			ALT_reads > 0 & REF_reads > 0 ~ log2(ALT_reads / REF_reads),
			ALT_reads == 0 & REF_reads > 0 ~ -Inf,
			ALT_reads > 0 & REF_reads == 0 ~  Inf,
			TRUE ~ NA_real_
		)
	)

# keep only those ALT cells that actually exist in the Seurat object
alt_ai <- alt_ai %>% dplyr::filter(.data$CB %in% colnames(obj))

# -----------------------------
# 3) Pull PRRC2A expression for these ALT cells
# -----------------------------
# case-insensitive rescue for feature name
feat <- rownames(obj[["RNA"]])
if (!(gene %in% feat)) {
	cand <- feat[toupper(feat) == toupper(gene)]
	if (length(cand) == 1) gene <- cand else stop("Feature ", gene, " not found in RNA assay.")
}

expr_df <- FetchData(obj, vars = gene, cells = alt_ai$CB, slot = "data")
colnames(expr_df) <- "expr"
alt_ai <- dplyr::bind_cols(alt_ai, expr_df[alt_ai$CB, , drop = FALSE])

alt_ai$cell_type <- meta.st$Cell_type[ match(alt_ai$CB, meta.st$Index) ]

# -----------------------------
# AF
# -----------------------------
#alt_ai <- subset(alt_ai, REF_reads != 0)

suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

# alt_ai expected to have: CB, REF_reads, ALT_reads, total_reads, expr (PRRC2A)
stopifnot(all(c("ALT_reads","REF_reads","total_reads","expr") %in% names(alt_ai)))

# Pseudocount (Jeffreys 0.5 is standard; you can try 1.0 if depth is tiny)
eps <- 0.5

alt_ai <- alt_ai %>%
	mutate(
		AF        = ifelse(total_reads > 0, ALT_reads / total_reads, NA_real_),
		AI        = abs(AF - 0.5),
		log2_AR   = log2((ALT_reads + eps) / (REF_reads + eps)),
		logit_AF  = {
			p <- (ALT_reads + eps) / (total_reads + 2*eps)     # smoothed AF in (0,1)
			log(p / (1 - p))
		}
	)

# -----------------------------------------------------------------------------
# Expression vs. AF and total reads
# Last Modified: 20/09/25
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

# Spearman tests
cor_AF    <- suppressWarnings(stats::cor.test(alt_ai$expr, alt_ai$AF,          method = "spearman"))
cor_depth <- suppressWarnings(stats::cor.test(alt_ai$expr, alt_ai$total_reads, method = "spearman"))

rho_AF <- as.numeric(cor_AF$estimate)
p_AF   <- cor_AF$p.value
lab_AF <- paste0("Spearman \u03C1 = ", signif(rho_AF, 3), "\nP = ", signif(p_AF, 3))

rho_D  <- as.numeric(cor_depth$estimate)
p_D    <- cor_depth$p.value
lab_D  <- paste0("Spearman \u03C1 = ", signif(rho_D, 3), "\nP = ", signif(p_D, 3))

# Plots with top-right annotation
p1 <- ggplot(alt_ai, aes(x = AF, y = expr)) +
	geom_point(size = 3, alpha = 0.85, color = "red") +
	geom_smooth(method = "lm", se = FALSE, linetype = 2, color = "black") +
	labs(x = "ALT allele fraction (AF)", y = "CBL expression") +
	annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5, label = lab_AF, size = 5) +
	theme_classic(base_size = 14)

p2 <- ggplot(alt_ai, aes(x = total_reads, y = expr)) +
	geom_point(size = 3, alpha = 0.85, color = "blue") +
	geom_smooth(method = "lm", se = FALSE, linetype = 2, color = "black") +
	labs(x = "Total reads (REF + ALT)", y = "CBL expression") +
	annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5, label = lab_D, size = 5) +
	theme_classic(base_size = 14)

print(p1); print(p2)

# (optional) save
ggsave(file.path(wd.nf.data, "CBL_expr_vs_ALT_AF.png"),  p1, width = 6, height = 5, dpi = 300)
ggsave(file.path(wd.nf.data, "CBL_expr_vs_depth.png"), p2, width = 6, height = 5, dpi = 300)

write.table(
	alt_ai,
	file = file.path(wd.nf.data, "CBL_expr_vs_depth.txt"),
	sep = "\t",
	quote = FALSE,
	row.names = FALSE,
	col.names = TRUE
)






# -----------------------------------------------------------------------------
# For each cell type
# Last Modified: 06/10/25
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

for (ct in sort(unique(alt_ai$cell_type))) {
	df  <- dplyr::filter(alt_ai, cell_type == ct)
	AF_col <- if ("AF" %in% names(df)) "AF" else "AF_alt"
	safe   <- gsub("[^A-Za-z0-9._-]+", "_", ct)
	
	# keep only finite x–y pairs for each test
	df_AF <- dplyr::filter(df, is.finite(expr), is.finite(.data[[AF_col]]))
	df_D  <- dplyr::filter(df, is.finite(expr), is.finite(total_reads))
	
	n_AF <- nrow(df_AF); n_D <- nrow(df_D)
	
	if (n_AF >= 3 && stats::sd(df_AF$expr) > 0 && stats::sd(df_AF[[AF_col]]) > 0) {
		cor_AF <- suppressWarnings(stats::cor.test(df_AF$expr, df_AF[[AF_col]], method = "spearman"))
		lab_AF <- paste0("Spearman \u03C1 = ", signif(unname(cor_AF$estimate), 3),
							  "\nP = ", signif(cor_AF$p.value, 3))
	} else {
		lab_AF <- paste0("Spearman \u03C1 = NA\nP = NA (n=", n_AF, ")")
	}
	
	if (n_D >= 3 && stats::sd(df_D$expr) > 0 && stats::sd(df_D$total_reads) > 0) {
		cor_depth <- suppressWarnings(stats::cor.test(df_D$expr, df_D$total_reads, method = "spearman"))
		lab_D <- paste0("Spearman \u03C1 = ", signif(unname(cor_depth$estimate), 3),
							 "\nP = ", signif(cor_depth$p.value, 3))
	} else {
		lab_D <- paste0("Spearman \u03C1 = NA\nP = NA (n=", n_D, ")")
	}
	
	p1 <- ggplot(df, aes_string(x = AF_col, y = "expr")) +
		geom_point(size = 3, alpha = 0.85, color = "red") +
		labs(title = ct,
			  x = if (AF_col == "AF") "ALT allele fraction (AF)" else "ALT allele fraction (AF_alt)",
			  y = "CBL expression") +
		annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5, label = lab_AF, size = 5) +
		theme_classic(base_size = 14)
	
	p2 <- ggplot(df, aes(x = total_reads, y = expr)) +
		geom_point(size = 3, alpha = 0.85, color = "blue") +
		labs(title = ct, x = "Total reads (REF + ALT)", y = "CBL expression") +
		annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5, label = lab_D, size = 5) +
		theme_classic(base_size = 14)
	
	print(p1); print(p2)
	
	ggsave(file.path(wd.nf.data, paste0("CBL_", safe, "_expr_vs_AF.png")), p1, width = 6, height = 5, dpi = 300)
	ggsave(file.path(wd.nf.data, paste0("CBL_", safe, "_expr_vs_depth.png")), p2, width = 6, height = 5, dpi = 300)
	write.table(df, file = file.path(wd.nf.data, paste0("CBL_", safe, "_expr_table.txt")),
					sep = "\t", quote = FALSE, row.names = FALSE)
}










# -----------------------------------------------------------------------------
# 
# Last Modified: 18/09/25
# -----------------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(readr)

# -----------------------------
# 0) Inputs
# -----------------------------
obj <- so.integrated.st           # Seurat with RNA/SCT/integrated
gene <- "PRRC2A"                  # target gene (case-insensitive handled)
# barcodes.ref / barcodes.alt come from your snippet

# Helper to extract a barcode vector from your data frames
get_cb <- function(df){
	if ("CB" %in% colnames(df)) return(as.character(df$CB))
	if ("Index" %in% colnames(df)) return(as.character(df$Index))
	if ("cell" %in% colnames(df)) return(as.character(df$cell))
	if (!is.null(rownames(df)))    return(rownames(df))
	stop("Couldn't find a barcode column in the input (looked for CB/Index/cell).")
}

# 1) Take only the meta cells and keep only those with coordinates
meta_cells_only <- spatial_data %>%
	dplyr::filter(.data$meta_flag == "yes") %>%
	dplyr::filter(complete.cases(.data$plot_x, .data$plot_y))

# 2) Desired biological order for cell types
ct_order_master <- c(
	"Stage 0","Stage 0A","Stage 0B","Stage 1","Stage 2","Stage 3",
	"Leptotene","Zygotene","Pachytene","Diplotene","Meiotic division",
	"Early spermatid","Late spermatid"
)

# 2) Keep ONLY the types present in the 10 meta cells, in biological order
ct_levels <- ct_order_master[ct_order_master %in% meta_cells_only$cell_type]  # preserves order
# (equivalent to: base::intersect(ct_order_master, unique(meta_cells_only$cell_type)))

# Relevel & drop unused types
meta_cells_only$cell_type <- factor(meta_cells_only$cell_type, levels = ct_levels)

# 3) Palette restricted to these 7 (or however many) types
pal_meta <- setNames(scales::hue_pal()(length(ct_levels)), ct_levels)

# 4) Plot only these meta cells, colored by their original cell_type
p_meta <- ggplot(meta_cells_only, aes(x = plot_x, y = plot_y, color = cell_type)) +
	geom_point(shape = 19, size = 4, alpha = 1) +
	scale_color_manual(values = pal_meta, breaks = ct_levels, name = NULL) +  # no legend title
	coord_fixed() +
	labs(x = "X (µm)", y = "Y (µm)") +
	theme_classic() +
	theme(
		axis.title = element_text(size = 16),
		axis.text  = element_text(size = 14),
		legend.text = element_text(size = 14),
		legend.position = c(0.05, 0.05),
		legend.justification = c("left", "bottom"),
		legend.background = element_rect(fill = "white", color = "grey80"),
		legend.margin = margin(3, 6, 3, 6)
	)

print(p_meta)

# Save
outfile <- file.path(wd.de.plots, "Spatial_ST_META_ONLY_7types_origcolors.png")
ggsave(outfile, plot = p_meta, width = 6, height = 6, dpi = 300)
cat("Saved:", outfile, "\n")

# -----------------------------------------------------------------------------
# Differential expression analysis
# Last Modified: 19/09/25
# -----------------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(readr)

# -----------------------------
# 0) Inputs
# -----------------------------
obj <- so.integrated.st           # Seurat with RNA/SCT/integrated
gene <- "PRRC2A"                  # target gene (case-insensitive handled)

barcodes.ref <- subset(subset(barcodes, CHROM == "chr6"), Base_observed == "A")
barcodes.alt <- subset(subset(barcodes, CHROM == "chr6"), Base_observed == "G")

## Keep only cells from "Zygotene","Pachytene","Diplotene","Meiotic division"
## Keep only cells from these four types
keep_types <- c("Zygotene","Pachytene","Diplotene","Meiotic division")

# Helper to extract a barcode column without changing suffixes
get_cb <- function(df){
	if ("CB" %in% colnames(df))    return(as.character(df$CB))
	if ("Index" %in% colnames(df)) return(as.character(df$Index))
	if ("cell" %in% colnames(df))  return(as.character(df$cell))
	if (!is.null(rownames(df)))    return(rownames(df))
	stop("Couldn't find a barcode column (looked for CB/Index/cell).")
}

# Exact cell IDs present in the Seurat object
cells_obj <- colnames(obj)

# Allowed cells by cell_type (exact IDs, no suffix stripping)
ct_vec <- obj$cell_type
names(ct_vec) <- cells_obj
allowed_cells <- cells_obj[!is.na(ct_vec) & ct_vec %in% keep_types]

# Filter barcodes.ref / barcodes.alt *in place* (preserve variable names)
ref_ids <- get_cb(barcodes.ref)
alt_ids <- get_cb(barcodes.alt)

barcodes.ref <- barcodes.ref[ref_ids %in% allowed_cells, , drop = FALSE]
barcodes.alt <- barcodes.alt[alt_ids %in% allowed_cells, , drop = FALSE]

# -----------------------------
# 1) REF / ALT cell IDs (EXACT MATCH; NO SUFFIX CHANGES)
# -----------------------------
ref_cells_raw <- unique(get_cb(barcodes.ref))
alt_cells_raw <- unique(get_cb(barcodes.alt))

cells_obj <- colnames(obj)

# exact intersection only — disambiguate to base::intersect
ref_cells <- base::intersect(ref_cells_raw, cells_obj)
alt_cells <- base::intersect(alt_cells_raw, cells_obj)

message(sprintf("Matched cells — REF: %d / %d, ALT: %d / %d",
					 length(ref_cells), length(ref_cells_raw),
					 length(alt_cells), length(alt_cells_raw)))
stopifnot(length(ref_cells) > 0, length(alt_cells) > 0)

# -----------------------------
# 2) Fetch PRRC2A from RNA assay (log-normalized 'data' slot)
# -----------------------------
DefaultAssay(obj) <- "RNA"
gene <- "PRRC2A"
rna_feats <- rownames(obj[["RNA"]])
if (!(gene %in% rna_feats)) {
	cand <- rna_feats[toupper(rna_feats) == toupper(gene)]
	if (length(cand) == 1) gene <- cand else stop("Gene 'PRRC2A' not found in RNA assay.")
}

expr_ref <- FetchData(obj, vars = gene, cells = ref_cells, slot = "data")[,1]
expr_alt <- FetchData(obj, vars = gene, cells = alt_cells, slot = "data")[,1]

df <- dplyr::bind_rows(
	data.frame(cell = ref_cells, expr = expr_ref, group = "REF (A)"),
	data.frame(cell = alt_cells, expr = expr_alt, group = "ALT (G)")
) |>
	dplyr::filter(is.finite(expr)) |>
	dplyr::mutate(group = factor(group, levels = c("REF (A)", "ALT (G)")))  # REF left

# -----------------------------
# 3) Wilcoxon + plot + save
# -----------------------------
w <- wilcox.test(expr ~ group, data = df, exact = FALSE)
pval <- w$p.value
p_label <- if (pval < 1e-4) "p < 1e-4" else paste0("p = ", signif(pval, 3))

stats_tbl <- df |>
	dplyr::group_by(group) |>
	dplyr::summarise(n = dplyr::n(),
						  median = median(expr, na.rm = TRUE),
						  mean = mean(expr, na.rm = TRUE),
						  .groups = "drop") |>
	dplyr::mutate(gene = gene, assay = "RNA", slot = "data", p_value = pval)

p <- ggplot(df, aes(x = group, y = expr, fill = group)) +
	geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.65) +
	geom_jitter(width = 0.15, size = 1.8, alpha = 0.75, aes(color = group), show.legend = FALSE) +
	scale_fill_manual(values = c("REF (A)" = "#9ecae1", "ALT (G)" = "#fb6a4a")) +
	scale_color_manual(values = c("REF (A)" = "#3182bd", "ALT (G)" = "#cb181d")) +
	labs(x = NULL, y = paste0(gene, " expression")) +
	annotate("text",
				x = 1.5,
				y = max(df$expr, na.rm = TRUE) + diff(range(df$expr, na.rm = TRUE))*0.08,
				label = p_label, size = 5) +
	theme_classic() +
	theme(legend.position = "none",
			axis.text = element_text(size = 12),
			axis.title = element_text(size = 13),
			plot.margin = margin(10,10,10,10))

print(p)

readr::write_csv(stats_tbl, file.path(wd.de.plots, "PRRC2A_REF_ALT_wilcox_stats_RNA_exactIDs_G2M.csv"))
ggsave(file.path(wd.de.plots, "PRRC2A_REF_vs_ALT_boxplot_RNA_exactIDs_G2M.png"),
		 plot = p, width = 4.8, height = 4.8, dpi = 300)

# -----------------------------------------------------------------------------
# MAST
# Last Modified: 19/09/25
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
	library(Seurat)
	library(dplyr)
	library(ggplot2)
	library(readr)
})

# -----------------------------
# Inputs
# -----------------------------
obj  <- so.integrated.st
gene <- "DHX9"

# Your variant calls (exact cell IDs; no suffix changes here)
barcodes.ref <- subset(subset(barcodes, CHROM == "chr1"), Base_observed == "T")
barcodes.alt <- subset(subset(barcodes, CHROM == "chr1"), Base_observed == "C")

# Helper to extract a barcode column robustly
get_cb <- function(df){
	if ("CB" %in% colnames(df))    return(as.character(df$CB))
	if ("Index" %in% colnames(df)) return(as.character(df$Index))
	if ("cell" %in% colnames(df))  return(as.character(df$cell))
	if (!is.null(rownames(df)))    return(rownames(df))
	stop("Couldn't find a barcode column (looked for CB/Index/cell).")
}

# REF/ALT (exact matches only)
ref_cells_use <- unique(get_cb(barcodes.ref))
alt_cells_use <- unique(get_cb(barcodes.alt))

# Keep only those that actually exist in the object
ref_cells_use <- base::intersect(ref_cells_use, colnames(obj))
alt_cells_use <- base::intersect(alt_cells_use, colnames(obj))

message(sprintf("Using REF=%d ALT=%d cells for MAST.", length(ref_cells_use), length(alt_cells_use)))
stopifnot(length(ref_cells_use) > 0, length(alt_cells_use) > 0)

# -----------------------------
# Subset the object to REF+ALT only (recommended)
# -----------------------------
obj_sub <- subset(obj, cells = c(ref_cells_use, alt_cells_use))
DefaultAssay(obj_sub) <- "RNA"

# Grouping factor
obj_sub$genotype_grp <- NA_character_
obj_sub$genotype_grp[colnames(obj_sub) %in% ref_cells_use] <- "REF"
obj_sub$genotype_grp[colnames(obj_sub) %in% alt_cells_use] <- "ALT"
Idents(obj_sub) <- factor(obj_sub$genotype_grp, levels = c("REF","ALT"))



# -----------------------------
# Covariates present (use only what exists)
# -----------------------------
cand_covars <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "pct.mt", "percent_mt", "cell_type")
latent_covars <- base::intersect(cand_covars, colnames(obj_sub@meta.data))

# Optional: if groups are tiny, you may want to drop 'cell_type' to avoid overfitting:
# latent_covars <- setdiff(latent_covars, "cell_type")

# -----------------------------
# Run MAST (normalized 'data' slot)
# -----------------------------
if (!requireNamespace("MAST", quietly = TRUE)) {
	stop("Package 'MAST' is required (install.packages('MAST')).")
}

mast_res <- FindMarkers(
	object          = obj_sub,
	ident.1         = "REF",
	ident.2         = "ALT",
	test.use        = "MAST",
	assay           = "RNA",
	slot            = "data",
	logfc.threshold = 0,
	min.pct         = 0,
	latent.vars     = latent_covars
)

mast_res$gene <- rownames(mast_res)

# Pull PRRC2A row (case-insensitive)
prrc2a_row <- mast_res %>% dplyr::filter(toupper(.data$gene) == toupper(gene))
if (nrow(prrc2a_row) == 0) stop("PRRC2A not found in MAST results (check feature naming).")

print(prrc2a_row[, c("gene","avg_log2FC","p_val","p_val_adj")])

# -----------------------------
# Boxplot of PRRC2A expression (REF left, ALT right)
# -----------------------------
# Fetch expression from the same subset used in MAST
expr_ref <- FetchData(obj_sub, vars = gene, cells = ref_cells_use, slot = "data")[,1]
expr_alt <- FetchData(obj_sub, vars = gene, cells = alt_cells_use, slot = "data")[,1]

df <- dplyr::bind_rows(
	data.frame(cell = ref_cells_use, expr = expr_ref, group = "REF (A)"),
	data.frame(cell = alt_cells_use, expr = expr_alt, group = "ALT (G)")
)
df$group <- factor(df$group, levels = c("REF (A)","ALT (G)"))

# p-value label from MAST (two-sided)
pval <- prrc2a_row$p_val[1]
p_label <- ifelse(pval < 1e-4, "MAST p < 1e-4", paste0("MAST p = ", signif(pval, 3)))

p <- ggplot(df, aes(x = group, y = expr, fill = group)) +
	geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.6) +
	geom_jitter(width = 0.15, size = 1.8, alpha = 0.7, aes(color = group), show.legend = FALSE) +
	scale_fill_manual(values = c("REF (A)" = "#9ecae1", "ALT (G)" = "#fb6a4a")) +
	scale_color_manual(values = c("REF (A)" = "#3182bd", "ALT (G)" = "#cb181d")) +
	labs(x = NULL, y = paste0(gene, " expression (normalized)")) +
	annotate("text", x = 1.5, y = max(df$expr, na.rm = TRUE) * 1.05, label = p_label, size = 5) +
	theme_classic() +
	theme(
		legend.position = "none",
		axis.text = element_text(size = 12),
		axis.title = element_text(size = 13),
		plot.margin = margin(10, 10, 10, 10)
	)

print(p)

# -----------------------------
# Save outputs
# -----------------------------
out_plot <- file.path(wd.de.plots, "PRRC2A_REF_vs_ALT_boxplot_MAST.png")
out_csv  <- file.path(wd.de.plots, "PRRC2A_MAST_results_full.csv")
ggsave(out_plot, plot = p, width = 4.8, height = 4.8, dpi = 300)
readr::write_csv(mast_res, out_csv)
cat("Saved plot:", out_plot, "\nSaved MAST table:", out_csv, "\n")













# -----------------------------------------------------------------------------
# Plot cell-type-by-cell-type - Each cell type highlighted separately in a loop
# Last Modified: 15/09/25
# -----------------------------------------------------------------------------
# --- mark meta indices
meta_idx <- unique(meta.st$Index)
spatial_data$meta_flag    <- ifelse(rownames(spatial_data) %in% meta_idx, "yes", "no")
spatial_data$variant_flag <- ifelse(spatial_data$meta_flag == "yes", "6_31631400_A_G", NA)

# palette for cell types
all_ct <- sort(unique(na.omit(spatial_data$cell_type)))
pal <- setNames(scales::hue_pal()(length(all_ct)), all_ct)

for (ct in all_ct) {
 target_cells <- dplyr::filter(spatial_data, .data$cell_type == .env$ct)
 other_cells  <- dplyr::filter(spatial_data, .data$cell_type != .env$ct | is.na(.data$cell_type))
 
 # only meta cells for THIS cell type (+ valid coords)
 meta_cells <- dplyr::filter(
  spatial_data,
  .data$meta_flag == "yes",
  .data$cell_type == .env$ct
 ) %>% dplyr::filter(complete.cases(.data[[spatial_coords[1]]], .data[[spatial_coords[2]]]))
 
 have_variant <- nrow(meta_cells) > 0
 
 # legend bits that must match the number of legend items
 legend_breaks <- if (have_variant) c(ct, "6_31631400_A_G") else ct
 legend_values <- c(setNames(pal[ct], ct), `6_31631400_A_G` = "red")
 override_list <- list(
  size   = rep(3, length(legend_breaks)),
  shape  = if (have_variant) c(16, 21) else 16,  # dot for ct, ring for variant
  fill   = rep(NA, length(legend_breaks)),
  stroke = if (have_variant) c(0, 1.2) else 0,
  alpha  = 1
 )
 
 dim_plot <- ggplot() +
  # background (others): grey hollow
  geom_point(
   data = other_cells,
   aes(x = .data[[spatial_coords[1]]], y = .data[[spatial_coords[2]]]),
   color = "#D3D3D3", fill = "white", size = 3, alpha = 1.0, shape = 21
  ) +
  # highlighted target cells
  geom_point(
   data = target_cells,
   aes(x = .data[[spatial_coords[1]]], y = .data[[spatial_coords[2]]], color = cell_type),
   size = 3, alpha = 1.0
  ) +
  # meta ring ONLY for this cell type (no fill)
  geom_point(
   data = meta_cells,
   aes(x = .data[[spatial_coords[1]]], y = .data[[spatial_coords[2]]], color = variant_flag),
   shape = 21, fill = NA, size = 3, stroke = 1.0, alpha = 1
  ) +
  geom_text_repel(
   data = cluster_centers,
   aes(x = coord1, y = coord2, label = seurat_clusters),
   size = 4, color = "black"
  ) +
  # merged legend
  scale_color_manual(
   values = legend_values,
   breaks = legend_breaks,
   na.translate = FALSE
  ) +
  coord_fixed() +
  labs(x = "X", y = "Y") +
  theme_classic() +
  theme(
   axis.title = element_text(size = 18),
   axis.text  = element_text(size = 16),
   legend.text  = element_text(size = 16),
   legend.title = element_text(size = 16),
   legend.position = c(0.98, 0.98),
   legend.justification = c("right", "top"),
   legend.background = element_rect(fill = "white", color = "grey80"),
   legend.margin = margin(3, 6, 3, 6)
  ) +
  guides(color = guide_legend(override.aes = override_list))
 
 safe_cell_type <- gsub("[^A-Za-z0-9]", "_", ct)
 filename <- paste0("Spatial_ST_", safe_cell_type, "_highlighted_dbscan.png")
 ggsave(file.path(wd.de.plots, filename), plot = dim_plot, width = 6, height = 6, dpi = 300)
 cat("Saved plot for", ct, "\n")
}
