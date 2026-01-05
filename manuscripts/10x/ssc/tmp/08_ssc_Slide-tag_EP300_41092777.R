# =============================================================================
# Name         : 
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 15/09/25
# =============================================================================
CHR  <- "chr22"
POS  <- 41092777          # if your table has Start column
REF  <- "T"
ALT  <- "C"
GENE <- "EP300"
MUT <- paste0(gsub("chr", "", CHR), "_", POS, "_", REF, "_", ALT)

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
wd.de.plots <- file.path(wd.de, "plots", GENE, MUT)
#dir.create(file.path(wd.de.plots, "origUM"), showWarnings = FALSE)
dest <- file.path(wd.de.plots, "origUM")
if (!dir.exists(dest)) {
	ok <- base::dir.create(dest, recursive = TRUE, showWarnings = FALSE, mode = "0775")
	if (!ok && !dir.exists(dest)) stop("Failed to create: ", dest)
}

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

barcodes.st <- subset(subset(subset(barcodes, CHROM == CHR), Start == POS), Base_observed == ALT)
meta.st <- subset(meta, Index %in% barcodes.st$CB)

# -----------------------------------------------------------------------------
# Standard Seurat re-processing workflow
# 01_QC
# https://satijalab.org/seurat/articles/pbmc3k_tutorial
# -----------------------------------------------------------------------------
library(Seurat)
library(qs)

so <- qread(file.path(wd.rna, "ngs/slide-tag/pipeline_outputs/PD53626B/Positions/seurat.qs"))
so.integrated.st <- subset(so.integrated, sample.id == "PD53626b_ST1")

# Explore
Reductions(so)
#> "pca" "umap" "dbscan" "kde" "spatial"

# Plot spatial embedding
p <- DimPlot(
   so, reduction = "spatial",
   group.by = "orig.ident",   # or "seurat_clusters", "Phase", etc.
   pt.size = 0.2
) +
   coord_fixed() +            # keep aspect ratio 1:1
   labs(x = "X", y = "Y")

# Many spatial pipelines have Y increasing downward; if upside-down, flip Y:
# p <- p + scale_y_reverse()

# -----------------------------------------------------------------------------
# Get so.st
# Last Modified: 17/09/25
# -----------------------------------------------------------------------------
# Step 1: Add helper cols (no need to rename colnames)
so.integrated.st$orig_cellname <- colnames(so.integrated.st)
so.integrated.st$temp_cellname <- sub("_[0-9]+$", "", colnames(so.integrated.st))
colnames(so.integrated.st) <- so.integrated.st$temp_cellname

# Use base::intersect to avoid the conflict
overlapping_cells <- base::intersect(colnames(so), colnames(so.integrated.st))
so.st <- subset(so, cells = overlapping_cells)

# Step 5: Add cell_type information from so.integrated.st
cell_type_mapping <- so.integrated.st$cell_type[match(colnames(so.st), colnames(so.integrated.st))]
so.st$cell_type <- cell_type_mapping

# Step 6: Add original suffix back to so.st colnames using orig_cellname
# Get the original names with suffix for the overlapping cells
orig_names_for_overlapping <- so.integrated.st$orig_cellname[match(colnames(so.st), colnames(so.integrated.st))]
colnames(so.st) <- orig_names_for_overlapping

# Step 7: Restore original colnames to so.integrated.st and remove temp column
colnames(so.integrated.st) <- so.integrated.st$orig_cellname

# -----------------------------------------------------------------------------
# Build plotting frame using ORIGINAL slide-tag coordinates (µm), not embeddings
# Last Modified: 2025-09-18
# -----------------------------------------------------------------------------
meta_idx <- unique(meta.st$Index)

suppressPackageStartupMessages({
	library(Seurat)
	library(dplyr)
	library(readr)
	# Optional but handy if you have Bioc packages loaded that mask verbs:
	if (requireNamespace("conflicted", quietly = TRUE)) {
		conflicted::conflict_prefer("select", "dplyr", quiet = TRUE)
		conflicted::conflict_prefer("filter", "dplyr", quiet = TRUE)
		conflicted::conflict_prefer("mutate", "dplyr", quiet = TRUE)
		conflicted::conflict_prefer("lag",    "dplyr", quiet = TRUE)
	}
})

# --- inputs -------------------------------------------------------------------
path_matrix    <- file.path(wd.rna, "ngs/slide-tag/pipeline_outputs/PD53626B/Positions/matrix.csv.gz")
path_whitelist <- file.path(wd.rna, "ngs/slide-tag/pipeline_outputs/PD53626B/Positions/cb_whitelist.txt")
flip_y         <- FALSE  # set TRUE if Y is inverted in image space

# --- helpers ------------------------------------------------------------------
canonicalize_barcode <- function(x) {
	x <- toupper(x)
	x <- sub("_[0-9]+$", "", x)  # drop slice suffix like _19
	x <- sub("-1$", "", x)       # drop trailing "-1" if present
	trimws(x)
}
rng <- function(v) paste0(round(range(v, na.rm = TRUE), 2), collapse = " .. ")

# --- load matrix --------------------------------------------------------------
df <- readr::read_csv(path_matrix, show_col_types = FALSE)
stopifnot(all(c("cb_index","x_um","y_um","umi") %in% names(df)))

# --- aggregate per-cell (UMI-weighted) coordinates ---------------------------
agg <- df %>%
	dplyr::group_by(cb_index) %>%
	dplyr::summarise(
		plot_x    = sum(x_um * umi, na.rm = TRUE) / sum(umi, na.rm = TRUE),
		plot_y    = sum(y_um * umi, na.rm = TRUE) / sum(umi, na.rm = TRUE),
		umi_total = sum(umi,          na.rm = TRUE),
		.groups   = "drop"
	)

if (isTRUE(flip_y)) agg$plot_y <- -agg$plot_y

# --- load whitelist; add both 1-based and 0-based indices ---------------------
# If cb_whitelist.txt is one barcode per line:
# Option A: build then mutate (uses dplyr::row_number)
wl <- tibble::tibble(barcode = readr::read_lines(path_whitelist)) %>%
	dplyr::mutate(
		cb_index_1 = dplyr::row_number(),
		cb_index_0 = cb_index_1 - 1L
	)

# --- build Seurat mapping table (string-only) ---------------------------------
stopifnot(exists("so.st"), inherits(so.st, "Seurat"))
map_seu <- tibble::tibble(
	cell         = colnames(so.st),
	base_barcode = canonicalize_barcode(colnames(so.st))
)

# --- try both index bases, keep the one with more whitelist matches -----------
cand1 <- agg %>%
	dplyr::inner_join(wl, by = c("cb_index" = "cb_index_1")) %>%
	dplyr::mutate(index_base = 1L)

cand0 <- agg %>%
	dplyr::inner_join(wl, by = c("cb_index" = "cb_index_0")) %>%
	dplyr::mutate(index_base = 0L)

agg_wl <- if (nrow(cand0) > nrow(cand1)) cand0 else cand1
message(sprintf(
	"Whitelist join candidates: 1-based=%d, 0-based=%d. Using %s-based.",
	nrow(cand1), nrow(cand0),
	if (identical(agg_wl, cand0)) "0" else "1"
))

# --- map to Seurat cells (canonicalize only; coordinates remain in µm) --------
coords_from_file <- agg_wl %>%
	dplyr::mutate(base_barcode = canonicalize_barcode(.data$barcode)) %>%
	dplyr::inner_join(map_seu, by = "base_barcode") %>%
	dplyr::select(dplyr::any_of(c("cell","plot_x","plot_y","umi_total","index_base")))

if (nrow(coords_from_file) == 0) {
	cat("No matches after canonicalization.\nExamples:\n",
		 "  whitelist: ", paste(utils::head(wl$barcode, 3), collapse = ", "), "\n",
		 "  Seurat:    ", paste(utils::head(colnames(so.st), 3), collapse = ", "), "\n", sep = "")
	stop("Check sample/slice suffixes and 10x '-1' formatting.")
}

# --- if duplicates per cell, keep highest UMI ---------------------------------
coords_from_file <- coords_from_file %>%
	dplyr::group_by(.data$cell) %>%
	dplyr::arrange(dplyr::desc(.data$umi_total), .by_group = TRUE) %>%
	dplyr::slice(1L) %>%
	dplyr::ungroup()

message(sprintf("Matched %d barcodes to Seurat cells.", nrow(coords_from_file)))

# --- attach to so.st meta (µm) and create a µm-space reduction "orig_um" ------
so.st$orig_x_um      <- NA_real_
so.st$orig_y_um      <- NA_real_
so.st$orig_umi_total <- NA_real_

m   <- match(colnames(so.st), coords_from_file$cell)
hit <- !is.na(m)

so.st$orig_x_um[hit]      <- coords_from_file$plot_x[m[hit]]
so.st$orig_y_um[hit]      <- coords_from_file$plot_y[m[hit]]
so.st$orig_umi_total[hit] <- coords_from_file$umi_total[m[hit]]

message(sprintf("Annotated %d / %d cells with ORIGINAL µm coords.",
					 sum(hit), ncol(so.st)))

# Build the embeddings matrix (as you did)
orig_um_mat <- matrix(NA_real_, nrow = ncol(so.st), ncol = 2,
							 dimnames = list(colnames(so.st), c("orig_x_um","orig_y_um")))
orig_um_mat[hit, 1] <- so.st$orig_x_um[hit]
orig_um_mat[hit, 2] <- so.st$orig_y_um[hit]

# Use an alphanumeric key and set column names to <key><dimension>
key <- "origUM_"                                # must be [A-Za-z0-9]+_
colnames(orig_um_mat) <- paste0(key, 1:ncol(orig_um_mat))

so.st[["orig_um"]] <- CreateDimReducObject(
	embeddings = orig_um_mat,
	key   = key,
	assay = DefaultAssay(so.st)
)

# --- quick sanity (still in micrometers) --------------------------------------
cat("orig_x_um range: ", rng(so.st$orig_x_um), "\n")
cat("orig_y_um range: ", rng(so.st$orig_y_um), "\n")

# -----------------------------------------------------------------------------
# Plot cell-type-by-cell-type - Each cell type highlighted separately in a loop
# Last Modified: 18/09/25
# -----------------------------------------------------------------------------
if (requireNamespace("conflicted", quietly = TRUE)) {
	conflicted::conflicts_prefer(base::as.factor)
}

# Use the µm-space reduction we added earlier
orig_emb <- Embeddings(so.st, "orig_um")  # columns: orig_x_um, orig_y_um
stopifnot(ncol(orig_emb) >= 2)

spatial_data <- as.data.frame(orig_emb)
# Make sure the column names are exactly what we’ll use below
colnames(spatial_data) <- c("plot_x", "plot_y")
spatial_data$cell <- rownames(spatial_data)

# Attach metadata you need for plotting
spatial_data$cell_type <- so.st$cell_type[spatial_data$cell]
spatial_data$seurat_clusters <- as.factor(Idents(so.st)[spatial_data$cell])

# Mark meta cells (match on full cell IDs, e.g., AAAC...-1_19)
spatial_data$meta_flag    <- ifelse(spatial_data$cell %in% meta_idx, "yes", "no")
spatial_data$variant_flag <- ifelse(spatial_data$meta_flag == "yes", MUT, NA)

# For your loop:
spatial_coords <- c("plot_x","plot_y")

# Sanity: how many meta cells total and with usable coords?
cat("Total meta cells flagged:", sum(spatial_data$meta_flag == "yes"), "\n")
cat("Meta with coords:", sum(spatial_data$meta_flag == "yes" &
									  	complete.cases(spatial_data[, spatial_coords])), "\n")

# -----------------------------------------------------------------------------
# Builds a plotting frame using the original embeddings
# Last Modified: 18/09/25
# -----------------------------------------------------------------------------
library(ggplot2)
library(Seurat)
library(dplyr)
library(ggrepel)

# Use the actual column names (likely spatial_1 and spatial_2)
spatial_coords <- colnames(spatial_data)[1:2]  # First two columns should be coordinates
cluster_centers <- spatial_data %>%
	group_by(seurat_clusters) %>%
	summarise(coord1 = median(.data[[spatial_coords[1]]]), 
				 coord2 = median(.data[[spatial_coords[2]]]))

# --- mark meta indices
meta_idx <- unique(meta.st$Index)
spatial_data$meta_flag    <- ifelse(rownames(spatial_data) %in% meta_idx, "yes", "no")
spatial_data$variant_flag <- ifelse(spatial_data$meta_flag == "yes", MUT, NA)

# --- sizing controls (tweak these if you want)
point_size   <- 3     # size of plotted points (others + target cells)
point_alpha  <- 1.0
ring_size    <- 3     # size of the red ring for meta cells
ring_stroke  <- 1.0   # ring thickness

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
	
	# legend bits (sizes for cell-type dot == ring_size)
	legend_breaks <- if (have_variant) c(ct, MUT) else ct
	legend_values <- c(
		setNames(pal[ct], ct),
		setNames("red", MUT)   # <- use the value of MUT as the name
	)
	
	override_list <- list(
		size   = if (have_variant) rep(ring_size, 2) else ring_size,
		shape  = if (have_variant) c(16, 21) else 16,   # dot for ct, ring for variant
		fill   = if (have_variant) c(NA, NA) else NA,
		stroke = if (have_variant) c(0, ring_stroke) else 0,
		alpha  = 1
	)
	
	dim_plot <- ggplot() +
		# background (others): grey hollow
		geom_point(
			data = other_cells,
			aes(x = .data[[spatial_coords[1]]], y = .data[[spatial_coords[2]]]),
			color = "#D3D3D3", fill = "white", size = point_size, alpha = point_alpha, shape = 21
		) +
		# highlighted target cells
		geom_point(
			data = target_cells,
			aes(x = .data[[spatial_coords[1]]], y = .data[[spatial_coords[2]]], color = cell_type),
			size = point_size, alpha = point_alpha
		) +
		# meta ring ONLY for this cell type (no fill; keeps underlying color visible)
		geom_point(
			data = meta_cells,
			aes(x = .data[[spatial_coords[1]]], y = .data[[spatial_coords[2]]], color = variant_flag),
			shape = 21, fill = NA, size = ring_size, stroke = ring_stroke, alpha = 1
		) +
		geom_text_repel(
			data = cluster_centers,
			aes(x = coord1, y = coord2, label = seurat_clusters),
			size = 4, color = "black"
		) +
		# merged legend
		scale_color_manual(
			name = NULL,
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
			legend.title     = element_blank(),   
			legend.position = c(0.98, 0.98),
			legend.justification = c("right", "top"),
			legend.background = element_rect(fill = "white", color = "grey80"),
			legend.margin = margin(3, 6, 3, 6),
			legend.key.size     = unit(6, "mm")        # optional: bigger legend keys
		) +
		guides(color = guide_legend(override.aes = override_list))
	
	safe_cell_type <- gsub("[^A-Za-z0-9]", "_", ct)
	filename <- paste0("Spatial_ST_", safe_cell_type, "_highlighted_origUM_", MUT, ".png")
	ggsave(file.path(wd.de.plots, "origUM", filename), plot = dim_plot, width = 6, height = 6, dpi = 300)
	cat("Saved plot for", ct, "\n")
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

# map long names -> short legend labels
short_labels <- c(
	"Stage 0"         = "S0",
	"Stage 0A"        = "S0A",
	"Stage 0B"        = "S0B",
	"Stage 1"         = "S1",
	"Stage 2"         = "S2",
	"Stage 3"         = "S3",
	"Leptotene"       = "L",
	"Zygotene"        = "Z",
	"Pachytene"       = "P",
	"Diplotene"       = "D",
	"Meiotic division"= "MD",
	"Early spermatid" = "ES",
	"Late spermatid"  = "LS"
)

# 2) Keep ONLY the types present in the 10 meta cells, in biological order
ct_levels <- ct_order_master[ct_order_master %in% meta_cells_only$cell_type]  # preserves order
# (equivalent to: base::intersect(ct_order_master, unique(meta_cells_only$cell_type)))

# build a labels vector in the *same order* as ct_levels, with fallbacks
legend_labels <- ifelse(ct_levels %in% names(short_labels),
								short_labels[ct_levels],
								ct_levels)

# Relevel & drop unused types
meta_cells_only$cell_type <- factor(meta_cells_only$cell_type, levels = ct_levels)
meta_cells_only <- meta_cells_only[!is.na(meta_cells_only$cell_type),]

# 3) Palette restricted to these 7 (or however many) types
pal_meta <- setNames(scales::hue_pal()(length(ct_levels)), ct_levels)

# 4) Plot only these meta cells, colored by their original cell_type
p_meta <- ggplot(meta_cells_only, aes(x = plot_x, y = plot_y, color = cell_type)) +
	geom_point(shape = 19, size = 4, alpha = 1) +
	scale_color_manual(values = pal_meta, breaks = ct_levels, name = NULL, labels = legend_labels) +  # no legend title
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
outfile <- file.path(wd.de.plots, "origUM", "Spatial_ST_META_ONLY_origcolors_short_labels.png")
ggsave(outfile, plot = p_meta, width = 6, height = 6, dpi = 300)
cat("Saved:", outfile, "\n")

# -----------------------------------------------------------------------------
# 
# Last Modified: 07/10/25
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_23_res=0.5_-C0_ordered_annotated_monocle3+phase.RData")))

# Extract cell type and phase information
cell_data <- data.frame(
	cell_type = so.integrated$seurat_clusters,
	phase = so.integrated$Phase
)

# Calculate proportions of each phase within each cell type
phase_proportions <- cell_data %>%
	group_by(cell_type, phase) %>%
	summarise(count = n(), .groups = "drop") %>%
	group_by(cell_type) %>%
	mutate(
		total = sum(count),
		proportion = count / total
	) %>%
	ungroup()

# Reorder levels so "Stage 0" appears at the top of the y-axis
phase_proportions$cell_type <- factor(
	phase_proportions$cell_type,
	levels = rev(levels(factor(phase_proportions$cell_type)))
)

# Reorder phase levels to show G1, S, G2/M from left to right
# Need to specify in reverse order for stacked bars
phase_proportions$phase <- factor(
	phase_proportions$phase,
	levels = c("G2M", "S", "G1")
)

# -----------------------------------------------------------------------------
# To compute average Euclidean distances between cell types
# Last Modified: 22/09/25
# -----------------------------------------------------------------------------
library(dplyr)

# --- define order
cell_types <- c("Stage 0","Stage 0A","Stage 0B","Stage 1","Stage 2","Stage 3","Leptotene","Zygotene","Pachytene","Diplotene","Meiotic division","Early spermatid","Late spermatid")
cell_type_order <- base::intersect(cell_types, as.vector(unique(meta_cells_only$cell_type)))

meta_cells_only <- meta_cells_only %>%
	dplyr::mutate(cell_type = factor(cell_type, levels = cell_type_order)) %>%
	dplyr::arrange(cell_type)

# --- helper: average Euclidean distance FROM all cells in `from` TO all cells in `to`
avg_dist <- function(df, from, to) {
	df_from <- dplyr::filter(df, .data$cell_type == from) |> dplyr::select(plot_x, plot_y)
	df_to   <- dplyr::filter(df, .data$cell_type == to)   |> dplyr::select(plot_x, plot_y)
	if (nrow(df_from) == 0L || nrow(df_to) == 0L) return(NA_real_)
	
	combo <- base::rbind(as.matrix(df_from), as.matrix(df_to))     # explicit base::
	dmat  <- as.matrix(stats::dist(combo, method = "euclidean"))   # explicit stats::
	
	n_from <- nrow(df_from); n_to <- nrow(df_to)
	submat <- dmat[seq_len(n_from), n_from + seq_len(n_to), drop = FALSE]
	mean(submat)
}

# --- compute pairwise averages along the chain
dist_map <- data.frame(
	from = head(cell_type_order, -1),
	to   = tail(cell_type_order, -1),
	stringsAsFactors = FALSE
)

dist_map$dist_to_next <- mapply(function(f, t) avg_dist(meta_cells_only, f, t),
										  dist_map$from, dist_map$to)

# --- join back (no rename needed)
meta_cells_only <- meta_cells_only |>
	dplyr::left_join(dplyr::select(dist_map, from, dist_to_next),
						  by = c("cell_type" = "from"))

save(meta_idx, spatial_data, spatial_coords, cluster_centers, meta_cells_only, phase_proportions, file=file.path(wd.de.plots, "origUM", paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_23_res=0.5_-C0_ordered_annotated_meta_cells_only_", GENE, ".RData")))

# ensure dir exists
if (!dir.exists(wd.de.plots)) dir.create(wd.de.plots, recursive = TRUE)

# move rownames (cell barcodes) into a proper column, then write
meta_cells_only_out <- tibble::rownames_to_column(as.data.frame(meta_cells_only), var = "cell_barcode")
outfile_tsv <- file.path(wd.de.plots, "origUM", "meta_cells_only.tsv")
utils::write.table(meta_cells_only_out, file = outfile_tsv, sep = "\t",
						 row.names = FALSE, col.names = TRUE, quote = FALSE)
message("Wrote: ", outfile_tsv)

# -----------------------------------------------------------------------------
# To test correlation between distances (dist_to_next) and cell cycle phase proportions (G1/S/G2M)
# Last Modified: 22/09/25
# -----------------------------------------------------------------------------
library(dplyr)
library(ggplot2)

# --- aggregate distances to per-cell_type (unique value for each "from")
dist_summary <- meta_cells_only %>%
	dplyr::select(cell_type, dist_to_next) %>%
	dplyr::distinct()

# --- join with phase proportions
phase_dist <- phase_proportions %>%
	dplyr::left_join(dist_summary, by = "cell_type")

cor_results <- phase_dist %>%
	group_by(phase) %>%
	summarise(
		rho = cor(dist_to_next, proportion, use = "complete.obs"),
		pval = cor.test(dist_to_next, proportion)$p.value,
		.groups = "drop"
	)

print(cor_results)

# Ensure your palette exists
phase_colors <- c("G1" = "#619CFF", "S" = "#F8766D", "G2M" = "#00BA38")

# Distances per cell_type (unique per "from")
dist_summary <- meta_cells_only %>%
	dplyr::select(cell_type, dist_to_next) %>%
	dplyr::distinct()

# Join with phase proportions
phase_dist <- phase_proportions %>%
	dplyr::left_join(dist_summary, by = "cell_type")

# Ensure output dir exists
if (!dir.exists(wd.de.plots)) dir.create(wd.de.plots, recursive = TRUE)

# Helper to make & save one phase plot to a specific directory with a clear filename
save_phase_plot <- function(df, phase_text, phase_label, color_hex, outdir,
									 prefix = paste0("Spatial_ST_META_", GENE)) {
	sub <- df %>% dplyr::filter(.data$phase == phase_text, is.finite(.data$dist_to_next))
	if (nrow(sub) < 3 || all(is.na(sub$dist_to_next))) {
		warning(sprintf("Not enough data to plot/correlate for phase %s", phase_text))
		return(invisible(NULL))
	}
	
	ct  <- stats::cor.test(sub$dist_to_next, sub$proportion)
	rho <- base::unname(ct$estimate)
	p   <- ct$p.value
	lab <- sprintf("rho = %.3f\np = %.3g", rho, p)
	
	x_pos <- max(sub$dist_to_next, na.rm = TRUE)
	y_pos <- max(sub$proportion,   na.rm = TRUE)
	
	p_out <- ggplot2::ggplot(sub, ggplot2::aes(x = .data$dist_to_next, y = .data$proportion)) +
		ggplot2::geom_point(size = 4, alpha = 0.8, color = color_hex) +
		ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, color = color_hex) +
		ggplot2::annotate("text", x = x_pos, y = y_pos, label = lab,
								hjust = 1.05, vjust = 1.4, size = 6) +
		ggplot2::labs(
			title = "",
			x = "Spatial distance to next stage",
			y = paste0(phase_label, " proportion")
		) +
		ggplot2::theme_bw(base_size = 12) +
	   ggplot2::theme(
			legend.position = "none",
			plot.title = ggplot2::element_text(face = "bold", size = 20, hjust = 0.5),
			axis.title = ggplot2::element_text(size = 18),
			axis.text  = ggplot2::element_text(size = 16)
	   )
	
	# Safe filename per phase
	fname <- file.path(outdir,
							 sprintf("%s_%s.png",
							 		  prefix, gsub("[^A-Za-z0-9]+", "_", phase_text)))
	
	ggplot2::ggsave(filename = fname, plot = p_out,
						 width = 5, height = 5, units = "in", dpi = 300)
	message("Saved: ", fname)
	invisible(fname)
}

# Save each phase plot with different file names in the same directory
save_phase_plot(phase_dist, "G1",  "G0/G1", phase_colors[["G1"]],  file.path(wd.de.plots, "origUM"))
save_phase_plot(phase_dist, "S",   "S",     phase_colors[["S"]],   file.path(wd.de.plots, "origUM"))
save_phase_plot(phase_dist, "G2M", "G2/M",  phase_colors[["G2M"]], file.path(wd.de.plots, "origUM"))

write.table(
	dist_summary,
	file = file.path(wd.de.plots, "origUM", paste0("Spatial_ST_META_", GENE, ".txt")),
	sep = "\t",
	quote = FALSE,
	row.names = FALSE,
	col.names = TRUE
)
