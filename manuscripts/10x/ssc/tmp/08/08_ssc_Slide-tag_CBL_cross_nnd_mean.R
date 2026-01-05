# =============================================================================
# Name         : 
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 15/09/25
# =============================================================================
CHR  <- "chr11"
POS  <- 119284966          # if your table has Start column
REF  <- "T"
ALT  <- "A"
GENE <- "CBL"
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
dest <- file.path(wd.de.plots, "origUM", "dist_avg_to_next")
if (!dir.exists(dest)) {
	ok <- base::dir.create(dest, recursive = TRUE, showWarnings = FALSE, mode = "0775")
	if (!ok && !dir.exists(dest)) stop("Failed to create: ", dest)
}

# -----------------------------------------------------------------------------
# After running SComatic genotyping
# # Last Modified: 17/09/25
# -----------------------------------------------------------------------------
library(Seurat)

load(file.path(wd.de.data, paste0("ssc_filtered_normalised_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_ordered_annotated.RData")))
so.integrated$cell_type <- Idents(so.integrated)

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
library(qs)
library(ggplot2) 

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
		plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
		axis.title = element_text(size = 16),
		axis.text  = element_text(size = 14),
		legend.text = element_text(size = 14),
		legend.position = c(0.05, 0.05),
		legend.justification = c("left", "bottom"),
		legend.background = element_rect(fill = "white", color = "grey80"),
		legend.margin = margin(3, 6, 3, 6)
	) +
	ggtitle("CBL (11:119284966 T>A, splice-region)\n")

print(p_meta)

# Save
outfile <- file.path(wd.de.plots, "origUM", "Spatial_ST_META_ONLY_origcolors_short_labels_title_test.png")
ggsave(outfile, plot = p_meta, width = 6, height = 6, dpi = 300)
cat("Saved:", outfile, "\n")

# -----------------------------------------------------------------------------
# 
# Last Modified: 07/10/25
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_ordered_annotated_monocle3+phase_RNA.RData")))

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
# Compute cross–nearest-neighbor mean distance BETWEEN adjacent cell types
# (directed: FROM all cells in `from` TO their nearest neighbor in `to`)
# Last Modified: 26/10/25
# -----------------------------------------------------------------------------
library(dplyr)

# --- define order
cell_types <- c("Stage 0","Stage 0A","Stage 0B","Stage 1","Stage 2","Stage 3",
					 "Leptotene","Zygotene","Pachytene","Diplotene",
					 "Meiotic division","Early spermatid","Late spermatid")
cell_type_order <- base::intersect(cell_types, as.vector(unique(meta_cells_only$cell_type)))

meta_cells_only <- meta_cells_only %>%
	dplyr::mutate(cell_type = factor(cell_type, levels = cell_type_order)) %>%
	dplyr::arrange(cell_type)

# --- metrics ---------------------------------------------------------------
# 1) Average all-pairs Euclidean distance (global separation)
avg_dist <- function(df, from, to) {
	df_from <- dplyr::filter(df, .data$cell_type == from) |> dplyr::select(plot_x, plot_y)
	df_to   <- dplyr::filter(df, .data$cell_type == to)   |> dplyr::select(plot_x, plot_y)
	if (nrow(df_from) == 0L || nrow(df_to) == 0L) return(NA_real_)
	
	combo <- base::rbind(as.matrix(df_from), as.matrix(df_to))
	dmat  <- as.matrix(stats::dist(combo, method = "euclidean"))
	
	n_from <- nrow(df_from); n_to <- nrow(df_to)
	submat <- dmat[seq_len(n_from), n_from + seq_len(n_to), drop = FALSE]
	mean(submat)
}

# 2) Directional cross–nearest-neighbour median (local interface A→B)
#    (median is robust; switch to mean if you prefer)
cross_nnd_dir <- function(df, from, to, fun = median) {
	df_from <- dplyr::filter(df, .data$cell_type == from) |> dplyr::select(plot_x, plot_y)
	df_to   <- dplyr::filter(df, .data$cell_type == to)   |> dplyr::select(plot_x, plot_y)
	if (nrow(df_from) < 1L || nrow(df_to) < 1L) return(NA_real_)
	nn <- FNN::get.knnx(as.matrix(df_to), as.matrix(df_from), k = 1)$nn.dist[,1]
	fun(nn)
}

# --- adjacent pairs along biological order --------------------------------
dist_map <- data.frame(
	from = head(cell_type_order, -1),
	to   = tail(cell_type_order, -1),
	stringsAsFactors = FALSE
)

dist_map$dist_avg_to_next <- mapply(function(f,t) avg_dist(meta_cells_only, f, t),
												dist_map$from, dist_map$to)

dist_map$crossNND_med_to_next <- mapply(function(f,t) cross_nnd_dir(meta_cells_only, f, t, median),
													 dist_map$from, dist_map$to)

# --- join back (retain both metrics) ---------------------------------------
meta_cells_only <- meta_cells_only %>%
	dplyr::left_join(dplyr::select(dist_map, from, dist_avg_to_next, crossNND_med_to_next),
						  by = c("cell_type" = "from"))

# --- collapse to one row per "from" stage ---------------------------------
dist_summary <- meta_cells_only %>%
	dplyr::select(cell_type, dist_avg_to_next, crossNND_med_to_next) %>%
	dplyr::distinct()

# --- join with phase proportions (G1/S/G2M) --------------------------------
phase_dist <- phase_proportions %>%
	dplyr::left_join(dist_summary, by = "cell_type")

# --- helper: do a Spearman and make a plot --------------------------------
plot_phase_vs_metric <- function(df, phase_text, ylab, xcol, color_hex, outfile_prefix){
	sub <- df %>% dplyr::filter(.data$phase == phase_text,
										 is.finite(.data[[xcol]]))
	if (nrow(sub) < 3 || all(is.na(sub[[xcol]]))) return(invisible(NULL))
	ct  <- suppressWarnings(cor.test(sub[[xcol]], sub$proportion, method="spearman", exact=FALSE))
	rho <- unname(ct$estimate); p <- ct$p.value
	lab <- sprintf("Spearman rho = %.3f\np = %.3g", rho, p)
	
	xp <- max(sub[[xcol]], na.rm = TRUE)
	yp <- max(sub$proportion,   na.rm = TRUE)
	
	p_out <- ggplot2::ggplot(sub, aes(x = .data[[xcol]], y = proportion)) +
		ggplot2::geom_point(size = 4, alpha = 0.8, color = color_hex) +
		ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, color = color_hex) +
		ggplot2::annotate("text", x = xp, y = yp, label = lab,
								hjust = 1.05, vjust = 1.4, size = 6) +
		ggplot2::labs(x = "Spatial distance (µm)", y = ylab) +
		ggplot2::theme_bw(base_size = 12) +
		ggplot2::theme(legend.position = "none",
							axis.title = element_text(size = 18),
							axis.text  = element_text(size = 16))
	
	fn <- file.path(wd.de.plots, "origUM", "cross_nnd_mean",
						 sprintf("%s_%s_%s.png",
						 		  outfile_prefix, gsub("[^A-Za-z0-9]+","_", phase_text), xcol))
	ggplot2::ggsave(fn, p_out, width = 5, height = 5, dpi = 300)
	message("Saved: ", fn)
	invisible(
		tibble::tibble(metric = xcol, phase = phase_text, rho = rho, p = p, n = nrow(sub))
	)
}

# --- colors for phases ------------------------------------------------------
phase_colors <- c("G1"="#619CFF","S"="#F8766D","G2M"="#00BA38")

# --- run six plots: (avg vs G1/S/G2M) + (crossNND vs G1/S/G2M) ---------------
res_stats <- dplyr::bind_rows(
	plot_phase_vs_metric(phase_dist, "G1",  "% G0/G1", "dist_avg_to_next",      phase_colors[["G1"]],  "phase_vs_dist"),
	plot_phase_vs_metric(phase_dist, "S",   "% S",     "dist_avg_to_next",      phase_colors[["S"]],   "phase_vs_dist"),
	plot_phase_vs_metric(phase_dist, "G2M", "% G2/M",  "dist_avg_to_next",      phase_colors[["G2M"]], "phase_vs_dist"),
	plot_phase_vs_metric(phase_dist, "G1",  "% G0/G1", "crossNND_med_to_next",  phase_colors[["G1"]],  "phase_vs_crossNND"),
	plot_phase_vs_metric(phase_dist, "S",   "% S",     "crossNND_med_to_next",  phase_colors[["S"]],   "phase_vs_crossNND"),
	plot_phase_vs_metric(phase_dist, "G2M", "% G2/M",  "crossNND_med_to_next",  phase_colors[["G2M"]], "phase_vs_crossNND")
)

# --- save a small summary table --------------------------------------------
readr::write_tsv(res_stats,
					  file.path(wd.de.plots, "origUM", "cross_nnd_mean", "phase_vs_distance_spearman.tsv"))

# -----------------------------------------------------------------------------
# To check whether the S-phase association is just riding on stage size (total cells) or mutant burden
# Last Modified: 28/10/25
# -----------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(broom)

# --- 1) One row per stage with %G1/%S/%G2M and the metrics ---------------
phase_wide <- phase_dist %>%
	select(cell_type, phase, proportion, dist_avg_to_next, crossNND_med_to_next) %>%
	distinct() %>%
	mutate(phase = factor(phase, levels = c("G1","S","G2M"))) %>%
	pivot_wider(names_from = phase, values_from = proportion, names_prefix = "pct_") %>%
	group_by(cell_type) %>%
	summarise(
		pct_G1  = dplyr::first(pct_G1),
		pct_S   = dplyr::first(pct_S),
		pct_G2M = dplyr::first(pct_G2M),
		dist_avg_to_next     = dplyr::first(dist_avg_to_next),
		crossNND_med_to_next = dplyr::first(crossNND_med_to_next),
		.groups = "drop"
	)

# --- 2) Stage counts: n_total from phase_dist$total; n_mut from meta_cells_only ----
n_total_tbl <- phase_dist %>%
	group_by(cell_type) %>%
	summarise(n_total = dplyr::first(total), .groups = "drop")

if (exists("meta_cells_only")) {
	n_mut_tbl <- meta_cells_only %>%
		group_by(cell_type = as.character(cell_type)) %>%
		summarise(n_mut = dplyr::n(), .groups = "drop")
} else {
	# Fallback if you don't have meta_cells_only in memory
	n_mut_tbl <- tibble(cell_type = n_total_tbl$cell_type, n_mut = 0L)
}

stage_counts <- n_total_tbl %>%
	left_join(n_mut_tbl, by = "cell_type") %>%
	mutate(n_mut = ifelse(is.na(n_mut), 0L, n_mut))

# --- 3) Merge features + counts -------------------------------------------
df_anal <- phase_wide %>%
	inner_join(stage_counts, by = "cell_type") %>%
	mutate(log_n_total = log1p(n_total))

# -----------------------------------------------------------------------------
# Merge n_total and n_mut to phase_dist
# -----------------------------------------------------------------------------
# 1) Per-stage totals from phase_dist (total is constant within stage)
n_total_tbl <- phase_dist %>%
	group_by(cell_type) %>%
	summarise(n_total = dplyr::first(total), .groups = "drop")

# 2) Per-stage mutant counts from meta_cells_only (fallback to 0 if not present)
if (exists("meta_cells_only")) {
	n_mut_tbl <- meta_cells_only %>%
		mutate(cell_type = as.character(cell_type)) %>%
		dplyr::count(cell_type, name = "n_mut")
} else {
	n_mut_tbl <- tibble(cell_type = unique(phase_dist$cell_type), n_mut = 0L)
}

# 3) Merge counts and compute log scale
stage_counts <- n_total_tbl %>%
	left_join(n_mut_tbl, by = "cell_type") %>%
	mutate(
		n_mut = tidyr::replace_na(n_mut, 0L),
		log_n_total = log1p(n_total)
	)

# 4) Add these columns onto every row of phase_dist
phase_dist <- phase_dist %>%
	left_join(stage_counts, by = "cell_type")

# (optional) quick check
# phase_dist %>% distinct(cell_type, n_total, n_mut, log_n_total) %>% print(n = 99)

# --- helper: do a Spearman and make a plot --------------------------------
plot_phase_vs_metric <- function(df, phase_text, ylab, xcol, color_hex, outfile_prefix, title){
	sub <- df %>% dplyr::filter(.data$phase == phase_text,
										 is.finite(.data[[xcol]]))
	if (nrow(sub) < 3 || all(is.na(sub[[xcol]]))) return(invisible(NULL))
	ct  <- suppressWarnings(cor.test(sub[[xcol]], sub$proportion, method="spearman", exact=FALSE))
	rho <- unname(ct$estimate); p <- ct$p.value
	lab <- sprintf("rho = %.3f; p = %.3g", rho, p)
	
	xp <- max(sub[[xcol]], na.rm = TRUE)
	yp <- max(sub$proportion,   na.rm = TRUE)
	
	p_out <- ggplot2::ggplot(sub, aes(x = .data[[xcol]], y = proportion)) +
		ggplot2::geom_point(size = 4, alpha = 0.8, color = color_hex) +
		ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, color = color_hex) +
		ggplot2::annotate("text", x = xp, y = yp, label = lab,
								hjust = 1.05, vjust = 1.4, size = 6) +
		ggplot2::labs(x = "# cells", y = ylab, title = title) +
		ggplot2::theme_bw(base_size = 12) +
		ggplot2::theme(legend.position = "none",
							plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
							axis.title = element_text(size = 16),
							axis.text  = element_text(size = 14))
	
	fn <- file.path(wd.de.plots, "origUM", "cross_nnd_mean",
						 sprintf("%s_%s_%s.png",
						 		  outfile_prefix, gsub("[^A-Za-z0-9]+","_", phase_text), xcol))
	ggplot2::ggsave(fn, p_out, width = 5, height = 5, dpi = 300)
	message("Saved: ", fn)
	invisible(
		tibble::tibble(metric = xcol, phase = phase_text, rho = rho, p = p, n = nrow(sub))
	)
}

# --- colors for phases ------------------------------------------------------
phase_colors <- c("G1"="#619CFF","S"="#F8766D","G2M"="#00BA38")

# --- run six plots: (avg vs G1/S/G2M) + (crossNND vs G1/S/G2M) ---------------
res_stats_n_mut_dis <- dplyr::bind_rows(
	plot_phase_vs_metric(phase_dist, "G1",  "% G0/G1", "n_total",      phase_colors[["G1"]],  "dis_vs_n_total", "Stage size"),
	plot_phase_vs_metric(phase_dist, "S",   "% S",     "n_total",      phase_colors[["S"]],   "dis_vs_n_total", "Stage size"),
	plot_phase_vs_metric(phase_dist, "G2M", "% G2/M",  "n_total",      phase_colors[["G2M"]], "dis_vs_n_total", "Stage size"),
	plot_phase_vs_metric(phase_dist, "G1",  "% G0/G1", "log_n_total",      phase_colors[["G1"]],  "dis_vs_log_n_total", "Stage size (log)"),
	plot_phase_vs_metric(phase_dist, "S",   "% S",     "log_n_total",      phase_colors[["S"]],   "dis_vs_log_n_total", "Stage size (log)"),
	plot_phase_vs_metric(phase_dist, "G2M", "% G2/M",  "log_n_total",      phase_colors[["G2M"]], "dis_vs_log_n_total", "Stage size (log)"),
	plot_phase_vs_metric(phase_dist, "G1",  "% G0/G1", "n_mut",  phase_colors[["G1"]],  "dis_vs_n_mut", "Mutant size"),
	plot_phase_vs_metric(phase_dist, "S",   "% S",     "n_mut",  phase_colors[["S"]],   "dis_vs_n_mut", "Mutant size"),
	plot_phase_vs_metric(phase_dist, "G2M", "% G2/M",  "n_mut",  phase_colors[["G2M"]], "dis_vs_n_mut", "Mutant size")
)

# --- save a small summary table --------------------------------------------
readr::write_tsv(res_stats_n_mut,
					  file.path(wd.de.plots, "origUM", "cross_nnd_mean", "phase_vs_distance_spearman_n_mut.tsv"))

# =============================================================================
# 
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })

# --- Collapse to stage-level (one row per cell_type) -----------------------
stage_df <- phase_dist %>%
	select(cell_type, dist_avg_to_next, n_total, log_n_total, n_mut) %>%
	distinct(cell_type, .keep_all = TRUE) %>%
	drop_na(dist_avg_to_next)

# --- Helper: Spearman between X = xcol and Y = crossNND_med_to_next --------
plot_cov_vs_crossNND <- function(df, xcol, xlab, color_hex = "black",
											outfile_prefix = "cov_vs_crossNND",
											outdir = file.path(wd.de.plots, "origUM", "cross_nnd_mean")) {
	stopifnot(xcol %in% names(df))
	sub <- df %>%
		select(cell_type, dist_avg_to_next, !!rlang::sym(xcol)) %>%
		drop_na(dist_avg_to_next, !!rlang::sym(xcol))
	
	if (nrow(sub) < 3) return(invisible(NULL))
	
	ct  <- suppressWarnings(cor.test(sub[[xcol]], sub$dist_avg_to_next,
												method = "spearman", exact = FALSE))
	rho <- unname(ct$estimate); p <- ct$p.value
	lab <- sprintf("rho = %.3f; p = %.3g", rho, p, nrow(sub))
	
	xp <- max(sub[[xcol]], na.rm = TRUE)
	yp <- max(sub$dist_avg_to_next, na.rm = TRUE)
	
	p_out <- ggplot(sub, aes(x = .data[[xcol]], y = dist_avg_to_next)) +
		geom_point(size = 4, alpha = 0.85, color = color_hex) +
		geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, color = color_hex) +
		annotate("text", x = xp, y = yp, label = lab, hjust = 1.05, vjust = 1.4, size = 5.5) +
		labs(x = "# cells", y = "dist_avg_to_next (µm)", title = xlab) +
		theme_bw(base_size = 12) +
		theme(legend.position = "none",
				plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
				axis.title = element_text(size = 16),
				axis.text  = element_text(size = 14))
	
	dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
	fn <- file.path(outdir, sprintf("%s_X=%s.png", outfile_prefix, xcol))
	ggsave(fn, p_out, width = 5, height = 5, dpi = 300)
	message("Saved: ", fn)
	
	invisible(tibble::tibble(x = xcol, y = "dist_avg_to_next",
									 rho = rho, p = p, n = nrow(sub)))
}

# --- Run the three tests (phase-free, flipped axes) ------------------------
res_stats_dist_avg_vs_cov <- dplyr::bind_rows(
	plot_cov_vs_crossNND(stage_df, "n_total",     "Stage size",      "black", "dist_avg_vs_size"),
	plot_cov_vs_crossNND(stage_df, "log_n_total", "Stage size (log)",  "black", "dist_avg_vs_size_log"),
	plot_cov_vs_crossNND(stage_df, "n_mut",       "Mutant burden",   "black", "dist_avg_vs_nmut")
)

print(res_stats_dist_avg_vs_cov)

# --- Save summary table -----------------------------------------------------
readr::write_tsv(res_stats_dist_avg_vs_cov,
					  file.path(wd.de.plots, "origUM", "cross_nnd_mean", "dist_avg_vs_size_mut_spearman.tsv"))

















# -----------------------------------------------------------------------------
# Partial Spearman that controls only for n_mut
# -----------------------------------------------------------------------------
# df_anal must have: dist_avg_to_next, crossNND_med_to_next, pct_S, n_mut

# Helper: partial Spearman of X ~ Y | n_mut (via residuals)
partial_spearman_nmut <- function(df, x, y) {
	stopifnot(all(c(x, y, "n_mut") %in% names(df)))
	d <- df %>% select(all_of(c(x, y, "n_mut"))) %>% na.omit()
	# residualize X and Y on n_mut (linear works fine; Spearman is on ranks afterward)
	rx <- lm(reformulate("n_mut", response = x), data = d)$residuals
	ry <- lm(reformulate("n_mut", response = y), data = d)$residuals
	out <- cor.test(rx, ry, method = "spearman", exact = FALSE)
	tibble::tibble(var_x = x, var_y = y, rho = unname(out$estimate), p = out$p.value, n = nrow(d))
}

# Run the two tests you care about (S vs dist, S vs crossNND), controlling for n_mut
res_partial <- bind_rows(
	partial_spearman_nmut(df_anal, "dist_avg_to_next",     "pct_S"),
	partial_spearman_nmut(df_anal, "crossNND_med_to_next", "pct_S")
)
print(res_partial)

# --- helper: partial Spearman (control for n_mut) + make a plot ------------
plot_phase_vs_metric_control <- function(df, phase_text, ylab, xcol, color_hex, outfile_prefix){
	sub <- df %>%
		dplyr::filter(.data$phase == phase_text,
						  is.finite(.data[[xcol]]),
						  is.finite(proportion)) %>%
		dplyr::select(cell_type, phase, proportion, !!rlang::sym(xcol),
						  dplyr::any_of("n_mut")) %>%
		tidyr::drop_na(proportion, !!rlang::sym(xcol))
	
	if (nrow(sub) < 3 || all(is.na(sub[[xcol]]))) return(invisible(NULL))
	
	# can we do partial | n_mut?
	can_partial <- ("n_mut" %in% names(sub)) &&
		all(is.finite(sub$n_mut)) &&
		stats::var(sub$n_mut) > 0
	
	if (can_partial) {
		# residualize X and Y on n_mut; then Spearman on residuals
		rx <- lm(reformulate("n_mut", response = xcol), data = sub)$residuals
		ry <- lm(reformulate("n_mut", response = "proportion"), data = sub)$residuals
		ct  <- suppressWarnings(cor.test(rx, ry, method = "spearman", exact = FALSE))
		rho <- unname(ct$estimate); p <- ct$p.value
		lab <- sprintf("Partial Spearman | n_mut\nrho = %.3f\np = %.3g", rho, p)
	} else {
		ct  <- suppressWarnings(cor.test(sub[[xcol]], sub$proportion, method = "spearman", exact = FALSE))
		rho <- unname(ct$estimate); p <- ct$p.value
		lab <- sprintf("Spearman rho = %.3f\np = %.3g", rho, p)
		message(sprintf("[WARN] %s: cannot compute partial | n_mut (missing/constant n_mut). Using simple Spearman.",
							 phase_text))
	}
	
	xp <- max(sub[[xcol]], na.rm = TRUE)
	yp <- max(sub$proportion,   na.rm = TRUE)
	
	p_out <- ggplot2::ggplot(sub, aes(x = .data[[xcol]], y = proportion)) +
		ggplot2::geom_point(size = 4, alpha = 0.8, color = color_hex) +
		ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, color = color_hex) +
		ggplot2::annotate("text", x = xp, y = yp, label = lab,
								hjust = 1.05, vjust = 1.4, size = 6) +
		ggplot2::labs(x = "Spatial distance (µm)", y = ylab) +
		ggplot2::theme_bw(base_size = 12) +
		ggplot2::theme(legend.position = "none",
							plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
							axis.title = element_text(size = 16),
							axis.text  = element_text(size = 14))
	
	fn <- file.path(wd.de.plots, "origUM", "cross_nnd_mean",
						 sprintf("%s_%s_%s_partial_nmut.png",
						 		  outfile_prefix, gsub("[^A-Za-z0-9]+","_", phase_text), xcol))
	ggplot2::ggsave(fn, p_out, width = 5, height = 5, dpi = 300)
	message("Saved: ", fn)
	
	invisible(
		tibble::tibble(metric = xcol, phase = phase_text, rho = rho, p = p,
							n = nrow(sub), partial_nmut = can_partial)
	)
}

# colors
phase_colors <- c("G1"="#619CFF","S"="#F8766D","G2M"="#00BA38")

# run the six analyses (now partial | n_mut)
res_stats_control <- dplyr::bind_rows(
	plot_phase_vs_metric_control(phase_dist, "G1",  "% G0/G1", "dist_avg_to_next",      phase_colors[["G1"]],  "phase_vs_dist"),
	plot_phase_vs_metric_control(phase_dist, "S",   "% S",     "dist_avg_to_next",      phase_colors[["S"]],   "phase_vs_dist"),
	plot_phase_vs_metric_control(phase_dist, "G2M", "% G2/M",  "dist_avg_to_next",      phase_colors[["G2M"]], "phase_vs_dist"),
	plot_phase_vs_metric_control(phase_dist, "G1",  "% G0/G1", "crossNND_med_to_next",  phase_colors[["G1"]],  "phase_vs_crossNND"),
	plot_phase_vs_metric_control(phase_dist, "S",   "% S",     "crossNND_med_to_next",  phase_colors[["S"]],   "phase_vs_crossNND"),
	plot_phase_vs_metric_control(phase_dist, "G2M", "% G2/M",  "crossNND_med_to_next",  phase_colors[["G2M"]], "phase_vs_crossNND")
)

readr::write_tsv(res_stats_control,
					  file.path(wd.de.plots, "origUM", "cross_nnd_mean", "phase_vs_distance_partial_spearman_nmut_control.tsv"))

















suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })

# --- Collapse to stage-level (one row per cell_type) -----------------------
stage_df <- phase_dist %>%
	select(cell_type, crossNND_med_to_next, n_total, log_n_total, n_mut) %>%
	distinct(cell_type, .keep_all = TRUE) %>%
	drop_na(crossNND_med_to_next)

# --- Helper: Spearman between X=dist_avg_to_next and a Y column ------------
plot_dist_vs <- function(df, ycol, ylab, color_hex = "#5555FF",
								 outfile_prefix = "dist_vs_cov",
								 outdir = file.path(wd.de.plots, "origUM", "cross_nnd_mean")) {
	stopifnot(ycol %in% names(df))
	sub <- df %>% select(cell_type, crossNND_med_to_next, !!rlang::sym(ycol)) %>%
		drop_na(crossNND_med_to_next, !!rlang::sym(ycol))
	
	if (nrow(sub) < 3) return(invisible(NULL))
	
	ct  <- suppressWarnings(cor.test(sub$crossNND_med_to_next, sub[[ycol]],
												method = "spearman", exact = FALSE))
	rho <- unname(ct$estimate); p <- ct$p.value
	lab <- sprintf("rho = %.3f; p = %.3g", rho, p, nrow(sub))
	
	xp <- max(sub$crossNND_med_to_next, na.rm = TRUE)
	yp <- max(sub[[ycol]],          na.rm = TRUE)
	
	p_out <- ggplot(sub, aes(x = crossNND_med_to_next, y = .data[[ycol]])) +
		geom_point(size = 4, alpha = 0.85, color = color_hex) +
		geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, color = color_hex) +
		annotate("text", x = xp, y = yp, label = lab, hjust = 1.05, vjust = 1.4, size = 5.5) +
		labs(x = "crossNND_med_to_next (µm)", y = "# cells", title = ylab) +
		theme_bw(base_size = 12) +
		theme(legend.position = "none",
				plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
				axis.title = element_text(size = 16),
				axis.text  = element_text(size = 14))
	
	dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
	fn <- file.path(outdir, sprintf("%s_Y=%s.png", outfile_prefix, ycol))
	ggsave(fn, p_out, width = 5, height = 5, dpi = 300)
	message("Saved: ", fn)
	
	invisible(tibble::tibble(x = "crossNND_med_to_next", y = ycol, rho = rho, p = p, n = nrow(sub)))
}

# --- Run the three tests (no phase anywhere) --------------------------------
res_stats_crossNND_vs_cov <- dplyr::bind_rows(
	plot_dist_vs(stage_df, "n_total",      "Stage size",      "black", "crossNND_med_to_next_vs_size"),
	plot_dist_vs(stage_df, "log_n_total",  "Stage size (log)",  "black", "crossNND_med_to_next_vs_size_log"),
	plot_dist_vs(stage_df, "n_mut",        "Mutant burden",   "black", "crossNND_med_to_next_vs_nmut")
)

print(res_stats_dist_vs_cov)

# --- Save summary table -----------------------------------------------------
readr::write_tsv(res_stats_crossNND_vs_cov,
					  file.path(wd.de.plots, "origUM", "cross_nnd_mean", "crossNND_med_to_next_vs_size_mut_spearman_phasefree.tsv"))












# =============================================================================
# WITHIN-STAGE dispersion vs G1/S/G2M fractions
# Using ppp_from_df() for within-stage NND (spatstat-based)
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(readr)})

# Keep only valid cells with coords
stage_order_full <- c("Stage 0","Stage 0A","Stage 0B","Stage 1","Stage 2","Stage 3",
							 "Leptotene","Zygotene","Pachytene","Diplotene",
							 "Meiotic division","Early spermatid","Late spermatid")

all_cells_xy <- spatial_data %>%
	dplyr::filter(cell_type %in% stage_order_full,
					  is.finite(plot_x), is.finite(plot_y)) %>%
	dplyr::mutate(cell_type = factor(cell_type, levels = stage_order_full))

# --- helpers -----------------------------------------------------------------
# Mean of all pairwise distances within a stage (global spread)
intra_avg_dist <- function(df_xy) {
	if (nrow(df_xy) < 2) return(NA_real_)
	d <- stats::dist(as.matrix(df_xy[, c("plot_x","plot_y")]), method = "euclidean")
	mean(as.numeric(d))
}

# -------------------------- Spatial helpers ----------------------------------
ppp_from_df <- function(df_stage) {
	stopifnot(nrow(df_stage) >= 3)
	W <- owin(poly = spatstat.geom::convexhull.xy(df_stage$x, df_stage$y)$bdry[[1]])
	ppp(df_stage$x, df_stage$y, window = W,
		 marks = data.frame(mutant = df_stage$mutant,
		 						 S_prob = df_stage$S_prob,
		 						 CB = df_stage$CB))
}

# Median nearest-neighbour distance using spatstat (local packing)
# Requires that ppp_from_df(df_xy) returns a spatstat.geom::ppp object
intra_nnd_ppp <- function(df_xy) {
	if (nrow(df_xy) < 2) return(NA_real_)
	p <- ppp_from_df(df_xy)
	nn <- spatstat.geom::nndist(p)
	median(nn, na.rm = TRUE)
}

# --- compute per-stage metrics ----------------------------------------------
intra_metrics <- all_cells_xy %>%
	dplyr::group_by(cell_type) %>%
	dplyr::summarise(
		n_cells        = n(),
		intra_avg_dist = intra_avg_dist(cur_data_all()),
		intra_nnd_ppp  = intra_nnd_ppp(cur_data_all()),
		.groups = "drop"
	) %>%
	dplyr::rename(stage = cell_type)

# --- join with phase proportions --------------------------------------------
phase_dist_within <- phase_proportions %>%
	dplyr::filter(cell_type %in% stage_order_full,
					  phase %in% c("G1","S","G2M")) %>%
	dplyr::rename(stage = cell_type) %>%
	dplyr::left_join(intra_metrics, by = "stage")

# --- helper: Spearman correlation + plot ------------------------------------
plot_phase_vs_intra <- function(df, phase_text, metric_col, ylab, color_hex, outprefix){
	sub <- df %>% dplyr::filter(phase == phase_text, is.finite(.data[[metric_col]]))
	if (nrow(sub) < 3 || all(is.na(sub[[metric_col]]))) return(invisible(NULL))
	
	ct <- suppressWarnings(cor.test(sub[[metric_col]], sub$proportion,
											  method = "spearman", exact = FALSE))
	rho <- unname(ct$estimate); p <- ct$p.value
	lab <- sprintf("Spearman rho = %.3f\np = %.3g", rho, p)
	
	xp <- max(sub[[metric_col]], na.rm = TRUE)
	yp <- max(sub$proportion,   na.rm = TRUE)
	
	gp <- ggplot2::ggplot(sub, aes(x = .data[[metric_col]], y = proportion)) +
		ggplot2::geom_point(size = 4, alpha = 0.8, color = color_hex) +
		ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, color = color_hex) +
		ggplot2::annotate("text", x = xp, y = yp, label = lab,
								hjust = 1.05, vjust = 1.4, size = 6) +
		ggplot2::labs(x = paste0(metric_col, " (µm)"), y = ylab) +
		ggplot2::theme_bw(base_size = 12) +
		ggplot2::theme(legend.position = "none",
							axis.title = element_text(size = 18),
							axis.text  = element_text(size = 16))
	
	fn <- file.path(wd.de.plots, "origUM", "cross_nnd_mean",
						 sprintf("%s_%s_%s.png",
						 		  outprefix, gsub("[^A-Za-z0-9]+","_", phase_text), metric_col))
	ggplot2::ggsave(fn, gp, width = 5, height = 5, dpi = 300)
	message("Saved: ", fn)
	
	tibble::tibble(scope = "within-stage",
						metric = metric_col, phase = phase_text,
						rho = rho, p = p, n_stages = nrow(sub))
}

phase_colors <- c("G1"="#619CFF","S"="#F8766D","G2M"="#00BA38")

# --- run 6 plots (2 metrics × 3 phases) -------------------------------------
stats_within <- dplyr::bind_rows(
	plot_phase_vs_intra(phase_dist_within, "G1",  "intra_avg_dist", "% G0/G1", phase_colors[["G1"]],  "within_phase_vs_avgdist"),
	plot_phase_vs_intra(phase_dist_within, "S",   "intra_avg_dist", "% S",     phase_colors[["S"]],   "within_phase_vs_avgdist"),
	plot_phase_vs_intra(phase_dist_within, "G2M", "intra_avg_dist", "% G2/M",  phase_colors[["G2M"]], "within_phase_vs_avgdist"),
	plot_phase_vs_intra(phase_dist_within, "G1",  "intra_nnd_ppp",  "% G0/G1", phase_colors[["G1"]],  "within_phase_vs_pppNND"),
	plot_phase_vs_intra(phase_dist_within, "S",   "intra_nnd_ppp",  "% S",     phase_colors[["S"]],   "within_phase_vs_pppNND"),
	plot_phase_vs_intra(phase_dist_within, "G2M", "intra_nnd_ppp",  "% G2/M",  phase_colors[["G2M"]], "within_phase_vs_pppNND")
)

readr::write_tsv(stats_within,
					  file.path(wd.de.plots, "origUM", "cross_nnd_mean", "within_stage_phase_vs_dispersion_spearman.tsv"))

print(stats_within)


































# --- helper: CROSS NND mean FROM `from` TO `to` ------------------------------
cross_nnd_mean <- function(df, from, to) {
	df_from <- dplyr::filter(df, .data$cell_type == from) |> dplyr::select(plot_x, plot_y)
	df_to   <- dplyr::filter(df, .data$cell_type == to)   |> dplyr::select(plot_x, plot_y)
	if (nrow(df_from) == 0L || nrow(df_to) == 0L) return(NA_real_)
	# Distance matrix between groups
	combo <- base::rbind(as.matrix(df_from), as.matrix(df_to))
	dmat  <- as.matrix(stats::dist(combo, method = "euclidean"))
	n_from <- nrow(df_from)
	n_to   <- nrow(df_to)
	submat <- dmat[seq_len(n_from), n_from + seq_len(n_to), drop = FALSE]
	mean(apply(submat, 1, min))  # mean of each FROM point’s nearest neighbor in TO
}

# If you ever want the symmetric version, average both directions:
# cross_nnd_mean_sym <- function(df, a, b) {
#   mean(c(cross_nnd_mean(df, a, b), cross_nnd_mean(df, b, a)), na.rm = TRUE)
# }

# --- compute cross-NND (directed) along the chain -----------------------------
dist_map <- data.frame(
	from = head(cell_type_order, -1),
	to   = tail(cell_type_order, -1),
	stringsAsFactors = FALSE
)

dist_map$dist_to_next <- mapply(function(f, t) cross_nnd_mean(meta_cells_only, f, t),
										  dist_map$from, dist_map$to)

# --- join back (no rename needed) ---------------------------------------------
meta_cells_only <- meta_cells_only %>%
	dplyr::left_join(dplyr::select(dist_map, from, dist_to_next),
						  by = c("cell_type" = "from"))

# Save objects (note: save the new helper name, not avg_dist)
save(meta_idx, spatial_data, spatial_coords, cluster_centers,
	  cross_nnd_mean, dist_map, meta_cells_only, phase_proportions,
	  file = file.path(wd.de.plots, "origUM",
	  					  paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_23_res=0.5_-C0_ordered_annotated_meta_cells_only_",
	  					  		 GENE, "_RNA.RData")))

# ensure dir exists
if (!dir.exists(wd.de.plots)) dir.create(wd.de.plots, recursive = TRUE)

# move rownames (cell barcodes) into a proper column, then write
meta_cells_only_out <- tibble::rownames_to_column(as.data.frame(meta_cells_only), var = "cell_barcode")
outfile_tsv <- file.path(wd.de.plots, "origUM", "meta_cells_only.tsv")
utils::write.table(meta_cells_only_out, file = outfile_tsv, sep = "\t",
						 row.names = FALSE, col.names = TRUE, quote = FALSE)
message("Wrote: ", outfile_tsv)

# -----------------------------------------------------------------------------
# Correlate cross-NND (dist_to_next) with cell-cycle phase proportions (G1/S/G2M)
# -----------------------------------------------------------------------------
library(ggplot2)

# Distances per cell_type (unique per "from")
dist_summary <- meta_cells_only %>%
	dplyr::select(cell_type, dist_to_next) %>%
	dplyr::distinct()

# Join with phase proportions
phase_dist <- phase_proportions %>%
	dplyr::left_join(dist_summary, by = "cell_type")

cor_results_obs <- phase_dist %>%
	dplyr::group_by(phase) %>%
	dplyr::summarise(
		rho  = suppressWarnings(cor(dist_to_next, proportion, use = "complete.obs", method = "spearman")),
		pval = suppressWarnings(cor.test(dist_to_next, proportion, method = "spearman"))$p.value,
		.groups = "drop"
	)
print(cor_results_obs)

# Colors
phase_colors <- c("G1" = "#619CFF", "S" = "#F8766D", "G2M" = "#00BA38")

# Helper to plot one phase
save_phase_plot <- function(df, phase_text, phase_label, color_hex, outdir,
									 prefix = "Spatial_ST_META_CBL_XNND") {
	sub <- df %>% dplyr::filter(.data$phase == phase_text, is.finite(.data$dist_to_next))
	if (nrow(sub) < 3 || all(is.na(sub$dist_to_next))) {
		warning(sprintf("Not enough data to plot/correlate for phase %s", phase_text))
		return(invisible(NULL))
	}
	ct  <- suppressWarnings(cor.test(sub$dist_to_next, sub$proportion, method = "spearman"))
	rho <- unname(ct$estimate)
	p   <- ct$p.value
	lab <- sprintf("Spearman rho = %.3f\np = %.3g", rho, p)
	
	x_pos <- max(sub$dist_to_next, na.rm = TRUE)
	y_pos <- max(sub$proportion,   na.rm = TRUE)
	
	p_out <- ggplot2::ggplot(sub, ggplot2::aes(x = .data$dist_to_next, y = .data$proportion)) +
		ggplot2::geom_point(size = 4, alpha = 0.85, color = color_hex) +
		ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.7, color = color_hex) +
		ggplot2::annotate("text", x = x_pos, y = y_pos, label = lab,
								hjust = 1.05, vjust = 1.35, size = 6) +
		ggplot2::labs(
			title = "",
			x = "Cross–nearest-neighbor to next stage (µm)",
			y = paste0("% ", phase_label)
		) +
		ggplot2::theme_bw(base_size = 12) +
		ggplot2::theme(
			legend.position = "none",
			plot.title = ggplot2::element_text(face = "bold", size = 20, hjust = 0.5),
			axis.title = ggplot2::element_text(size = 18),
			axis.text  = ggplot2::element_text(size = 16)
		)
	
	fname <- file.path(outdir,
							 sprintf("%s_%s.png",
							 		  prefix, gsub("[^A-Za-z0-9]+", "_", phase_text)))
	ggplot2::ggsave(filename = fname, plot = p_out,
						 width = 5, height = 5, units = "in", dpi = 300)
	message("Saved: ", fname)
	invisible(fname)
}

# Save the three phase plots
save_phase_plot(phase_dist, "G1",  "G0/G1", phase_colors[["G1"]],  file.path(wd.de.plots, "origUM"))
save_phase_plot(phase_dist, "S",   "S",     phase_colors[["S"]],   file.path(wd.de.plots, "origUM"))
save_phase_plot(phase_dist, "G2M", "G2/M",  phase_colors[["G2M"]], file.path(wd.de.plots, "origUM"))

# Also write the per-stage distance summary used for plotting
write.table(
	dist_summary,
	file = file.path(wd.de.plots, "origUM", paste0("Spatial_ST_META_", GENE, "_XNND.txt")),
	sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
)

save(so.st, meta_idx, spatial_data, spatial_coords, cluster_centers, dist_map, meta_cells_only, phase_proportions, dist_summary, phase_dist, cor_results_obs, ct_order_master, short_labels, ct_levels, phase_colors, file=file.path(wd.de.plots, "origUM", paste0("ssc_filtered_normalised_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_ordered_annotated_meta_cells_only_", GENE, "_RNA_so.st.RData")))

# -----------------------------------------------------------------------------
# After running 08_ssc_Slide-tag_CBL.R
# # Last Modified: 09/10/25
# -----------------------------------------------------------------------------
library(dplyr)

load(file=file.path(wd.de.plots, "origUM", paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_23_res=0.5_-C0_ordered_annotated_meta_cells_only_", GENE, "_RNA.RData")))

cor_results_obs <- phase_dist %>%
	dplyr::group_by(phase) %>%
	dplyr::summarise(
		rho  = cor(dist_to_next, proportion, use = "complete.obs"),
		pval = cor.test(dist_to_next, proportion)$p.value,
		.groups = "drop"
	) %>%
	dplyr::mutate(b = 0L) %>%
	dplyr::select(b, phase, rho, pval)

# -----------------------------------------------
# Permutation loop
# -----------------------------------------------
B <- 10000L  # number of permutations
res_list <- vector("list", B)  # each element will be a tibble with 3 rows (phases)

stage_labels <- meta_cells_only$cell_type  # original labels; this vector encodes stage counts

for (b in seq_len(B)) {
	meta_cells_only_perm <- meta_cells_only
	meta_cells_only_perm$cell_type <- sample(stage_labels, length(stage_labels), replace = FALSE)
	
	# --- define order
	cell_types <- c("Stage 0","Stage 0A","Stage 0B","Stage 1","Stage 2","Stage 3","Leptotene","Zygotene","Pachytene","Diplotene","Meiotic division","Early spermatid","Late spermatid")
	cell_type_order <- base::intersect(cell_types, as.vector(unique(meta_cells_only_perm$cell_type)))
	
	meta_cells_only_perm <- meta_cells_only_perm %>%
		dplyr::mutate(cell_type = factor(cell_type, levels = cell_type_order)) %>%
		dplyr::arrange(cell_type)
	
	dist_map$dist_to_next <- mapply(function(f, t) avg_dist(meta_cells_only_perm, f, t),
											  dist_map$from, dist_map$to)
	
	# --- join back (no rename needed)
	meta_cells_only_perm <- meta_cells_only_perm |>
		dplyr::left_join(dplyr::select(dist_map, from, dist_to_next),
							  by = c("cell_type" = "from"))
	
	# --- aggregate distances to per-cell_type (unique value for each "from")
	dist_summary_perm <- meta_cells_only_perm %>%
		dplyr::select(cell_type, dist_to_next.y) %>%
		dplyr::distinct()
	
	# --- join with phase proportions
	phase_dist_perm <- phase_proportions %>%
		dplyr::left_join(dist_summary_perm, by = "cell_type")
	
	cor_results <- phase_dist_perm %>%
		group_by(phase) %>%
		summarise(
			rho = cor(dist_to_next.y, proportion, use = "complete.obs"),
			pval = cor.test(dist_to_next.y, proportion)$p.value,
			.groups = "drop"
		)
	
	res_list[[b]] <- cor_results
}

# -----------------------------------------------
# Combine & (optionally) save
# -----------------------------------------------
perm_cor_long <- dplyr::bind_rows(res_list)

perm_cor_long <- dplyr::bind_rows(
	Map(function(tbl, i) dplyr::mutate(tbl, b = as.integer(i)), res_list, seq_along(res_list))
) %>% dplyr::mutate(phase = as.character(phase)) %>%
	dplyr::select(b, phase, dplyr::everything())

# prepend the observed (b = 0) if you want both together
perm_cor_long_with_obs <- dplyr::bind_rows(cor_results_obs, perm_cor_long)

# wide convenience tables
perm_cor_rho_wide <- tidyr::pivot_wider(perm_cor_long, id_cols = b, names_from = phase, values_from = rho)
perm_cor_p_wide   <- tidyr::pivot_wider(perm_cor_long, id_cols = b, names_from = phase, values_from = pval)

# write to disk
readr::write_tsv(perm_cor_long, file.path(wd.de.plots, "origUM", "perm_cor_by_phase_long_RNA.tsv"))
readr::write_tsv(perm_cor_rho_wide, file.path(wd.de.plots, "origUM", "perm_cor_rho_wide_RNA.tsv"))
readr::write_tsv(perm_cor_p_wide,   file.path(wd.de.plots, "origUM", "perm_cor_p_wide_RNA.tsv"))

# -----------------------------------------------
# Example: empirical p
# -----------------------------------------------
# Directional rule per phase
phase_direction <- c(G1 = "less", S = "greater", G2M = "two.sided")

emp_pvals_one_phase_dir <- function(phase_name, direction = c("less","greater","two.sided")) {
	direction <- match.arg(direction)
	# observed
	rho_obs <- perm_cor_long_with_obs %>% dplyr::filter(b == 0, phase == phase_name) %>% pull(rho)
	p_obs   <- perm_cor_long_with_obs %>% dplyr::filter(b == 0, phase == phase_name) %>% pull(pval)
	# permuted
	rho_perm <- perm_cor_long_with_obs %>% dplyr::filter(b > 0, phase == phase_name) %>% pull(rho)
	p_perm   <- perm_cor_long_with_obs %>% dplyr::filter(b > 0, phase == phase_name) %>% pull(pval)
	
	# directional permutation p on rho
	p_perm_r <- switch(direction,
							 greater   = (sum(rho_perm >= rho_obs, na.rm = TRUE) + 1) / (length(rho_perm) + 1),
							 less      = (sum(rho_perm <= rho_obs, na.rm = TRUE) + 1) / (length(rho_perm) + 1),
							 two.sided = (sum(abs(rho_perm) >= abs(rho_obs), na.rm = TRUE) + 1) / (length(rho_perm) + 1)
	)
	
	# p-of-p version (kept as <= observed p, regardless of direction)
	p_perm_p <- (sum(p_perm <= p_obs, na.rm = TRUE) + 1) / (length(p_perm) + 1)
	
	dplyr::tibble(
		phase = phase_name,
		rho_obs = rho_obs,
		p_obs   = p_obs,
		perm_tail = direction,
		p_perm_r = p_perm_r,
		p_perm_p = p_perm_p
	)
}

phases <- c("G1","S","G2M")
emp_dir <- bind_rows(lapply(phases, function(ph) {
	emp_pvals_one_phase_dir(ph, direction = phase_direction[[ph]])
})) %>%
	mutate(
		q_perm_r = p.adjust(p_perm_r, method = "BH"),
		q_perm_p = p.adjust(p_perm_p, method = "BH")
	)

print(emp_dir)
readr::write_tsv(emp_dir, file.path(wd.de.plots, "origUM", "phase_empirical_pvalues_directional_RNA.tsv"))

save(phase_colors, stage_labels, res_list, perm_cor_long, perm_cor_rho_wide, perm_cor_p_wide, phase_direction, emp_dir, file=file.path(wd.de.plots, "origUM", paste0("ssc_filtered_normalised_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_ordered_annotated_meta_cells_only_", GENE, "_RNA_perm.RData")))

# -----------------------------------------------------------------------------
# To test correlation between distances (dist_to_next) and cell cycle phase proportions (G1/S/G2M)
# Last Modified: 09/10/25
# -----------------------------------------------------------------------------
# truncate to k decimals (works for +/- numbers)
trunc_dec <- function(x, k) trunc(x * 10^k) / 10^k

# Helper to make & save one phase plot to a specific directory with a clear filename
save_phase_plot_perm <- function(df, phase_text, phase_label, color_hex, outdir, emp_dir,
											prefix = "Spatial_ST_META_CBL_permutation_trunc_dec_RNA_test") {
	sub <- df %>% dplyr::filter(.data$phase == phase_text, is.finite(.data$dist_to_next))
	if (nrow(sub) < 3 || all(is.na(sub$dist_to_next))) {
		warning(sprintf("Not enough data to plot/correlate for phase %s", phase_text))
		return(invisible(NULL))
	}
	
	ct  <- stats::cor.test(sub$dist_to_next, sub$proportion)
	rho <- base::unname(ct$estimate)
	p   <- ct$p.value
	lab <- paste0("rho = ", round(rho, 2), "; p = ", trunc_dec(p, 4))
	lab <- paste0(lab, "\nempirical p = ", trunc_dec(subset(emp_dir, phase == phase_text)$p_perm_r, 4))
	
	x_pos <- max(sub$dist_to_next, na.rm = TRUE)
	y_pos <- max(sub$proportion,   na.rm = TRUE)
	
	p_out <- ggplot2::ggplot(sub, ggplot2::aes(x = .data$dist_to_next, y = .data$proportion)) +
		ggplot2::geom_point(size = 4, alpha = 0.8, color = color_hex) +
		ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, color = color_hex) +
		ggplot2::annotate("text", x = x_pos, y = y_pos, label = lab,
								hjust = 1.05, vjust = 1.4, size = 6) +
		ggplot2::labs(
			title = "",
			x = "Spatial distance (µm)",
			y = paste0("% ", phase_label)
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
save_phase_plot_perm(phase_dist, "G1",  "G0/G1", phase_colors[["G1"]],  file.path(wd.de.plots, "origUM"), emp_dir)
save_phase_plot_perm(phase_dist, "S",   "S",     phase_colors[["S"]],   file.path(wd.de.plots, "origUM"), emp_dir)
save_phase_plot_perm(phase_dist, "G2M", "G2/M",  phase_colors[["G2M"]], file.path(wd.de.plots, "origUM"), emp_dir)











# =============================================================================
# Table: per-stage callable counts and WT availability in space (MUT header)
# Columns: Cell type | Total callable cells | Expected | Spatial available | Missing
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(tibble); library(stringr)})

# --- helpers / joins ----------------------------------------------------------
# Make sure spatial_data has a CB key to join on
spatial_data_cb <- spatial_data
if (!"CB" %in% names(spatial_data_cb)) spatial_data_cb <- spatial_data_cb %>% mutate(CB = cell)

# Mutant barcodes to exclude from WT pool
mut_CBs <- unique(if ("CB" %in% names(meta_cells_only)) meta_cells_only$CB else meta_cells_only$cell)

# 1) Total callable cells per stage (from calls2)
tot_callable <- calls2 %>%
	filter(callable, stage %in% stage_order) %>%
	count(stage, name = "Total callable cells") %>%
	mutate(stage = factor(as.character(stage), levels = stage_order))

# 2) Expected WT per stage = n - k (from your callable-only k/n table)
expected_wt <- k_n_by_stage %>%
	transmute(stage = factor(as.character(stage), levels = stage_order),
				 Expected = n - k)

# 3) Spatial-available callable WT per stage
#    (callable cells, NOT mutant CBs, that have coords in spatial_data)
callable_wt <- calls2 %>%
	filter(callable, stage %in% stage_order) %>%
	transmute(CB, stage = factor(as.character(stage), levels = stage_order)) %>%
	filter(!(CB %in% mut_CBs))

join_probe <- callable_wt %>%
	left_join(spatial_data_cb %>% select(CB, plot_x, plot_y), by = "CB") %>%
	mutate(has_coords = is.finite(plot_x) & is.finite(plot_y))

spatial_avail <- join_probe %>%
	filter(has_coords) %>%
	count(stage, name = "Spatial available")

# 4) Assemble table, fill zeros for missing stages, compute Missing, add Total row
tbl_core <- full_join(tot_callable, expected_wt, by = "stage") %>%
	full_join(spatial_avail, by = "stage") %>%
	arrange(stage) %>%
	mutate(
		`Total callable cells` = coalesce(`Total callable cells`, 0L),
		Expected               = coalesce(Expected,               0L),
		`Spatial available`    = coalesce(`Spatial available`,    0L),
		Missing                = pmax(Expected - `Spatial available`, 0L)
	) %>%
	# keep only your 12 stages, in order
	filter(!is.na(stage)) %>%
	mutate(`Cell type` = as.character(stage)) %>%
	select(`Cell type`, `Total callable cells`, Expected, `Spatial available`, Missing)

# Totals row
tot_row <- tibble(
	`Cell type`            = "Total",
	`Total callable cells` = sum(tbl_core$`Total callable cells`, na.rm = TRUE),
	Expected               = sum(tbl_core$Expected,               na.rm = TRUE),
	`Spatial available`    = sum(tbl_core$`Spatial available`,    na.rm = TRUE),
	Missing                = sum(tbl_core$Missing,                na.rm = TRUE)
)

tbl_out <- bind_rows(tbl_core, tot_row)

# Print and (optionally) save
print(tbl_out, n = nrow(tbl_out))
readr::write_tsv(tbl_out, file.path(wd.de.plots, "origUM", "callable_vs_spatial_WT_table.tsv"))






























# -----------------------------------------------------------------------------
# Cell cycle analysis (PD53626b_ST1)
# Last Modified: 08/10/25
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# To calculate the cell cycle phase proportions (G1/S/G2M) in this sample
# Last Modified: 08/10/25
# -----------------------------------------------------------------------------
so.integrated <- so.integrated.st
#so.integrated <- so.st

# Identify the clusters you want to keep
clusters_to_remove <- c("Sertoli", "Fibrotic PMC", "PMC", "Leydig", "Endothelial and Macrophage", "4", "20")
clusters_to_keep <- base::setdiff(unique(Idents(so.integrated)), clusters_to_remove)

so.integrated$idents <- Idents(so.integrated) 
so.integrated <- subset(so.integrated, idents = clusters_to_keep)

# Rename specific cluster identity
Idents(so.integrated) <- plyr::mapvalues(
	Idents(so.integrated),
	from = "Meiotic division",
	to   = "Meiotic\ndivision"
)
Idents(so.integrated) <- plyr::mapvalues(
	Idents(so.integrated),
	from = "Early spermatid",
	to   = "Early\nspermatid"
)

dim_plot <- DimPlot(so.integrated, label = TRUE) + NoLegend() +
	labs(x = "UMAP 1", y = "UMAP 2") +  # Set x-axis and y-axis labels
	theme(
		plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
		axis.title = element_text(size = 18),  # Increase axis title size
		axis.text = element_text(size = 16),  # Increase axis tick label size
	) +
	ggtitle("PD53626b_ST1")
ggsave(file.path(wd.de.plots, paste0("DimPlot_SSC_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_6*6_SSC_PD53626b_ST1.png")), plot = dim_plot, width = 6, height = 6, dpi = 300)

Idents(so.integrated) <- plyr::mapvalues(
	Idents(so.integrated),
	from = "Meiotic\ndivision",
	to   = "Meiotic division"
)
Idents(so.integrated) <- plyr::mapvalues(
	Idents(so.integrated),
	from = "Early\nspermatid",
	to   = "Early spermatid"
)

# -----------------------------------------------------------------------------
# Cell cycle analysis
# Last Modified: 10/02/25
# -----------------------------------------------------------------------------
cc.genes <- Seurat::cc.genes.updated.2019  # Updated gene lists
s.genes <- cc.genes$s.genes  # S phase markers
g2m.genes <- cc.genes$g2m.genes  # G2/M phase markers

so.integrated <- CellCycleScoring(so.integrated, 
											 s.features = s.genes, 
											 g2m.features = g2m.genes, 
											 set.ident = F)  # Updates the cell identities (Idents(seurat_obj)) with the assigned cell cycle phases if set.ident = TRUE.

# Ensure Phase is a factor and ordered correctly
so.integrated@meta.data$Phase <- factor(so.integrated@meta.data$Phase, levels = c("G1", "S", "G2M"))

##
library(ggplot2)
library(Seurat)
library(dplyr)
library(ggrepel)

# Ensure `seurat_clusters` exists
so.integrated$seurat_clusters <- Idents(so.integrated)

# Define consistent Phase colors
phase_colors <- c("G1" = "#619CFF", "S" = "#F8766D", "G2M" = "#00BA38")

# Extract UMAP embeddings
umap_data <- as.data.frame(Embeddings(so.integrated, reduction = "umap"))
umap_data$Phase <- so.integrated$Phase  # Add this line!
umap_data$seurat_clusters <- as.factor(so.integrated$seurat_clusters)  # Ensure clusters are factors

# Compute cluster centroids
cluster_centers <- umap_data %>%
	group_by(seurat_clusters) %>%
	summarise(UMAP_1 = median(umap_1), UMAP_2 = median(umap_2))  # Use median to find label position

# Create UMAP plot with manual ggplot (same style as your loop code)
dim_plot <- ggplot() +
	# Plot all cells with their respective phase colors
	geom_point(data = umap_data, 
				  aes(x = umap_1, y = umap_2, color = Phase), 
				  size = 0.5, alpha = 0.8) +
	# Add cluster labels
	geom_text_repel(data = cluster_centers, 
						 aes(x = UMAP_1, y = UMAP_2, label = seurat_clusters), 
						 size = 4, color = "black") +
	# Manual color scale
	scale_color_manual(values = phase_colors, 
							 breaks = c("G1", "S", "G2M"), 
							 labels = c("G0/G1", "S", "G2/M"),
							 name = "Phase") +
	labs(x = "UMAP 1", y = "UMAP 2") +
	theme_classic() +
	theme(
		plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
		axis.title = element_text(size = 18),  
		axis.text = element_text(size = 16),  
		legend.text = element_text(size = 16),  
		legend.title = element_text(size = 16),
		legend.position = c(0.95, 0.02),
		legend.justification = c("right", "bottom")
	) +
	ggtitle("PD53626b_ST1") +
	guides(color = guide_legend(override.aes = list(size = 3)))

# Save the final plot
ggsave(file.path(wd.de.plots, paste0("Cycle_SSC_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_6x6_so.st.png")), 
		 plot = dim_plot, width = 6, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# Bar chart: Cell cycle phase proportions by cell type
# Cell types on y-axis, proportions on x-axis
# Phase order: G1, S, G2/M from left to right
# No percentage labels on bars
# Last Modified: 15/09/25
# -----------------------------------------------------------------------------
library(ggplot2)
library(Seurat)
library(dplyr)

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

# Define phase colors
phase_colors <- c("G1" = "#619CFF", "S" = "#F8766D", "G2M" = "#00BA38")

# Create stacked bar chart without percentage labels
bar_plot <- ggplot(phase_proportions, aes(x = proportion, y = cell_type, fill = phase)) +
	geom_bar(stat = "identity", width = 0.7) +
	scale_fill_manual(values = phase_colors, 
							breaks = c("G1", "S", "G2M"), 
							labels = c("G0/G1", "S", "G2/M")) +
	scale_x_continuous(labels = scales::percent_format(), 
							 expand = c(0, 0)) +
	labs(
		x = "Proportion", 
		y = "Cell type", 
		fill = "Phase"
	) +
	theme_classic() +
	theme(
		plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
		axis.title = element_text(size = 18),
		axis.text = element_text(size = 16),
		legend.text = element_text(size = 16),
		legend.title = element_text(size = 16),
		legend.position = "right",
		panel.grid.major.x = element_line(color = "grey90", size = 0.5),
		panel.grid.minor.x = element_line(color = "grey95", size = 0.25)
	) +
	ggtitle("PD53626b_ST1")

# Save the plot
ggsave(file.path(wd.de.plots, "cell_cycle_proportions_by_cell_typ_6x14_PD53626b_ST1.png"), 
		 plot = bar_plot, width = 14, height = 6, dpi = 300)

# Print the proportions table
print(phase_proportions)

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
									 prefix = "Spatial_ST_META_CBL_PD53626b_ST1") {
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
			x = "Spatial distance (µm)",
			y = paste0("% ", phase_label)
		) +
		ggplot2::theme_bw(base_size = 12) +
		ggplot2::theme(
			legend.position = "none",
			plot.title = ggplot2::element_text(face = "bold", size = 20, hjust = 0.5),
			axis.title = ggplot2::element_text(size = 18),
			axis.text  = ggplot2::element_text(size = 16)
		) +
		ggtitle("PD53626b_ST1")
	
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
	file = file.path(wd.de.plots, "origUM", paste0("Spatial_ST_META_", GENE, "_PD53626b_ST1.txt")),
	sep = "\t",
	quote = FALSE,
	row.names = FALSE,
	col.names = TRUE
)

# -----------------------------------------------------------------------------
# After running 08_ssc_Slide-tag_CBL.R
# # Last Modified: 09/10/25
# -----------------------------------------------------------------------------
cor_results_obs <- phase_dist %>%
	dplyr::group_by(phase) %>%
	dplyr::summarise(
		rho  = cor(dist_to_next, proportion, use = "complete.obs"),
		pval = cor.test(dist_to_next, proportion)$p.value,
		.groups = "drop"
	) %>%
	dplyr::mutate(b = 0L) %>%
	dplyr::select(b, phase, rho, pval)

# -----------------------------------------------
# Permutation loop
# -----------------------------------------------
B <- 10000L  # number of permutations
res_list <- vector("list", B)  # each element will be a tibble with 3 rows (phases)

for (b in seq_len(B)) {
	meta_cells_only_perm <- meta_cells_only
	meta_cells_only_perm$cell_type <- sample(stage_labels, length(stage_labels), replace = FALSE)
	
	# --- define order
	cell_types <- c("Stage 0","Stage 0A","Stage 0B","Stage 1","Stage 2","Stage 3","Leptotene","Zygotene","Pachytene","Diplotene","Meiotic division","Early spermatid","Late spermatid")
	cell_type_order <- base::intersect(cell_types, as.vector(unique(meta_cells_only_perm$cell_type)))
	
	meta_cells_only_perm <- meta_cells_only_perm %>%
		dplyr::mutate(cell_type = factor(cell_type, levels = cell_type_order)) %>%
		dplyr::arrange(cell_type)
	
	dist_map$dist_to_next <- mapply(function(f, t) avg_dist(meta_cells_only_perm, f, t),
											  dist_map$from, dist_map$to)
	
	# --- join back (no rename needed)
	meta_cells_only_perm <- meta_cells_only_perm |>
		dplyr::left_join(dplyr::select(dist_map, from, dist_to_next),
							  by = c("cell_type" = "from"))
	
	# --- aggregate distances to per-cell_type (unique value for each "from")
	dist_summary_perm <- meta_cells_only_perm %>%
		dplyr::select(cell_type, dist_to_next.y) %>%
		dplyr::distinct()
	
	# --- join with phase proportions
	phase_dist_perm <- phase_proportions %>%
		dplyr::left_join(dist_summary_perm, by = "cell_type")
	
	cor_results <- phase_dist_perm %>%
		group_by(phase) %>%
		summarise(
			rho = cor(dist_to_next.y, proportion, use = "complete.obs"),
			pval = cor.test(dist_to_next.y, proportion)$p.value,
			.groups = "drop"
		)
	
	res_list[[b]] <- cor_results
}

# -----------------------------------------------
# Combine & (optionally) save
# -----------------------------------------------
perm_cor_long <- dplyr::bind_rows(res_list)

perm_cor_long <- dplyr::bind_rows(
	Map(function(tbl, i) dplyr::mutate(tbl, b = as.integer(i)), res_list, seq_along(res_list))
) %>% dplyr::mutate(phase = as.character(phase)) %>%
	dplyr::select(b, phase, dplyr::everything())

# prepend the observed (b = 0) if you want both together
perm_cor_long_with_obs <- dplyr::bind_rows(cor_results_obs, perm_cor_long)

# wide convenience tables
perm_cor_rho_wide <- tidyr::pivot_wider(perm_cor_long, id_cols = b, names_from = phase, values_from = rho)
perm_cor_p_wide   <- tidyr::pivot_wider(perm_cor_long, id_cols = b, names_from = phase, values_from = pval)

# write to disk
readr::write_tsv(perm_cor_long, file.path(wd.de.plots, "origUM", "perm_cor_by_phase_long_PD53626b_ST1.tsv"))
readr::write_tsv(perm_cor_rho_wide, file.path(wd.de.plots, "origUM", "perm_cor_rho_wide_PD53626b_ST1.tsv"))
readr::write_tsv(perm_cor_p_wide,   file.path(wd.de.plots, "origUM", "perm_cor_p_wide_PD53626b_ST1.tsv"))

# -----------------------------------------------
# Example: empirical p
# -----------------------------------------------
suppressPackageStartupMessages({library(dplyr)})

# Directional rule per phase
phase_direction <- c(G1 = "less", S = "greater", G2M = "two.sided")

emp_pvals_one_phase_dir <- function(phase_name, direction = c("less","greater","two.sided")) {
	direction <- match.arg(direction)
	# observed
	rho_obs <- perm_cor_long_with_obs %>% filter(b == 0, phase == phase_name) %>% pull(rho)
	p_obs   <- perm_cor_long_with_obs %>% filter(b == 0, phase == phase_name) %>% pull(pval)
	# permuted
	rho_perm <- perm_cor_long_with_obs %>% filter(b > 0, phase == phase_name) %>% pull(rho)
	p_perm   <- perm_cor_long_with_obs %>% filter(b > 0, phase == phase_name) %>% pull(pval)
	
	# directional permutation p on rho
	p_perm_r <- switch(direction,
							 greater   = (sum(rho_perm >= rho_obs, na.rm = TRUE) + 1) / (length(rho_perm) + 1),
							 less      = (sum(rho_perm <= rho_obs, na.rm = TRUE) + 1) / (length(rho_perm) + 1),
							 two.sided = (sum(abs(rho_perm) >= abs(rho_obs), na.rm = TRUE) + 1) / (length(rho_perm) + 1)
	)
	
	# p-of-p version (kept as <= observed p, regardless of direction)
	p_perm_p <- (sum(p_perm <= p_obs, na.rm = TRUE) + 1) / (length(p_perm) + 1)
	
	tibble(
		phase = phase_name,
		rho_obs = rho_obs,
		p_obs   = p_obs,
		perm_tail = direction,
		p_perm_r = p_perm_r,
		p_perm_p = p_perm_p
	)
}

phases <- c("G1","S","G2M")
emp_dir <- bind_rows(lapply(phases, function(ph) {
	emp_pvals_one_phase_dir(ph, direction = phase_direction[[ph]])
})) %>%
	mutate(
		q_perm_r = p.adjust(p_perm_r, method = "BH"),
		q_perm_p = p.adjust(p_perm_p, method = "BH")
	)

print(emp_dir)
readr::write_tsv(emp_dir, file.path(wd.de.plots, "origUM", "phase_empirical_pvalues_directional_PD53626b_ST1.tsv"))

# -----------------------------------------------------------------------------
# To test correlation between distances (dist_to_next) and cell cycle phase proportions (G1/S/G2M)
# Last Modified: 09/10/25
# -----------------------------------------------------------------------------
# truncate to k decimals (works for +/- numbers)
trunc_dec <- function(x, k) trunc(x * 10^k) / 10^k

# Helper to make & save one phase plot to a specific directory with a clear filename
save_phase_plot_perm <- function(df, phase_text, phase_label, color_hex, outdir, emp_dir,
											prefix = "Spatial_ST_META_CBL_permutation_trunc_dec_PD53626b_ST1") {
	sub <- df %>% dplyr::filter(.data$phase == phase_text, is.finite(.data$dist_to_next))
	if (nrow(sub) < 3 || all(is.na(sub$dist_to_next))) {
		warning(sprintf("Not enough data to plot/correlate for phase %s", phase_text))
		return(invisible(NULL))
	}
	
	ct  <- stats::cor.test(sub$dist_to_next, sub$proportion)
	rho <- base::unname(ct$estimate)
	p   <- ct$p.value
	lab <- paste0("rho = ", round(rho, 2), "; p = ", trunc_dec(p, 4))
	lab <- paste0(lab, "\nEmpirical p = ", trunc_dec(subset(emp_dir, phase == phase_text)$p_perm_r, 4))
	
	x_pos <- max(sub$dist_to_next, na.rm = TRUE)
	y_pos <- max(sub$proportion,   na.rm = TRUE)
	
	p_out <- ggplot2::ggplot(sub, ggplot2::aes(x = .data$dist_to_next, y = .data$proportion)) +
		ggplot2::geom_point(size = 4, alpha = 0.8, color = color_hex) +
		ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, color = color_hex) +
		ggplot2::annotate("text", x = x_pos, y = y_pos, label = lab,
								hjust = 1.05, vjust = 1.4, size = 6) +
		ggplot2::labs(
			title = "PD53626b_ST1",
			x = "Spatial distance (µm)",
			y = paste0("% ", phase_label)
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
save_phase_plot_perm(phase_dist, "G1",  "G0/G1", phase_colors[["G1"]],  file.path(wd.de.plots, "origUM"), emp_dir)
save_phase_plot_perm(phase_dist, "S",   "S",     phase_colors[["S"]],   file.path(wd.de.plots, "origUM"), emp_dir)
save_phase_plot_perm(phase_dist, "G2M", "G2/M",  phase_colors[["G2M"]], file.path(wd.de.plots, "origUM"), emp_dir)
