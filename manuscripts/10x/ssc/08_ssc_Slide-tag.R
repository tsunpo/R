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
wd.de.plots <- file.path(wd.de, "plots/ARID1A")
dir.create(wd.de.plots, showWarnings = FALSE)

load(file.path(wd.de.data, paste0("ssc_filtered_normalised_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_ordered_annotated.RData")))

so.integrated$cell_type <- Idents(so.integrated)

# -----------------------------------------------------------------------------
# After running SComatic genotyping
# # Last Modified: 17/09/25
# -----------------------------------------------------------------------------
wd.nf.data  <- file.path(wd.rna, "/ngs/pacbio/nf-scomatic")
barcodes <- read.table(
   file = file.path(wd.nf.data, "scomatics_PD53626b.tsv"),
   sep = "\t",
   header = T,
   stringsAsFactors = FALSE
)

meta <- read.table(
   file = file.path(wd.nf.data, "celltypes.tsv"),
   sep = "\t",
   header = T,
   stringsAsFactors = FALSE
)

barcodes.st <- subset(subset(barcodes, CHROM == "chr1"), Base_observed == "A")
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
# 
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

# --- Gather embeddings
spat <- Embeddings(so.st, "spatial")
dbs  <- Embeddings(so.st, "dbscan")
kde  <- Embeddings(so.st, "kde")

# --- Align rows (cells) across reductions
cells_union <- Reduce(union, list(rownames(spat), rownames(dbs), rownames(kde)))
pad <- function(mat, cells) {
	out <- matrix(NA_real_, nrow = length(cells), ncol = ncol(mat),
					  dimnames = list(cells, colnames(mat)))
	hit <- intersect(cells, rownames(mat))
	out[hit, ] <- mat[hit, , drop = FALSE]
	out
}
spatA <- pad(spat, cells_union)
dbsA  <- pad(dbs,  cells_union)
kdeA  <- pad(kde,  cells_union)

# --- Row-wise NA flags (TRUE if any NA in that row for that reduction)
rowNA_spat <- apply(spatA, 1, function(x) any(is.na(x)))
rowNA_dbs  <- apply(dbsA,  1, function(x) any(is.na(x)))
rowNA_kde  <- apply(kdeA,  1, function(x) any(is.na(x)))

# --- Summary table per reduction
na_summary <- data.frame(
	reduction   = c("spatial","dbscan","kde"),
	NA_rows     = c(sum(rowNA_spat), sum(rowNA_dbs), sum(rowNA_kde)),
	complete_rows = c(sum(!rowNA_spat), sum(!rowNA_dbs), sum(!rowNA_kde)),
	total_rows  = length(cells_union)
)
na_summary$prop_NA <- round(na_summary$NA_rows / na_summary$total_rows, 4)
na_summary

# --- Intersections: which combos are NA?
flags <- data.frame(
	cell = cells_union,
	spatial_NA = rowNA_spat,
	dbscan_NA  = rowNA_dbs,
	kde_NA     = rowNA_kde,
	stringsAsFactors = FALSE
)

# Human-readable pattern label per cell
flags$pattern <- paste0(
	ifelse(flags$spatial_NA, "spatial=NA", "spatial=OK"), " | ",
	ifelse(flags$dbscan_NA,  "dbscan=NA",  "dbscan=OK"),  " | ",
	ifelse(flags$kde_NA,     "kde=NA",     "kde=OK")
)

# Counts by exact NA/OK pattern
pattern_counts <- as.data.frame(sort(table(flags$pattern), decreasing = TRUE))
colnames(pattern_counts) <- c("pattern", "n_cells")
pattern_counts

# --- Distribution: how many reductions are NA per cell? (0–3)
flags$n_reductions_NA <- rowSums(flags[, c("spatial_NA","dbscan_NA","kde_NA")])
na_count_distribution <- as.data.frame(table(flags$n_reductions_NA))
colnames(na_count_distribution) <- c("num_reductions_with_NA", "n_cells")
na_count_distribution

# (Optional) write to CSV
# write.csv(na_summary, "na_summary_by_reduction.csv", row.names = FALSE)
# write.csv(pattern_counts, "na_intersection_patterns.csv", row.names = FALSE)
# write.csv(flags, "na_flags_per_cell.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# NAs in spatial, dbscan, and kde
# Last Modified: 15/09/25
# -----------------------------------------------------------------------------
count_na_rows <- function(mat) {
	na_row <- rowSums(is.na(mat)) > 0
	data.frame(
		total_rows     = nrow(mat),
		na_rows        = sum(na_row),
		complete_rows  = nrow(mat) - sum(na_row),
		prop_na_rows   = round(mean(na_row), 4)
	)
}

spatial_mat <- Embeddings(so.st, "spatial")
dbscan_mat  <- Embeddings(so.st, "dbscan")
kde_mat     <- Embeddings(so.st, "kde")

na_by_reduction <- dplyr::bind_rows(
	spatial = count_na_rows(spatial_mat),
	dbscan  = count_na_rows(dbscan_mat),
	kde     = count_na_rows(kde_mat),
	.id = "reduction"
)

na_by_reduction
# reduction total_rows na_rows complete_rows prop_na_rows
# 1   spatial       7772    5997          1775       0.7716
# 2    dbscan       7772    3977          3795       0.5117
# 3       kde       7772    5434          2338       0.6992

# where to save
out_csv <- file.path(wd.de.plots, "NA_by_reduction.csv")

# tidy + write
na_by_reduction %>%
	mutate(prop_na_rows = round(prop_na_rows, 4)) %>%
	write_csv(out_csv)

cat("Saved:", out_csv, "\n")

# -----------------------------------------------------------------------------
# Plot cell-type-by-cell-type - Each cell type highlighted separately in a loop
# Last Modified: 15/09/25
# -----------------------------------------------------------------------------
library(ggplot2)
library(Seurat)
library(dplyr)
library(ggrepel)

# Ensure `seurat_clusters` exists
so.st$seurat_clusters <- Idents(so.st)

# Extract spatial information from so.st
spatial_data <- as.data.frame(Embeddings(so.st, reduction = "dbscan"))
spatial_data$cell_type <- so.st$cell_type
spatial_data$seurat_clusters <- as.factor(so.st$seurat_clusters)

# Compute cluster centroids (same for all plots)
# First check the column names of spatial coordinates
print("Spatial coordinate column names:")
print(colnames(spatial_data))

# Use the actual column names (likely spatial_1 and spatial_2)
spatial_coords <- colnames(spatial_data)[1:2]  # First two columns should be coordinates
cluster_centers <- spatial_data %>%
 group_by(seurat_clusters) %>%
 summarise(coord1 = median(.data[[spatial_coords[1]]]), 
           coord2 = median(.data[[spatial_coords[2]]]))

# Get unique cell types to plot
cell_types_to_plot <- unique(spatial_data$cell_type)
cell_types_to_plot <- cell_types_to_plot[!is.na(cell_types_to_plot)]  # Remove NA values

# Use the colors from so.integrated (same as your previous code)
cell_type_colors <- Idents(so.integrated)  # Extract the colors used for each cell type
# Extract the colors for each identity (cell type) programmatically
plot_colors <- scales::hue_pal()(length(unique(cell_type_colors)))
names(plot_colors) <- levels(cell_type_colors)

# Build a palette keyed to the cell types you actually plot
all_ct <- sort(unique(na.omit(spatial_data$cell_type)))
pal <- setNames(scales::hue_pal()(length(all_ct)), all_ct)

for (ct in all_ct) {
 target_cells <- dplyr::filter(spatial_data, .data$cell_type == .env$ct)
 other_cells  <- dplyr::filter(spatial_data, .data$cell_type != .env$ct | is.na(.data$cell_type))
 
 dim_plot <- ggplot() +
  # background (others) in grey
  geom_point(data = other_cells,
             aes(x = .data[[spatial_coords[1]]], y = .data[[spatial_coords[2]]]),
             color = "#D3D3D3", fill = "white", size = 3, alpha = 1, shape = 21) +
  # highlighted target cells
  geom_point(data = target_cells,
             aes(x = .data[[spatial_coords[1]]], y = .data[[spatial_coords[2]]], color = cell_type),
             size = 3, alpha = 1) +
  geom_text_repel(data = cluster_centers,
                  aes(x = coord1, y = coord2, label = seurat_clusters), size = 4, color = "black") +
  # only one legend entry: the current ct color
  scale_color_manual(values = pal[ct], breaks = ct, name = "Cell type") +
  coord_fixed() +
  labs(x = "X", y = "Y") +
  theme_classic() +
  theme(
   plot.title = element_blank(),
   axis.title = element_text(size = 18),
   axis.text  = element_text(size = 16),
   legend.text  = element_text(size = 16),
   legend.title = element_text(size = 16),
   legend.position = c(0.98, 0.98),          # bottom-left (x,y in [0,1])
   legend.justification = c("right", "top"),
   legend.background = element_rect(fill = "white", color = "grey80"),
   legend.margin = margin(3, 6, 3, 6)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))
 
 safe_cell_type <- gsub("[^A-Za-z0-9]", "_", ct)
 filename <- paste0("Spatial_ST_", safe_cell_type, "_highlighted_dbscan.png")
 ggsave(file.path(wd.de.plots, filename), plot = dim_plot, width = 6, height = 6, dpi = 300)
 cat("Saved plot for", ct, "\n")
}

# -----------------------------------------------------------------------------
# Builds a plotting frame using the dbscan embeddings
# Last Modified: 18/09/25
# -----------------------------------------------------------------------------
# --- mark meta indices
meta_idx <- unique(meta.st$Index)
spatial_data$meta_flag    <- ifelse(rownames(spatial_data) %in% meta_idx, "yes", "no")
spatial_data$variant_flag <- ifelse(spatial_data$meta_flag == "yes", "6_31631400_A_G", NA)

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
	legend_breaks <- if (have_variant) c(ct, "6_31631400_A_G") else ct
	legend_values <- c(setNames(pal[ct], ct), `6_31631400_A_G` = "red")
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
	filename <- paste0("Spatial_ST_", safe_cell_type, "_highlighted_origUM_6_31631400_A_G.png")
	ggsave(file.path(wd.de.plots, filename), plot = dim_plot, width = 6, height = 6, dpi = 300)
	cat("Saved plot for", ct, "\n")
}
# > meta_cells
# d_1      d_2        cell_type seurat_clusters meta_flag   variant_flag
# AGTGCCGGTAAGGCTG-1_19 3203.926 1758.197 Meiotic division               9       yes 6_31631400_A_G
# ATCCATTGTGCGGTAA-1_19 4773.835 3406.395        Diplotene               1       yes 6_31631400_A_G
# ATTTCTGCAAAGGATT-1_19 4564.617 1730.516         Zygotene               6       yes 6_31631400_A_G
# CATTCATCATCAGTGT-1_19 4058.873 2022.929        Diplotene               9       yes 6_31631400_A_G
# CTACCCAAGACGCAGT-1_19       NA       NA        Pachytene               1       yes 6_31631400_A_G
# CTAGGTACACATTCTT-1_19 2002.293 3405.891 Meiotic division               9       yes 6_31631400_A_G
# GAATCGTCAATCTAGC-1_19 3967.537 1659.792  Early spermatid              11       yes 6_31631400_A_G
# TCATGGACATTGAAAG-1_19       NA       NA          Stage 3               3       yes 6_31631400_A_G
# TCATGTTGTGATACTC-1_19       NA       NA  Early spermatid              11       yes 6_31631400_A_G
# TTGAACGCAGGTTCAT-1_19 2670.183 3023.895          Stage 1               0       yes 6_31631400_A_G

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
spatial_data$variant_flag <- ifelse(spatial_data$meta_flag == "yes", "1_26779496_G_A", NA)

# For your loop:
spatial_coords <- c("plot_x","plot_y")

# Sanity: how many meta cells total and with usable coords?
cat("Total meta cells flagged:", sum(spatial_data$meta_flag == "yes"), "\n")
cat("Meta with coords:", sum(spatial_data$meta_flag == "yes" &
									  	complete.cases(spatial_data[, spatial_coords])), "\n")

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
gene <- "ARID1A"                  # target gene (case-insensitive handled)
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
outfile <- file.path(wd.de.plots, "Spatial_ST_META_ONLY_origcolors.png")
ggsave(outfile, plot = p_meta, width = 6, height = 6, dpi = 300)
cat("Saved:", outfile, "\n")

# -----------------------------------------------------------------------------
# To compute average Euclidean distances between cell types
# Last Modified: 22/09/25
# -----------------------------------------------------------------------------
library(dplyr)

# --- define order
cell_type_order <- c("Stage 0A","Stage 3","Leptotene","Zygotene","Early spermatid")

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

save(meta_idx, spatial_data, meta_cells_only, phase_proportions, file=file.path(wd.de.plots, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_23_res=0.5_-C0_ordered_annotated_meta_cells_only_ARID1A.RData")))

# ensure dir exists
if (!dir.exists(wd.de.plots)) dir.create(wd.de.plots, recursive = TRUE)

# move rownames (cell barcodes) into a proper column, then write
meta_cells_only_out <- tibble::rownames_to_column(as.data.frame(meta_cells_only), var = "cell_barcode")
outfile_tsv <- file.path(wd.de.plots, "meta_cells_only.tsv")
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
save_phase_plot <- function(df, phase_label, color_hex, outdir,
									 prefix = "Spatial_ST_META_ONLY_5") {
	sub <- df %>% dplyr::filter(.data$phase == phase_label, is.finite(.data$dist_to_next))
	if (nrow(sub) < 3 || all(is.na(sub$dist_to_next))) {
		warning(sprintf("Not enough data to plot/correlate for phase %s", phase_label))
		return(invisible(NULL))
	}
	
	ct  <- stats::cor.test(sub$dist_to_next, sub$proportion)
	rho <- base::unname(ct$estimate)
	p   <- ct$p.value
	lab <- sprintf("rho = %.3f\np = %.3g", rho, p)
	
	x_pos <- max(sub$dist_to_next, na.rm = TRUE)
	y_pos <- max(sub$proportion,   na.rm = TRUE)
	
	p_out <- ggplot2::ggplot(sub, ggplot2::aes(x = .data$dist_to_next, y = .data$proportion)) +
		ggplot2::geom_point(size = 2, alpha = 0.8, color = color_hex) +
		ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, color = color_hex) +
		ggplot2::annotate("text", x = x_pos, y = y_pos, label = lab,
								hjust = 1.05, vjust = 1.4, size = 4.2) +
		ggplot2::labs(
			title = paste0("Distance vs ", phase_label, " proportion"),
			x = "Average distance to next stage",
			y = paste0(phase_label, " proportion")
		) +
		ggplot2::theme_bw(base_size = 12) +
		ggplot2::theme(legend.position = "none",
							plot.title = ggplot2::element_text(face = "bold"))
	
	# Safe filename per phase
	fname <- file.path(outdir,
							 sprintf("%s_%s.png",
							 		  prefix, gsub("[^A-Za-z0-9]+", "_", phase_label)))
	
	ggplot2::ggsave(filename = fname, plot = p_out,
						 width = 5, height = 5, units = "in", dpi = 300)
	message("Saved: ", fname)
	invisible(fname)
}

# Save each phase plot with different file names in the same directory
save_phase_plot(phase_dist, "G1",   phase_colors[["G1"]],  wd.de.plots)
save_phase_plot(phase_dist, "S",    phase_colors[["S"]],   wd.de.plots)
save_phase_plot(phase_dist, "G2M",  phase_colors[["G2M"]], wd.de.plots)

# -----------------------------------------------------------------------------
# 
# Last Modified: 20/09/25
# -----------------------------------------------------------------------------




# Choose k adaptively: ~30 for thousands of cells; scale with log10(n)
n_cells <- nrow(df0)
k <- max(10L, min(60L, round(10 * log10(n_cells + 10))))

# (Optional) drop any rows with non-finite coords
bad <- !is.finite(df0$plot_x) | !is.finite(df0$plot_y)
if (any(bad)) {
	df0 <- df0[!bad, ]
	coords <- as.matrix(df0[, c("plot_x","plot_y")])
}

# kNN graph
kn <- FNN::get.knn(coords, k = k)

# Local phase proportions
get_local_prop <- function(phase_vec, kn_idx, target) {
	n <- length(phase_vec); out <- numeric(n)
	for (i in seq_len(n)) {
		neigh <- c(i, kn_idx$nn.index[i,])   # include self
		out[i] <- mean(phase_vec[neigh] == target, na.rm = TRUE)
	}
	out
}

df0$prop_G1_local  <- get_local_prop(df0$Phase, kn, "G1")
df0$prop_S_local   <- get_local_prop(df0$Phase, kn, "S")
df0$prop_G2M_local <- get_local_prop(df0$Phase, kn, "G2M")












# -----------------------------------------------------------------------------
# To quantify per-cell allelic imbalance (AI) at chr6:31631400
# Last Modified: 20/09/25
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
	library(dplyr)
	library(tidyr)
	library(ggplot2)
	library(readr)
})

# -----------------------------
# CONFIG
# -----------------------------
locus_chr <- "chr6"
locus_pos <- 31631400           # position of interest
min_total_reads <- 2            # minimum reads to keep a cell for AI stats

# Optional: if your barcode column is 'Index' rather than 'CB', set it here:
barcode_col <- if ("CB" %in% names(barcodes)) "CB" else if ("Index" %in% names(barcodes)) "Index" else "CB"

# -----------------------------
# 1) Filter to locus & sanitize inputs
# -----------------------------
dat <- barcodes %>%
	filter(.data$CHROM == !!locus_chr, .data$Start == !!locus_pos) %>%
	mutate(
		CB = .data[[barcode_col]],
		Base_observed = toupper(Base_observed)
	)

stopifnot(nrow(dat) > 0)

# -----------------------------
# 2) Collapse per cell: counts of A and G
# -----------------------------
per_cell <- dat %>%
	mutate(is_A = (Base_observed == "A"),
			 is_G = (Base_observed == "G")) %>%
	group_by(CB) %>%
	summarise(
		A_reads = sum(Num_reads[is_A], na.rm = TRUE),
		G_reads = sum(Num_reads[is_G], na.rm = TRUE),
		total_reads = A_reads + G_reads,
		.groups = "drop"
	) %>%
	filter(total_reads >= min_total_reads) %>%              # depth filter
	mutate(
		AF_alt = ifelse(total_reads > 0, G_reads / total_reads, NA_real_),
		log2_A_over_G = ifelse(A_reads > 0 & G_reads > 0, log2(A_reads / G_reads),
									  ifelse(G_reads == 0 & A_reads > 0,  Inf,   # all-A
									  		 ifelse(A_reads == 0 & G_reads > 0, -Inf, NA_real_))),
		winner = case_when(
			G_reads > A_reads ~ "ALT",
			A_reads > G_reads ~ "REF",
			TRUE              ~ "TIE"
		)
	)

# -----------------------------
# 3) Binomial test per cell vs p0 = 0.5 (allelic balance)
# -----------------------------
bt <- per_cell %>%
	rowwise() %>%
	mutate(
		p_binom = ifelse(total_reads > 0,
							  binom.test(G_reads, total_reads, p = 0.5, alternative = "two.sided")$p.value,
							  NA_real_)
	) %>%
	ungroup() %>%
	mutate(p_adj = p.adjust(p_binom, method = "BH"))

# -----------------------------
# 4) Join Seurat metadata (optional)
# -----------------------------
if (exists("so.st")) {
	meta_cols <- c("cell_type", "seurat_clusters", "nCount_RNA", "nFeature_RNA", "percent.mt")
	meta_cols <- base::intersect(meta_cols, colnames(so.st@meta.data))
	bt <- bt %>%
		mutate(cell = CB) %>%                   # your Seurat colnames look like "AAAA...-1_19"
		filter(cell %in% colnames(so.st)) %>%   # keep only cells present in object
		bind_cols(so.st@meta.data[.$cell, meta_cols, drop = FALSE])
}

# If you want to highlight your 10 twist-captured “meta” indices:
if (exists("meta.st") && "Index" %in% colnames(meta.st)) {
	bt$meta_flag <- ifelse(bt$CB %in% meta.st$Index, "yes", "no")
}

# -----------------------------
# 5) Save tidy results
# -----------------------------
out_tbl <- file.path(wd.de.plots, sprintf("AI_chr6_%d_per_cell.csv", locus_pos))
write_csv(bt, out_tbl)
cat("Saved per-cell AI table to:", out_tbl, "\n")

# -----------------------------
# 6) Plot 1: AF by winner group (REF vs ALT), REF left, ALT right
# -----------------------------
bt$winner_f <- factor(bt$winner, levels = c("REF","ALT","TIE"))

p1 <- ggplot(bt %>% filter(winner %in% c("REF","ALT")),
				 aes(x = winner_f, y = AF_alt, color = winner_f, size = total_reads)) +
	geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.25, color = "grey40") +
	geom_jitter(width = 0.12, alpha = 0.9) +
	scale_color_manual(values = c(REF = "#3182bd", ALT = "#cb181d", TIE = "grey50")) +
	scale_size_continuous(range = c(1.2, 4), breaks = c(1,2,3,4,5)) +
	geom_hline(yintercept = 0.5, linetype = 2, color = "grey50") +
	coord_cartesian(ylim = c(0,1)) +
	labs(x = NULL, y = sprintf("ALT allele fraction at %s:%d (G / (A+G))", locus_chr, locus_pos),
		  color = "Winner", size = "Depth") +
	theme_classic() +
	theme(legend.position = "right")

out_png1 <- file.path(wd.de.plots, sprintf("AI_chr6_%d_AF_by_winner.png", locus_pos))
ggsave(out_png1, plot = p1, width = 6, height = 4.5, dpi = 300)
cat("Saved:", out_png1, "\n")

# -----------------------------
# 7) Plot 2: Histogram of AF (all cells)
# -----------------------------
p2 <- ggplot(bt, aes(x = AF_alt, fill = winner_f)) +
	geom_histogram(bins = 20, alpha = 0.7, position = "identity", color = "white") +
	geom_vline(xintercept = 0.5, linetype = 2, color = "grey40") +
	scale_fill_manual(values = c(REF = "#6baed6", ALT = "#fc9272", TIE = "grey70"), name = "Winner") +
	labs(x = sprintf("ALT allele fraction at %s:%d", locus_chr, locus_pos), y = "Cells") +
	theme_classic()

out_png2 <- file.path(wd.de.plots, sprintf("AI_chr6_%d_AF_hist.png", locus_pos))
ggsave(out_png2, plot = p2, width = 6, height = 4.0, dpi = 300)
cat("Saved:", out_png2, "\n")

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
locus_chr <- "chr6"
locus_pos <- 31631400          # if your table has Start column
gene      <- "PRRC2A"

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
		REF_reads   = sum(Num_reads[Base_observed == "A"], na.rm = TRUE),
		ALT_reads   = sum(Num_reads[Base_observed == "G"], na.rm = TRUE),
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

# -----------------------------
# 4) Correlate AI with PRRC2A (Spearman) — conflict-proof
# -----------------------------

# manual NA mask (avoids complete.cases namespace issues)
keep_idx <- !is.na(alt_ai$AI) & !is.na(alt_ai$expr)

cor_res <- if (sum(keep_idx) >= 3) {
	suppressWarnings(stats::cor.test(alt_ai$AI[keep_idx],
												alt_ai$expr[keep_idx],
												method = "spearman"))
} else {
	NULL
}

rho  <- if (!is.null(cor_res)) base::unname(cor_res$estimate) else NA_real_
pval <- if (!is.null(cor_res)) cor_res$p.value else NA_real_

p_lab <- if (is.na(pval)) {
	"Too few points"
} else if (pval < 1e-4) {
	paste0("Spearman \u03C1 = ", signif(rho, 3), "\nP < 1e-4")
} else {
	paste0("Spearman \u03C1 = ", signif(rho, 3), "\nP = ", signif(pval, 3))
}

# (optional) Weighted linear fit by depth, if you have total_reads:
# fit <- stats::lm(expr ~ AI, data = alt_ai[keep_idx, ], weights = total_reads[keep_idx])
# summary(fit)

# -----------------------------
# 5) Plot scatter + save
# -----------------------------
p <- ggplot(alt_ai, aes(x = AI, y = expr)) +
	geom_point(size = 3, alpha = 0.85, color = "red") +
	geom_smooth(method = "lm", se = FALSE, linetype = 2, color = "black") +
	labs(
		x = "Allelic imbalance (|ALT/(ALT+REF) − 0.5|)",
		y = paste0(gene, " expression (normalized)"),
		title = "ALT cells: AI vs PRRC2A"
	) +
	annotate("text",
				x = Inf, y = Inf, hjust = 1.05, vjust = 1.5,
				label = if (is.na(pval)) "Too few points" else paste0("Spearman ρ = ", round(rho, 2),
																						"\nP = ", signif(pval, 3)),
				size = 5) +
	theme_classic(base_size = 14)

print(p)
ggsave(file.path(wd.de.plots, "ALT_cells_AI_vs_PRRC2A.png"), p, width = 6, height = 5, dpi = 300)

# -----------------------------
# 6) Save table
# -----------------------------
out_csv <- file.path(wd.de.plots, "ALT_cells_AI_PRRC2A.csv")
readr::write_csv(alt_ai, out_csv)
cat("Saved:", out_csv, "\n")




# -----------------------------
# AF
# -----------------------------
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

# Spearman correlations
spearman_AF      <- suppressWarnings(stats::cor.test(alt_ai$AF,       alt_ai$expr, method = "spearman"))
spearman_logitAF <- suppressWarnings(stats::cor.test(alt_ai$logit_AF, alt_ai$expr, method = "spearman"))
spearman_log2AR  <- suppressWarnings(stats::cor.test(alt_ai$log2_AR,  alt_ai$expr, method = "spearman"))

summ_cor <- data.frame(
	metric = c("AF", "logitAF", "log2AR"),
	rho    = c(as.numeric(spearman_AF$estimate),
				  as.numeric(spearman_logitAF$estimate),
				  as.numeric(spearman_log2AR$estimate)),
	pval   = c(spearman_AF$p.value,
				  spearman_logitAF$p.value,
				  spearman_log2AR$p.value)
)

# pretty print
print(transform(summ_cor,
					 rho  = signif(rho,  3),
					 pval = signif(pval, 3)
))

# save if you like
# readr::write_csv(summ_cor, file.path(wd.de.plots, "ALT_AI_PRRC2A_spearman_summary.csv"))

# Weighted linear fits (optional), weighting by depth
fit_logit <- stats::lm(expr ~ logit_AF, data = alt_ai, weights = total_reads)
fit_log2  <- stats::lm(expr ~ log2_AR,  data = alt_ai, weights = total_reads)

summary(fit_logit)  # slope & p-value on logit(AF)
summary(fit_log2)   # slope & p-value on log2(ALT/REF)

# Plot AF (human-friendly) with Spearman label based on logit(AF)
lab_rho <- paste0("Spearman \u03C1 (logitAF) = ",
						signif(unname(spearman_logitAF$estimate), 3),
						"\nP = ", signif(spearman_logitAF$p.value, 3))

p <- ggplot(alt_ai, aes(x = AF, y = expr)) +
	geom_point(size = 3, alpha = 0.85, color = "red") +
	geom_smooth(method = "lm", se = FALSE, linetype = 2, color = "black") +
	labs(x = "ALT allele fraction (AF = ALT / (ALT+REF))",
		  y = "PRRC2A expression (normalized)",
		  title = "ALT cells: AI vs PRRC2A") +
	annotate("text", x = Inf, y = Inf, label = lab_rho,
				hjust = 1.05, vjust = 1.5, size = 5) +
	theme_classic(base_size = 14)
print(p)

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
	labs(x = "ALT allele fraction (AF)", y = "PRRC2A expression") +
	annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5, label = lab_AF, size = 5) +
	theme_classic(base_size = 14)

p2 <- ggplot(alt_ai, aes(x = total_reads, y = expr)) +
	geom_point(size = 3, alpha = 0.85, color = "blue") +
	geom_smooth(method = "lm", se = FALSE, linetype = 2, color = "black") +
	labs(x = "Total reads (REF + ALT)", y = "PRRC2A expression") +
	annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5, label = lab_D, size = 5) +
	theme_classic(base_size = 14)

print(p1); print(p2)

# (optional) save
# ggsave(file.path(wd.de.plots, "ALT_AF_vs_expr.png"),  p1, width = 6, height = 5, dpi = 300)
# ggsave(file.path(wd.de.plots, "ALT_depth_vs_expr.png"), p2, width = 6, height = 5, dpi = 300)













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
gene <- "PRRC2A"

# Your variant calls (exact cell IDs; no suffix changes here)
barcodes.ref <- subset(subset(barcodes, CHROM == "chr6"), Base_observed == "A")
barcodes.alt <- subset(subset(barcodes, CHROM == "chr6"), Base_observed == "G")

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
