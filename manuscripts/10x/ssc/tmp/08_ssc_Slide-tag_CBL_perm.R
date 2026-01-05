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
dest <- file.path(wd.de.plots, "origUM")
if (!dir.exists(dest)) {
	ok <- base::dir.create(dest, recursive = TRUE, showWarnings = FALSE, mode = "0775")
	if (!ok && !dir.exists(dest)) stop("Failed to create: ", dest)
}

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
									 prefix = "Spatial_ST_META_CBL_permutation_trunc_dec_RNA") {
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
















# -----------------------------------------------------------------------------
# WT-matched resampling (N times)
# Last Modified: 09/10/25
# -----------------------------------------------------------------------------

# --- 0) What we assume you already have -------------------------------
# spatial_data: all cells with columns at least: cell, cell_type, plot_x, plot_y
# meta_cells_only: the 53 mutant cells with the same columns
# ct_order_master: the biological stage order

# --- 1) Stages covered by mutants (respect your biological order) -----
stages_mut <- ct_order_master[ct_order_master %in% unique(as.character(meta_cells_only$cell_type))]

# --- 2) Count mutants per stage (targets to match) --------------------
mut_counts <- meta_cells_only %>%
	mutate(cell_type = as.character(cell_type)) %>%
	dplyr::count(cell_type, name = "n_target") %>%   # <-- explicitly call dplyr::count
	right_join(tibble(cell_type = stages_mut), by = "cell_type") %>%
	mutate(n_target = ifelse(is.na(n_target), 0L, n_target))

# --- 3) Build a WT pool: all cells minus mutants, only these stages ---
# Prefer matching by 'cell' id; fall back to coords+stage if needed
if ("cell" %in% names(spatial_data) && "cell" %in% names(meta_cells_only)) {
	wt_pool <- dplyr::anti_join(
		spatial_data,
		meta_cells_only %>% dplyr::select(cell),
		by = "cell"
	)
} else {
	wt_pool <- dplyr::anti_join(
		spatial_data,
		meta_cells_only %>% dplyr::select(plot_x, plot_y, cell_type) %>% dplyr::distinct(),
		by = c("plot_x","plot_y","cell_type")
	)
}

wt_pool <- wt_pool %>%
	mutate(cell_type = as.character(cell_type)) %>%
	dplyr::filter(cell_type %in% stages_mut) %>%
	dplyr::filter(is.finite(plot_x), is.finite(plot_y))

# --- 4) Attach targets and sample per stage ---------------------------
wt_pool_targeted <- wt_pool %>%
	left_join(mut_counts, by = "cell_type")

# -----------------------------------------------
# Resampling loop
# -----------------------------------------------
B <- 1000L  # number of permutations
res_list_resampling <- vector("list", B)  # each element will be a tibble with 3 rows (phases)

for (b in seq_len(B)) {
	# 5) Sample per stage using group_modify (per-group scalar n)
	meta_cells_only_resampling <- wt_pool_targeted %>%
		dplyr::group_by(cell_type) %>%
		dplyr::group_modify(function(.x, .key) {
			n_take <- min(.x$n_target[1], nrow(.x))
			dplyr::slice_sample(.x, n = n_take, replace = FALSE)
		}) %>%
		dplyr::ungroup() %>%
		dplyr::select(dplyr::any_of(c("cell","plot_x","plot_y","cell_type","seurat_clusters"))) %>%
		dplyr::mutate(cell_type = factor(cell_type, levels = stages_mut))
	
	# --- define order
	cell_types <- c("Stage 0","Stage 0A","Stage 0B","Stage 1","Stage 2","Stage 3","Leptotene","Zygotene","Pachytene","Diplotene","Meiotic division","Early spermatid","Late spermatid")
	cell_type_order <- base::intersect(cell_types, as.vector(unique(meta_cells_only_resampling$cell_type)))
	
	meta_cells_only_resampling <- meta_cells_only_resampling %>%
		dplyr::mutate(cell_type = factor(cell_type, levels = cell_type_order)) %>%
		dplyr::arrange(cell_type)
	
	dist_map$dist_to_next <- mapply(function(f, t) avg_dist(meta_cells_only_resampling, f, t),
											  dist_map$from, dist_map$to)
	
	# --- join back (no rename needed)
	meta_cells_only_resampling <- meta_cells_only_resampling |>
		dplyr::left_join(dplyr::select(dist_map, from, dist_to_next),
							  by = c("cell_type" = "from"))
	
	# --- aggregate distances to per-cell_type (unique value for each "from")
	dist_summary_resampling <- meta_cells_only_resampling %>%
		dplyr::select(cell_type, dist_to_next) %>%
		dplyr::distinct()
	
	# --- join with phase proportions
	phase_dist_resampling <- phase_proportions %>%
		dplyr::left_join(dist_summary_resampling, by = "cell_type")
	
	cor_results <- phase_dist_resampling %>%
		group_by(phase) %>%
		summarise(
			rho = cor(dist_to_next, proportion, use = "complete.obs"),
			pval = cor.test(dist_to_next, proportion)$p.value,
			.groups = "drop"
		)
	
	res_list_resampling[[b]] <- cor_results
}

save(res_list_resampling, file=file.path(wd.de.plots, "origUM", paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_23_res=0.5_-C0_ordered_annotated_meta_cells_only_", GENE, "_perm_resampling.RData")))

# =============================================================================
# Compare WT-matched resampling with observed mutant result (S phase)
# =============================================================================

# 1) Combine the 1,000 resamples into one long tibble with permutation index b
perm_cor_long_resampling <- dplyr::bind_rows(res_list_resampling, .id = "b") %>%
	dplyr::mutate(b = as.integer(b),
					  phase = as.character(phase))

save_hist_plot <- function(df, phase_text, phase_label, color_hex, outdir)  {
	# 2) Extract the WT null for S-phase (ρ per resample)
	rhos_wt_S <- df %>%
		dplyr::filter(phase == phase_text) %>%
		dplyr::pull(rho)
	
	# HERE!
	rho_mut_S <- subset(cor_results_obs, phase == phase_text)$rho
	
	# 4) Empirical p-value and 95% CI of the WT null
	N <- length(rhos_wt_S)
	emp_p_wt_vs_mut <- (sum(rhos_wt_S >= rho_mut_S, na.rm = TRUE) + 1) / (N + 1)
	ci_wt <- quantile(rhos_wt_S, c(0.025, 0.975), na.rm = TRUE)
	
	summary_tbl <- tibble::tibble(
		metric        = "WT-matched_Sphase_pairwise",
		rho_mut       = rho_mut_S,
		emp_p         = emp_p_wt_vs_mut,
		wt_null_ci_lo = ci_wt[1],
		wt_null_ci_hi = ci_wt[2],
		N_resamples   = N
	)
	print(summary_tbl)
	readr::write_tsv(summary_tbl, file.path(wd.de.plots, "origUM", paste0("WTmatched_", phase_text, "rho_summary.tsv")))
	
	# 5) Save the WT null as TSV (optional)
	readr::write_tsv(tibble::tibble(rho = rhos_wt_S),
						  file.path(wd.de.plots, "origUM", paste0("WTmatched_null_rho_", phase_text, ".tsv")))
	
	# 6) Plot histogram of WT null with a vertical line at mutant ρ and empirical p in subtitle
	out_png <- file.path(wd.de.plots, "origUM", paste0("WTmatched_null_hist_", phase_text, "rho_pairwise.png"))
	p_hist <- ggplot2::ggplot(data.frame(rho = rhos_wt_S), aes(x = rho)) +
		geom_histogram(bins = 40) +
		geom_vline(xintercept = rho_mut_S, linewidth = 1) +
		labs(
			title = "WT-matched resamples",
			subtitle = sprintf("empirical p = %.4f", emp_p_wt_vs_mut),
			x = "Spearman ρ (S-phase)",
			y = "Count"
		) +
		theme_bw(base_size = 12) +
		ggplot2::theme(
			legend.position = "none",
			plot.title = ggplot2::element_text(face = "bold", size = 20, hjust = 0.5),
			axis.title = ggplot2::element_text(size = 18),
			axis.text  = ggplot2::element_text(size = 16)
		)
	ggsave(out_png, p_hist, width = 6, height = 4.5, dpi = 300)
	message("Saved: ", out_png)	
}

# Save each phase plot with different file names in the same directory
save_hist_plot(perm_cor_long_resampling, "G1",  "G0/G1", phase_colors[["G1"]])
save_hist_plot(perm_cor_long_resampling, "S",   "S",     phase_colors[["S"]])
save_hist_plot(perm_cor_long_resampling, "G2M", "G2/M",  phase_colors[["G2M"]])

