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
wd.de.plots <- file.path(wd.de.plots, "origUM")
# ensure dir exists
if (!dir.exists(file.path(wd.de.plots)))                           dir.create(file.path(wd.de.plots), recursive = TRUE)
if (!dir.exists(file.path(wd.de.plots, "avg_dist", "resampling"))) dir.create(file.path(wd.de.plots, "avg_dist", "resampling"), recursive = TRUE)

# -----------------------------------------------------------------------------
# After running SComatic genotyping
# # Last Modified: 17/09/25
# -----------------------------------------------------------------------------
#library(Seurat)
#
#load(file.path(wd.de.data, paste0("ssc_filtered_normalised_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_ordered_annotated.RData")))
#so.integrated$cell_type <- Idents(so.integrated)

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

# -----------------------------------------------------------------------------
# Clonal expansion across 12 stages — CALLABLE-ONLY tallies
# Prevalence p_i with misclassification (s_i, fpr_i) – inputs k_i, n_i
# Last Modified: 13/10/25
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
	library(dplyr); library(tidyr); library(stringr); library(tibble)
})

# ----- parameters you can tweak ----------------------------------------------
MIN_TOTAL_UMI <- 1L   # callable threshold: total (REF+ALT) UMIs at the locus
MIN_ALT_UMI   <- 1L   # strict ALT-positive threshold (set 1L if you want lenient)
REQUIRE_BOTH_STRANDS <- FALSE  # set TRUE if you have per-cell fwd/rev info
# Column names used below; adjust if your SComatic export differs:
COL_CB   <- "CB"
COL_READS<- "Num_reads"   # per-row read/UMI count
COL_ALT  <- "ALT_UMI"     # per-row ALT UMI count (if not present, we infer)
COL_REF  <- "REF_UMI"     # per-row REF UMI count (if not present, we infer)
# ------------------------------------------------------------------------------

stage_order <- c(
	"Stage 0","Stage 0A","Stage 0B","Stage 1","Stage 2","Stage 3",
	"Leptotene","Zygotene","Pachytene","Diplotene",
	"Meiotic division","Early spermatid","Late spermatid"
)

# --- 0) Build per-cell coverage at the locus (union of REF & ALT) -------------
# (Do NOT restrict to germ here; we’ll keep only your 12 germ stages after joining labels)

# helper to safely pull a column if present
getcol <- function(df, nm) if (nm %in% names(df)) df[[nm]] else NULL

ref_calls <- barcodes %>%
	dplyr::filter(.data$CHROM == CHR, .data$Start == POS, .data$Base_observed == REF) %>%
	dplyr::select(all_of(COL_CB), any_of(c(COL_REF, COL_READS))) %>%
	dplyr::group_by(.data[[COL_CB]]) %>%
	dplyr::summarise(
		ref_umis = if (!is.null(getcol(cur_data_all(), COL_REF)))
			sum(.data[[COL_REF]], na.rm = TRUE) else
				sum(.data[[COL_READS]], na.rm = TRUE),
		.groups = "drop"
	)

alt_calls <- barcodes %>%
	dplyr::filter(.data$CHROM == CHR, .data$Start == POS, .data$Base_observed == ALT) %>%
	dplyr::select(all_of(COL_CB), any_of(c(COL_ALT, COL_READS))) %>%
	dplyr::group_by(.data[[COL_CB]]) %>%
	dplyr::summarise(
		alt_umis = if (!is.null(getcol(cur_data_all(), COL_ALT)))
			sum(.data[[COL_ALT]], na.rm = TRUE) else
				sum(.data[[COL_READS]], na.rm = TRUE),
		.groups = "drop"
	)

# These three mutant cells contain REF reads
both <- intersect(alt_calls$CB, ref_calls$CB)
# Remove these three cells from ref_calls
ref_calls <- ref_calls %>%
	filter(!CB %in% both)

# For using IGV manually
readr::write_lines(ref_calls$CB,
					  file.path(wd.de.plots, "avg_dist", "resampling", "wildtype_cell_barcodes.txt"))
readr::write_lines(alt_calls$CB,
					  file.path(wd.de.plots, "avg_dist", "resampling", "mutant_cell_barcodes.txt"))

# merge and compute totals
cov_tbl <- dplyr::full_join(ref_calls, alt_calls, by = COL_CB) %>%
	dplyr::mutate(
		ref_umis = dplyr::coalesce(.data$ref_umis, 0L),
		alt_umis = dplyr::coalesce(.data$alt_umis, 0L),
		total_umis = ref_umis + alt_umis
	)

# (optional) if you have per-cell strand/start fields in barcodes, join & use them:
# e.g., unique starts, forward/reverse counts -> set REQUIRE_BOTH_STRANDS = TRUE
if (REQUIRE_BOTH_STRANDS && all(c("BCf","BCr") %in% strsplit(getOption("SComatic_fields",""), ",")[[1]])) {
	warning("REQUIRE_BOTH_STRANDS requested, but per-cell strand info not wired in; set to FALSE or join those columns here.")
}

# --- 1) Stage labels from Seurat (PD53626b_ST1) with meta as fallback ----------
seu_stage <- Seurat::FetchData(so.st, vars = "cell_type") %>%
	tibble::rownames_to_column(COL_CB) %>%
	dplyr::rename(stage_seu = cell_type)

meta_stage <- meta %>%
	dplyr::select(Index, Cell_type) %>%
	dplyr::rename(!!COL_CB := Index, stage_meta = Cell_type)

# join labels; keep only your 12-stage scheme
calls2 <- cov_tbl %>%
	dplyr::left_join(seu_stage, by = COL_CB) %>%
	dplyr::left_join(meta_stage, by = COL_CB) %>%
	dplyr::mutate(
		stage_raw = dplyr::if_else(!is.na(stage_seu), stage_seu, stage_meta),
		stage_raw = stringr::str_squish(as.character(stage_raw)),
		stage     = dplyr::if_else(stage_raw %in% stage_order, stage_raw, NA_character_)
	) %>%
	dplyr::filter(!is.na(stage)) %>%
	dplyr::mutate(stage = factor(stage, levels = stage_order))

# --- 2) Define CALLABLE & ALT-positive per cell -------------------------------
calls2 <- calls2 %>%
	dplyr::mutate(
		callable   = total_umis >= MIN_TOTAL_UMI,
		alt_pos    = alt_umis >= MIN_ALT_UMI    # tighten/loosen as needed
		# if you add strand/start filters, include them here (e.g., alt_pos & both_strands & ≥2 starts)
	)

# ensure both levels exist, then widen with a prefix
callable_stage_counts <- calls2 %>%
	dplyr::count(stage, callable, name = "n") %>%
	tidyr::complete(stage, callable = c(FALSE, TRUE), fill = list(n = 0)) %>%
	tidyr::pivot_wider(names_from = callable, values_from = n, values_fill = 0) %>%
	dplyr::rename(`non_callable` = `FALSE`, `callable` = `TRUE`) %>%
	dplyr::arrange(stage)

print(callable_stage_counts)

# --- 3) Tally k (ALT+) and n (CALLABLE) per stage -----------------------------
k_n_by_stage <- calls2 %>%
	dplyr::filter(callable) %>%
	dplyr::group_by(stage) %>%
	dplyr::summarise(
		k = sum(alt_pos, na.rm = TRUE),
		n = dplyr::n(),            # number of callable cells at the locus
		.groups = "drop"
	) %>%
	dplyr::arrange(stage) %>%
	tidyr::complete(stage = factor(stage_order, levels = stage_order),
						 fill = list(k = 0L, n = 0L))

print(k_n_by_stage)

# 4) Write to disk
readr::write_tsv(k_n_by_stage,
					  file.path(wd.de.plots, "avg_dist", "resampling", "k_n_by_stage.tsv"))

# --- 4) Export vectors for prevalence model -----------------------------------
k_vec <- as.integer(k_n_by_stage$k)
n_vec <- as.integer(k_n_by_stage$n)
names(k_vec) <- names(n_vec) <- as.character(k_n_by_stage$stage)

cat("\n# k (ALT-positive, callable-only) per stage:\n"); print(k_vec)
cat("\n# n (callable cells) per stage:\n");           print(n_vec)

# =============================================================================
# To generate BAM files
# Last Modified: 12/11/25
# =============================================================================
# --- Split into mutant and reference calls -------------------------------
calls2.alt <- calls2 %>%
	filter(alt_pos | alt_umis > 0)

calls2.ref <- calls2 %>%
	filter(!alt_pos & alt_umis == 0)

# Move CAATGACGTGTTACTG-1_19 from calls2.ref to calls2.alt
cell_to_move <- "CAATGACGTGTTACTG-1_19"

# 1. Extract the row from ref and modify alt_pos
cell_row <- calls2.ref %>%
	filter(CB == cell_to_move) %>%
	mutate(alt_pos = TRUE)  # ensure at least 1 ALT UMI

# 2. Remove it from ref
calls2.ref <- calls2.ref %>%
	filter(CB != cell_to_move)

# 3. Add it to alt
calls2.alt <- bind_rows(calls2.alt, cell_row)

# For using IGV manually
readr::write_lines(calls2.ref$CB,
						 file.path(wd.de.plots, "BAM", "wildtype_cell_barcodes_99.txt"))
readr::write_lines(calls2.alt$CB,
						 file.path(wd.de.plots, "BAM", "mutant_cell_barcodes_54.txt"))

# =============================================================================
# WT-matched resampling (callable-only at 11:119284966)
# Last Modified: 20/10/25
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr)})

load(file=file.path(wd.de.plots, paste0("ssc_filtered_normalised_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_ordered_annotated_meta_cells_only_", GENE, "_RNA_so.st.RData")))

# --- Helper: ensure we can join by CB -----------------------------------------
spatial_data_cb <- spatial_data
if (!"CB" %in% names(spatial_data_cb)) {
	# assume 'cell' is your barcode; adjust if you keep suffixes elsewhere
	spatial_data_cb <- spatial_data_cb %>% dplyr::mutate(CB = .data$cell)
}
meta_cells_only_cb <- meta_cells_only
if (!"CB" %in% names(meta_cells_only_cb)) {
	meta_cells_only_cb <- meta_cells_only_cb %>% dplyr::mutate(CB = .data$cell)
}

# --- 0) Stages covered by mutants (respect biological order) ------------------
stages_mut <- ct_order_master[ct_order_master %in% unique(as.character(meta_cells_only_cb$cell_type))]

# --- 1) Mutant per-stage counts to match --------------------------------------
mut_counts <- meta_cells_only_cb %>%
	dplyr::mutate(cell_type = as.character(cell_type)) %>%
	dplyr::count(cell_type, name = "n_target") %>%
	dplyr::right_join(tibble::tibble(cell_type = stages_mut), by = "cell_type") %>%
	dplyr::mutate(n_target = ifelse(is.na(n_target), 0L, n_target))

# --- 2) Callable-only WT pool at the locus ------------------------------------
# calls2 came from your "callable-only tallies" step:
#   columns: CB, ref_umis, alt_umis, total_umis, stage (12 germ stages), callable (TRUE/FALSE)
# Keep only callable cells, exclude mutants, keep only stages you use
mut_CBs <- unique(meta_cells_only_cb$CB)

wt_pool_callable <- calls2 %>%
	dplyr::filter(callable) %>%
	dplyr::select(CB, stage) %>%
	dplyr::rename(cell_type = stage) %>%
	dplyr::mutate(cell_type = as.character(cell_type)) %>%
	dplyr::filter(!(CB %in% mut_CBs), cell_type %in% stages_mut) %>%
	dplyr::inner_join(
		spatial_data_cb %>% dplyr::select(CB, cell, plot_x, plot_y, cell_type_spatial = cell_type, seurat_clusters),
		by = "CB"
	) %>%
	# trust the stage label from calls2 (callable table) for matching
	dplyr::select(CB, cell, plot_x, plot_y, cell_type, seurat_clusters)

# sanity: check availability vs targets
avail_vs_target <- wt_pool_callable %>%
	dplyr::count(cell_type, name = "available") %>%
	dplyr::full_join(mut_counts, by = "cell_type") %>%
	dplyr::mutate(
		available = ifelse(is.na(available), 0L, available),
		n_target  = ifelse(is.na(n_target), 0L, n_target),
		shortfall = pmax(n_target - available, 0L)
	) %>% dplyr::arrange(match(cell_type, stages_mut))
print(avail_vs_target, n = 20)
if (any(avail_vs_target$shortfall > 0)) {
	message("Warning: some stages have fewer callable WT cells than mutant targets; ",
			  "we will sample as many as available for those stages.")
}

# --- 3) Attach targets and prepare sampling frame -----------------------------
wt_pool_targeted <- wt_pool_callable %>%
	dplyr::left_join(mut_counts, by = "cell_type")

# --- 4) Resampling loop (callable-only WT) ------------------------------------
B <- 1000L
res_list_resampling <- vector("list", B)

for (b in seq_len(B)) {
	# per-stage matched sample (min(target, available))
	meta_cells_only_resampling <- wt_pool_targeted %>%
		dplyr::group_by(cell_type) %>%
		dplyr::group_modify(function(.x, .key) {
			n_take <- min(.x$n_target[1], nrow(.x))
			dplyr::slice_sample(.x, n = n_take, replace = FALSE)
		}) %>%
		dplyr::ungroup() %>%
		dplyr::mutate(cell_type = factor(cell_type, levels = stages_mut)) %>%
		dplyr::arrange(cell_type)
	
	# recompute adjacent-stage distances on this resample
	# (dist_map has 'from' and 'to' stages already defined in your environment)
	dist_map$dist_to_next <- mapply(function(f, t) avg_dist(meta_cells_only_resampling, f, t),
											  dist_map$from, dist_map$to)
	
	# per-stage unique distance rows
	dist_summary_resampling <- dist_map %>%
		dplyr::select(from, dist_to_next) %>%
		dplyr::rename(cell_type = from) %>%
		dplyr::distinct()
	
	# join with phase proportions and compute phase-wise correlations
	phase_dist_resampling <- phase_proportions %>%
		dplyr::left_join(dist_summary_resampling, by = c("cell_type"))
	
	cor_results <- phase_dist_resampling %>%
		dplyr::group_by(phase) %>%
		dplyr::summarise(
			rho = unname(cor.test(dist_to_next, proportion, method = "spearman", exact = FALSE)$estimate),
			pval = cor.test(dist_to_next, proportion, method = "spearman", exact = FALSE)$p.value,
			.groups = "drop"
		)
	
	res_list_resampling[[b]] <- cor_results
}

save(res_list_resampling,
	  file = file.path(wd.de.plots,
	  					  paste0("ssc_filtered_normalised_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_ordered_annotated_meta_cells_only_",
	  					  		 GENE, "_perm_resampling_callable_only_spearman.RData")))

# -----------------------------------------------
# Combine callable-only WT resamples across phases
# and compare to observed (b = 0)
# -----------------------------------------------
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr)})

# 1) Bind resamples with permutation index b
perm_cor_long <- dplyr::bind_rows(
	Map(function(tbl, i) dplyr::mutate(tbl, b = as.integer(i)),
		 res_list_resampling, seq_along(res_list_resampling))
) %>%
	dplyr::mutate(phase = as.character(phase)) %>%
	dplyr::select(b, phase, dplyr::everything())

# 2) Prepend the observed (b = 0)
stopifnot(all(c("phase","rho","pval") %in% names(cor_results_obs)))
cor_results_obs <- cor_results_obs %>%
	dplyr::mutate(b = 0L, phase = as.character(phase)) %>%
	dplyr::select(b, phase, rho, pval)

perm_cor_long_with_obs <- dplyr::bind_rows(cor_results_obs, perm_cor_long)

# 3) Wide convenience tables
perm_cor_rho_wide <- tidyr::pivot_wider(
	perm_cor_long, id_cols = "b", names_from = "phase", values_from = "rho"
)
perm_cor_p_wide <- tidyr::pivot_wider(
	perm_cor_long, id_cols = "b", names_from = "phase", values_from = "pval"
)

# 4) Write to disk
readr::write_tsv(perm_cor_long,
					  file.path(wd.de.plots, "avg_dist", "resampling", "perm_cor_by_phase_long_callable_spearman.tsv"))
readr::write_tsv(perm_cor_rho_wide,
					  file.path(wd.de.plots, "avg_dist", "resampling", "perm_cor_rho_wide_callable_spearman.tsv"))
readr::write_tsv(perm_cor_p_wide,
					  file.path(wd.de.plots, "avg_dist", "resampling", "perm_cor_p_wide_callable_spearman.tsv"))

# -----------------------------------------------
# Directional empirical p-values per phase
# -----------------------------------------------
phase_direction <- c(G1 = "less", S = "greater", G2M = "two.sided")

resampling_pvals_one_phase_dir <- function(phase_name,
												direction = c("less","greater","two.sided")) {
	direction <- match.arg(direction)
	
	# observed
	rho_obs <- perm_cor_long_with_obs %>%
		dplyr::filter(b == 0, phase == phase_name) %>% dplyr::pull(rho)
	p_obs   <- perm_cor_long_with_obs %>%
		dplyr::filter(b == 0, phase == phase_name) %>% dplyr::pull(pval)
	
	# permuted (WT-matched callable-only)
	rho_perm <- perm_cor_long_with_obs %>%
		dplyr::filter(b > 0, phase == phase_name) %>% dplyr::pull(rho)
	p_perm   <- perm_cor_long_with_obs %>%
		dplyr::filter(b > 0, phase == phase_name) %>% dplyr::pull(pval)
	
	# directional empirical p on rho
	p_perm_r <- switch(direction,
							 greater   = (sum(rho_perm >=  rho_obs, na.rm = TRUE) + 1) / (length(rho_perm) + 1),
							 less      = (sum(rho_perm <=  rho_obs, na.rm = TRUE) + 1) / (length(rho_perm) + 1),
							 two.sided = (sum(abs(rho_perm) >= abs(rho_obs), na.rm = TRUE) + 1) / (length(rho_perm) + 1)
	)
	
	# "p-of-p" version (how often permuted p ≤ observed p)
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
resampling_dir <- dplyr::bind_rows(lapply(phases, function(ph) {
	resampling_pvals_one_phase_dir(ph, direction = phase_direction[[ph]])
})) %>%
	dplyr::mutate(
		q_perm_r = p.adjust(p_perm_r, method = "BH"),
		q_perm_p = p.adjust(p_perm_p, method = "BH")
	)

print(resampling_dir)
readr::write_tsv(resampling_dir,
					  file.path(wd.de.plots, "avg_dist", "resampling", "phase_empirical_pvalues_callable_spearman.tsv"))

# -----------------------------------------------
# Optional: extract WT null vectors per phase for plotting
# -----------------------------------------------
#rhos_wt_G1 <- perm_cor_long %>% dplyr::filter(phase == "G1")  %>% dplyr::pull(rho)
#rhos_wt_S  <- perm_cor_long %>% dplyr::filter(phase == "S")   %>% dplyr::pull(rho)
#rhos_wt_G2M<- perm_cor_long %>% dplyr::filter(phase == "G2M") %>% dplyr::pull(rho)

# -----------------------------------------------------------------------------
# To test correlation between distances (dist_to_next) and cell cycle phase proportions (G1/S/G2M)
# Last Modified: 09/10/25
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.plots, paste0("ssc_filtered_normalised_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_ordered_annotated_meta_cells_only_", GENE, "_RNA_perm_spearman.RData")))

# truncate to k decimals (works for +/- numbers)
trunc_dec <- function(x, k) trunc(x * 10^k) / 10^k

# Helper to make & save one phase plot to a specific directory with a clear filename
save_phase_plot_perm_resampling <- function(df, phase_text, phase_label, color_hex, outdir, emp_dir, resampling_dir,
											prefix = "Spatial_ST_META_CBL_permutation_trunc_dec_RNA_callable_spearman") {
	sub <- df %>% dplyr::filter(.data$phase == phase_text, is.finite(.data$dist_to_next))
	if (nrow(sub) < 3 || all(is.na(sub$dist_to_next))) {
		warning(sprintf("Not enough data to plot/correlate for phase %s", phase_text))
		return(invisible(NULL))
	}
	
	ct  <- cor.test(sub$dist_to_next, sub$proportion, method = "spearman", exact = FALSE)
	rho <- base::unname(ct$estimate)
	p   <- ct$p.value
	lab <- paste0("rho = ", round(rho, 2), "; p = ", trunc_dec(p, 4))
	lab <- paste0(lab, "\nempirical p = ", trunc_dec(subset(emp_dir, phase == phase_text)$p_perm_r, 4))
	lab <- paste0(lab, "\nWT-matched p = ", trunc_dec(subset(resampling_dir, phase == phase_text)$p_perm_r, 4))
	
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
save_phase_plot_perm_resampling(phase_dist, "G1",  "G0/G1", phase_colors[["G1"]],  file.path(wd.de.plots, "avg_dist", "resampling"), emp_dir, resampling_dir)
save_phase_plot_perm_resampling(phase_dist, "S",   "S",     phase_colors[["S"]],   file.path(wd.de.plots, "avg_dist", "resampling"), emp_dir, resampling_dir)
save_phase_plot_perm_resampling(phase_dist, "G2M", "G2/M",  phase_colors[["G2M"]], file.path(wd.de.plots, "avg_dist", "resampling"), emp_dir, resampling_dir)
