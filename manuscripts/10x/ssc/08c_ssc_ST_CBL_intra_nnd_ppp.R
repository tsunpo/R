# =============================================================================
# Name         :
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 24/10/25
# Purpose      : Callable-only within-stage spatial tests for CBL mutant,
#                using per-cell S/G2M from so.integrated (S.Score, G2M.Score).
# =============================================================================

# -------------------------- Config -------------------------------------------
CHR  <- "chr11"
POS  <- 119284966
REF  <- "T"
ALT  <- "A"
GENE <- "CBL"
MUT  <- paste0(gsub("chr", "", CHR), "_", POS, "_", REF, "_", ALT)

wd <- "/lustre/scratch125/casm/staging/team294/ty2"   ## farm22
BASE <- "SSC"
base <- tolower(BASE)
wd.rna <- file.path(wd, BASE)
wd.de       <- file.path(wd.rna, "analysis", paste0(base, ""))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots", GENE, MUT)
wd.de.plots <- file.path(wd.de.plots, "origUM")
# ensure dir exists
dir.create(file.path(wd.de.plots), recursive = TRUE)
dir.create(file.path(wd.de.plots, "intra_nnd_ppp"), recursive = TRUE)

# -------------------------- Libraries ----------------------------------------
suppressPackageStartupMessages({
	library(dplyr); library(tidyr); library(stringr); library(tibble); library(readr); library(purrr)
	library(Seurat)
	library(spatstat.geom); library(spatstat.explore); library(spatstat.core)
	library(dbscan); library(FNN); library(spdep); library(scales)
})

# -------------------------- Data load ----------------------------------------
wd.nf.data  <- file.path(wd.rna, "ngs/pacbio/SComatic/results/Step4_VariantCalling/new_beta")

barcodes <- read.table(file = file.path(wd.nf.data, "scomatics_PD53626b.tsv"),
							  sep = "\t", header = TRUE, stringsAsFactors = FALSE) |>
	subset(Cell_type_observed == "germ")

meta <- read.table(file = file.path(wd.nf.data, "celltypes.tsv"),
						 sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Load Seurat objects
# so.st is your cell subset for this sample/gene; so.integrated holds S/G2M scores
load(file=file.path(
	wd.de.plots,
	paste0("ssc_filtered_normalised_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_ordered_annotated_meta_cells_only_", GENE, "_RNA_so.st.RData")
))
# If not already loaded in your session, load so.integrated from your data dir
if (!exists("so.integrated")) {
	# Adjust filename if yours differs
	int_rdata <- file.path(wd.de.data, "ssc_filtered_normalised_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_ordered_annotated_monocle3+phase_RNA.RData")
	if (file.exists(int_rdata)) load(int_rdata)
}

# -------------------------- Stage order --------------------------------------
stage_order <- c(
	"Stage 0","Stage 0A","Stage 0B","Stage 1","Stage 2","Stage 3",
	"Leptotene","Zygotene","Pachytene","Diplotene",
	"Meiotic division","Early spermatid","Late spermatid"
)

# -------------------------- Callable coverage at locus -----------------------
MIN_TOTAL_UMI <- 1L
MIN_ALT_UMI   <- 1L
COL_READS <- "Num_reads"

getcol <- function(df, nm) if (nm %in% names(df)) df[[nm]] else NULL

ref_calls <- barcodes |>
	dplyr::filter(.data$CHROM == CHR, .data$Start == POS, .data$Base_observed == REF) |>
	dplyr::select(CB, any_of(c("REF_UMI", COL_READS))) |>
	dplyr::group_by(CB) |>
	dplyr::summarise(
		ref_umis = if (!is.null(getcol(cur_data_all(), "REF_UMI")))
			sum(.data[["REF_UMI"]], na.rm = TRUE) else sum(.data[[COL_READS]], na.rm = TRUE),
		.groups = "drop"
	)

alt_calls <- barcodes |>
	dplyr::filter(.data$CHROM == CHR, .data$Start == POS, .data$Base_observed == ALT) |>
	dplyr::select(CB, any_of(c("ALT_UMI", COL_READS))) |>
	dplyr::group_by(CB) |>
	dplyr::summarise(
		alt_umis = if (!is.null(getcol(cur_data_all(), "ALT_UMI")))
			sum(.data[["ALT_UMI"]], na.rm = TRUE) else sum(.data[[COL_READS]], na.rm = TRUE),
		.groups = "drop"
	)

cov_tbl <- dplyr::full_join(ref_calls, alt_calls, by = "CB") |>
	dplyr::mutate(
		ref_umis   = dplyr::coalesce(ref_umis, 0L),
		alt_umis   = dplyr::coalesce(alt_umis, 0L),
		total_umis = ref_umis + alt_umis
	)

# Stage labels (prefer Seurat cell_type on so.st)
seu_stage <- Seurat::FetchData(so.st, vars = "cell_type") |>
	tibble::rownames_to_column("CB") |>
	dplyr::rename(stage_seu = cell_type)

meta_stage <- meta |>
	dplyr::select(Index, Cell_type) |>
	dplyr::rename(CB = Index, stage_meta = Cell_type)

calls2 <- cov_tbl |>
	dplyr::left_join(seu_stage, by = "CB") |>
	dplyr::left_join(meta_stage, by = "CB") |>
	dplyr::mutate(
		stage_raw = dplyr::if_else(!is.na(stage_seu), stage_seu, stage_meta),
		stage_raw = stringr::str_squish(as.character(stage_raw)),
		stage     = dplyr::if_else(stage_raw %in% stage_order, stage_raw, NA_character_)
	) |>
	dplyr::filter(!is.na(stage)) |>
	dplyr::mutate(stage = factor(stage, levels = stage_order),
					  callable = total_umis >= MIN_TOTAL_UMI,
					  alt_pos  = alt_umis   >= MIN_ALT_UMI)

# Save callable tallies
k_n_by_stage <- calls2 |>
	dplyr::filter(callable) |>
	dplyr::group_by(stage) |>
	dplyr::summarise(k = sum(alt_pos, na.rm = TRUE), n = dplyr::n(), .groups = "drop") |>
	dplyr::arrange(stage) |>
	tidyr::complete(stage = factor(stage_order, levels = stage_order), fill = list(k = 0L, n = 0L))
print(k_n_by_stage)
readr::write_tsv(k_n_by_stage, file.path(wd.de.plots, "intra_nnd_ppp", "k_n_by_stage.tsv"))

# -------------------------- Coordinates (FORCE original µm) -------------------
meta_df <- so.st@meta.data %>% tibble::rownames_to_column("CB")

if (!all(c("orig_x_um","orig_y_um") %in% names(meta_df))) {
	stop("orig_x_um/orig_y_um not found in so.st@meta.data")
}

coords_df <- meta_df %>%
	transmute(CB,
				 x = as.numeric(orig_x_um),
				 y = as.numeric(orig_y_um)) %>%
	distinct(CB, .keep_all = TRUE) %>%
	filter(is.finite(x), is.finite(y))

message("# coords_df (original µm) built: n=", nrow(coords_df),
		  ", dupCB=", sum(duplicated(coords_df$CB)))

# Optional sanity: how many callable cells will get coords?
message("# callable cells w/ coords (µm): ",
		  sum(calls2$callable & calls2$CB %in% coords_df$CB), " / ",
		  sum(calls2$callable))

# -------------------------- Per-cell callable frame (rebuild) -----------------
# (unchanged except it now joins the forced-µm coords above)
per_cell <- calls2 %>%
	dplyr::filter(callable) %>%
	dplyr::transmute(CB,
						  stage  = factor(stage, levels = stage_order),
						  mutant = as.integer(alt_pos)) %>%
	dplyr::inner_join(coords_df, by = "CB") %>%
	dplyr::inner_join(sg_probs,  by = "CB") %>%
	dplyr::filter(is.finite(x), is.finite(y))

message("# Callable (joined) per_cell using ORIGINAL µm coords: n=", nrow(per_cell))
print(
	per_cell %>%
		dplyr::count(stage, mutant) %>%
		tidyr::pivot_wider(names_from = mutant, values_from = n, values_fill = 0) %>%
		dplyr::rename(WT = `0`, MUT = `1`) %>%
		dplyr::mutate(total = WT + MUT)
)

# -------------------------- Per-cell S/G2M probabilities ---------------------
# Pull S.Score & G2M.Score from so.integrated for the same CBs
meta_int <- so.integrated@meta.data |>
	tibble::rownames_to_column("CB")

if (!all(c("S.Score","G2M.Score") %in% colnames(meta_int))) {
	stop("so.integrated@meta.data must contain S.Score and G2M.Score")
}

# Stable 3-way softmax to convert (G1=0, S=S.Score, G2M=G2M.Score) -> probabilities
softmax3 <- function(s, g2m) {
	# logits: (0, s, g2m)
	m <- pmax(0, pmax(s, g2m))  # stability shift
	e0   <- exp(0 - m)
	eS   <- exp(s - m)
	eG2M <- exp(g2m - m)
	Z <- e0 + eS + eG2M
	cbind(G1_prob = e0/Z, S_prob = eS/Z, G2M_prob = eG2M/Z)
}

sg_tbl <- meta_int |>
	dplyr::select(CB, S.Score, G2M.Score) |>
	dplyr::mutate(
		S.Score   = as.numeric(S.Score),
		G2M.Score = as.numeric(G2M.Score)
	)

probs <- softmax3(s = sg_tbl$S.Score, g2m = sg_tbl$G2M.Score)
sg_probs <- cbind(sg_tbl["CB"], as.data.frame(probs)) |>
	tibble::as_tibble()

# -------------------------- Per-cell callable frame --------------------------
per_cell <- calls2 |>
	dplyr::filter(callable) |>
	dplyr::transmute(CB,
						  stage = factor(stage, levels = stage_order),
						  mutant = as.integer(alt_pos)) |>
	dplyr::inner_join(coords_df, by = "CB") |>
	dplyr::inner_join(sg_probs, by = "CB") |>
	dplyr::filter(is.finite(x), is.finite(y))  # coordinates must be finite

message("# Callable (joined) per_cell: n=", nrow(per_cell))
print(
	per_cell %>%
		dplyr::count(stage, mutant) %>%
		tidyr::pivot_wider(names_from = mutant, values_from = n, values_fill = 0) %>%
		dplyr::rename(WT = `0`, MUT = `1`) %>%
		dplyr::mutate(total = WT + MUT)
)

# -------------------------- Spatial helpers ----------------------------------
ppp_from_df <- function(df_stage) {
	stopifnot(nrow(df_stage) >= 3)
	W <- owin(poly = spatstat.geom::convexhull.xy(df_stage$x, df_stage$y)$bdry[[1]])
	ppp(df_stage$x, df_stage$y, window = W,
		 marks = data.frame(mutant = df_stage$mutant,
		 						 S_prob = df_stage$S_prob,
		 						 CB = df_stage$CB))
}

nnd_mutant_perm <- function(ppp_all, B = 1999, seed = 1) {
	set.seed(seed)
	m <- ppp_all$marks$mutant
	n_mut <- sum(m == 1)
	if (n_mut < 2) return(tibble(nnd_mm = NA_real_, p_nnd = NA_real_))
	obs <- mean(nndist(ppp_all[m == 1]))
	N <- npoints(ppp_all)
	stats <- replicate(B, {
		lab <- rep(0L, N); lab[sample.int(N, n_mut)] <- 1L
		pts <- ppp_all[lab == 1L]
		if (npoints(pts) < 2) return(NA_real_)
		mean(nndist(pts))
	})
	stats <- stats[is.finite(stats)]
	p <- (1 + sum(stats <= obs)) / (1 + length(stats))  # one-sided, clustering => smaller NND
	tibble(nnd_mm = obs, p_nnd = p)
}

kinhom_summary <- function(ppp_all) {
	if (sum(ppp_all$marks$mutant == 1) < 2) return(tibble(LmR_kinhom_medr = NA_real_))
	lam_all <- density.ppp(ppp_all, at = "points", edge = TRUE)
	kin <- Kinhom(ppp_all[ppp_all$marks$mutant == 1],
					  lambda = lam_all[ppp_all$marks$mutant == 1])
	LmR <- sqrt(kin$iso/pi) - kin$r
	tibble(LmR_kinhom_medr = LmR[which.min(abs(kin$r - median(kin$r, na.rm = TRUE)))])
}

# Bivariate K: mutants (X) vs S-high (Y) within stage, inhomogeneous, random S relabeling
bivarK_Shigh_summary <- function(ppp_all, q = 0.75, B = 999, seed = 1) {
	set.seed(seed)
	S <- ppp_all$marks$S_prob
	fin <- is.finite(S)
	if (sum(fin) < 3) return(tibble(bivK_outside_prop = NA_real_, Sthr = NA_real_))
	thr <- as.numeric(quantile(S[fin], q, na.rm = TRUE))
	type <- ifelse(S >= thr, "Shigh", "Slow")
	
	X <- ppp_all[ppp_all$marks$mutant == 1]
	Y <- ppp_all[type == "Shigh"]
	if (npoints(X) < 2 || npoints(Y) < 2) return(tibble(bivK_outside_prop = NA_real_, Sthr = thr))
	
	XU <- spatstat.geom::unmark(X); YU <- spatstat.geom::unmark(Y)
	XY <- spatstat.geom::superimpose(X = XU, Y = YU, W = ppp_all$window)
	spatstat.geom::marks(XY) <- factor(spatstat.geom::marks(XY), levels = c("X","Y"))
	
	lambdaI <- density.ppp(XU, at = "points")
	lambdaJ <- density.ppp(YU, at = "points")
	
	Kxy <- Kcross.inhom(XY, i = "X", j = "Y", lambdaI = lambdaI, lambdaJ = lambdaJ)
	Lo  <- sqrt(Kxy$iso/pi) - Kxy$r
	
	sims <- replicate(B, {
		typ_perm <- sample(type, replace = FALSE)
		Yperm <- ppp_all[typ_perm == "Shigh"]
		if (npoints(Yperm) < 2) return(rep(NA_real_, length(Kxy$r)))
		YUP <- spatstat.geom::unmark(Yperm)
		XYP <- spatstat.geom::superimpose(X = XU, Y = YUP, W = ppp_all$window)
		spatstat.geom::marks(XYP) <- factor(spatstat.geom::marks(XYP), levels = c("X","Y"))
		Kp <- Kcross.inhom(XYP, i = "X", j = "Y",
								 lambdaI = lambdaI,
								 lambdaJ = density.ppp(YUP, at = "points"))
		sqrt(Kp$iso/pi) - Kp$r
	})
	sims <- sims[, colSums(!is.finite(sims)) == 0, drop = FALSE]
	if (ncol(sims) == 0) return(tibble(bivK_outside_prop = NA_real_, Sthr = thr))
	
	env_lo <- apply(sims, 1, quantile, 0.025, na.rm = TRUE)
	env_hi <- apply(sims, 1, quantile, 0.975, na.rm = TRUE)
	
	tibble(bivK_outside_prop = mean(Lo > env_hi | Lo < env_lo, na.rm = TRUE),
			 Sthr = thr)
}

# Moran-type: correlation(mutant, spatial lag of S_prob)
moran_bivariate <- function(df_stage, k = 10, B = 1999, seed = 1) {
	n <- nrow(df_stage)
	if (n < 3) return(tibble(moran_bivar = NA_real_, p_moran = NA_real_))
	
	xy <- as.matrix(df_stage[, c("x","y")])
	k  <- min(k, max(1, n - 1), max(1, floor(n/3)))  # tame for small n
	
	knn <- try(spdep::knearneigh(xy, k = k), silent = TRUE)
	if (inherits(knn, "try-error")) return(tibble(moran_bivar = NA_real_, p_moran = NA_real_))
	nb <- try(spdep::knn2nb(knn), silent = TRUE)
	if (inherits(nb, "try-error")) return(tibble(moran_bivar = NA_real_, p_moran = NA_real_))
	lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
	
	y <- df_stage$S_prob
	if (!any(is.finite(y))) return(tibble(moran_bivar = NA_real_, p_moran = NA_real_))
	y_lag <- spdep::lag.listw(lw, y, zero.policy = TRUE)
	if (!any(is.finite(y_lag)) || sd(y_lag, na.rm = TRUE) == 0) {
		return(tibble(moran_bivar = NA_real_, p_moran = NA_real_))
	}
	
	x <- df_stage$mutant
	Iobs <- suppressWarnings(cor(x, y_lag, use = "pairwise.complete.obs"))
	
	set.seed(seed)
	stats <- replicate(B, {
		xp <- sample(x, replace = FALSE)
		suppressWarnings(cor(xp, y_lag, use = "pairwise.complete.obs"))
	})
	p <- (1 + sum(abs(stats) >= abs(Iobs))) / (1 + B)
	
	tibble(moran_bivar = Iobs, p_moran = p)
}

clusters_dbscan <- function(df_stage, eps = NULL, minPts = 3) {
	dm <- df_stage |> dplyr::filter(mutant == 1)
	if (nrow(dm) < minPts) return(tibble(n_clusters = 0, largest = 0, total_hull_area = 0, eps = NA_real_))
	if (is.null(eps)) {
		nnd <- FNN::get.knn(as.matrix(dm[,c("x","y")]), k = 1)$nn.dist[,1]
		eps <- as.numeric(quantile(nnd, 0.10, na.rm = TRUE))
		if (!is.finite(eps)) eps <- max(nnd[is.finite(nnd)], na.rm = TRUE)
	}
	cl <- dbscan(as.matrix(dm[,c("x","y")]), eps = eps, minPts = minPts)
	dm$cl <- cl$cluster
	ncl <- max(0, dm$cl, na.rm = TRUE)
	if (ncl == 0) return(tibble(n_clusters = 0, largest = 0, total_hull_area = 0, eps = eps))
	areas <- dm |> dplyr::filter(cl > 0) |>
		dplyr::group_by(cl) |>
		dplyr::reframe(area = { if (n() < 3) 0 else { ch <- chull(x, y); sum(spatstat.geom::ewpolyarea(x[ch], y[ch])) } }) |>
		dplyr::pull(area)
	tibble(n_clusters = ncl,
			 largest = max(table(dm$cl[dm$cl > 0])),
			 total_hull_area = sum(areas),
			 eps = eps)
}

kde_ratio_q95 <- function(df_stage, grid_n = 40) {
	if (nrow(df_stage) < 5) return(tibble(kde_ratio_q95 = NA_real_))
	ppp_all <- ppp_from_df(df_stage)
	kd_all <- density.ppp(ppp_all, at = "pixels", dimyx = c(grid_n, grid_n))$v
	kd_mut <- density.ppp(ppp_all[ppp_all$marks$mutant == 1], at = "pixels", dimyx = c(grid_n, grid_n))$v
	ratio <- (kd_mut + 1e-9) / (kd_all + 1e-9)
	tibble(kde_ratio_q95 = as.numeric(quantile(ratio[is.finite(ratio)], 0.95, na.rm = TRUE)))
}

# -------------------------- Per-stage analysis --------------------------------
analyze_stage <- function(df, stage_name, B_perm = 1999, min_mutants_heavy = 5) {
	d0 <- df |> dplyr::filter(stage == stage_name)
	out <- tibble(
		stage = stage_name,
		n_callable = nrow(d0),
		n_mut = sum(d0$mutant),
		n_wt  = nrow(d0) - sum(d0$mutant),
		S_prob_mean = mean(d0$S_prob, na.rm = TRUE),
		G1_prob_mean = mean(d0$G1_prob, na.rm = TRUE),
		G2M_prob_mean = mean(d0$G2M_prob, na.rm = TRUE)
	)
	if (nrow(d0) < 5 || sum(d0$mutant) < 1) return(dplyr::bind_cols(out, tibble(note = "Too few callable or no mutants")))
	ppp_all <- try(ppp_from_df(d0), silent = TRUE)
	if (inherits(ppp_all, "try-error")) return(dplyr::bind_cols(out, tibble(note = "Degenerate window")))
	# Light metrics
	nnd <- nnd_mutant_perm(ppp_all, B = B_perm)
	kin <- kinhom_summary(ppp_all)
	out <- dplyr::bind_cols(out, nnd, kin)
	# S-based association (now per-cell)
	kb  <- bivarK_Shigh_summary(ppp_all, q = 0.75, B = ceiling(B_perm/2))
	mor <- moran_bivariate(d0, k = 10, B = B_perm)
	out <- dplyr::bind_cols(out, kb, mor)
	# Heavy metrics only when powered
	if (sum(d0$mutant) < min_mutants_heavy) {
		out$note <- "Limited power; heavy clustering metrics skipped (n_mut < 5)"
		return(out)
	}
	cl  <- clusters_dbscan(d0, minPts = 3)
	kdq <- kde_ratio_q95(d0, grid_n = 40)
	dplyr::bind_cols(out, cl, kdq)
}

run_all_stages <- function(df, stages = NULL, B_perm = 1999) {
	if (is.null(stages)) {
		stages <- df |> dplyr::group_by(stage) |>
			dplyr::summarise(n_mut = sum(mutant), .groups="drop") |>
			dplyr::filter(n_mut >= 1) |> dplyr::pull(stage) |> as.character()
	}
	purrr::map_dfr(stages, function(stg) {
		out <- try(analyze_stage(df, stg, B_perm = B_perm), silent = TRUE)
		if (inherits(out, "try-error")) {
			tibble::tibble(stage = stg, note = paste0("Stage failed: ", as.character(out)))
		} else out
	})
}

# -------------------------- Run per-stage -------------------------------------
stage_spatial <- run_all_stages(per_cell, B_perm = 1999) |>
	dplyr::arrange(factor(stage, levels = stage_order))

message("\n# Within-stage spatial summary (callable-only, with per-cell S/G2M):")
print(stage_spatial, n = 100)
readr::write_tsv(stage_spatial, file.path(wd.de.plots, "intra_nnd_ppp", "stage_spatial_summary_with_S.tsv"))

# -------------------------- Between-stage correlations ------------------------
# Use per-stage mean S_prob (from callable cells) vs clustering metrics
stage_summary <- stage_spatial |>
	dplyr::mutate(stage = factor(stage, levels = stage_order)) |>
	dplyr::arrange(stage)

btw_df <- stage_summary |>
	dplyr::transmute(stage,
						  S_prob = S_prob_mean,
						  neg_nnd = -nnd_mm,
						  LmR_med = LmR_kinhom_medr,
						  bivK_outside_prop = bivK_outside_prop,
						  moran_bivar = moran_bivar,
						  n_clusters = n_clusters,
						  largest = largest,
						  kde_ratio_q95 = kde_ratio_q95)

spearman_safe <- function(x, y) {
	ok <- is.finite(x) & is.finite(y)
	if (sum(ok) < 3) return(c(rho = NA_real_, p = NA_real_, n = sum(ok)))
	ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
	c(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}

metrics_to_check <- c("neg_nnd","LmR_med","bivK_outside_prop","moran_bivar","n_clusters","largest","kde_ratio_q95")
btw_results <- tibble::tibble(metric = metrics_to_check) |>
	dplyr::rowwise() |>
	dplyr::mutate(
		rho = spearman_safe(btw_df$S_prob, btw_df[[metric]])["rho"],
		p   = spearman_safe(btw_df$S_prob, btw_df[[metric]])["p"],
		n   = as.integer(spearman_safe(btw_df$S_prob, btw_df[[metric]])["n"])
	) |>
	dplyr::ungroup()

message("\n# Between-stage: mean per-cell S_prob vs clustering metrics (Spearman)")
print(btw_results)
readr::write_tsv(btw_results, file.path(wd.de.plots, "intra_nnd_ppp", "between_stage_correlations_with_S.tsv"))

# =============================================================================
# NND difference via label-permutation within stage (handles n_wt != n_mut)
# Computes:
#   NND_M = mean NND among mutants
#   NND_W = mean NND among WT
#   delta = NND_M - NND_W
#   p_perm = two-sided permutation p with labels shuffled preserving counts
# Also returns z_perm = (delta - mean_perm)/sd_perm for an effect size
# =============================================================================

mean_nnd_pp <- function(xy) {
	if (nrow(xy) < 2) return(NA_real_)
	W  <- suppressWarnings(spatstat.geom::owin(
		poly = spatstat.geom::convexhull.xy(xy[,1], xy[,2])$bdry[[1]]
	))
	pp <- spatstat.geom::ppp(xy[,1], xy[,2], window = W)
	mean(spatstat.geom::nndist(pp))
}

stage_nnd_delta_perm <- function(df_stage, B = 4999, seed = 1) {
	stg <- as.character(df_stage$stage[1])
	n_m <- sum(df_stage$mutant == 1)
	n_w <- sum(df_stage$mutant == 0)
	n   <- nrow(df_stage)
	
	out <- tibble::tibble(
		stage = stg, n_mut = n_m, n_wt = n_w, n_tot = n,
		NND_M = NA_real_, NND_W = NA_real_,
		delta = NA_real_, p_perm = NA_real_, z_perm = NA_real_
	)
	
	# need at least 2 per group to form NNDs
	if (n_m < 2 || n_w < 2) return(out)
	
	XY <- as.matrix(df_stage[, c("x","y")])
	lab <- as.integer(df_stage$mutant)  # 1 = mutant, 0 = WT
	
	NND_M <- mean_nnd_pp(XY[lab==1,,drop=FALSE])
	NND_W <- mean_nnd_pp(XY[lab==0,,drop=FALSE])
	if (!is.finite(NND_M) || !is.finite(NND_W)) return(out)
	
	delta_obs <- NND_M - NND_W
	
	set.seed(seed)
	stats <- replicate(B, {
		labp <- sample(lab, replace = FALSE)  # preserve counts
		m <- mean_nnd_pp(XY[labp==1,,drop=FALSE])
		w <- mean_nnd_pp(XY[labp==0,,drop=FALSE])
		if (!is.finite(m) || !is.finite(w)) return(NA_real_)
		m - w
	})
	stats <- stats[is.finite(stats)]
	if (!length(stats)) return(out)
	
	p <- (1 + sum(abs(stats) >= abs(delta_obs))) / (1 + length(stats))
	z <- (delta_obs - mean(stats)) / sd(stats)
	
	out$NND_M <- NND_M; out$NND_W <- NND_W
	out$delta <- delta_obs; out$p_perm <- p; out$z_perm <- z
	out
}

# ---- run per stage
nnd_delta_by_stage <- per_cell %>%
	dplyr::group_split(stage) %>%
	purrr::map_dfr(~ stage_nnd_delta_perm(.x, B = 4999, seed = 1)) %>%
	dplyr::arrange(factor(stage, levels = stage_order))

message("\n# NND difference (mutant − WT) via label-permutation within stage:")
print(nnd_delta_by_stage, n = 50)
readr::write_tsv(nnd_delta_by_stage,
					  file.path(wd.de.plots, "intra_nnd_ppp", "nnd_delta_by_stage.tsv"))

# -----------------------------------------------------------------------------
# Define the pooled order you want in the output
# -----------------------------------------------------------------------------
pool_levels <- c("Stage 0 pooled", as.character(stage_order))

nnd_delta_pooled <- per_cell_pooled %>%
	dplyr::group_split(stage_pool, .keep = TRUE) %>%
	purrr::map_dfr(function(d) {
		# overwrite 'stage' with the pooled label so the helper uses/returns it
		pooled_lab <- as.character(unique(d$stage_pool))
		d2 <- d %>% dplyr::mutate(stage = pooled_lab)
		stage_nnd_delta_perm(d2, B = 4999, seed = 1)
	}) %>%
	dplyr::arrange(factor(stage, levels = pool_levels))

message("\n# NND difference (mutant − WT) via label-permutation — pooled Stage 0/0A/0B:")
print(nnd_delta_pooled, n = 50)
readr::write_tsv(nnd_delta_pooled,
					  file.path(wd.de.plots, "origUM", "nnd_delta_by_stage_pooled.tsv"))

# -----------------------------------------------------------------------------
# ---- Permutation CI for delta NND (reuses your stage_nnd_delta_perm logic) ----
# -----------------------------------------------------------------------------
stage_nnd_delta_perm_ci <- function(df_stage, B = 4999, seed = 1) {
	stg <- as.character(df_stage$stage[1])
	n_m <- sum(df_stage$mutant == 1); n_w <- sum(df_stage$mutant == 0)
	n   <- nrow(df_stage)
	res <- tibble::tibble(stage = stg, n_mut = n_m, n_wt = n_w, n_tot = n,
								 NND_M = NA_real_, NND_W = NA_real_,
								 delta = NA_real_, p_perm = NA_real_,
								 sd_perm = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_)
	if (n_m < 2 || n_w < 2) return(res)
	XY  <- as.matrix(df_stage[, c("x","y")])
	lab <- as.integer(df_stage$mutant)
	
	mean_nnd_pp <- function(xy) {
		if (nrow(xy) < 2) return(NA_real_)
		W  <- suppressWarnings(spatstat.geom::owin(
			poly = spatstat.geom::convexhull.xy(xy[,1], xy[,2])$bdry[[1]]
		))
		pp <- spatstat.geom::ppp(xy[,1], xy[,2], window = W)
		mean(spatstat.geom::nndist(pp))
	}
	
	NND_M <- mean_nnd_pp(XY[lab==1,,drop=FALSE])
	NND_W <- mean_nnd_pp(XY[lab==0,,drop=FALSE])
	if (!is.finite(NND_M) || !is.finite(NND_W)) return(res)
	delta_obs <- NND_M - NND_W
	
	set.seed(seed)
	stats <- replicate(B, {
		labp <- sample(lab, replace = FALSE)
		m <- mean_nnd_pp(XY[labp==1,,drop=FALSE])
		w <- mean_nnd_pp(XY[labp==0,,drop=FALSE])
		if (!is.finite(m) || !is.finite(w)) return(NA_real_)
		m - w
	})
	stats <- stats[is.finite(stats)]
	if (!length(stats)) return(res)
	
	p  <- (1 + sum(abs(stats) >= abs(delta_obs))) / (1 + length(stats))
	sd <- stats::sd(stats)
	ci <- delta_obs + stats::quantile(stats - mean(stats), probs = c(0.025, 0.975), na.rm = TRUE)
	
	res$NND_M  <- NND_M; res$NND_W <- NND_W
	res$delta  <- delta_obs; res$p_perm <- p
	res$sd_perm <- sd; res$ci_lo <- unname(ci[1]); res$ci_hi <- unname(ci[2])
	res
}

# Run on powered stages (your pooled object or original per_cell)
nnd_delta_pooled_ci <- per_cell_pooled %>%
	dplyr::group_split(stage_pool) %>%
	purrr::map_dfr(function(d){
		lab <- as.character(unique(d$stage_pool))
		stage_nnd_delta_perm_ci(d %>% dplyr::mutate(stage = lab), B = 4999, seed = 1)
	})

message("\n# ΔNND with permutation SD and 95% CI:")
print(nnd_delta_pooled_ci, n = 50)
readr::write_tsv(nnd_delta_pooled_ci,
					  file.path(wd.de.plots, "origUM", "nnd_delta_by_stage_pooled_CI.tsv"))

# -----------------------------------------------------------------------------
# ---- Small-r pair-correlation difference: mean_{r<=r_max}[ gMM(r) - gWW(r) ] ----
# -----------------------------------------------------------------------------
pcf_delta_smallr <- function(df_stage, r_max = 50, B = 1999, seed = 1) {
	stg <- as.character(df_stage$stage[1])
	n_m <- sum(df_stage$mutant == 1); n_w <- sum(df_stage$mutant == 0)
	out <- tibble::tibble(stage = stg, n_mut = n_m, n_wt = n_w,
								 pcf_delta = NA_real_, r_max = r_max,
								 p_perm = NA_real_, z_perm = NA_real_)
	if (n_m < 3 || n_w < 3) return(out)  # need a few points for pcf stability
	
	W  <- suppressWarnings(spatstat.geom::owin(
		poly = spatstat.geom::convexhull.xy(df_stage$x, df_stage$y)$bdry[[1]]
	))
	P  <- spatstat.geom::ppp(df_stage$x, df_stage$y, window = W,
									 marks = factor(ifelse(df_stage$mutant==1,"M","W")))
	
	# Inhom pair correlation for each group separately
	g_group <- function(pp, grp) {
		sub <- pp[pp$marks==grp]
		lam <- spatstat.core::density.ppp(sub, at = "points")
		# pcfinhom returns g(r) on a grid of r
		spatstat.core::pcfinhom(sub, lambda = lam, stoyan = 0.15, correction = "Ripley")
	}
	
	Gm <- g_group(P, "M"); Gw <- g_group(P, "W")
	r  <- Gm$r
	keep <- which(r <= r_max & is.finite(Gm$iso) & is.finite(Gw$iso))
	if (!length(keep)) return(out)
	delta_obs <- mean(Gm$iso[keep] - Gw$iso[keep], na.rm = TRUE)
	
	set.seed(seed)
	stats <- replicate(B, {
		# permute labels preserving counts
		labp <- sample(P$marks, replace = FALSE)
		Pp <- P; spatstat.geom::marks(Pp) <- labp
		Gmp <- g_group(Pp, "M"); Gwp <- g_group(Pp, "W")
		if (!all(length(Gmp$iso) == length(r), length(Gwp$iso) == length(r))) return(NA_real_)
		kp <- which(r <= r_max & is.finite(Gmp$iso) & is.finite(Gwp$iso))
		if (!length(kp)) return(NA_real_)
		mean(Gmp$iso[kp] - Gwp$iso[kp], na.rm = TRUE)
	})
	stats <- stats[is.finite(stats)]
	if (!length(stats)) return(out)
	
	p <- (1 + sum(abs(stats) >= abs(delta_obs))) / (1 + length(stats))
	z <- (delta_obs - mean(stats)) / stats::sd(stats)
	
	out$pcf_delta <- delta_obs; out$p_perm <- p; out$z_perm <- z
	out
}

# Run on your pooled object
pcf_smallr_pooled <- per_cell_pooled %>%
	dplyr::group_split(stage_pool) %>%
	purrr::map_dfr(function(d){
		lab <- as.character(unique(d$stage_pool))
		pcf_delta_smallr(d %>% dplyr::mutate(stage = lab), r_max = 50, B = 1999, seed = 1)
	}) %>%
	dplyr::arrange(match(stage, c("Stage 0 pooled", as.character(stage_order))))

message("\n# Small-r pair-correlation difference (gMM - gWW), r<=50 µm:")
print(pcf_smallr_pooled, n = 50)
readr::write_tsv(pcf_smallr_pooled,
					  file.path(wd.de.plots, "origUM", "pcf_smallr_by_stage_pooled.tsv"))

# =============================================================================
# Pool stages into 4 composite cell-cycle bins
# =============================================================================
per_cell_pooled4 <- per_cell %>%
	mutate(stage_pool4 = case_when(
		stage %in% c("Stage 0","Stage 0A","Stage 0B","Stage 1") ~ "G0/G1 (pre)",
		stage %in% c("Stage 2","Stage 3","Leptotene") ~ "S",
		stage %in% c("Zygotene","Pachytene","Diplotene","Meiotic division") ~ "G2/M",
		stage %in% c("Early spermatid","Late spermatid") ~ "G0 (post)",
		TRUE ~ as.character(stage)
	)) %>%
	mutate(stage_pool4 = factor(stage_pool4,
										 levels = c("G0/G1 (pre)", "S", "G2/M", "G0 (post)")))

# Quick sanity check (force dplyr:: to avoid plyr/other masks)
counts_pooled4 <- per_cell_pooled4 %>%
	dplyr::count(stage_pool4, mutant, name = "n") %>%
	tidyr::pivot_wider(names_from = mutant, values_from = n, values_fill = 0) %>%
	dplyr::rename(WT = `0`, MUT = `1`) %>%
	dplyr::mutate(total = WT + MUT) %>%
	dplyr::arrange(stage_pool4)

print(counts_pooled4, n = 10)

# =============================================================================
# ΔNND (mutant − WT) with permutation CI per pooled stage
# =============================================================================
nnd_delta_pooled4 <- per_cell_pooled4 %>%
	dplyr::group_split(stage_pool4, .keep = TRUE) %>%
	purrr::map_dfr(function(d){
		pool_lab <- as.character(unique(d$stage_pool4))
		d2 <- d %>% dplyr::mutate(stage = pool_lab)
		stage_nnd_delta_perm_ci(d2, B = 20000, seed = 1)
	}) %>%
	dplyr::arrange(factor(stage, levels = c("G0/G1 (pre)", "S", "G2/M", "G0 (post)")))

message("\n# ΔNND (mutant − WT) by pooled 4-phase groups:")
print(nnd_delta_pooled4, n = 10)

readr::write_tsv(nnd_delta_pooled4,
					  file.path(wd.de.plots, "origUM", "nnd_delta_by_4phase_pooled_B20000.tsv"))

# =============================================================================
# Optional: small-r pair-correlation (gMM − gWW, r ≤ 50 µm)
# =============================================================================
pcf_smallr_pooled4 <- per_cell_pooled4 %>%
	dplyr::group_split(stage_pool4, .keep = TRUE) %>%
	purrr::map_dfr(function(d){
		pool_lab <- as.character(unique(d$stage_pool4))
		pcf_delta_smallr(d %>% dplyr::mutate(stage = pool_lab), r_max = 50, B = 1999, seed = 1)
	}) %>%
	dplyr::arrange(factor(stage, levels = c("G0/G1 (pre)", "S", "G2/M", "G0 (post)")))

message("\n# Small-r pair-correlation difference (gMM − gWW), r≤50 µm by 4-phase pools:")
print(pcf_smallr_pooled4, n = 10)

readr::write_tsv(pcf_smallr_pooled4,
					  file.path(wd.de.plots, "origUM", "pcf_smallr_by_4phase_pooled.tsv"))

# =============================================================================
# Build global PC1 axis once (if not already)
# =============================================================================
XY_all <- per_cell %>% dplyr::select(x, y) %>% as.matrix()
pc <- prcomp(XY_all, center = TRUE, scale. = FALSE)
u1 <- pc$rotation[,1, drop = FALSE]
s_all <- as.numeric(scale(XY_all, center = pc$center, scale = FALSE) %*% u1)

per_cell_axis <- per_cell %>%
	dplyr::mutate(s_axis = s_all) %>%
	dplyr::select(CB, stage, mutant, s_axis, x, y, S_prob, G1_prob, G2M_prob)

# Pool to 4 groups
per_cell_axis4 <- per_cell_axis %>%
	dplyr::mutate(stage_pool4 = dplyr::case_when(
		stage %in% c("Stage 0","Stage 0A","Stage 0B","Stage 1") ~ "G0/G1 (pre)",
		stage %in% c("Stage 2","Stage 3","Leptotene") ~ "S",
		stage %in% c("Zygotene","Pachytene","Diplotene","Meiotic division") ~ "G2/M",
		stage %in% c("Early spermatid","Late spermatid") ~ "G0 (post)",
		TRUE ~ as.character(stage)
	)) %>%
	dplyr::mutate(stage_pool4 = factor(stage_pool4,
												  levels = c("G0/G1 (pre)", "S", "G2/M", "G0 (post)")))

# Δsd (mutants − matched-WT) along axis with resampling
axis_spread_pooled4 <- per_cell_axis4 %>%
	dplyr::group_split(stage_pool4, .keep = TRUE) %>%
	purrr::map_dfr(function(d){
		pool_lab <- as.character(unique(d$stage_pool4))
		n_mut <- sum(d$mutant==1); n_wt <- sum(d$mutant==0)
		if (n_mut < 2 || n_wt < n_mut) {
			return(tibble::tibble(stage = pool_lab, n_mut=n_mut, n_wt=n_wt,
										 sd_mut=NA_real_, sd_wt_mean=NA_real_, sd_wt_sd=NA_real_,
										 delta_sd=NA_real_, p_delta_sd=NA_real_))
		}
		sd_mut <- stats::sd(d$s_axis[d$mutant==1])
		s_wt   <- d$s_axis[d$mutant==0]
		R <- 9999
		sd_wt <- replicate(R, { stats::sd(sample(s_wt, n_mut, replace = FALSE)) })
		sd_wt <- sd_wt[is.finite(sd_wt)]
		delta_sd <- sd_mut - mean(sd_wt)
		p_delta_sd <- (1 + sum(abs(sd_wt - mean(sd_wt)) >= abs(delta_sd))) / (1 + length(sd_wt))
		tibble::tibble(stage = pool_lab, n_mut=n_mut, n_wt=n_wt,
							sd_mut=sd_mut, sd_wt_mean=mean(sd_wt), sd_wt_sd=sd(sd_wt),
							delta_sd=delta_sd, p_delta_sd=p_delta_sd)
	}) %>%
	dplyr::arrange(factor(stage, levels = c("G0/G1 (pre)", "S", "G2/M", "G0 (post)")))

message("\n# Axis Δsd (mutant − WT) along tubule PC1 by 4 pools:")
print(axis_spread_pooled4, n = 10)
readr::write_tsv(axis_spread_pooled4,
					  file.path(wd.de.plots, "origUM", "axis_spread_by_4phase_pooled.tsv"))

library(ggplot2)
library(dplyr)

plot_df <- nnd_delta_pooled4 %>%
	mutate(stage = factor(stage, levels = c("G0/G1 (pre)", "S", "G2/M", "G0 (post)")))

p1 <- ggplot(plot_df, aes(x = stage, y = delta)) +
	geom_hline(yintercept = 0, linetype = 2) +
	geom_point(size = 3) +
	geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.15) +
	labs(x = NULL, y = expression(Delta*"NND (µm), mutant − WT"),
		  title = "Nearest-neighbor difference by pooled cycle phase") +
	theme_classic() +
	# add per-row p-value labels at the top of each CI
	geom_text(
		data = dplyr::filter(plot_df, is.finite(ci_hi), is.finite(p_perm)),
		aes(y = ci_hi, label = paste0("p=", signif(p_perm, 3))),
		vjust = -0.6, size = 3, na.rm = TRUE
	) +
	# give a little headroom so labels aren’t clipped
	expand_limits(y = max(plot_df$ci_hi, na.rm = TRUE) * 1.12) +
	coord_cartesian(clip = "off")

print(p1)

library(ggplot2)
library(dplyr)

plot_df2 <- axis_spread_pooled4 %>%
	mutate(stage = factor(stage, levels = c("G0/G1 (pre)", "S", "G2/M", "G0 (post)")))

p2 <- ggplot(plot_df2, aes(x = stage, y = delta_sd)) +
	geom_hline(yintercept = 0, linetype = 2) +
	geom_point(size = 3) +
	labs(x = NULL, y = expression(Delta*sd["axis"]~"(µm), mutant − WT"),
		  title = "Axis spread along tubule (PC1) by pooled phase") +
	theme_classic() +
	geom_text(
		data = dplyr::filter(plot_df2, is.finite(delta_sd), is.finite(p_delta_sd)),
		aes(y = delta_sd, label = paste0("p=", signif(p_delta_sd, 3))),
		vjust = -1, size = 3, na.rm = TRUE
	) +
	expand_limits(y = max(plot_df2$delta_sd, na.rm = TRUE) * 1.15) +
	coord_cartesian(clip = "off")

print(p2)

ggsave(file.path(wd.de.plots, "origUM", "pooled4_deltaNND_B20000.pdf"), p1, width = 6, height = 3.2)
ggsave(file.path(wd.de.plots, "origUM", "pooled4_deltaSDaxis.pdf"), p2, width = 6, height = 3.2)

