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
wd.de.plots <- file.path(wd.de.plots, "CCF", "anchor")
if (!dir.exists(wd.de.plots)) {
	ok <- base::dir.create(wd.de.plots, recursive = TRUE, showWarnings = FALSE, mode = "0775")
	if (!ok && !dir.exists(wd.de.plots)) stop("Failed to create: ", dest)
}

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

#barcodes.st <- subset(subset(subset(barcodes, CHROM == CHR), Start == POS), Base_observed == ALT)
#meta.st <- subset(meta, Index %in% barcodes.st$CB)

# -----------------------------------------------------------------------------
# Standard Seurat re-processing workflow
# 01_QC
# https://satijalab.org/seurat/articles/pbmc3k_tutorial
# -----------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(stringr)
library(Seurat)
library(ggplot2)
library(scales)

#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_23_res=0.5_-C0_ordered_annotated_monocle3+phase.RData")))
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_ordered_annotated_monocle3+phase_RNA.RData")))
#so.integrated.st <- subset(so.integrated, sample.id == "PD53626b_ST1")

# -----------------------------------------------------------------------------
# Get so.st
# Last Modified: 17/09/25
# -----------------------------------------------------------------------------
load(file=file.path(wd.de, "plots", GENE, MUT, "origUM", paste0("ssc_filtered_normalised_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_ordered_annotated_meta_cells_only_", GENE, "_RNA_so.st.RData")))

# -----------------------------------------------------------------------------
# Clonal expansion across 12 stages
# Prevalence p_i with misclassification (s_i, fpr_i)
# Posterior by grid, trend tests, ribbon plot
# 
# # Last Modified: 13/10/25
# -----------------------------------------------------------------------------
library(dplyr)
if (requireNamespace("conflicted", quietly = TRUE)) {
	conflicted::conflict_prefer("select",   "dplyr", quiet = TRUE)
	conflicted::conflict_prefer("filter",   "dplyr", quiet = TRUE)
	conflicted::conflict_prefer("summarise","dplyr", quiet = TRUE)
	conflicted::conflict_prefer("mutate",   "dplyr", quiet = TRUE)
	conflicted::conflict_prefer("rename",   "dplyr", quiet = TRUE)
}

stage_order <- c(
	"Stage 0", "Stage 0A", "Stage 0B", "Stage 1", "Stage 2", "Stage 3", "Leptotene", "Zygotene", "Pachytene", "Diplotene", "Meiotic division", "Early spermatid","Late spermatid"
)

# --- 0) Build per-cell call table (union of REF & ALT) ------------------------
# barcodes already restricted to Cell_type_observed == "germ"
ref_calls <- subset(subset(barcodes, CHROM == CHR), Base_observed == REF) %>%
	dplyr::select(CB, Num_reads) %>%
	group_by(CB) %>%
	summarise(ref_reads = sum(Num_reads), .groups = "drop")

alt_calls <- subset(subset(barcodes, CHROM == CHR), Base_observed == ALT) %>%
	dplyr::select(CB, Num_reads) %>%
	group_by(CB) %>%
	summarise(alt_reads = sum(Num_reads), .groups = "drop")

calls <- full_join(ref_calls, alt_calls, by = "CB") %>%
	mutate(
		call = dplyr::case_when(
			!is.na(alt_reads) ~ "ALT",     # prefer ALT if seen at least once
			!is.na(ref_reads) ~ "REF",
			TRUE ~ NA_character_
		)
	)

# 1) Get stage labels from Seurat (PD53626b_ST1) and meta as fallback
seu_stage <- FetchData(so.st, vars = "cell_type") %>%
	tibble::rownames_to_column("CB") %>%
	dplyr::rename(stage_seu = cell_type)

meta_stage <- meta %>%
	dplyr::select(Index, Cell_type) %>%
	dplyr::rename(CB = Index, stage_meta = Cell_type)

# 2) Join labels (prefer Seurat), clean whitespace
calls2 <- calls %>%
	left_join(seu_stage, by = "CB") %>%
	left_join(meta_stage, by = "CB") %>%
	mutate(stage_raw = if_else(!is.na(stage_seu), stage_seu, stage_meta),
			 stage_raw = str_squish(as.character(stage_raw)))

# 3) Keep only your 12 germ stages EXACTLY as given
calls2 <- calls2 %>%
	mutate(stage = if_else(stage_raw %in% stage_order, stage_raw, NA_character_)) %>%
	filter(!is.na(stage)) %>%
	mutate(stage = factor(stage, levels = stage_order))

# (Optional) see what got dropped (helpful sanity check)
dropped <- calls %>%
	anti_join(calls2 %>% dplyr::select(CB), by = "CB") %>%
	left_join(meta_stage, by = "CB") %>%
	dplyr::count(stage_meta, sort = TRUE)
if (nrow(dropped)) {
	message("Dropped (labels not in your 12-state scheme):")
	print(dropped, n = 20)
}

# 4) Tally k (ALT+ cells) and n (evaluable REF ∪ ALT) per stage
k_n_by_stage <- calls2 %>%
	group_by(stage) %>%
	summarise(
		k = sum(call == "ALT", na.rm = TRUE),
		n = n_distinct(CB),
		.groups = "drop"
	) %>%
	arrange(stage) %>%
	complete(stage = factor(stage_order, levels = stage_order),
				fill = list(k = 0L, n = 0L))

print(k_n_by_stage)

# 5) Export vectors for the prevalence model
k_vec <- as.integer(k_n_by_stage$k)
n_vec <- as.integer(k_n_by_stage$n)
names(k_vec) <- names(n_vec) <- as.character(k_n_by_stage$stage)

cat("\n# k (ALT-positive) per stage:\n")
print(k_vec)
cat("\n# n (evaluable) per stage:\n")
print(n_vec)

# -----------------------------------------------------------------------------
# Estimate stage-wise mutant prevalence (pᵢ) correcting for stage-specific sensitivity (sᵢ) 
# and false-positive rate (fprᵢ), then test for a monotonic increase across your 12 ordered stages 
# and plot a ribbon with 95% credible intervals
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
	library(dplyr); library(tidyr); library(stringr)
	library(Seurat); library(ggplot2); library(scales)
})

# ---- 1) Derive stage-specific sensitivity s_i from coverage among ALL cells ----
# Build per-cell metadata (stage + CB list)
all_cells <- FetchData(so.st, vars = "cell_type") %>%
	tibble::rownames_to_column("CB") %>%
	transmute(CB, stage = factor(str_squish(as.character(cell_type)), levels = stage_order)) %>%
	filter(!is.na(stage))

# Coverage at locus per CB from calls (any ref/alt reads => covered)
cov_tbl <- calls %>%
	mutate(total_reads = dplyr::coalesce(ref_reads,0) + dplyr::coalesce(alt_reads,0),
			 has_cov    = total_reads > 0) %>%
	select(CB, has_cov, total_reads)

# Coverage fraction among ALL cells of each stage
cov_stage <- all_cells %>%
	left_join(cov_tbl, by = "CB") %>%
	mutate(has_cov = tidyr::replace_na(has_cov, FALSE),
			 total_reads = tidyr::replace_na(total_reads, 0)) %>%
	group_by(stage) %>%
	summarise(n_cells = n(),
				 cov_frac = mean(has_cov),
				 reads_mu = ifelse(any(has_cov), mean(total_reads[has_cov]), 0),
				 .groups = "drop") %>%
	arrange(stage)

s_floor <- 0.20
s_max   <- 0.99

# pick Pachytene if present and reasonably powered; otherwise best-covered stage
anchor <- cov_stage |>
	dplyr::filter(stage == "Pachytene", n_cells >= 20)
if (nrow(anchor) == 0) {
	anchor <- cov_stage |>
		dplyr::filter(n_cells >= 20) |>
		dplyr::arrange(dplyr::desc(cov_frac)) |>
		dplyr::slice(1)
}
cov_ref <- anchor$cov_frac[1]
if (!is.finite(cov_ref) || cov_ref <= 0) {
	cov_ref <- stats::median(cov_stage$cov_frac, na.rm = TRUE)
}

# linear-by-coverage mapping anchored to s_max at cov_ref
s_by_stage <- cov_stage |>
	dplyr::mutate(
		s = pmin(s_max, pmax(s_floor, s_max * cov_frac / cov_ref))
	) |>
	dplyr::select(stage, s)

# ---- 2) False positive rate (stage-constant here; replace if you have neg controls) ----
fpr_default <- 1e-3
fpr_by_stage <- tibble(stage = factor(stage_order, levels = stage_order),
							  fpr = fpr_default)

# ---- 3) Assemble prevalence inputs per stage (k, n, s_i, fpr_i) ---------------
df_prev <- k_n_by_stage %>%
	mutate(stage = factor(as.character(stage), levels = stage_order)) %>%
	left_join(s_by_stage,   by = "stage") %>%
	left_join(fpr_by_stage, by = "stage")

# ---- 4) Posterior p_i on a grid with misclassification correction -------------
posterior_p_grid <- function(k, n, s, fpr, a=1, b=1, grid_len=5001) {
	p <- seq(0, 1, length.out = grid_len)
	q <- p * s + (1 - p) * fpr
	loglik   <- dbinom(k, size = n, prob = q, log = TRUE)
	logprior <- dbeta(p, shape1 = a, shape2 = b, log = TRUE)
	logpost  <- loglik + logprior
	w <- exp(logpost - max(logpost)); w <- w/sum(w)
	cdf <- cumsum(w)
	tibble(
		p_mean = sum(p*w),
		p_low  = p[which.max(cdf >= 0.025)],
		p_high = p[which.max(cdf >= 0.975)],
		p_mode = p[which.max(w)]
	)
}

# ------------------------ Overall adjusted CCF (pooled) ----------------------
summ <- df_prev %>%
	rowwise() %>%
	mutate(res = list(posterior_p_grid(k, n, s, fpr))) %>%
	tidyr::unnest(res) %>%
	ungroup() %>%
	mutate(stage_idx = as.integer(stage))

# Add an observed fraction column
summ <- df_prev %>%
	rowwise() %>%
	mutate(res = list(posterior_p_grid(k, n, s, fpr))) %>%
	tidyr::unnest(res) %>%
	ungroup() %>%
	mutate(stage_idx = as.integer(stage)) %>%
	mutate(p_obs = ifelse(n > 0, k / n, NA_real_))   # <-- NEW: observed detection rate

overall_posterior <- function(df, a=1, b=1, grid_len=5001) {
	p <- seq(0, 1, length.out = grid_len)
	logprior <- dbeta(p, a, b, log=TRUE)
	loglik <- Reduce(`+`, lapply(seq_len(nrow(df)), function(i) {
		q <- df$s[i] * p + (1 - p) * df$fpr[i]
		dbinom(df$k[i], df$n[i], q, log=TRUE)
	}))
	logpost <- logprior + loglik
	w <- exp(logpost - max(logpost)); w <- w/sum(w)
	cdf <- cumsum(w)
	list(p_grid = p, w = w,
		  mean = sum(p*w),
		  low  = p[which.max(cdf >= 0.025)],
		  high = p[which.max(cdf >= 0.975)])
}
ov <- overall_posterior(df_prev)
p_overall_adj <- ov$mean
ci_overall    <- c(ov$low, ov$high)

# ---- 5) Monotonic increase tests across ordered stages ------------------------
sp <- suppressWarnings(cor.test(summ$p_mean, summ$stage_idx, method = "spearman", exact = FALSE))
kd <- suppressWarnings(cor.test(summ$p_mean, summ$stage_idx, method = "kendall"))
fit <- lm(qlogis(pmin(pmax(p_mean,1e-6),1-1e-6)) ~ stage_idx, data = summ)

cat("\n--- Trend tests (misclassification-corrected p_i vs stage) ---\n")
cat(sprintf("Spearman rho = %.3f, p = %.3g\n", sp$estimate, sp$p.value))
cat(sprintf("Kendall tau  = %.3f, p = %.3g\n", kd$estimate, kd$p.value))
cat(sprintf("Logit slope per stage = %.3f (p = %.3g)\n",
				coef(summary(fit))['stage_idx','Estimate'],
				coef(summary(fit))['stage_idx','Pr(>|t|)']))

# ---- 6) Plot: ribbon (95% CrI) of p_i across stages --------------------------
lab_pct <- function(x) scales::percent(x, accuracy = 1)

p1 <- ggplot(summ, aes(x = stage_idx, y = p_mean, group = 1)) +
	geom_ribbon(aes(ymin = p_low, ymax = p_high), alpha = 0.20) +
	geom_line(size = 0.9) +
	geom_point(size = 2) +
	scale_x_continuous(breaks = summ$stage_idx, labels = levels(summ$stage)) +
	scale_y_continuous(labels = lab_pct, limits = c(0, NA)) +
	labs(x = "Cell type",
		  y = "Cancer cell fraction (CCF)",
		  title = "CBL (11:119284966 T>A)",
		  subtitle = "") +
	theme_classic(base_size = 12) +
	theme(axis.text.x = element_text(angle = 35, hjust = 1),
			plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5),
			axis.title = ggplot2::element_text(size = 14),
			axis.text  = ggplot2::element_text(size = 12))
print(p1)

# add observed detection rate (k/n) as a blue line + points
p1 <- p1 +
	geom_line(data = summ, aes(x = stage_idx, y = p_obs), linewidth = 0.9, color = "blue") +
	geom_point(data = summ, aes(x = stage_idx, y = p_obs), size = 2, color = "blue")

print(p1)

# Save
outfile <- file.path(wd.de.plots, paste0(GENE, "_prevalence_by_stage_misclass_corrected_0.20.png"))
ggsave(outfile, plot = p1, width = 6, height = 6, dpi = 300)
cat("Saved:", outfile, "\n")

p1 <- p1 +
	geom_hline(yintercept = p_overall_adj, linetype = "dashed", linewidth = 0.6) +
	annotate("text", x = 1, y = p_overall_adj, vjust = -0.6, hjust = 0,
				label = paste0("Overall (adj) = ",
									scales::percent(p_overall_adj, 1),
									" [", scales::percent(ci_overall[1], 1), ", ",
									scales::percent(ci_overall[2], 1), "]"),
				size = 3.5)
# Save
outfile <- file.path(wd.de.plots, paste0(GENE, "_prevalence_by_stage_misclass_corrected_adj_0.20.png"))
ggsave(outfile, plot = p1, width = 6, height = 6, dpi = 300)
cat("Saved:", outfile, "\n")

# ---- add legend for black vs blue lines --------------------------------------
# Tell ggplot what each line represents (so legend works)
p1 <- p1 +
	geom_line(data = summ, aes(x = stage_idx, y = p_mean, color = "Adjusted CCF"), linewidth = 0.9) +
	geom_point(data = summ, aes(x = stage_idx, y = p_mean, color = "Adjusted CCF"), size = 2) +
	geom_line(data = summ, aes(x = stage_idx, y = p_obs, color = "Observed detection rate"), linewidth = 0.9) +
	geom_point(data = summ, aes(x = stage_idx, y = p_obs, color = "Observed detection rate"), size = 2)

p1 <- p1 +
	scale_color_manual(
		name = NULL,
		values = c("Adjusted CCF" = "black", "Observed detection rate" = "blue")
	) +
	guides(
		color = guide_legend(override.aes = list(linewidth = 1.2, size = 3)),
		fill  = guide_legend(override.aes = list(alpha = 0.3))
	) +
	theme(legend.position = "top")

# Save
outfile <- file.path(wd.de.plots, paste0(GENE, "_prevalence_by_stage_misclass_corrected_adj_0.20_kn.png"))
ggsave(outfile, plot = p1, width = 6, height = 6, dpi = 300)
cat("Saved:", outfile, "\n")

#Mutant prevalence pᵢ (posterior mean)
#CBL (11:119284966 T>A) prevalence across germ cell stages
#Misclassification-corrected posterior mean ± 95% credible interval

# ---- 7) Export tidy results table --------------------------------------------
out <- summ %>%
	select(stage, k, n, s, fpr, p_mean, p_low, p_high) %>%
	mutate(across(c(p_mean,p_low,p_high), ~ round(., 4)))
print(out, n = Inf)
readr::write_csv(out, file.path(wd.de.plots, paste0(GENE, "_prevalence_by_stage_misclass_corrected.csv")))

# -----------------------------------------------------------------------------
# Bin schemes
# -----------------------------------------------------------------------------
# 1) 3 bins: pre-meiotic, meiotic, post-meiotic
bins_3 <- list(
	`Pre-meiotic`   = c("Stage 0","Stage 0A","Stage 0B","Stage 1"),
	`Meiotic`       = c("Stage 2","Stage 3","Leptotene","Zygotene","Pachytene","Diplotene","Meiotic division"),
	`Post-meiotic`  = c("Early spermatid","Late spermatid")
)

# 2) 4 bins by dominant phase (no overlaps)
bins_4 <- list(
	`G0/G1 (pre)`   = c("Stage 0","Stage 0A","Stage 0B","Stage 1"),
	`S`             = c("Stage 2","Stage 3","Leptotene"),
	`G2/M`          = c("Zygotene","Pachytene","Diplotene","Meiotic division"),
	`G0 (post)`     = c("Early spermatid","Late spermatid")
)

# ---------- Stage-specific sensitivity from coverage ----------
# all cells with stage
all_cells <- FetchData(so.st, vars = "cell_type") %>%
	tibble::rownames_to_column("CB") %>%
	transmute(CB, stage = factor(str_squish(as.character(cell_type)), levels = stage_order)) %>%
	filter(!is.na(stage))

# coverage at locus
cov_tbl <- calls %>%
	mutate(total_reads = coalesce(ref_reads,0) + coalesce(alt_reads,0),
			 has_cov    = total_reads > 0) %>%
	select(CB, has_cov, total_reads)

cov_stage <- all_cells %>%
	left_join(cov_tbl, by = "CB") %>%
	mutate(has_cov = replace_na(has_cov, FALSE),
			 total_reads = replace_na(total_reads, 0)) %>%
	group_by(stage) %>%
	summarise(n_cells=n(),
				 cov_frac=mean(has_cov),
				 reads_mu=ifelse(any(has_cov), mean(total_reads[has_cov]), 0),
				 .groups="drop")

# Pick a robust anchor (prefer Pachytene with many cells)
anchor <- cov_stage |>
	dplyr::filter(stage == "Pachytene", n_cells >= 20) |>
	dplyr::slice(1)
if (nrow(anchor) == 0) {
	anchor <- cov_stage |>
		dplyr::filter(n_cells >= 20) |>
		dplyr::arrange(dplyr::desc(cov_frac)) |>
		dplyr::slice(1)
}
cov_ref <- anchor$cov_frac[1]
s_floor <- 0.20; s_max <- 0.99

s_by_stage <- cov_stage |>
	dplyr::mutate(s = pmin(s_max, pmax(s_floor, s_max * cov_frac / cov_ref))) |>
	dplyr::select(stage, s)

fpr_default <- 1e-3
fpr_by_stage <- tibble(stage=factor(stage_order,levels=stage_order), fpr=fpr_default)

# ---------- Helper: fit a bin scheme ----------
posterior_p_grid <- function(k, n, s, fpr, a=1, b=1, grid_len=5001) {
	p <- seq(0, 1, length.out = grid_len)
	q <- p * s + (1 - p) * fpr
	loglik   <- dbinom(k, size=n, prob=q, log=TRUE)
	logprior <- dbeta(p, shape1=a, shape2=b, log=TRUE)
	w <- exp(loglik + logprior - max(loglik + logprior))
	w <- w/sum(w)
	cdf <- cumsum(w)
	tibble(
		p_mean = sum(p*w),
		p_low  = p[which.max(cdf >= 0.025)],
		p_high = p[which.max(cdf >= 0.975)]
	)
}

# --- Robust fitter: avoid Hessian inversion; borrow CIs from bin posteriors ---
fit_bin_scheme <- function(k_n_by_stage, s_by_stage, fpr_by_stage, bins_named_list) {
	stage_order <- levels(k_n_by_stage$stage)
	
	# stage -> bin map
	map <- tibble(bin = names(bins_named_list),
					  stages = I(unname(bins_named_list))) |>
		tidyr::unnest_longer(stages, values_to = "stage") |>
		mutate(stage = factor(stage, levels = stage_order),
				 bin   = factor(bin, levels = names(bins_named_list)))
	
	df <- k_n_by_stage |>
		mutate(stage = factor(as.character(stage), levels = stage_order)) |>
		left_join(map, by = "stage") |>
		left_join(s_by_stage, by = "stage") |>
		left_join(fpr_by_stage, by = "stage") |>
		arrange(bin, stage)
	
	stopifnot(!any(is.na(df$bin)))
	
	# Neg log-likelihood for bin logits theta
	bins <- levels(df$bin)
	npar <- length(bins)
	nll <- function(theta) {
		p <- plogis(theta)  # length npar
		q <- p[as.integer(df$bin)] * df$s + (1 - p[as.integer(df$bin)]) * df$fpr
		-sum(dbinom(df$k, size = df$n, prob = q, log = TRUE))
	}
	
	opt <- optim(par = rep(qlogis(0.2), npar), fn = nll, method = "BFGS",
					 control = list(maxit = 1000))
	
	loglik <- -opt$value
	aic    <- 2*npar - 2*loglik
	
	# Posterior-by-bin for CIs (aggregating k,n and s,fpr with n-weights)
	posterior_p_grid <- function(k, n, s, fpr, a=1, b=1, grid_len=5001) {
		p <- seq(0, 1, length.out = grid_len)
		q <- p * s + (1 - p) * fpr
		loglik   <- dbinom(k, size = n, prob = q, log = TRUE)
		logprior <- dbeta(p, shape1 = a, shape2 = b, log = TRUE)
		w <- exp(loglik + logprior - max(loglik + logprior)); w <- w/sum(w)
		cdf <- cumsum(w)
		tibble(
			p_mean = sum(p*w),
			p_low  = p[which.max(cdf >= 0.025)],
			p_high = p[which.max(cdf >= 0.975)]
		)
	}
	
	# Don’t re-anchor in bin aggregation
	agg <- df |>
		group_by(bin) |>
		summarise(k_bin   = sum(k),
					 n_bin   = sum(n),
					 # DO NOT re-scale s within each bin; just average the already-anchored stage s
					 s_bin   = sum(n * s) / sum(n),          # or use sum(n*s)/sum(n) as you had
					 # If you want FPR to reflect information, weight it by n*s (optional, safer)
					 fpr_bin = sum(n * s * fpr) / sum(n * s),
					 .groups = "drop") |>
		rowwise() |>
		mutate(post = list(posterior_p_grid(k_bin, n_bin, s_bin, fpr_bin))) |>
		tidyr::unnest(post) |>
		ungroup()
	
	# Report MLE p plus posterior CIs (no Hessian)
	p_hat <- tibble(bin = bins, p = plogis(opt$par)) |>
		left_join(agg |> select(bin, p_low, p_high), by = "bin")
	
	list(df_stage = df, p_hat_mle = p_hat, agg_posterior = agg,
		  logLik = loglik, AIC = aic)
}

# (Re)define bins (as you had)
bins_1  <- list(`All stages` = stage_order)
bins_12 <- as.list(stage_order); names(bins_12) <- stage_order

fit1  <- fit_bin_scheme(k_n_by_stage, s_by_stage, fpr_by_stage, bins_1)
fit3  <- fit_bin_scheme(k_n_by_stage, s_by_stage, fpr_by_stage, bins_3)
fit4  <- fit_bin_scheme(k_n_by_stage, s_by_stage, fpr_by_stage, bins_4)
fit12 <- fit_bin_scheme(k_n_by_stage, s_by_stage, fpr_by_stage, bins_12)

res <- dplyr::tibble(
	model  = c("1-bin (null)", "3-bin (pre/meiotic/post)", "4-bin (G0→S→G2M→G0)", "12-bin (per-stage)"),
	k      = c(1, 3, 4, 12),
	logLik = c(fit1$logLik, fit3$logLik, fit4$logLik, fit12$logLik),
	AIC    = c(fit1$AIC,    fit3$AIC,    fit4$AIC,    fit12$AIC)
) |> arrange(AIC)
print(res)

# Monotone trend across bins (posterior means)
mono_test <- function(agg, order_labels) {
	x <- seq_along(order_labels); m <- match(agg$bin, order_labels)
	sp <- suppressWarnings(cor.test(agg$p_mean, x[m], method = "spearman", exact = FALSE))
	kd <- suppressWarnings(cor.test(agg$p_mean, x[m], method = "kendall"))
	c(Spearman_rho = unname(sp$estimate), Spearman_p = sp$p.value,
	  Kendall_tau  = unname(kd$estimate), Kendall_p  = kd$p.value)
}
mt3 <- mono_test(fit3$agg_posterior, names(bins_3))
mt4 <- mono_test(fit4$agg_posterior, names(bins_4))
cat(sprintf("\n3-bin trend (Pre→Meiotic→Post): rho=%.3f (p=%.3g)\n", mt3["Spearman_rho"], mt3["Spearman_p"]))
cat(sprintf("4-bin trend (G0/G1→S→G2/M→G0):  rho=%.3f (p=%.3g)\n", mt4["Spearman_rho"], mt4["Spearman_p"]))

# Plots
plot_bins <- function(agg, title) {
	agg <- agg |> mutate(bin_idx = as.integer(factor(bin, levels = unique(bin))))
	ggplot(agg, aes(x = bin, y = p_mean, group = 1)) +
		geom_ribbon(aes(ymin = p_low, ymax = p_high), alpha = 0.20) +
		geom_line(aes(x = bin_idx), linewidth = 0.9) +
		geom_point(size = 2) +
		scale_y_continuous("Cancer cell fraction (CCF)", labels = scales::percent_format(accuracy = 1),
								 limits = c(0, NA)) +
		labs(x = "Cell cycle phase", title = title,
			  subtitle = "") +
		theme_classic(base_size = 12) + 
		theme(axis.text.x = element_text(angle = 35, hjust = 1),
				plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5),
				axis.title = ggplot2::element_text(size = 14),
				axis.text  = ggplot2::element_text(size = 12))
}
#Posterior mean ± 95% credible interval (misclassification-corrected)

p3 <- print(plot_bins(fit3$agg_posterior, "CBL (11:119284966 T>A)"))
#3-bin: Pre-meiotic → Meiotic → Post-meiotic
outfile <- file.path(wd.de.plots, paste0(GENE, "_prevalence_by_stage_misclass_corrected_3-bin.png"))
ggsave(outfile, plot = p3, width = 6, height = 6, dpi = 300)
cat("Saved:", outfile, "\n")

p3 <- p3 +
	geom_hline(yintercept = p_overall_adj, linetype = "dashed", linewidth = 0.6) +
	annotate("text", x = 1, y = p_overall_adj, vjust = -0.6, hjust = 0,
				label = paste0("Overall (adj) = ",
									scales::percent(p_overall_adj, 1),
									" [", scales::percent(ci_overall[1], 1), ", ",
									scales::percent(ci_overall[2], 1), "]"),
				size = 3.5)
outfile <- file.path(wd.de.plots, paste0(GENE, "_prevalence_by_stage_misclass_corrected_3-bin_adj.png"))
ggsave(outfile, plot = p3, width = 6, height = 6, dpi = 300)
cat("Saved:", outfile, "\n")

p4 <- print(plot_bins(fit4$agg_posterior, "CBL (11:119284966 T>A)"))
#4-bin: G0/G1 (pre) → S → G2/M → G0 (post)
outfile <- file.path(wd.de.plots, paste0(GENE, "_prevalence_by_stage_misclass_corrected_4-bin.png"))
ggsave(outfile, plot = p4, width = 6, height = 6, dpi = 300)
cat("Saved:", outfile, "\n")

p4 <- p4 +
	geom_hline(yintercept = p_overall_adj, linetype = "dashed", linewidth = 0.6) +
	annotate("text", x = 1, y = p_overall_adj, vjust = -0.6, hjust = 0,
				label = paste0("Overall (adj) = ",
									scales::percent(p_overall_adj, 1),
									" [", scales::percent(ci_overall[1], 1), ", ",
									scales::percent(ci_overall[2], 1), "]"),
				size = 3.5)
outfile <- file.path(wd.de.plots, paste0(GENE, "_prevalence_by_stage_misclass_corrected_4-bin_adj.png"))
ggsave(outfile, plot = p4, width = 6, height = 6, dpi = 300)
cat("Saved:", outfile, "\n")

# Useful tables
cat("\nBin prevalence (posterior) — 3-bin:\n")
print(fit3$agg_posterior |> mutate(across(starts_with("p_"), ~round(., 4))))
readr::write_csv(fit3$agg_posterior, file.path(wd.de.plots, paste0(GENE, "_prevalence_by_stage_misclass_corrected_3-bin.csv")))

cat("\nBin prevalence (posterior) — 4-bin:\n")
print(fit4$agg_posterior |> mutate(across(starts_with("p_"), ~round(., 4))))
readr::write_csv(fit4$agg_posterior, file.path(wd.de.plots, paste0(GENE, "_prevalence_by_stage_misclass_corrected_4-bin.csv")))

save(calls, calls2, k_n_by_stage, fit1, fit3, fit4, fit12, res, file=file.path(wd.de.plots, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_23_res=0.5_-C0_ordered_annotated_meta_cells_only_", GENE, "_CCF.RData")))













# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
library(dplyr)

# 1) Build the 12-type CCF vector from your `out`
ccf_one <- out %>%
	transmute(cell_type = as.character(stage),
				 CCF = p_mean, CCF_lo = p_low, CCF_hi = p_high)

allowed_ct <- ccf_one$cell_type   # the 12 stages with CCFs

# Predict composition only for those 12 types (re-normalize to sum=1)
mean_age <- mean(df$age, na.rm=TRUE)
sd_age   <- sd(df$age,  na.rm=TRUE)
sc       <- function(a) (a - mean_age) / sd_age
invlogit <- function(x) 1/(1+exp(-x))

pred_comp12 <- function(a_years) {
	pr <- sapply(allowed_ct, function(ct) {
		m <- mods[[ct]]; if (is.null(m)) return(NA_real_)
		predict(m, newdata = data.frame(age_c = sc(a_years)),
				  type = "response", re.form = NA)
	})
	pr <- pr / sum(pr, na.rm = TRUE)     # normalize across the 12 only
	tibble(cell_type = allowed_ct, prop = as.numeric(pr))
}

# Project overall CCF at age 60 → 70 using this donor’s within-type CCFs
proj_CCF12 <- function(a_years) {
	pred_comp12(a_years) %>%
		inner_join(ccf_one, by = "cell_type") %>%
		summarise(projected_CCF = sum(prop * CCF, na.rm = TRUE), .groups="drop") %>%
		pull(projected_CCF)
}

proj60 <- proj_CCF12(60)
proj70 <- proj_CCF12(70)
delta_per_decade <- proj70 - proj60
delta_per_decade

# -----------------------------------------------------------------------------
# A) Per-stage contributions (deterministic point estimates)
# -----------------------------------------------------------------------------
library(dplyr)
# 1) Get predicted compositions at 60 and 70 (population-average)
comp60 <- pred_comp12(60)  # tibble: cell_type, prop
comp70 <- pred_comp12(70)

# 2) Join and compute per-stage contributions
drivers <- comp70 %>%
	rename(prop70 = prop) %>%
	inner_join(rename(comp60, prop60 = prop), by = "cell_type") %>%
	inner_join(ccf_one, by = "cell_type") %>%
	mutate(
		delta_prop      = prop70 - prop60,             # change in stage fraction (absolute)
		contribution    = delta_prop * CCF,            # adds up to ΔCCF/decade
		delta_prop_pp   = 100 * delta_prop,            # percentage points of composition
		contribution_pp = 100 * contribution           # percentage points of CCF
	) %>%
	arrange(desc(contribution))

# 3) Check that the sum equals your ΔCCF/decade
sum_contrib <- sum(drivers$contribution, na.rm = TRUE)
# Should match your delta_per_decade ~ 0.003880074

# 4) Add each stage’s share of the total change
drivers <- drivers %>%
	mutate(share_pct = 100 * contribution / sum_contrib)

drivers

# -----------------------------------------------------------------------------
# Quick bar plot (who pushes the increase?)
# -----------------------------------------------------------------------------
library(ggplot2)

ggplot(drivers, aes(x = reorder(cell_type, contribution_pp), y = contribution_pp)) +
	geom_hline(yintercept = 0, linetype = 2) +
	geom_col() +
	coord_flip() +
	labs(x = NULL, y = "Contribution to ΔCCF per decade (pp)",
		  title = "Stage drivers of composition-only ΔCCF (age 60→70)") +
	theme_classic()

# -----------------------------------------------------------------------------
# the projected overall CCF at 60 and 70,
# the relative % change per decade,
# and a small table of the top 5 positive/negative stage contributions (in pp).
# -----------------------------------------------------------------------------
library(dplyr)

## 1) Project overall CCF at 60 and 70
comp60 <- pred_comp12(60)   # tibble: cell_type, prop
comp70 <- pred_comp12(70)

tbl60 <- comp60 %>% inner_join(ccf_one, by = "cell_type")
tbl70 <- comp70 %>% inner_join(ccf_one, by = "cell_type")

proj60 <- sum(tbl60$prop * tbl60$CCF, na.rm = TRUE)
proj70 <- sum(tbl70$prop * tbl70$CCF, na.rm = TRUE)

delta          <- proj70 - proj60                   # absolute change (fraction units)
abs_pp         <- 100 * delta                       # percentage points per decade
rel_pct        <- 100 * delta / proj60              # relative % per decade

cat(sprintf("Projected overall CCF at 60: %.4f  (%.2f%%)\n", proj60, 100*proj60))
cat(sprintf("Projected overall CCF at 70: %.4f  (%.2f%%)\n", proj70, 100*proj70))
cat(sprintf("Change per decade: %+0.4f ( %+0.2f pp ); relative: %+0.2f%% per decade\n\n",
				delta, abs_pp, rel_pct))

## 2) Per-stage driver decomposition (sums to delta)
drivers <- comp70 %>%
	rename(prop70 = prop) %>%
	inner_join(rename(comp60, prop60 = prop), by = "cell_type") %>%
	inner_join(ccf_one, by = "cell_type") %>%
	mutate(
		delta_prop      = prop70 - prop60,      # change in composition (absolute)
		contribution    = delta_prop * CCF,     # change in overall CCF due to this stage
		delta_prop_pp   = 100 * delta_prop,     # pp change in composition
		contribution_pp = 100 * contribution,   # pp contribution to overall CCF
		share_pct       = 100 * contribution / sum(contribution, na.rm = TRUE)
	) %>%
	arrange(desc(contribution_pp))

# Take top 5 positive and top 5 negative contributors
pos5 <- drivers %>% slice_max(contribution_pp, n = 5, with_ties = FALSE)
neg5 <- drivers %>% slice_min(contribution_pp, n = 5, with_ties = FALSE)

report <- bind_rows(pos5, neg5) %>%
	transmute(
		cell_type,
		`Δ composition (pp/dec)` = round(delta_prop_pp, 3),
		`CCF (donor)`            = round(CCF, 3),
		`Contribution (pp/dec)`  = round(contribution_pp, 3),
		`Share of total (%)`     = round(share_pct, 1)
	) %>%
	arrange(desc(`Contribution (pp/dec)`))

print(report, n = nrow(report))

# Sanity check: contributions sum to the overall change
cat(sprintf("\nCheck: sum of contributions = %+0.4f ( %+0.2f pp )\n",
				sum(drivers$contribution, na.rm = TRUE),
				sum(drivers$contribution_pp, na.rm = TRUE)))












# ---------- 2) Build S-fraction per stage from your phase_proportions ----------
phases_wide <- phase_proportions %>%
	mutate(cell_type = str_squish(as.character(cell_type)),
			 phase     = factor(as.character(phase), levels = c("G1","S","G2M"))) %>%
	filter(cell_type %in% stage_order, !is.na(phase)) %>%
	select(cell_type, phase, proportion) %>%
	pivot_wider(names_from = phase, values_from = proportion) %>%
	mutate(across(c(G1,S,G2M), ~ replace_na(., 0))) %>%
	rowwise() %>%
	mutate(.sum = sum(c_across(c(G1,S,G2M))),
			 G1 = ifelse(.sum > 0, G1/.sum, NA_real_),
			 S  = ifelse(.sum > 0, S /.sum, NA_real_),
			 G2M= ifelse(.sum > 0, G2M/.sum, NA_real_)) %>%
	ungroup() %>%
	select(-.sum) %>%
	mutate(Cycling = S + G2M,
			 stage = factor(cell_type, levels = stage_order)) %>%
	arrange(stage)

# ---------- 3) Misclassification-corrected prevalence (build 'summ') ----------
s_default   <- 0.65   # sensitivity Pr(call | truly mutant) — replace with stage-specific if you have them
fpr_default <- 1e-3   # false positive rate Pr(call | WT)

posterior_p_grid <- function(k, n, s, fpr, a=1, b=1, grid_len=5001) {
	p <- seq(0, 1, length.out = grid_len)
	q <- p * s + (1 - p) * fpr
	loglik <- dbinom(k, n, q, log = TRUE)
	logpost <- loglik + dbeta(p, a, b, log = TRUE)
	m <- max(logpost); w <- exp(logpost - m); w <- w / sum(w)
	cdf <- cumsum(w)
	tibble(
		p_mean = sum(p*w),
		p_low  = p[which.max(cdf >= 0.025)],
		p_high = p[which.max(cdf >= 0.975)],
		p_mode = p[which.max(w)]
	)
}

df_prev <- k_n_by_stage %>%
	mutate(s = s_default, fpr = fpr_default)

summ <- df_prev %>%
	rowwise() %>%
	mutate(res = list(posterior_p_grid(k, n, s, fpr))) %>%
	unnest(res) %>%
	ungroup() %>%
	mutate(stage_idx = as.integer(stage),
			 logit_p   = qlogis(pmin(pmax(p_mean, 1e-6), 1-1e-6)))

# ---------- 4) Join prevalence with G1/S/G2M and run correlations ----------
dat <- summ %>%
	select(stage, stage_idx, p_mean, p_low, p_high) %>%
	inner_join(phases_wide %>% select(stage, G1, S, G2M, Cycling), by = "stage") %>%
	arrange(stage)

# Spearman correlations (raw)
vars <- c("G1","S","G2M","Cycling")
corrs <- lapply(vars, function(v) {
	ct <- suppressWarnings(cor.test(dat$p_mean, dat[[v]], method = "spearman", exact = FALSE))
	tibble(var = v, rho = unname(ct$estimate), p = ct$p.value, n = nrow(dat))
}) %>% bind_rows() %>% mutate(q = p.adjust(p, "BH"))
print(corrs)

# Residualized (beyond stage order trend)
residual_cor <- function(y, x, idx) {
	ry <- lm(y ~ idx)$residuals
	rx <- lm(x ~ idx)$residuals
	suppressWarnings(cor.test(ry, rx, method = "spearman", exact = FALSE))
}
pcorrs <- lapply(vars, function(v) {
	ct <- residual_cor(dat$p_mean, dat[[v]], dat$stage_idx)
	tibble(var = v, rho_resid = unname(ct$estimate), p_resid = ct$p.value)
}) %>% bind_rows() %>% mutate(q_resid = p.adjust(p_resid, "BH"))
print(pcorrs)

# ---------- 5) Quick plot ----------
lab_pct <- function(x) percent(x, accuracy = 1)
long <- dat %>%
	select(stage, p_mean, G1, S, G2M, Cycling) %>%
	pivot_longer(cols = c(G1, S, G2M, Cycling), names_to = "phase_metric", values_to = "prop")

gg1 <- ggplot(long, aes(x = prop, y = p_mean)) +
	geom_point(size = 2) +
	geom_smooth(method = "lm", se = FALSE, linetype = 2) +
	facet_wrap(~ phase_metric, nrow = 1, scales = "free_x") +
	scale_x_continuous(labels = lab_pct, limits = c(0, 1)) +
	scale_y_continuous(labels = lab_pct, limits = c(0, NA)) +
	labs(x = "Phase proportion", y = "Mutant prevalence p\u1D62",
		  title = "CBL mutant prevalence vs phase composition across stages",
		  subtitle = "Spearman correlations shown in console; Cycling = S + G2M") +
	theme_classic(base_size = 12)
print(gg1)

# -----------------------------------------------------------------------------
# 1) Is CBL expression lower in S/G2M (by stage × phase)?
# -----------------------------------------------------------------------------
library(dplyr); library(tidyr); library(stringr); library(Seurat)

DefaultAssay(so.integrated) <- "RNA"  # or "RNA" if that's your working assay
# Per-cell CBL expression
cb_vec <- tryCatch(GetAssayData(so.integrated, slot="data")["CBL", ], error=function(e) NULL)
if (is.null(cb_vec)) stop("Gene 'CBL' not found in current assay; check feature name/case or switch assay.")

md <- FetchData(so.integrated, vars = c("cell_type","Phase")) %>%
	tibble::rownames_to_column("CB") %>%
	mutate(cell_type = str_squish(cell_type),
			 Phase = factor(Phase, levels=c("G1","S","G2M"))) %>%
	mutate(CBL = as.numeric(cb_vec[CB]),
			 is_expr = CBL > 0)

# summarize by stage × phase
cb_stage_phase <- md %>%
	filter(cell_type %in% stage_order, !is.na(Phase)) %>%
	group_by(cell_type, Phase) %>%
	summarise(
		n = n(),
		det_frac = mean(is_expr),
		mean_expr = mean(CBL),       # log-normalized mean
		med_expr  = median(CBL),
		.groups = "drop"
	) %>%
	mutate(stage = factor(cell_type, levels = stage_order)) %>%
	arrange(stage, Phase)

print(cb_stage_phase, n=36)

# -----------------------------------------------------------------------------
# tests detection bias (coverage + mutant call) ~ CBL expression, total UMIs, Phase,
# recomputes local density and next-stage distance,
# fits an adjusted regression: distance ~ %S + density + UMIs + CBL + slide,
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
	library(dplyr); library(tidyr); library(stringr); library(purrr)
	library(Seurat); library(FNN); library(lme4); library(lmerTest); library(broom)
})

# -------- 1) Per-cell metadata table ------------------------------------------
# Fetch stage, phase, counts and coordinates
vars <- c("cell_type","Phase","nCount_RNA","sample.id")
md <- FetchData(so.integrated, vars = vars) %>%
	tibble::rownames_to_column("CB") %>%
	mutate(cell_type = str_squish(as.character(cell_type)),
			 Phase     = factor(as.character(Phase), levels = c("G1","S","G2M")),
			 stage     = factor(cell_type, levels = stage_order))

# Try UMAP; if missing, fall back to plot_x/plot_y
umap_ok <- "umap" %in% names(so.integrated@reductions)
if (umap_ok) {
	md <- md %>% bind_cols(
		FetchData(so.st, vars = c("UMAP_1","UMAP_2"))
	)
} else {
	# If you stored custom coords:
	if (all(c("plot_x","plot_y") %in% colnames(so.integrated@meta.data))) {
		xy <- FetchData(so.integrated, vars = c("plot_x","plot_y"))
		colnames(xy) <- c("UMAP_1","UMAP_2")
		md <- md %>% bind_cols(xy)
	} else {
		stop("No UMAP_1/UMAP_2 or plot_x/plot_y found.")
	}
}

# CBL expression (log-normalised) from active assay
DefaultAssay(so.integrated) <- DefaultAssay(so.integrated) # keep current
cb_vec <- tryCatch(GetAssayData(so.integrated, slot = "data")["CBL", ], error=function(e) NULL)
if (is.null(cb_vec)) stop("Feature 'CBL' not found in current assay; switch assay or check name.")
md$CBL_expr <- as.numeric(cb_vec[md$CB])

# Total UMIs (adjust if you prefer nCount_SCT)
md$log_umis <- log1p(md$nCount_RNA)

# -------- 2) Join mutation calls & coverage -----------------------------------
# calls has union of REF/ALT per CB; build coverage and mutant flags
stopifnot(all(c("CB","ref_reads","alt_reads") %in% colnames(calls)))
calls_cov <- calls %>%
	mutate(total_reads = coalesce(ref_reads,0) + coalesce(alt_reads,0),
			 has_cov    = total_reads > 0,
			 is_mut1    = !is.na(alt_reads) & alt_reads >= 1,  # ≥1 ALT read
			 is_mut2    = !is.na(alt_reads) & alt_reads >= 2)  # ≥2 ALT reads

dat <- md %>%
	left_join(calls_cov %>% select(CB, has_cov, total_reads, is_mut1, is_mut2), by = "CB") %>%
	mutate(across(c(has_cov, is_mut1, is_mut2), ~ replace_na(., FALSE)),
			 total_reads = replace_na(total_reads, 0)) %>%
	filter(!is.na(stage), !is.na(Phase))

# -------- 3) CHECK 1: Detection bias models -----------------------------------
# (A) Coverage presence ~ CBL + UMIs + Phase + Stage
m_cov <- glm(has_cov ~ scale(CBL_expr) + scale(log_umis) + Phase + stage,
				 data = dat, family = binomial())
cat("\n--- Coverage model (has_cov) ---\n")
print(broom::tidy(m_cov) %>% mutate(or = exp(estimate)))

# (B) Mutation call among covered cells ~ CBL + UMIs + Phase + Stage
covered <- dat %>% filter(has_cov)
m_mut1 <- glm(is_mut1 ~ scale(CBL_expr) + scale(log_umis) + Phase + stage,
				  data = covered, family = binomial())
m_mut2 <- glm(is_mut2 ~ scale(CBL_expr) + scale(log_umis) + Phase + stage,
				  data = covered, family = binomial())
cat("\n--- Mutation call model (≥1 ALT) ---\n")
print(broom::tidy(m_mut1) %>% mutate(or = exp(estimate)))
cat("\n--- Mutation call model (≥2 ALT) ---\n")
print(broom::tidy(m_mut2) %>% mutate(or = exp(estimate)))

# Interpretation guide:
#  - If PhaseS or PhaseG2M are positive/NS after adjusting for CBL_expr and UMIs,
#    phase-by-coverage bias is not driving calls; if negative, it suggests under-calling in S/G2M.

# -------- 4) Stage %S and other stage covariates ------------------------------
# %S by stage (from your table)
S_by_stage <- phase_proportions %>%
	filter(phase %in% c("G1","S","G2M")) %>%
	group_by(cell_type) %>%
	mutate(prop = proportion / sum(proportion)) %>%  # re-normalize defensively
	ungroup() %>%
	filter(phase == "S") %>%
	transmute(stage = factor(cell_type, levels = stage_order),
				 S_frac = prop) %>%
	distinct()

# Stage sizes
n_by_stage <- md %>% count(stage, name = "n_stage")

stage_covars <- S_by_stage %>%
	full_join(n_by_stage, by = "stage")

# -------- 5) Local density and next-stage distance ----------------------------
# kNN density (inverse mean distance to k nearest neighbours among all germ stages)
k <- 20
coords <- dat %>% select(CB, stage, UMAP_1, UMAP_2)
mat <- as.matrix(coords %>% select(UMAP_1, UMAP_2))
nn <- FNN::get.knn(mat, k = min(k, nrow(mat)-1))

inv_mean_d <- 1 / rowMeans(nn$nn.dist)
dens_tbl <- tibble(CB = coords$CB, local_density = as.numeric(scale(inv_mean_d)))

# Next-stage mapping
stage_next <- setNames(c(stage_order[-1], NA), stage_order)
dat2 <- dat %>%
	left_join(dens_tbl, by = "CB") %>%
	mutate(next_stage = stage_next[as.character(stage)])

# Fast nearest distance to any cell in the next stage
get_next_dist <- function(df) {
	out <- rep(NA_real_, nrow(df))
	for (s in stage_order) {
		idx <- which(df$stage == s)
		ns  <- stage_next[s]
		if (is.na(ns)) next
		jdx <- which(df$stage == ns)
		if (length(idx) == 0 || length(jdx) == 0) next
		m1 <- as.matrix(df[idx, c("UMAP_1","UMAP_2")])
		m2 <- as.matrix(df[jdx, c("UMAP_1","UMAP_2")])
		d  <- FNN::get.knnx(m2, m1, k = 1)$nn.dist[,1]
		out[idx] <- d
	}
	out
}
dat2$next_dist <- get_next_dist(dat2)

# Bring in stage-level covariates
dat2 <- dat2 %>%
	left_join(stage_covars, by = "stage") %>%
	# n_next: size of the next stage
	left_join(n_by_stage %>% rename(next_stage = stage, n_next = n_stage),
				 by = "next_stage") %>%
	mutate(log_n_stage = log1p(n_stage),
			 log_n_next  = log1p(n_next))

# -------- 6) CHECK 2+3: Adjusted regression of distance -----------------------
# distance ~ %S + density + UMIs + CBL + (1|stage) + (1|sample.id)
dist_df <- dat2 %>%
	filter(!is.na(next_dist), is.finite(next_dist)) %>%
	mutate(S_frac = replace_na(S_frac, 0))

m_dist <- lmer(
	next_dist ~ S_frac + local_density + scale(log_umis) + scale(CBL_expr) +
		log_n_stage + log_n_next + (1|stage) + (1|sample.id),
	data = dist_df, REML = FALSE
)
cat("\n--- Adjusted distance model ---\n")
print(summary(m_dist))

# Optional: also test S vs G2M separately at cell level (encode cell-cycle as dummies)
# Here we add cell-level Phase to ensure %S isn’t just reflecting micro-phase composition
m_dist2 <- lmer(
	next_dist ~ S_frac + local_density + scale(log_umis) + scale(CBL_expr) +
		log_n_stage + log_n_next + Phase + (1|stage) + (1|sample.id),
	data = dist_df, REML = FALSE
)
cat("\n--- Adjusted distance model (with cell-level Phase) ---\n")
print(summary(m_dist2))

# Quick tidy tables for reporting
tab_cov <- broom.mixed::tidy(m_cov)
tab_mut1 <- broom.mixed::tidy(m_mut1)
tab_mut2 <- broom.mixed::tidy(m_mut2)
tab_dist <- broom.mixed::tidy(m_dist, effects = "fixed")
tab_dist2<- broom.mixed::tidy(m_dist2, effects = "fixed")

list(
	coverage_model = tab_cov,
	mutation_model_ge1 = tab_mut1,
	mutation_model_ge2 = tab_mut2,
	distance_model = tab_dist,
	distance_model_with_phase = tab_dist2
) -> results_tables

# Example: write to CSVs
# readr::write_csv(tab_cov,  file.path(wd.de.plots, "model_coverage.csv"))
# readr::write_csv(tab_mut1, file.path(wd.de.plots, "model_mut_ge1.csv"))
# readr::write_csv(tab_mut2, file.path(wd.de.plots, "model_mut_ge2.csv"))
# readr::write_csv(tab_dist, file.path(wd.de.plots, "model_distance.csv"))
# readr::write_csv(tab_dist2,file.path(wd.de.plots, "model_distance_with_phase.csv"))














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
outfile <- file.path(wd.de.plots, "origUM", "Spatial_ST_META_ONLY_origcolors_short_labels_title.png")
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

save(meta_idx, spatial_data, spatial_coords, cluster_centers, avg_dist, dist_map, meta_cells_only, phase_proportions, file=file.path(wd.de.plots, "origUM", paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_23_res=0.5_-C0_ordered_annotated_meta_cells_only_", GENE, ".RData")))

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

cor_results_obs <- phase_dist %>%
	group_by(phase) %>%
	summarise(
		rho = cor(dist_to_next, proportion, use = "complete.obs"),
		pval = cor.test(dist_to_next, proportion)$p.value,
		.groups = "drop"
	)

print(cor_results_obs)
save(meta_idx, spatial_data, spatial_coords, cluster_centers, dist_map, meta_cells_only, phase_proportions, dist_summary, phase_dist, cor_results_obs, ct_order_master, short_labels, ct_levels, phase_colors, file=file.path(wd.de.plots, "origUM", paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_23_res=0.5_-C0_ordered_annotated_meta_cells_only_", GENE, ".RData")))

# Ensure your palette exists
phase_colors <- c("G1" = "#619CFF", "S" = "#F8766D", "G2M" = "#00BA38")

# Helper to make & save one phase plot to a specific directory with a clear filename
save_phase_plot <- function(df, phase_text, phase_label, color_hex, outdir,
									 prefix = "Spatial_ST_META_CBL") {
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

# -----------------------------------------------------------------------------
# Cell cycle analysis (PD53626b_ST1)
# Last Modified: 08/10/25
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# To calculate the cell cycle phase proportions (G1/S/G2M) in this sample
# Last Modified: 08/10/25
# -----------------------------------------------------------------------------
so.integrated <- so.integrated.st

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
ggsave(file.path(wd.de.plots, paste0("Cycle_SSC_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_6x6_PD53626b_ST1.png")), 
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
