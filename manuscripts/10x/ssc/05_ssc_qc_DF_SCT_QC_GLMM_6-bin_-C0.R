# =============================================================================
# Name         : 
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 11/09/25
# =============================================================================
wd.src <- "/nfs/users/nfs_t/ty2/dev/R"            ## ty2@farm
#wd.src <- "/Users/ty2/Work/dev/R"                ## ty2@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Commons.R", "SingleCellTranscriptomics.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
wd <- "/lustre/scratch125/casm/staging/team294/ty2"   ## ty2@farm22
#wd <- "/Users/ty2/Work/sanger/ty2"                   ## ty2@localhost
BASE <- "SSC"
base <- tolower(BASE)

wd.rna <- file.path(wd, BASE)
wd.rna.raw <- file.path(wd.rna, "ngs", "10x")

wd.de    <- file.path(wd.rna, "analysis", paste0(base, ""))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

load(file.path(wd.de.data, paste0("ssc_filtered_normalised_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_ordered_annotated_monocle3+phase.RData")))

so.integrated$orig.ident <- ifelse(
   so.integrated$orig.ident == "SeuratProject",
   "VL00297",
   so.integrated$orig.ident
)
unique(so.integrated$orig.ident)
so.integrated$patient_id <- str_extract(so.integrated$orig.ident, "^[^_]+")

so.integrated$cell_type <- so.integrated$idents

so.integrated@meta.data <- so.integrated@meta.data %>%
   mutate(age_group = cut(
      age,
      breaks = c(20, 30, 40, 50, 60, 70, 80),
      labels = c("20-30", "30-40", "40-50", "50-60", "60-70", "70-80"),
      right = FALSE
   ))

save(cds, so.integrated, phase_colors, plot_colors, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_23_res=0.5_-C0_ordered_annotated_monocle3+phase.RData")))

# -----------------------------------------------------------------------------
# Proportion of cell type ~ age + (1 | patient_id), per cell type
# -----------------------------------------------------------------------------
library(lme4)
library(dplyr)
library(stringr)
library(purrr)
library(ggplot2)

# 1. Total cells per sample
total_cells <- so.integrated@meta.data %>%
	group_by(patient_id) %>%
	summarise(total_cells = n(), .groups = "drop")

# 2. Cell counts per sample × cell type
cell_counts <- so.integrated@meta.data %>%
	group_by(patient_id, cell_type = cell_type, age) %>%
	summarise(n_cells = n(), .groups = "drop")

# 3. Merge total cells into cell_counts
cell_counts <- left_join(cell_counts, total_cells, by = "patient_id")

# 4. Get list of unique cell types
cell_types <- unique(cell_counts$cell_type)

# -----------------------------------------------------------------------------
# Binomial GLMM (better than Poisson here)
# -----------------------------------------------------------------------------
library(lme4)
library(broom.mixed)
library(dplyr)

# per-sample counts (you already have cell_counts with n_cells and total_cells)
df <- cell_counts %>%
	mutate(success   = n_cells,
			 failure   = pmax(total_cells - n_cells, 0),
			 age_c     = scale(age, center = TRUE, scale = TRUE))  # helps convergence

fit_binom <- function(ct) {
	dd <- filter(df, cell_type == ct)
	if (n_distinct(dd$success) < 2) return(NULL)
	glmer(cbind(success, failure) ~ age_c + (1|patient_id),
			data = dd,
			family = binomial(link = "logit"),
			control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
}

mods <- setNames(lapply(unique(df$cell_type), fit_binom), unique(df$cell_type))
res  <- bind_rows(lapply(names(mods), function(ct) {
	m <- mods[[ct]]; if (is.null(m)) return(NULL)
	tt <- tidy(m, effects = "fixed") %>% filter(term == "age_c")
	tt %>% transmute(cell_type = ct,
						  logit_per_SD = estimate, se = std.error, p = p.value)
}))

# FDR across cell types
res <- res %>% mutate(q = p.adjust(p, method = "BH"))

# Report odds ratios per year (or per decade) instead of per SD (helps biology readers)
res <- res %>%
	mutate(OR_per_SD = exp(logit_per_SD),
			 CI_lo = exp(logit_per_SD - 1.96*se),
			 CI_hi = exp(logit_per_SD + 1.96*se))

sd_age <- sd(df$age, na.rm = TRUE)
res <- res %>%
	mutate(
		OR_per_year   = exp(logit_per_SD / sd_age),
		OR_per_decade = exp(logit_per_SD * (10 / sd_age))
	)

# Order stages biologically and make a quick forest plot (effects ±95% CI):
stage_order <- c("Stage 0","Stage 0A","Stage 0B","Stage 1","Stage 2","Stage 3",
					  "Leptotene","Zygotene","Pachytene","Diplotene",
					  "Meiotic division","Early spermatid","Late spermatid")

res$cell_type <- factor(res$cell_type, levels = rev(stage_order))

library(ggplot2)
ggplot(res, aes(x = cell_type, y = OR_per_SD)) +
	geom_hline(yintercept = 1, linetype = 2) +
	geom_point() +
	geom_errorbar(aes(ymin = CI_lo, ymax = CI_hi), width = 0.15) +
	coord_flip() +
	labs(y = "Odds ratio per 1 SD age", x = NULL,
		  title = "Age effect on cell-type proportion (GLMM, per SD age)") +
	theme_classic()




# -----------------------------------------------------------------------------
# Proportion of cell type ~ age + (1 | patient_id), per cell type
# -----------------------------------------------------------------------------
library(lme4)
library(dplyr)
library(stringr)
library(purrr)
library(ggplot2)

# 1. Total cells per sample
total_cells <- so.integrated@meta.data %>%
   group_by(patient_id) %>%
   summarise(total_cells = n(), .groups = "drop")

# 2. Cell counts per sample × cell type
cell_counts <- so.integrated@meta.data %>%
   group_by(patient_id, cell_type = cell_type, age) %>%
   summarise(n_cells = n(), .groups = "drop")

# 3. Merge total cells into cell_counts
cell_counts <- left_join(cell_counts, total_cells, by = "patient_id")

# 4. Get list of unique cell types
cell_types <- unique(cell_counts$cell_type)

# 5. Fit a binomial GLMM per cell type
models <- map(cell_types, function(ct) {
   df_ct <- filter(cell_counts, cell_type == ct)
 
   if (n_distinct(df_ct$n_cells) > 1) {
      tryCatch({
         model <- glmer(n_cells ~ age + (1 | patient_id), offset = log(total_cells),
                        data = df_ct,
                        family = poisson())
         summary(model)
      }, error = function(e) {
         message(paste("Model failed for", ct, ":", e$message))
         NULL
      })
   } else {
      message(paste("Skipped", ct, "due to lack of variation."))
      NULL
   }
})

# 6. Name the models
names(models) <- cell_types

# 7. Extract coefficients
results_df <- map_dfr(cell_types, function(ct) {
   model <- models[[ct]]
   if (!is.null(model)) {
      coef_summary <- coef(summary(model))
      data.frame(
         cell_type = ct,
         estimate = coef_summary["age", "Estimate"],
         std_error = coef_summary["age", "Std. Error"],
         p_value = coef_summary["age", "Pr(>|z|)"]
      )
   }
})

# 8. Inspect results
print(results_df)

# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
# Add significance column based on p-value
results_df <- results_df %>%
   mutate(
      significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ ""
      ) 
   )

# Reorder cell types biologically
cell_order <- c(
   "Stage 0", "Stage 0A", "Stage 0B", "Stage 1", "Stage 2", "Stage 3",
   "Leptotene", "Zygotene", "Pachytene", "Diplotene", "Meiotic division",
   "Early spermatid", "Late spermatid"
)
results_df$cell_type <- factor(results_df$cell_type, levels = rev(cell_order))

# Create the plot
plot <- ggplot(results_df, aes(x = cell_type, y = estimate, fill = estimate)) +
   geom_bar(stat = "identity") +
   coord_flip() +
   scale_fill_gradient2(
      low = blue,
      mid = "white",
      high = red,
      midpoint = 0,
      name = "Age effect"
   ) +
   geom_text(aes(label = significance, 
                 hjust = ifelse(estimate < 0, 1.5, -0.5)),
             vjust = 0.5, color = "black", size = 6) +
   labs(
      title = "",
      x = "",
      y = "Age effect (GLMM estimate)"
   ) +
   theme_minimal() +
   theme(
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      axis.title.x = element_text(size = 12, margin = margin(t = 10)),
      plot.title = element_text(hjust = 0.5),
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 12)
   ) +
   scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

# Save to PDF
pdf(file = file.path(wd.de.plots, "barplot_age_effect_by_cell_type_GLMM_cell_counts_merged_age+total_cells_poisson.pdf"), width = 6, height = 6)
print(plot)
dev.off()

# -----------------------------------------------------------------------------
# Methods: Plot cell counts by patient ID
# Last Modified: 15/01/25
# -----------------------------------------------------------------------------
# Add patient_label and order by age
patient_age_map <- cell_counts %>%
   distinct(patient_id, age) %>%
   mutate(patient_label = paste0(patient_id, " (", age, ")"))

cell_counts <- left_join(cell_counts, patient_age_map, by = c("patient_id", "age"))
cell_counts$patient_label <- factor(cell_counts$patient_label, levels = patient_age_map %>% arrange(age) %>% pull(patient_label))

# Merge n_cells across samples with the same patient_id + cell_type
library(dplyr)
cell_counts_merged <- cell_counts %>%
   group_by(patient_id, cell_type, age, patient_label) %>%
   summarise(n_cells = sum(n_cells), .groups = "drop")

cell_counts_merged <- cell_counts_merged %>%
   group_by(patient_id) %>%
   mutate(total_cells = sum(n_cells)) %>%
   ungroup()

# Plot: One stacked bar per patient, colored by cell type
plot <- ggplot(cell_counts_merged, aes(x = patient_label, y = n_cells, fill = cell_type)) +
   geom_bar(stat = "identity", position = "stack", color = "black") +
   scale_fill_manual(values = plot_colors, name = "Cell type") +
   labs(title = "", x = "Patient (Age)", y = "Number of cells") +
   theme_minimal(base_size = 14) +
   theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11, margin = margin(t = -15)),
      axis.text.y = element_text(size = 11),
      axis.title.x = element_text(size = 13),
      axis.title.y = element_text(size = 13),
      plot.title = element_text(hjust = 0.5),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12)
   ) +
   guides(fill = guide_legend(ncol = 1))

# Save to PNG
png(file = file.path(wd.de.plots, "barplot_patient_cell_type_composition.png"),
    width = 5.3, height = 7, units = "in", res = 300)
print(plot)
dev.off()

# -----------------------------------------------------------------------------
# Proportion of cell type ~ age + (1 | patient_id), per cell type
# -----------------------------------------------------------------------------
totals <- cell_counts_merged %>%
   group_by(patient_id) %>%
   summarise(total_cells = sum(n_cells), .groups = "drop")

# Merge back into cell_counts_merged
cell_counts_binomial <- cell_counts_merged %>%
   left_join(totals, by = "patient_id") %>%
   mutate(other_cells = total_cells - n_cells)

cell_types <- unique(cell_counts_binomial$cell_type)

binom_models <- map(cell_types, function(ct) {
   df_ct <- filter(cell_counts_binomial, cell_type == ct)
 
   # Need at least 3 patients and variation in n_cells
   if (n_distinct(df_ct$patient_id) >= 3 && var(df_ct$n_cells) > 0) {
      tryCatch({
         model <- glmer(cbind(n_cells, other_cells) ~ age + (1 | patient_id),
                        data = df_ct,
                        family = binomial())
         model
      }, error = function(e) {
         message(paste("Model failed for", ct, ":", e$message))
         NULL
      })
   } else {
      message(paste("Skipped", ct, "due to insufficient variation."))
      NULL
   }
})
names(binom_models) <- cell_types

results_binom <- map_dfr(cell_types, function(ct) {
   model <- binom_models[[ct]]
   if (!is.null(model)) {
      coef_summary <- coef(summary(model))
      tibble(
         cell_type = ct,
         estimate = coef_summary["age", "Estimate"],
         std_error = coef_summary["age", "Std. Error"],
         p_value = coef_summary["age", "Pr(>|z|)"]
      )
   }
})

results_binom
results_df <- results_binom
