# =============================================================================
# Name         : 
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 19/08/25
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

load(file=file.path(wd.de.data, "ssc_filtered_normalised_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0.RData"))

# -----------------------------------------------------------------------------
# Di Persio et al; DotPlot (resolution = 0.25)
# -----------------------------------------------------------------------------
genes_of_interest <- c("DDX4", "MAGEA4", "DMRT1", "SOX4", "ID4", "FGFR3", "TCF3", "GFRA1", "NANOS2", "KIT", "MKI67", "NANOS3", "STRA8", "SYCP1", "SYCP3", "MLH3", "SPO11", "MEIOB", "SCML1", "TEX19", "DPH7", "DMC1", "LY6K", "SELENOT", "TDRG1", "PIWIL1", "POU5F2", "OVOL2", "CCDC112", "AURKA", "CCNA1", "C9orf116", "SLC26A3", "SIRPG", "TEX29", "TNP1", "PRM2", "PRM1", "VIM", "CITED1", "SOX9", "FATE1", "HSD17B3", "STAR", "INSL3", "CLEC3B", "CFD", "MYH11", "ACTA2", "PECAM1", "VWF", "CD68", "LYZ", "C1QA", "CD14")

DefaultAssay(so.integrated) <- "SCT"
dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels

pdf(file = file.path(wd.de.plots, paste0("Di Persio_DotPlot_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0.pdf")), width = 14, height = 6)
print(dot_plot)
dev.off()

# Apply this custom order
# Note: Only change levels if 'Idents' are factor. If they are characters, convert them to factors first.
new_order <- c("5", "8", "19", "1", "16", "15", "3", "6", "2", "10", "12", "18", "14", "9", "13", "11", "20", "4", "23", "22", "7", "21", "0", "17")
#new_order <- rev(new_order)
Idents(so.integrated) <- factor(Idents(so.integrated), levels = new_order)

dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
   scale_y_discrete(limits = new_order) # Ensure the new order is used in plotting

pdf(file = file.path(wd.de.plots, paste0("Di Persio_DotPlot_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_ordered.pdf")), width = 14, height = 6)
print(dot_plot)
dev.off()

##
dim_plot <- DimPlot(so.integrated, label = TRUE) + NoLegend() +
	  labs(x = "UMAP 1", y = "UMAP 2") +  # Set x-axis and y-axis labels
	  theme(
		    plot.title = element_blank(),  # Remove title
		    axis.title = element_text(size = 18),  # Increase axis title size
		    axis.text = element_text(size = 16),  # Increase axis tick label size
	  )
ggsave(file.path(wd.de.plots, paste0("DimPlot_SSC_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_new-color_6*6.png")), plot = dim_plot, width = 6, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# Di Persio et al; DotPlot (resolution = 0.4; Six SPG states)
# -----------------------------------------------------------------------------
genes_of_interest <- c("TAF6", "ST3GAL4", "SH2B2", "MSL3", "PHGDH", "C19orf84", "LIN7B", "FSD1", "TSPAN33", "EGR4", "PIWIL4", "CELF4", "UTF1", "FGFR3", "A2M", "ENO3", "SERPINE2", "SRRT", "BAG6", "DND1", "PELP1", "NANOS2", "C1QBP", "NANOS3", "GFRA2", "GFRA1", "ID2", "ASB9", "L1TD1", "ID4", "MKI67", "PDPN", "KIT", "DMRT1", "DNMT1", "CALR", "SYCP3", "STRA8")

DefaultAssay(so.integrated) <- "SCT"
dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels

pdf(file = file.path(wd.de.plots, paste0("Di Persio_SSC_DotPlot_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0.pdf")), width = 14, height = 6)
print(dot_plot)
dev.off()

# Apply this custom order
# Note: Only change levels if 'Idents' are factor. If they are characters, convert them to factors first.
#new_order <- c("8", "2", "21", "20", "11", "25", "19", "4", "10", "22", "15", "23", "17", "14", "3", "6", "9", "12", "13", "24", "7", "18", "0", "16", "5", "1")
#new_order <- rev(new_order)
Idents(so.integrated) <- factor(Idents(so.integrated), levels = new_order)

dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
scale_y_discrete(limits = new_order) # Ensure the new order is used in plotting

pdf(file = file.path(wd.de.plots, paste0("Di Persio_SSC_DotPlot_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_ordered.pdf")), width = 14, height = 6)
print(dot_plot)
dev.off()

# Create a mapping from cluster ID to cell type name
# This mapping should be adjusted according to your specific dataset and clustering results
table(Idents(so.integrated))
write.table(table(Idents(so.integrated)), file=file.path(wd.de.data, paste0("Di Persio_SSC_DotPlot_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_ordered_annotated.txt")), row.names=T, col.names=T, quote=F, sep='\t')

cluster_to_celltype <- c('5' = 'Stage 0', '8' = 'Stage 0A', '19' = 'Stage 0B', '1' = 'Stage 1',
																									'16' = 'Stage 2', '15' = 'Stage 3',
																									'3' = 'Leptotene',
																									'6' = 'Zygotene', '2' = 'Zygotene',
																								 '10' = 'Pachytene', '12' = 'Pachytene',
																									'18' = 'Diplotene',
																									'14' = 'Meiotic division',
																									'9' = 'Early spermatid',
																									'13' = 'Late spermatid', '11' = 'Late spermatid',
																									'23' = 'Endothelial and Macrophage',
																									'22' = 'PMC',
																									'7' = 'Fibrotic PMC',
																									'21' = 'Leydig',
																									'0' = 'Sertoli', '17' = 'Sertoli')

# Update the identities using this mapping
#new_order <- c("16", "14", "13", "11", "10", "7", "5", "6", "8", "12", "15", "9", "1", "2", "4", "0", "3")
Idents(so.integrated) <- factor(Idents(so.integrated), levels = new_order)

Idents(so.integrated) <- plyr::mapvalues(x = Idents(so.integrated), from = names(cluster_to_celltype), to = cluster_to_celltype)

DefaultAssay(so.integrated) <- "SCT"
dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
scale_y_discrete(limits = new_order) # Ensure the new order is used in plotting

pdf(file = file.path(wd.de.plots, paste0("Di Persio_DotPlot_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_ordered_annotated.pdf")), width = 14, height = 6)
print(dot_plot)
dev.off()

##
# Rename specific cluster identity
# Rename specific cluster identity
Idents(so.integrated) <- plyr::mapvalues(
   Idents(so.integrated),
   from = "Endothelial and Macrophage",
   to   = "Endothelial\nand Macrophage"
)
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
      plot.title = element_blank(),  # Remove title
      axis.title = element_text(size = 18),  # Increase axis title size
      axis.text = element_text(size = 16),  # Increase axis tick label size
   )
ggsave(file.path(wd.de.plots, paste0("DimPlot_SSC_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_7*7.png")), plot = dim_plot, width = 7, height = 7, dpi = 300)

Idents(so.integrated) <- plyr::mapvalues(
   Idents(so.integrated),
   from = "Endothelial\nand Macrophage",
   to   = "Endothelial and Macrophage"
)
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

save(prin_comp, so.integrated, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_ordered_annotated.RData")))
