# Date: Dec 10 2025
# Author: Caolán Maguire
# Script 04: Visualization (PCA and Heatmap)

cat("\n=== VISUALISATION ===\n")

# Load libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# Load preprocessed data
cat("Loading preprocessed data... \n")
dds <- readRDS("results/dds_preprocessed.rds")

# Get VST transformed data
cat("Performing variance stabilizing transformation...\n")
vst_data <- vst(dds, blind = FALSE)

# - - - - - PCA PLOT - - - - -
cat("\n--- Creating PCA Plot ---\n")

# Run PCA
pca_data <- plotPCA(vst_data, intgroup = "ERBB2_status", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"), 1)

cat("PC1 explains", percent_var[1], "% variance\n")
cat("PC2 explains", percent_var[2], "% variance\n")

# Create PCA plot
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = ERBB2_status)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(
    values = c("Not_Amplified" = "blue", "Amplified" = "red"),
    labels = c("Not_Amplified" = "Non-HER2+", "Amplified" = "HER2+")
  ) +
  labs(
    title = "Principal Component Analysis",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance"),
    color = "ERBB2 Status"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("figures/PCA_plot.pdf", pca_plot, width = 7, height = 6)
cat("Saved: figures/PCA_plot.pdf\n")

# - - - - - HEATMAP - - - - -
cat("\n--- Creating Heatmap ---\n")

# Load DE results to get top genes
de_results <- read.csv("results/DE_results_significant.csv")
de_results_sorted <- de_results[order(de_results$padj), ]
top10_genes <- de_results_sorted$gene[1:10]

cat("Using top 50 genes by adjusted p-value\n")

# Extract VST values for top 50 genes
vst_matrix <- assay(vst_data)
top10_matrix <- vst_matrix[top10_genes, ]

# Z-score normalization (by row)
top10_scaled <- t(scale(t(top10_matrix)))

# Create annotation for samples
annotation_col <- data.frame(
  ERBB2_Status = colData(vst_data)$ERBB2_status,
  row.names = colnames(vst_data)
)

# Define colors
annotation_colors <- list(
  ERBB2_Status = c(Not_Amplified = "blue", Amplified = "red")
)

# Create heatmap
pdf("figures/heatmap_top10_genes.pdf", width = 10, height = 12)
pheatmap(
  top10_scaled,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  main = "Top 50 Differentially Expressed Genes",
  fontsize_row = 6
)
dev.off()
cat("Saved: figures/heatmap_top50_genes.pdf\n")

cat("\n=== VISUALIZATION COMPLETE ===\n")