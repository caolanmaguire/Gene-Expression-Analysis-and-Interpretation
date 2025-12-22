# Date: Dec 02 2025
# Author: Caolán Maguire
# Assignment: Assignment 2

library(DESeq2)

# Load the dds object saved in pre-processing script
dds <- readRDS("results/dds_preprocessed.rds")

# apply DESeq
dds <- DESeq(dds)

# Extract results comparing Amplified vs Not_Amplified
res <- results(dds, contrast = c("ERBB2_status", "Amplified", "Not_Amplified"))

# summary
summary(res)

# save results to dataframe
res_df <- as.data.frame(res)

# transform gene names from rownaes to columns
res_df$gene <- rownames(res_df)

# save all results now
write.csv(res_df, "results/DE_results_all_genes.csv", row.names = FALSE)

# Filter for genes (padj < 0.05 AND |log2FC| > 1)
res_sig <- res_df[!is.na(res_df$padj) &
                    res_df$padj < 0.05 &
                    abs(res_df$log2FoldChange) > 1, ]

#print(res_sig)

# save these significant genes filtered above
write.csv(res_sig, "results/DE_results_significant.csv", row.names = FALSE)

# summary
cat("Total significant genes:", nrow(res_sig), "\n")


# Get Top 10 Genes (5 up, 5 down)


cat("\n--- Getting Top 10 Genes ---\n")

# Sort by absolute fold change
res_sig_sorted <- res_sig[order(abs(res_sig$log2FoldChange), decreasing = TRUE), ]

# Get top 5 upregulated
top5_up <- res_sig_sorted[res_sig_sorted$log2FoldChange > 0, ][1:5, ]

# Get top 5 downregulated
top5_down <- res_sig_sorted[res_sig_sorted$log2FoldChange < 0, ][1:5, ]

print(top5_up)
print(top5_down)

# Combine
top10 <- rbind(top5_up, top5_down)

# Save
write.csv(top10, "results/top10_DE_genes.csv", row.names = FALSE)

cat("Top 10 genes saved\n")
print(top10[, c("gene", "log2FoldChange", "pvalue", "padj")])

# - - - - - - - - - - - - - - -
# Create Volcano Plot
# - - - - - - - - - - - - - - -

cat("\n--- Creating Volcano Plot ---\n")

library(ggplot2)
library(ggrepel)

# Add significance column for coloring
res_df$significant <- ifelse(
  !is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,
  "Significant",
  "Not Significant"
)

# Add labels for top genes and ERBB2
res_df$label <- ""
res_df$label[res_df$gene %in% top10$gene] <- res_df$gene[res_df$gene %in% top10$gene]
if ("ERBB2" %in% res_df$gene) {
  res_df$label[res_df$gene == "ERBB2"] <- "ERBB2"
}

# Create plot
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = significant), alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 20) +
  labs(
    title = "Differential Gene Expression: HER2+ vs Non-HER2+",
    x = "Log2 Fold Change",
    y = "-Log10 P-value",
    color = ""
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

# Create figures directory if needed
if (!dir.exists("figures")) {
  dir.create("figures")
}

# Save plot
ggsave("figures/volcano_plot.pdf", volcano_plot, width = 8, height = 6)

cat("Volcano plot saved: figures/volcano_plot.pdf\n")
cat("\n=== Script 02 Complete ===\n")