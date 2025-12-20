# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# SET CORRECT WORKING DIRECTORY
#setwd("C:/Users/caola/OneDrive/Documents/GitHub/Gene-Expression-Analysis-and-Interpretation")

# Verify it worked
#cat("\nWorking directory:", getwd(), "\n")

# Load your significant DE genes
de_sig <- read.csv("results/DE_results_significant.csv")
cat("Total DE genes: ", nrow(de_sig), "\n")

# Separate by fold change direction
up_genes <- de_sig$gene[de_sig$log2FoldChange > 0]
down_genes <- de_sig$gene[de_sig$log2FoldChange < 0]

# Check split
cat("Upregulated genes: ", length(up_genes), "\n")
cat("Downregulated genes: ", length(down_genes), "\n")

# enrichment - for up genes
go_up <- enrichGO(
  gene = up_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)
cat("Enriched pathways for UP genes:", nrow(go_up), "\n")

# encrichment - for down genes
cat("\n--- GO Enrichment: Downregulated Genes ---\n")
go_down <- enrichGO(
  gene = down_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)
cat("Enriched pathways for DOWN genes:", nrow(go_down), "\n")

# save to csv
write.csv(as.data.frame(go_up), "results/GO_enrichment_upregulated.csv", row.names = FALSE)
write.csv(as.data.frame(go_down), "results/GO_enrichment_downregulated.csv", row.names = FALSE)
cat("\nResults saved to CSV files\n")

# check if directory exsists if not - create directory
if (!dir.exists("figures")) {
  dir.create("figures")
}

# creating plots
pdf("figures/GO_dotplot_upregulated.pdf", width = 10, height = 8)
print(dotplot(go_up, showCategory = 20, title = "GO Enrichment - Upregulated in HER2+"))
dev.off()

pdf("figures/GO_dotplot_downregulated.pdf", width = 10, height = 8)
print(dotplot(go_down, showCategory = 20, title = "GO Enrichment - Downregulated in HER2+"))
dev.off()

cat("Plots saved to figures/ folder\n")
cat("\n=== Pathway Enrichment Complete ===\n")