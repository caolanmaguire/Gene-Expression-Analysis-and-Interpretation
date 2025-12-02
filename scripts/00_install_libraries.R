# Libraries I hadn't installed on my personal machine
install.packages("BiocManager")
BiocManager::install("DESeq2")

BiocManager::install(c(
  "DESeq2",
  "clusterProfiler",
  "org.Hs.eg.db",
  "EnhancedVolcano",
  "ComplexHeatmap"
))


# Install all CRAN packages needed
install.packages("tidyverse")
install.packages("pheatmap")
install.packages("glmnet")
install.packages("survival")
install.packages("survminer")
install.packages("ggrepel")
install.packages("RColorBrewer")