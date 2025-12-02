# Date: Dec 02 2025
# Author: Caolán Maguire
# Assignment: Assignment 2

# set working dir
setwd("./")

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
install.packages(c(
  "tidyverse",
  "pheatmap",
  "glmnet",
  "survival",
  "survminer",
  "ggrepel",
  "RColorBrewer"
))


# Load libraries
library(tidyverse)
# library(DESeq2)

cat("Env setup\n")