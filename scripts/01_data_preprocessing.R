# Date: Dec 02 2025
# Author: Caolán Maguire
# Assignment: Assignment 2

# set working dir
setwd("C:/Users/caola/OneDrive/Documents/GitHub/Gene-Expression-Analysis-and-Interpretation/")

# Load libraries
library(tidyverse)
library(DESeq2)

cat("Env setup\n")

# - - - - - - - - - - - - - - - 
# Load Data Files
# - - - - - - - - - - - - - - - 

cat(getwd())
cat("\n")

cat("- - - Loading Files - - - \n")
rna_data <- read.delim("data/data_mrna_seq_v2_rsem.txt",
                       header = TRUE,
                       check.names = FALSE)

cat("RNA data loaded \n")

cat("✓ RNA data loaded (raw)\n")
cat("Dimensions:", nrow(rna_data), "rows x", ncol(rna_data), "columns\n")

# Check first few columns
cat("\nFirst few column names:\n")
print(head(colnames(rna_data)))
