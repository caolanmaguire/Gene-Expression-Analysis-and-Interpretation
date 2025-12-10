# Date: Dec 02 2025
# Author: Caolán Maguire
# Assignment: Assignment 2

# set working dir

if (require("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  project_dir <- rstudioapi::selectDirectory(
    caption = "Select Project Root Directory",
    label = "Select"
  )
}else{
  setwd("C:/Users/caola/OneDrive/Documents/GitHub/Gene-Expression-Analysis-and-Interpretation/")
}


# Load libraries
library(tidyverse)
library(DESeq2)

cat("Env setup\n")

# - - - - - - - - - - - - - - - 
# Load Data Files
# - - - - - - - - - - - - - - - 

# cat(getwd())
# cat("\n")

# Define files to load
data_files <- list(
  rna = "data/data_mrna_seq_v2_rsem.txt",
  clinical = "data/data_clinical_patient.txt",
  cna = "data/data_cna.txt"
)

# Create empty list to store data
datasets <- list()

# Loop through and load each file
for (name in names(data_files)) {
  file_path <- data_files[[name]]
  
  cat("Loading", name, "data from:", basename(file_path), "\n")
  
  if (file.exists(file_path)) {
    
    # Special handling for clinical data (has comment lines)
    if (name == "clinical") {
      datasets[[name]] <- read.delim(
        file_path,
        header = TRUE,
        comment.char = "#",
        na.strings = c("", "NA", "[Not Available]", "[Not Applicable]"),
        check.names = FALSE
      )
    } else {
      # RNA and CNA data
      datasets[[name]] <- read.delim(
        file_path,
        header = TRUE,
        check.names = FALSE
      )
    }
    
    cat("  ✓ Loaded:", nrow(datasets[[name]]), "rows x", 
        ncol(datasets[[name]]), "columns\n")
    
    cat("\nFirst few column names:\n")
    print(head(colnames(datasets[[name]])))
    
  } else {
    stop("ERROR: File not found - ", file_path)
  }
}

# cat("✓ RNA data loaded (raw)\n")
# cat("Dimensions:", nrow(rna_data), "rows x", ncol(rna_data), "columns\n")

# Check first few columns
# cat("\nFirst few column names:\n")
# print(head(colnames(rna_data)))

rna_data <- datasets$rna
clinical_data <- datasets$clinical
cna_data <- datasets$cna

# - - - - - Data Inspection - - - -\n")
cat("✓ Data extracted:\n")
cat("  - rna_data: ", nrow(rna_data), "x", ncol(rna_data), "\n")
cat("  - clinical_data: ", nrow(clinical_data), "x", ncol(clinical_data), "\n")
cat("  - cna_data: ", nrow(cna_data), "x", ncol(cna_data), "\n")


# - - - - - Getter for sample ID formats - - - - \n")
# RNA-sequence sample IDS (exlcude the first column)
rna_samples <- colnames(rna_data)[-1] # -1 to remove the first column
cat("RNA-seq ID (first):", rna_samples[1], "\n")
cat("RNA-seq ID length:", nchar(rna_samples[1]), "characters\n")

# CNA sample IDs (exclude Hugo_Symbol column)
cna_samples <- colnames(cna_data)[-1]  # Remove first column
cat("CNA ID (first):", cna_samples[1], "\n")
cat("CNA ID length:", nchar(cna_samples[1]), "characters\n")

# Clinical patient IDs
clinical_samples <- clinical_data$PATIENT_ID
cat("Clinical ID (first):", clinical_samples[1], "\n")
cat("Clinical ID length:", nchar(clinical_samples[1]), "characters\n")

# Check if they match format
cat("\n--- ID Format Check ---\n")
if (rna_samples[1] == cna_samples[1]) {
  cat("✓ RNA and CNA IDs already match!\n")
} else {
  cat("⚠ RNA and CNA IDs have different formats\n")
  cat("  RNA:", rna_samples[1], "\n")
  cat("  CNA:", cna_samples[1], "\n")
}
