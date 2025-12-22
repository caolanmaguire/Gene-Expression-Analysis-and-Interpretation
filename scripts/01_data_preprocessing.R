# Date: Dec 02 2025
# Author: Caolán Maguire
# Assignment: Assignment 2

# set working dir
if (require("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  project_dir <- rstudioapi::selectDirectory(
    caption = "Select Project Root Directory",
    label = "Select"
  )
  setwd(project_dir)
} else {
  setwd("C:/Users/caola/OneDrive/Documents/GitHub/Gene-Expression-Analysis-and-Interpretation/")
}

# Check working directory
cat("\n--- Working Directory Check ---\n")
cat("Current directory:", getwd(), "\n")

# Load libraries
library(tidyverse)
library(DESeq2)
cat("Env setup\n")

# - - - - - - - - - - - - - - - 
# Load Data Files
# - - - - - - - - - - - - - - - 

# Define files to load
data_files <- list(
  rna = "data/data_mrna_seq_v2_rsem.txt",
  clinical = "data/data_clinical_patient.txt",
  cna = "data/data_cna.txt"
)

# Check if all files exist first
cat("\n--- File Existence Check ---\n")
all_files_exist <- TRUE
for (name in names(data_files)) {
  file_path <- data_files[[name]]
  exists <- file.exists(file_path)
  
  if (exists) {
    cat("Found", name, ":", file_path, "\n")
  } else {
    cat("MISSING:", file_path, "\n")
    all_files_exist <- FALSE
  }
}

if (!all_files_exist) {
  stop("\nERROR: Some data files are missing. Please check:\n",
       "  1. Working directory is correct: ", getwd(), "\n",
       "  2. Data files are in the 'data/' folder\n",
       "  3. Files are named correctly\n")
}

cat("\nAll files found. Proceeding with data loading...\n\n")

# Create empty list to store data
datasets <- list()

# Loop through and load each file
for (name in names(data_files)) {
  file_path <- data_files[[name]]
  
  cat("Loading", name, "data from:", basename(file_path), "\n")
  
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
  
  cat("  Loaded:", nrow(datasets[[name]]), "rows x", 
      ncol(datasets[[name]]), "columns\n")
  
  cat("\nFirst few column names:\n")
  print(head(colnames(datasets[[name]])))
}

rna_data <- datasets$rna
clinical_data <- datasets$clinical
cna_data <- datasets$cna

# - - - - - Data Inspection - - - -
cat("\n--- Data Extraction ---\n")
cat("Data extracted:\n")
cat("  - rna_data: ", nrow(rna_data), "x", ncol(rna_data), "\n")
cat("  - clinical_data: ", nrow(clinical_data), "x", ncol(clinical_data), "\n")
cat("  - cna_data: ", nrow(cna_data), "x", ncol(cna_data), "\n")

# - - - - - Getter for sample ID formats - - - -
cat("\n--- Sample ID Format Check ---\n")

# RNA-sequence sample IDS (exclude the first column)
rna_samples <- colnames(rna_data)[-1]
cat("RNA-seq ID (first):", rna_samples[1], "\n")
cat("RNA-seq ID length:", nchar(rna_samples[1]), "characters\n")

# CNA sample IDs (exclude Hugo_Symbol column)
cna_samples <- colnames(cna_data)[-1]
cat("CNA ID (first):", cna_samples[1], "\n")
cat("CNA ID length:", nchar(cna_samples[1]), "characters\n")

# Clinical patient IDs
clinical_samples <- clinical_data$PATIENT_ID
cat("Clinical ID (first):", clinical_samples[1], "\n")
cat("Clinical ID length:", nchar(clinical_samples[1]), "characters\n")

# Check if they match format
cat("\n--- ID Format Comparison ---\n")
if (rna_samples[1] == cna_samples[1]) {
  cat("RNA and CNA IDs already match.\n")
} else {
  cat("RNA and CNA IDs have different formats\n")
  cat("  RNA:", rna_samples[1], "\n")
  cat("  CNA:", cna_samples[1], "\n")
}

cat("\nPreprocessing complete.\n")

# - - - - - - - - - - - - - - - 
# Match Sample IDs
# - - - - - - - - - - - - - - - 

cat("\n--- Matching Sample IDs ---\n")

# Get sample IDs from each dataset (skip first TWO columns: Hugo_Symbol and Entrez_Gene_Id)
rna_samples <- colnames(rna_data)[-(1:2)]
cna_samples <- colnames(cna_data)[-(1:2)]
clinical_samples <- clinical_data$PATIENT_ID

cat("RNA samples:", length(rna_samples), "\n")
cat("CNA samples:", length(cna_samples), "\n")
cat("Clinical samples:", length(clinical_samples), "\n")

# Extract patient IDs from RNA/CNA (remove -01 suffix to match clinical)
rna_patient_ids <- substr(rna_samples, 1, 12)  # Get first 12 characters
cna_patient_ids <- substr(cna_samples, 1, 12)

# Find common patient IDs
common_patient_ids <- Reduce(intersect, list(rna_patient_ids, cna_patient_ids, clinical_samples))

cat("Common patient IDs across all datasets:", length(common_patient_ids), "\n")

# Map back to original sample names (with -01 suffix)
rna_samples_matched <- rna_samples[rna_patient_ids %in% common_patient_ids]
cna_samples_matched <- cna_samples[cna_patient_ids %in% common_patient_ids]

cat("Matched RNA samples:", length(rna_samples_matched), "\n")
cat("Matched CNA samples:", length(cna_samples_matched), "\n")

# Subset datasets
rna_data_matched <- rna_data[, c("Hugo_Symbol", "Entrez_Gene_Id", rna_samples_matched)]
cna_data_matched <- cna_data[, c("Hugo_Symbol", "Entrez_Gene_Id", cna_samples_matched)]
clinical_data_matched <- clinical_data[clinical_data$PATIENT_ID %in% common_patient_ids, ]

cat("Matched datasets created.\n")

# - - - - - - - - - - - - - - - 
# Extract ERBB2 Status from CNA
# - - - - - - - - - - - - - - - 

cat("\n--- Extracting ERBB2 Status ---\n")

# Find ERBB2 row in CNA data
erbb2_idx <- which(cna_data_matched$Hugo_Symbol == "ERBB2")

if (length(erbb2_idx) == 0) {
  stop("ERROR: ERBB2 not found in CNA data")
}

# Extract ERBB2 CNA values for each sample (skip gene name columns)
erbb2_cna <- as.numeric(cna_data_matched[erbb2_idx, -(1:2)])
names(erbb2_cna) <- cna_samples_matched

# Define amplification status (CNA > 0 = amplified)
erbb2_status <- ifelse(erbb2_cna > 0, "Amplified", "Not_Amplified")

cat("ERBB2 amplified samples:", sum(erbb2_status == "Amplified"), "\n")
cat("ERBB2 not amplified samples:", sum(erbb2_status == "Not_Amplified"), "\n")

# - - - - - - - - - - - - - - - 
# Create Metadata
# - - - - - - - - - - - - - - - 

cat("\n--- Creating Metadata ---\n")

# Use RNA sample names (with -01 suffix) as rownames
metadata <- data.frame(
  sample_id = rna_samples_matched,
  patient_id = substr(rna_samples_matched, 1, 12),
  ERBB2_status = factor(erbb2_status[cna_samples_matched], levels = c("Not_Amplified", "Amplified")),
  row.names = rna_samples_matched
)

cat("Metadata created with", nrow(metadata), "samples\n")
print(table(metadata$ERBB2_status))


# - - - - - - - - - - - - - - - 
# Prepare Count Matrix
# - - - - - - - - - - - - - - - 

cat("\n--- Preparing Count Matrix ---\n")

# Remove empty gene names
rna_data_matched <- rna_data_matched[rna_data_matched$Hugo_Symbol != "", ]

# Remove duplicate gene names (keep first occurrence)
dup_idx <- duplicated(rna_data_matched$Hugo_Symbol)
rna_data_matched <- rna_data_matched[!dup_idx, ]

cat("Rows after removing empty/duplicate genes:", nrow(rna_data_matched), "\n")

# Extract count data (remove Hugo_Symbol and Entrez_Gene_Id columns)
count_matrix <- as.matrix(rna_data_matched[, -(1:2)])
rownames(count_matrix) <- rna_data_matched$Hugo_Symbol

# Round to integers
count_matrix <- round(count_matrix)

cat("Count matrix dimensions:", nrow(count_matrix), "genes x", ncol(count_matrix), "samples\n")


# - - - - - - - - - - - - - - - 
# Filter Low-Count Genes
# - - - - - - - - - - - - - - - 

cat("\n--- Filtering Low-Count Genes ---\n")

# Keep genes with at least 10 counts in at least 70% of samples
min_count <- 10
min_samples <- ceiling(0.7 * ncol(count_matrix))

keep_genes <- rowSums(count_matrix >= min_count) >= min_samples

count_matrix_filtered <- count_matrix[keep_genes, ]

cat("Original genes:", nrow(count_matrix), "\n")
cat("After filtering:", nrow(count_matrix_filtered), "\n")
cat("Genes removed:", nrow(count_matrix) - nrow(count_matrix_filtered), "\n")

# - - - - - - - - - - - - - - - 
# Create DESeq2 Dataset
# - - - - - - - - - - - - - - - 

cat("\n--- Creating DESeq2 Dataset ---\n")

# Check that sample order matches
if (!all(colnames(count_matrix_filtered) == rownames(metadata))) {
  stop("ERROR: Sample order mismatch between count matrix and metadata")
}

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_filtered,
  colData = metadata,
  design = ~ ERBB2_status
)

cat("DESeq2 dataset created successfully\n")
cat("Dimensions:", nrow(dds), "genes x", ncol(dds), "samples\n")

# Save Preprocessed Data
# Create results directory if it doesn't exist
if (!dir.exists("results")) {
  dir.create("results")
}

# Save DESeq2 object
saveRDS(dds, "results/dds_preprocessed.rds")
cat("Saved: results/dds_preprocessed.rds\n")

# Save metadata
write.csv(metadata, "results/metadata.csv", row.names = TRUE)
cat("Saved: results/metadata.csv\n")

cat("\n=== PREPROCESSING COMPLETE ===\n\n\n\n")
