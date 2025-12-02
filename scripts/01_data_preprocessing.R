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
    print(head(colnames(datasets[[]])))
    
  } else {
    stop("ERROR: File not found - ", file_path)
  }
}

# cat("✓ RNA data loaded (raw)\n")
# cat("Dimensions:", nrow(rna_data), "rows x", ncol(rna_data), "columns\n")

# Check first few columns
# cat("\nFirst few column names:\n")
# print(head(colnames(rna_data)))
