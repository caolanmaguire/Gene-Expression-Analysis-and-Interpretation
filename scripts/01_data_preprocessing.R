# Date: Dec 02 2025
# Author: Caolán Maguire
# Assignment: Assignment 2

# set working dir
setwd("./")

# Load libraries
library(tidyverse)
library(DESeq2)

cat("Env setup\n")

# - - - - - - - - - - - - - - - 
# Load Data Files
# - - - - - - - - - - - - - - - 

cat(getwd())

cat("- - - Loading Files - - - \n")
rna_data <- read.delim("data/data_mrna_seq_v2_rsem.txt",
                       header = TRUE,
                       row.names = 1,
                       check.names = FALSE)

cat("RNA data loaded \n")
