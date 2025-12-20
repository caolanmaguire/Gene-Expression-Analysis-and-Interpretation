# main_analysis.R | Run all analysis scripts in order



# Record start time
start_time <- Sys.time()

# Step 1: Data Preprocessing
cat("Step 1/5: Data Preprocessing...\n")
source("scripts/01_data_preprocessing.R")

# Step 2: Differential Expression
cat("\nStep 2/5: Differential Expression Analysis...\n")
source("scripts/02_differential_expression.R")

# Step 3: Pathway Enrichment
cat("\nStep 3/5: Pathway Enrichment Analysis...\n")
source("scripts/03_pathway_enrichment.R")

# Step 4: Visualisation
cat("\nStep 4/5: Visualization (PCA & Heatmap)...\n")
source("scripts/04_visualization.R")

# Step 5: Survival Analysis
cat("\nStep 5/5: Survival Analysis...\n")
source("scripts/05_survival_analysis.R")

# Pipeline complete
end_time <- Sys.time()
elapsed_time <- difftime(end_time, start_time, units = "mins")

cat("Total runtime:", round(elapsed_time, 2), "minutes\n")
cat("\nOutputs saved to:results and figures directories \n")