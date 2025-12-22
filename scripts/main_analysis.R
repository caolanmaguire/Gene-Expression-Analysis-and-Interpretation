# main_analysis.R | Run all analysis scripts in order

sink(file = "application_run_output.txt")


cat("
   __  ____________                        _           __              ____  _                     _            _       __         
  / / / / ____/ __ \\   ____  _________    (_)__  _____/ /_            / __ \\(_)___     ____  _____(_)___  _____(_)___  / /__  _____
 / / / / /   / / / /  / __ \\/ ___/ __ \\  / / _ \\/ ___/ __/  ______   / __  / / __ \\   / __ \\/ ___/ / __ \\/ ___/ / __ \\/ / _ \\/ ___/
/ /_/ / /___/ /_/ /  / /_/ / /  / /_/ / / /  __/ /__/ /_   /_____/  / /_/ / / /_/ /  / /_/ / /  / / / / / /__/ / /_/ / /  __(__  ) 
\\____/\\____/_____/  / .___/_/   \\____/_/ /\\___/\\___/\\__/           /_____/_/\\____/  / .___/_/  /_/_/ /_/\\___/_/ .___/_/\\___/____/  
                   /_/              /___/                                          /_/                       /_/
Caolán Maguire - UCD -2525669
\n")



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

sink(file = NULL)
