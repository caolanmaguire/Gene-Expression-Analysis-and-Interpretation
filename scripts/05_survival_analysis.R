# Load libraries
library(DESeq2)
library(survival)
library(glmnet)
library(survminer)

# Load 3 datasets I need for survival analysisGG

# Load preprocessed DESeq2 object
cat("Loading preprocessed data...\n")
dds <- readRDS("results/dds_preprocessed.rds")

# Load clinical data
cat("Loading clinical data...\n")
clinical <- read.delim("data/data_clinical_patient.txt", comment.char = "#")
warnings()

# Load DE genes
cat("Loading DE genes...\n")
de_sig <- read.csv("results/DE_results_significant.csv")

# Check what we loaded
cat("DDS object:", nrow(dds), "genes x", ncol(dds), "samples\n")
cat("Clinical data:", nrow(clinical), "patients\n")
cat("DE genes:", nrow(de_sig), "genes\n")

vst_data <- vst(dds, blind = FALSE)
vst_matrix <- assay(vst_data)
de_genes <- de_sig$gene
vst_de <- vst_matrix[de_genes, ]
cat("Expression matrix for DE genes:", nrow(vst_de), "genes x", ncol(vst_de), "samples\n")

# Transpose
vst_de_t <- t(vst_de)

sample_ids <- rownames(vst_de_t)
patient_ids <- substr(sample_ids, 1, 12)
cat("Sample ID example:", sample_ids[1], "\n")
cat("Patient ID example:", patient_ids[1], "\n")

matched_indices <- match(patient_ids, clinical$PATIENT_ID)
clinical_matched <- clinical[matched_indices, ]
cat("Samples with matching clinical data:", sum(!is.na(matched_indices)), "\n")

# preparing for survival data analysis
surv_data <- data.frame(
  time = clinical_matched$OS_MONTHS,
  clinical = clinical_matched$OS_STATUS,
  vst_de_t,
  row.names = rownames(vst_de_t)
)

cat("Combined data dimensions:", nrow(surv_data), "samples x", ncol(surv_data), "variables\n")

cat("\nChecking for missing data:\n")
cat("NAs in time:", sum(is.na(surv_data$time)), "\n")
cat("NAs in clinical:", sum(is.na(surv_data$clinical)), "\n")
cat("Total rows with ANY NA:", sum(!complete.cases(surv_data)), "\n")
cat("NAs in gene expression columns:", sum(is.na(surv_data[, -(1:2)])), "\n")
cat("Rows before filtering:", nrow(surv_data), "\n")

surv_data_complete <- surv_data
cat("Rows after filtering:", nrow(surv_data_complete), "\n")

cat("\nSurvival status breakdown:\n")
print(table(surv_data_complete$clinical))

# Remove patients who unfortunately died immediately (time = 0 or negative)
surv_data_complete <- surv_data_complete[surv_data_complete$time > 0, ]
cat("After removing time <= 0:", nrow(surv_data_complete), "samples\n")

# --- LASSO COX REGRESSION ---
cat("\n--- Running Lasso Cox Regression ---\n")

# Prepare data matrices
x <- as.matrix(surv_data_complete[, -(1:2)])  # Gene expression matrix
y <- Surv(surv_data_complete$time, surv_data_complete$clinical == "1:DECEASED")  # Survival object

# Run cross-validated Lasso (FIX: was glmnet, should be cv.glmnet)
cv_fit <- cv.glmnet(x, y, family="cox", alpha=1, nfolds=10)

# best_lambda <- cv_fit$lambda.1se
best_lambda <- cv_fit$lambda.min

cat("Best Lambda:", best_lambda, "\n")

# Fit final model with best lambda
final_model <- glmnet(x, y, family="cox", alpha=1, lambda=best_lambda)
coefs <- coef(final_model, s=best_lambda)
selected_genes <- rownames(coefs)[coefs[,1] != 0]

cat("Number of genes selected by Lasso:", length(selected_genes), "\n")

# --- HANDLE CASE WHERE LASSO SELECTS 0 GENES ---
if (length(selected_genes) == 0) {
  cat("\nLasso selected 0 genes (regularization too strong).\n")
  cat("Using alternative approach: top 20 most variable DE genes\n")
  
  # Select top 20 genes by variance across samples
  gene_vars <- apply(x, 2, var)
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:20]
  selected_genes <- top_genes
  
  cat("Selected genes:", length(selected_genes), "\n")
  print(selected_genes)
} else {
  cat("Selected genes:", selected_genes, "\n")
}

# --- CREATE SURVIVAL MODEL WITH SELECTED GENES ---
cat("\n--- Creating Survival Model ---\n")

# Extract expression data for selected genes only
x_selected <- x[, selected_genes, drop = FALSE]

# Fit Cox model with selected genes
surv_df <- data.frame(
  time = surv_data_complete$time,
  status = surv_data_complete$clinical,
  x_selected
)

# Build Cox formula
cox_formula <- as.formula(paste("Surv(time, status == '1:DECEASED') ~", 
                                paste(selected_genes, collapse = " + ")))
cox_model <- coxph(cox_formula, data = surv_df)

# Calculate risk scores for each patient
risk_score <- predict(cox_model, type = "risk")
risk_group <- ifelse(risk_score > median(risk_score), "High", "Low")

cat("High risk patients:", sum(risk_group == "High"), "\n")
cat("Low risk patients:", sum(risk_group == "Low"), "\n")

# Save selected genes with coefficients
prognostic_genes <- data.frame(
  gene = selected_genes,
  coefficient = coef(cox_model)[selected_genes]
)
write.csv(prognostic_genes, "results/lasso_selected_genes.csv", row.names = FALSE)
cat("Saved: results/lasso_selected_genes.csv\n")

# --- KAPLAN-MEIER SURVIVAL PLOT ---
cat("\n--- Creating Kaplan-Meier Plot ---\n")

# Create survival object
surv_obj <- Surv(surv_data_complete$time, surv_data_complete$clinical == "1:DECEASED")

# Fit survival curves for high vs low risk groups
fit <- survfit(surv_obj ~ risk_group)

# Create KM plot
km_plot <- ggsurvplot(
  fit,
  data = data.frame(risk_group = risk_group),
  pval = TRUE,                    # Show p-value from log-rank test
  risk.table = TRUE,              # Show risk table
  title = "Survival Analysis: High vs Low Risk Groups",
  xlab = "Time (months)",
  ylab = "Survival Probability",
  legend.title = "Risk Group",
  legend.labs = c("High Risk", "Low Risk")
)

# Save plot
if (!dir.exists("figures")) {
  dir.create("figures")
}

pdf("figures/kaplan_meier_risk_groups.pdf", width = 10, height = 8)
print(km_plot)
dev.off()

cat("Saved: figures/kaplan_meier_risk_groups.pdf\n")
cat("\n=== SURVIVAL ANALYSIS COMPLETE ===\n")