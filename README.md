# HER2+ Breast Cancer Gene Expression Analysis

[![R Version](https://img.shields.io/badge/R-4.3.0+-blue.svg)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Status](https://img.shields.io/badge/status-active-success.svg)]()

> Differential gene expression analysis of HER2-amplified vs non-amplified breast cancer tumours using TCGA data

---

## Table of Contents

- [Overview](#overview)
- [Key Findings](#key-findings)
- [Project Structure](#project-structure)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Data Download](#data-download)
- [Usage](#usage)
- [Outputs](#outputs)
- [Methods Summary](#methods-summary)
- [Results](#results)

---

## Overview

This project analyses RNA-seq data from The Cancer Genome Atlas (TCGA) breast cancer cohort to identify differentially expressed genes between **HER2-amplified (ERBB2+)** and **non-amplified** tumours.

### Objectives

- Identify differentially expressed genes in HER2+ breast cancer
- Perform pathway enrichment analysis to understand biological mechanisms
- Develop prognostic gene signatures using survival analysis
- Visualise expression patterns through dimensionality reduction

### Clinical Relevance

HER2+ breast cancer represents ~15-20% of cases and is associated with aggressive disease. Despite targeted therapies (trastuzumab, pertuzumab), response rates remain around 40%. Understanding the molecular landscape may reveal:

- Novel therapeutic targets
- Biomarkers for treatment response
- Mechanisms of resistance

---

## Key Findings

### TODO findings to be added
---

## Project Structure

### TODO project structure to be added

---

## Prerequisites

### Software Requirements

- **R** ≥ 4.3.0
- **RStudio** (recommended, optional)
- **Git** for version control

### R Packages Required

#### Bioconductor Packages
```r
DESeq2              # Differential expression analysis
clusterProfiler     # Pathway enrichment
org.Hs.eg.db        # Human gene annotations
EnhancedVolcano     # Volcano plots
ComplexHeatmap      # Advanced heatmaps
```

#### CRAN Packages
```r
tidyverse           # Data manipulation and visualisation
pheatmap            # Heatmaps
glmnet              # Lasso regression
survival            # Survival analysis
survminer           # Survival visualisation
ggrepel             # Label placement
RColorBrewer        # Colour palettes
```

---

## 🚀 Installation

### Step 1: Clone the Repository
```bash
git clone https://github.com/yourusername/her2-breast-cancer-analysis.git
cd her2-breast-cancer-analysis
```

### Step 2: Install R Packages

Open R or RStudio and run:
```r
source 00_install_libararies.R # - A install script that installs all required libararies
```

### Step 3: Verify Installation
```r
# Check package versions
packageVersion("DESeq2")
packageVersion("clusterProfiler")
```

---

## Data Download

### Direct Download from cBioPortal

1. Visit [TCGA Breast Cancer Study](https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018)

2. Click **"Download"** in the top right

3. Download these files:
   - `data_mrna_seq_v2_rsem.txt` (RNA-seq data)
   - `data_clinical_patient.txt` (Clinical annotations)
   - `data_cna.txt` (Copy number alterations)

4. Place files in the `data/` directory


### Data Structure

After download, your `data/` directory should contain:
```
data/
├── data_mrna_seq_v2_rsem.txt      (~20,000 genes × ~1,000 samples)
├── data_clinical_patient.txt       (Patient clinical data)
└── data_cna.txt                    (Copy number alteration data)
```

---

## Usage

### Quick Start: Run Complete Analysis
```r
# Set working directory
setwd("path/to/her2-breast-cancer-analysis")

# Run complete pipeline
source("scripts/main_analysis.R")
```

This will:
1. Load and preprocess data
2. Perform differential expression analysis
3. Run pathway enrichment
4. Generate all visualisations
5. Perform survival analysis
6. Save results to `results/` and `figures/`

**Estimated runtime**: ~10-15 minutes

---

### Step-by-Step Execution

If you prefer to run analyses individually:
```r
# Step 1: Data Preprocessing
source("scripts/01_data_preprocessing.R")
# Outputs: cleaned count matrix, metadata

# Step 2: Differential Expression
source("scripts/02_differential_expression.R")
# Outputs: DE results, volcano plot, MA plot

# Step 3: Pathway Enrichment
source("scripts/03_pathway_enrichment.R")
# Outputs: GO/KEGG enrichment results and plots

# Step 4: Visualisation
source("scripts/04_visualization.R")
# Outputs: PCA plot, heatmap

# Step 5: Survival Analysis
source("scripts/05_survival_analysis.R")
# Outputs: Prognostic genes, Kaplan-Meier curves
```

---

## Outputs

### Figures Generated

| File | Description |
|------|-------------|
| `volcano_plot.pdf` | Differential expression volcano plot |
| `MA_plot.pdf` | Mean-average plot showing DE distribution |
| `GO_dotplot_upregulated.pdf` | GO enrichment for upregulated genes |
| `GO_dotplot_downregulated.pdf` | GO enrichment for downregulated genes |
| `KEGG_barplot_upregulated.pdf` | KEGG pathway enrichment |
| `PCA_plot.pdf` | Principal component analysis |
| `heatmap_top50_genes.pdf` | Expression heatmap of top 50 DE genes |
| `kaplan_meier_risk_groups.pdf` | Survival curves for risk groups |
| `lasso_coefficients_barplot.pdf` | Prognostic gene coefficients |

### Results Files

| File | Description |
|------|-------------|
| `DE_results_all_genes.csv` | All genes with statistics |
| `DE_results_significant.csv` | Significantly DE genes only |
| `top10_DE_genes.csv` | Top 10 up/downregulated genes |
| `GO_enrichment_upregulated.csv` | GO terms enriched in upregulated genes |
| `GO_enrichment_downregulated.csv` | GO terms enriched in downregulated genes |
| `lasso_selected_genes.csv` | Prognostic genes from survival analysis |

---

## Methods Summary

### Differential Expression Analysis

- **Tool**: DESeq2 v1.40.0
- **Normalization**: Median of ratios method
- **Statistical test**: Wald test
- **Significance threshold**: padj < 0.05, |log2FC| > 1
- **Multiple testing correction**: Benjamini-Hochberg (FDR)

### Pathway Enrichment

- **Tool**: clusterProfiler v4.8.0
- **Databases**: 
  - GO Biological Process
  - KEGG pathways
- **Background**: All expressed genes in the dataset
- **Significance**: adjusted p-value < 0.05

### Dimensionality Reduction

- **Transformation**: Variance Stabilising Transformation (VST)
- **PCA**: Base R `prcomp()` function
- **Heatmap**: Top 50 DE genes, z-score normalised

### Survival Analysis

- **Method**: Lasso Cox proportional hazards regression
- **Tool**: glmnet v4.1-6
- **Cross-validation**: 10-fold CV
- **Lambda selection**: Lambda.min
- **Risk stratification**: Median split of risk scores

---

## Results

### Differential Expression

- **Total genes analyzed**: 20,530
- **Significant DE genes**: XXX (padj < 0.05, |log2FC| > 1)
  - Upregulated: XXX genes
  - Downregulated: XXX genes

### Top Upregulated Genes

| Gene | log2FC | padj | Function |
|------|--------|------|----------|
| ERBB2 | X.XX | <0.001 | Receptor tyrosine kinase |
| GRB7 | X.XX | <0.001 | Adaptor protein |
| ... | ... | ... | ... |

### Enriched Pathways

**Upregulated in HER2+:**
- Cell proliferation
- EGFR signalling pathway
- PI3K-Akt signaling

**Downregulated in HER2+:**
- Immune response
- ECM organisation
- Cell adhesion

### Survival Analysis

- **Prognostic genes selected**: XX
- **Risk groups**: High vs Low (median split)
- **Log-rank test**: p = X.XXX
- **Clinical implications**: [To be completed]

---
