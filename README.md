# HER2+ Breast Cancer Gene Expression Analysis

[![R Version](https://img.shields.io/badge/R-4.3.0+-blue.svg)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Status](https://img.shields.io/badge/status-active-success.svg)]()

> Differential gene expression analysis of HER2-amplified vs non-amplified breast cancer tumours using TCGA data

---

## 📋 Table of Contents

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
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)
- [Acknowledgments](#acknowledgments)

---

## 🔬 Overview

This project analyzes RNA-seq data from The Cancer Genome Atlas (TCGA) breast cancer cohort to identify differentially expressed genes between **HER2-amplified (ERBB2+)** and **non-amplified** tumours.

### Objectives

- Identify differentially expressed genes in HER2+ breast cancer
- Perform pathway enrichment analysis to understand biological mechanisms
- Develop prognostic gene signatures using survival analysis
- Visualize expression patterns through dimensionality reduction

### Clinical Relevance

HER2+ breast cancer represents ~15-20% of cases and is associated with aggressive disease. Despite targeted therapies (trastuzumab, pertuzumab), response rates remain around 40%. Understanding the molecular landscape may reveal:

- Novel therapeutic targets
- Biomarkers for treatment response
- Mechanisms of resistance

---

## 🎯 Key Findings

> **Note**: This section will be updated after analysis completion

- **XXX genes** significantly differentially expressed (padj < 0.05, |log2FC| > 1)
- Top pathways enriched: [To be completed]
- Prognostic signature identified: **XX genes** associated with survival
- Risk stratification achieved significant separation (log-rank p = X.XXX)

---

## 📁 Project Structure
```
her2-breast-cancer-analysis/
│
├── README.md                          # This file
├── .gitignore                         # Git ignore rules
├── requirements.txt                   # R package dependencies
│
├── scripts/                           # Analysis scripts
│   ├── 01_data_preprocessing.R        # Data loading and QC
│   ├── 02_differential_expression.R   # DESeq2 analysis
│   ├── 03_pathway_enrichment.R        # GO/KEGG enrichment
│   ├── 04_visualization.R             # PCA, heatmaps
│   ├── 05_survival_analysis.R         # Lasso Cox regression
│   └── main_analysis.R                # Run complete pipeline
│
├── functions/                         # Helper functions
│   └── helper_functions.R             # Reusable utility functions
│
├── data/                              # Raw data (not tracked by git)
│   ├── .gitkeep                       # Placeholder
│   └── README.txt                     # Data download instructions
│
├── figures/                           # Generated plots
│   ├── volcano_plot.pdf
│   ├── PCA_plot.pdf
│   ├── heatmap_top50_genes.pdf
│   └── kaplan_meier_risk_groups.pdf
│
└── results/                           # Analysis outputs
    ├── DE_results_all_genes.csv
    ├── DE_results_significant.csv
    ├── top10_DE_genes.csv
    ├── GO_enrichment_upregulated.csv
    └── lasso_selected_genes.csv
```

---

## ⚙️ Prerequisites

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
tidyverse           # Data manipulation and visualization
pheatmap            # Heatmaps
glmnet              # Lasso regression
survival            # Survival analysis
survminer           # Survival visualization
ggrepel             # Label placement
RColorBrewer        # Color palettes
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
# Install BiocManager if not already installed
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c(
  "DESeq2",
  "clusterProfiler",
  "org.Hs.eg.db",
  "EnhancedVolcano",
  "ComplexHeatmap"
))

# Install CRAN packages
install.packages(c(
  "tidyverse",
  "pheatmap",
  "glmnet",
  "survival",
  "survminer",
  "ggrepel",
  "RColorBrewer"
))
```

### Step 3: Verify Installation
```r
# Check package versions
packageVersion("DESeq2")
packageVersion("clusterProfiler")
```

---

## 📊 Data Download

### Option 1: Direct Download from cBioPortal

1. Visit [TCGA Breast Cancer Study](https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018)

2. Click **"Download"** in the top right

3. Download these files:
   - `data_mrna_seq_v2_rsem.txt` (RNA-seq data)
   - `data_clinical_patient.txt` (Clinical annotations)
   - `data_cna.txt` (Copy number alterations)

4. Place files in the `data/` directory

### Option 2: Using wget (Linux/Mac)
```bash
cd data/

# Download dataset (replace with actual URLs if available)
wget https://cbioportal-datahub.s3.amazonaws.com/brca_tcga_pan_can_atlas_2018.tar.gz

# Extract
tar -xzvf brca_tcga_pan_can_atlas_2018.tar.gz

# Move required files to data/
mv brca_tcga_pan_can_atlas_2018/data_*.txt .
```

### Data Structure

After download, your `data/` directory should contain:
```
data/
├── data_mrna_seq_v2_rsem.txt      (~20,000 genes × ~1,000 samples)
├── data_clinical_patient.txt       (Patient clinical data)
└── data_cna.txt                    (Copy number alteration data)
```

---

## 🎮 Usage

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
4. Generate all visualizations
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

# Step 4: Visualization
source("scripts/04_visualization.R")
# Outputs: PCA plot, heatmap

# Step 5: Survival Analysis
source("scripts/05_survival_analysis.R")
# Outputs: Prognostic genes, Kaplan-Meier curves
```

---

## 📈 Outputs

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

## 🔬 Methods Summary

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
- **Background**: All expressed genes in dataset
- **Significance**: adjusted p-value < 0.05

### Dimensionality Reduction

- **Transformation**: Variance Stabilizing Transformation (VST)
- **PCA**: Base R `prcomp()` function
- **Heatmap**: Top 50 DE genes, z-score normalized

### Survival Analysis

- **Method**: Lasso Cox proportional hazards regression
- **Tool**: glmnet v4.1-6
- **Cross-validation**: 10-fold CV
- **Lambda selection**: Lambda.min
- **Risk stratification**: Median split of risk scores

---

## 📊 Results

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
- EGFR signaling pathway
- PI3K-Akt signaling

**Downregulated in HER2+:**
- Immune response
- ECM organization
- Cell adhesion

### Survival Analysis

- **Prognostic genes selected**: XX
- **Risk groups**: High vs Low (median split)
- **Log-rank test**: p = X.XXX
- **Clinical implications**: [To be completed]

---

## 🤝 Contributing

This is an academic assignment project. However, suggestions and feedback are welcome!

If you find issues or have suggestions:
1. Open an issue describing the problem
2. Include relevant code/output
3. Suggest potential solutions

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 📧 Contact

**Caolán Maguire**

- GitHub: [@yourusername](https://github.com/yourusername)
- Email: your.email@ucd.ie
- Institution: University College Dublin

---

## 🙏 Acknowledgments

- **Data Source**: The Cancer Genome Atlas (TCGA) via cBioPortal
- **Course**: ANAT40040 - Bio Principles & Cellular Org (2025/26)
- **Instructor**: [Instructor Name]
- **Tools**: 
  - [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) - Love MI et al. (2014)
  - [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) - Yu G et al. (2012)
  - [glmnet](https://cran.r-project.org/web/packages/glmnet/index.html) - Friedman J et al. (2010)

---

## 📚 References

- Cancer Genome Atlas Network. (2012). Comprehensive molecular portraits of human breast tumours. *Nature*, 490(7418), 61-70.
- Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12), 550.
- Yu, G., Wang, L. G., Han, Y., & He, Q. Y. (2012). clusterProfiler: an R package for comparing biological themes among gene clusters. *OMICS*, 16(5), 284-287.

---

## 📝 Session Info
```r
R version 4.3.2 (2023-10-31)
Platform: x86_64-apple-darwin20 (64-bit)
Running under: macOS Sonoma 14.1

# Package versions will be added after analysis
```

---

## 🔄 Project Status

- [x] Project setup
- [x] Data download
- [x] Data preprocessing
- [ ] Differential expression analysis
- [ ] Pathway enrichment
- [ ] Visualization
- [ ] Survival analysis
- [ ] Report writing

**Last Updated**: December 2025

---

<div align="center">

Made with ❤️ for ANAT40040 Assignment 2

⭐ Star this repo if you found it useful!

</div>