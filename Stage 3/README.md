## COAD Biomarker Discovery: ML and Differential Expression Analysis
This repository contains the code and results for our project on integrating machine learning and differential expression analysis for biomarker discovery in Colon Adenocarcinoma (COAD).

**Project Overview**

**Objective:** Identify potential biomarkers for COAD using TCGA data

**Methods:** Differential Expression Analysis, K-Nearest Neighbor (KNN) model

**Data:** TCGA COAD dataset (RNA-seq and clinical data)

## Key Files

**R script:** [Stage3_script.R](https://github.com/MercyOFrancis/hackbio-cancer-internship/blob/main/Stage%203/Code/stage3_script.R)

**Figures:** [Figures&Visualization](https://github.com/MercyOFrancis/hackbio-cancer-internship/tree/a374ccce1e603791a8c4cab8ccc2a9d7e8a2c42f/Stage%203/Figures%20%26%20Visualization)


## Results

Identified 752 upregulated and 371 downregulated genes

KNN model showed limitations, potentially due to small sample size

Key genes identified: RPS4Y1, APOH, CYP1A1, GSTM1, PF4V1

## Analysis Chain

**Data Preprocessing**

Download and prepare TCGA COAD dataset

Normalize RNA-seq data and filter lowly expressed genes


**Differential Expression Analysis**

Identify DEGs between groups (adjusted p-value < 0.05, |log fold change| > 1)


**Enrichment Analysis**

Perform GO term enrichment on DEGs


**Machine Learning**

Apply KNN model for gender prediction and biomarker identification


**Visualization**

Generate heatmaps, volcano plots, and enrichment bar plots



**How to Use**

Clone this repository
Install required R packages
Run the R script to reproduce the analysis

For more detailed information, please refer to the full report in this repository.
