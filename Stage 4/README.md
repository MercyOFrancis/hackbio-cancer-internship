# Glioma Biomarker Discovery: Integrating Machine Learning and Differential Expression Analysis

## Project Overview
This project aims to discover biomarkers in glioma by integrating machine learning techniques with differential expression analysis. The focus is on distinguishing between IDH-mutant and IDH-wildtype gliomas using gene expression data from TCGA LGG and GBM datasets.

## Authors
- Mercy Francis (@Mercylee)
- Ojiaku Confidence Chisom (@Areta)
- Manal Agdada (@Manal)
- Rahma Nabil Sallam (@rahmanabil2002)
- Pascal Onaho (@PascalOnaho)
- Hagar Haitham Elazab (@HBONH33)
- Ariyo Adesokan (@Adesokan_ariyo1)

## Repository Contents
- R script: [stage4_script.R](https://github.com/MercyOFrancis/hackbio-cancer-internship/blob/0d4d132f6167f1a0306140432d7e4aee791d20c8/Stage%204/Code/stage4_script.R)
- Figures: [Figures directory](https://github.com/MercyOFrancis/hackbio-cancer-internship/tree/0d4d132f6167f1a0306140432d7e4aee791d20c8/Stage%204/Figures)

## Methodology
1. Data preprocessing using TCGAbiolinks package in R
2. Differential expression analysis using edgeR
3. K-means clustering to analyze gene expression patterns

## Key Findings
- 9899 differentially expressed genes identified
- K-means clustering (k=4) did not distinctly separate IDH-mutant and IDH-wildtype samples
- Gene expression data alone may not provide sufficient separation between IDH-mutant and IDH-wildtype gliomas

## Future Directions
- Explore additional features to enhance clustering performance
- Consider alternative clustering methods or parameter tuning
- Investigate the biological significance of overlapping expression profiles between IDH-mutant and IDH-wildtype samples

## References
1. Louis DN, et al. (2016). The 2016 World Health Organization Classification of Tumors of the Central Nervous System: a summary. Acta Neuropathol.
2. Ceccarelli M, et al. (2016). Molecular Profiling Reveals Biologically Discrete Subsets and Pathways of Progression in Diffuse Glioma. Cell.
3. Agarwal A, et al. (2022). A comparative analysis of machine learning classifiers for predicting protein-binding nucleotides in RNA sequences. Comput Struct Biotechnol J.

## How to Use
1. Clone the repository
2. Install required R packages (TCGAbiolinks, edgeR)
3. Run the R script to reproduce the analysis
4. Refer to the Figures directory for visual results

For any questions or collaborations, please contact the authors via their GitHub profiles.
