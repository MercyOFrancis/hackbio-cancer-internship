# Integrating Machine Learning and Differential Expression Analysis for Biomarker Discovery in Glioma 🧬💻🤖

 **Authors:** Manal Agdada (@Manal), Ojiaku Confidence Chisom (@Areta),  Rahma Nabil Sallam (@rahmanabil2002), Pascal Onaho (@PascalOnaho), Hagar Haitham Elazab (@HBONH33), Mercy Francis (@Mercylee), Ariyo Adesokan (@Adesokan_ariyo1)

**R script:** https://github.com/manal-agdada/Stage-4---HackBio-internship/blob/main/Code/stage4_script.R

**Figures:** https://github.com/manal-agdada/Stage-4---HackBio-internship/tree/main/Figures

## Introduction 
Gliomas are the most common malignant brain tumors, classified into low-grade gliomas (LGG) and glioblastoma (GBM) [1]. Isocitrate dehydrogenase (IDH) mutations are key prognostic markers, as IDH-mutant gliomas have better prognoses and distinct molecular characteristics compared to IDH-wildtype tumors [2]. IDH-mutant gliomas are associated with a hypermethylation phenotype and tend to have more favorable outcomes [2]. In contrast, IDH-wildtype gliomas exhibit more aggressive behavior [2]. Therefore, methylation patterns in gliomas can serve as a valuable tool for molecular and clinical classification of gliomas [2].

## Dataset Description and preprocessing steps
For this analysis, we used TCGA LGG (534 samples) and GBM (176 samples) datasets. RNA-seq data were processed in RStudio using the ‘TCGAbiolinks’ package. After filtering samples by IDH status(wt/mut), and normalizing expression data, differential expression analysis was performed using the edgeR protocol, identifying 9899 differentially expressed genes.

![heatmap](https://github.com/user-attachments/assets/a51c0122-917a-44dd-945e-eb8dabf76f47)
![caption1](https://github.com/user-attachments/assets/2d5293b1-2bf2-47de-afc1-cd80f48f2c7d)

![volcanoplot](https://github.com/user-attachments/assets/9a4e6128-cbdf-460e-b080-e0e4bc8840fe)
![caption2](https://github.com/user-attachments/assets/dc06f511-48ce-4a57-83e1-dfbeccab70ed)

## Methodology and Results
For the next part of the task, we employed the K-means clustering algorithm, an unsupervised machine learning technique [3], to analyze a dataset comprising 667 samples of gliomas (LGG = 513, GBM = 154). The K-means algorithm was implemented in RStudio, specifying 4 clusters to classify the samples based on gene expression data. After performing the clustering, we cross-referenced the resulting clusters with the corresponding metadata indicating the IDH status of the samples (Mutant and Wild Type).  
 ![Rplot06](https://github.com/user-attachments/assets/e20d88e1-e300-422f-b884-5a308dfd37e2)

*_Figure 3. K-means clusters_*    
  
  The K-means clustering resulted in the identification of four clusters. However, upon cross-referencing with the IDH status, both Mutant and Wild Type samples were found in all four clusters. This indicates that the model did not distinctly and exclusively classify the samples into Wild Type and Mutant subgroups. Additionally, two of the four clusters overlapped, suggesting that the boundaries between clusters were not well-defined. The clustering failed to separate Mutant from Wild Type samples as anticipated.
## Conclusion
  The clustering results suggest that gene expression data alone may not provide sufficient separation between IDH-mutant and IDH-wildtype gliomas. Possible reasons for this outcome include the possible inherent similarity in the expression profiles of some Mutant and Wild Type samples or the need for additional features to enhance the clustering performance. The overlapping clusters may also indicate that the chosen number of clusters (k = 4) did not align well with the natural structure of the data, potentially requiring further tuning of the model parameters or the use of alternative clustering methods.



## References
1. Louis DN, Perry A, Reifenberger G, von Deimling A, Figarella-Branger D, Cavenee WK, Ohgaki H, Wiestler OD, Kleihues P, Ellison DW. **The 2016 World Health Organization Classification of Tumors of the Central Nervous System: a summary**. Acta Neuropathol. 2016 Jun;131(6):803-20. 
2. Ceccarelli M, Barthel FP, Malta TM, Sabedot TS, Salama SR, Murray BA, Morozova O, Newton Y, Radenbaugh A, Pagnotta SM, Anjum S, Wang J, Manyam G, Zoppoli P, Ling S, Rao AA, Grifford M, Cherniack AD, Zhang H, Poisson L, Carlotti CG Jr, Tirapelli DP, Rao A, Mikkelsen T, Lau CC, Yung WK, Rabadan R, Huse J, Brat DJ, Lehman NL, Barnholtz-Sloan JS, Zheng S, Hess K, Rao G, Meyerson M, Beroukhim R, Cooper L, Akbani R, Wrensch M, Haussler D, Aldape KD, Laird PW, Gutmann DH; TCGA Research Network; Noushmehr H, Iavarone A, Verhaak RG. **Molecular Profiling Reveals Biologically Discrete Subsets and Pathways of Progression in Diffuse Glioma**. Cell. 2016 Jan 28;164(3):550-63.
3. Agarwal A, Singh K, Kant S, Bahadur RP. **A comparative analysis of machine learning classifiers for predicting protein-binding nucleotides in RNA sequences**. Comput Struct Biotechnol J. 2022 Jun 17;20:3195-3207..
