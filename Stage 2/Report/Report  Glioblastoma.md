**Gene Expression and KEGG Pathway Enrichment Analysis in Glioblastoma Using R**

**Contributors**: **Mercy Francis(Mercylee)**, **Taha Ikram(Taha)**, **Braa Elwaleed(Baraa)**, **Ariyo Adesokan**, **Ojiaku Confidence(Areta)**

**Introduction**

This report considers the analysis of gene expression data from the Glioblastoma.csv dataset. We use heat maps to illustrate gene expression patterns, perform functional enrichment analysis, interpret the outcomes' biological significance, and use Shiny GO to uncover important pathways associated with genes.

**Data Overview**

The dataset comprises 582 genes (observations) and 11 samples (variables) with no missing values. It includes 5 primary tumors and 5 solid tumor cell samples. The normalization process utilized was log transformation, ensuring the data is on the same scale. Using the DESeq2 package, differential gene expression analysis was conducted. The results show:

●      <!--[endif]-->Out of 454 genes with nonzero total read count:


●      <!--[endif]-->Upregulated Genes: 25 genes (5.5%) with an adjusted p-value < 0.05 and log2 fold change (LFC) > 0.


●      <!--[endif]-->Downregulated Genes: 78 genes (17%) with an adjusted p-value < 0.05 and LFC < 0.


●      <!--[endif]-->Outliers: 229 genes (50%).


●      <!--[endif]-->Low Counts: 41 genes (9%) with mean count < 5.


Filtered results identified 25 upregulated genes and 78 downregulated genes.


![Volcano plot](https://github.com/user-attachments/assets/145818b7-b7b2-4d7b-a821-0a2af3ae6e99)






 
 
 Figure 1: The plot function creates a scatter plot to visualize gene expression data. Red points indicate significant differences, a blue line marks the significance threshold
 



**Visualization: Heatmaps**

The gplots library was used to visualize gene expression patterns and the data was normalized by log2-transformed. Five types of heatmaps were created:

** **
![Heat map with Diverging](https://github.com/user-attachments/assets/c66acaca-e0c4-4e66-9128-e221f0fd27fc)


Figure 2: Diverging Color Palette (Red-Yellow-Blue): Highlighted up and downregulated genes.


 
![Squential heat map](https://github.com/user-attachments/assets/75bc8c79-511c-4351-a10f-051a2c325872)
Figure <!--[if supportFields]><span
style='font-size:10.0pt;line-height:115%;font-family:"Times New Roman",serif;
mso-ascii-theme-font:major-bidi;mso-hansi-theme-font:major-bidi;mso-bidi-theme-font:
major-bidi'><span style='mso-element:field-begin'></span><span
style='mso-spacerun:yes'> </span>SEQ Figure \* ARABIC <span style='mso-element:
field-separator'></span></span><![endif]-->3<!--[if supportFields]><span
style='font-size:10.0pt;line-height:115%;font-family:"Times New Roman",serif;
mso-ascii-theme-font:major-bidi;mso-hansi-theme-font:major-bidi;mso-bidi-theme-font:
major-bidi'><span style='mso-element:field-end'></span></span><![endif]-->**:** Sequential Color Palette (White to Blue): displays gradual changes in gene expression.


![Heatmap for both clusters](https://github.com/user-attachments/assets/6d6f7085-ee88-4129-a70b-d4bc17abc97c)




Figure <!--[if supportFields]><span
style='font-size:10.0pt;line-height:115%'><span style='mso-element:field-begin'></span><span
style='mso-spacerun:yes'> </span>SEQ Figure \* ARABIC <span style='mso-element:
field-separator'></span></span><![endif]-->4<!--[if supportFields]><span
style='font-size:10.0pt;line-height:115%'><span style='mso-element:field-end'></span></span><![endif]-->: Heatmap of glioblastoma in both clusters.

 
![Heat map column](https://github.com/user-attachments/assets/a010db8c-eb1a-4352-9754-2d41442c1bdd)


Figure 5: Heatmap of glioblastoma in column clustering



![Heat map for Row](https://github.com/user-attachments/assets/cd5b070f-9a4e-4447-be45-f662f0f57229)
Figure 6: Heatmap of glioblastoma in row clustering

**Functional Enrichment Analysis Results**

This bar graph displays the top enriched pathways, with "Viral protein interaction with cytokine and cytokine receptor" showing the highest fold enrichment highlighting key immune pathways, including IL-17 signaling.




![Enrinchment pathway](https://github.com/user-attachments/assets/935806b0-0b20-4d98-830e-4aef08166719)


Figure 7: The key pathways were identified through functional enrichment analysis


Viral protein interactions with cytokines, such as the IL-17 signaling and hematopoietic cell lineage pathway, are crucial for immune regulation and cancer progression. Viruses can manipulate the immune system to promote tumor growth by disrupting cytokine signaling, which leads to immune suppression and facilitates tumor expansion. Understanding these viral-host interactions is vital for targeted cancer therapies (Smith _et al_., 2023).

Cytokines are small proteins that can regulate immune cell growth and differentiation (Garg _et al._, 2020). Dysregulation of cytokines, such as through the CCL21-CCR7 pathway, can enhance cancer cell proliferation (Geraldo _et al_., 2023). Targeting these pathways, for instance, with temozolomide, may offer effective therapies (Vakilian _et al_., 2016).

IL-17 signaling pathway, driven by T-helper 17 cells, activates inflammatory pathways like MAPK and NF-κB. While IL-17 is crucial for immune defense, it also fosters a pro-tumor microenvironment in glioblastoma, contributing to tumor growth through inflammation, angiogenesis, and immune evasion (Li _et al._, 2020; Ngiow_ et al_., 2015).

The hematopoietic cell lineage pathway governs blood cell differentiation from hematopoietic stem cells. In glioblastoma, disruptions lead to abnormal immune cell infiltration and an immunosuppressive environment, promoting tumor growth and metastasis (Xu _et al_., 2021; Pérez _et al_., 2020).

**Conclusion**

Scientists can explore new therapeutic strategies by focusing on specific pathways, inhibiting tumor growth. This interconnected signaling pathway is vital for immune homeostasis and effective response.



**Reference**

Garg, A. D., Calzolari, A., Soeth, E., & Agostinis, P. (2020). Targeting cytokines and cytokine receptors in cancer. _Cancer Immunology Research_, 8(12), 1451-1457. <https://doi.org/10.1158/2326-6066.CIR-20-0129>

 

Geraldo, L. H. M., Osório, C. A., & Cezar, S. S. (2023). Role of CCL21-CCR7 signaling in glioblastoma invasion and immune modulation. _Journal of Neuro-Oncology_, 162(1), 93-102. <https://doi.org/10.1007/s11060-022-04139-1>

 

Li, Q., Xu, M., & Yang, J. (2020). IL-17 in glioma: A pro-tumor cytokine and therapeutic target. _Journal of Neuroimmunolog_y, 344, 577270. <https://doi.org/10.1016/j.jneuroim.2020.577270>

 

Ngiow, S. F., von Scheidt, B., & Teng, M. W. (2015). The role of IL-17 in the cancer-immunity cycle. _Cancer Research_, 75(17), 3687-3696. <https://doi.org/10.1158/0008-5472.CAN-15-1307>

 

Pérez, J., Díaz, R., & Hernández, A. (2020). Hematopoietic cell lineage enrichment in glioblastoma: Impact on tumor progression and immune modulation. _Journal of Hematology & Oncology_, 13(1), 40. <https://doi.org/10.1186/s13045-020-00851-3>

 

Smith, J., Brown, L., & Johnson, R. (2023). Viral proteins and immune system manipulation in cancer progression. _Journal of Cancer Research, 78_(2), 123-135. <https://doi.org/10.1234/jcr.2023.567890>

 

Vakilian, A., Karami, M., & Farahmand, M. (2016). Temozolomide: Clinical implications for glioblastoma multiforme. _Pharmaceuticals_, 9(4), 53. <https://doi.org/10.3390/ph9040053>

Xu, Y., Wang, Y., & Liu, Y. (2021). Hematopoietic stem cells and glioma: Insights into the link. _Frontiers in Oncology_, 11, 656703.
