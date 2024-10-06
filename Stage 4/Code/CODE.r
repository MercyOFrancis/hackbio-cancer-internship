######################################
#### HackBio internship - stage 4 #### 
######################################

# Biomarkers analysis #### 

## load necessary libraries 

library("TCGAbiolinks")
library(SummarizedExperiment)
library(edgeR)
library(gplots)
library(ggplot2)
library(biomaRt)
library(ggrepel)

# project information 
getProjectSummary("TCGA-LGG")
getProjectSummary("TCGA-GBM")

# download LGG dataset 
tcga_lgg <- GDCquery(project = "TCGA-LGG",
                     data.category = "Transcriptome Profiling",
                     data.type = "Gene Expression Quantification")
GDCdownload(tcga_lgg) 
lgg_data <- GDCprepare(tcga_lgg) 

# download GBM dataset
tcga_gbm <- GDCquery(project = "TCGA-GBM",
                     data.category = "Transcriptome Profiling",
                     data.type = "Gene Expression Quantification")
GDCdownload(tcga_gbm) 
gbm_data <- GDCprepare(tcga_gbm) 

# get raw counts from lgg_data
lgg_rawdata <- assays(lgg_data) 
dim(lgg_rawdata$unstranded) # 534 lgg samples

# get raw counts from gbm_data
gbm_rawdata <- assays(gbm_data) 
dim(gbm_rawdata$unstranded) # 176 gbm samples

# prepare metadata for lgg
lgg_data$barcode
lgg_data$paper_IDH.status
table(lgg_data$paper_IDH.status) # 94 wt and 419 mut

meta_lgg <- data.frame("barcode" = lgg_data$barcode,
                        "IDH status" = lgg_data$paper_IDH.status)
View(meta_lgg)
table(is.na(meta_lgg)) # 21 NAs -> eliminate samples with no IDH status
meta_lgg <- na.omit(meta_lgg)

# prepare metadata for gbm 
gbm_data$barcode
gbm_data$paper_IDH.status
table(gbm_data$paper_IDH.status) # 143 wt and 11 mut 

meta_gbm <- data.frame("barcode" = gbm_data$barcode,
                       "IDH status" = gbm_data$paper_IDH.status)
View(meta_gbm)
table(is.na(meta_gbm)) # 22 NAs -> eliminate samples with no IDH status
meta_gbm <- na.omit(meta_gbm)

# final dataset for lgg (513 samples)
lgg <- lgg_rawdata$unstranded
samples_lgg <- meta_lgg$barcode
lgg <- lgg[, colnames(lgg) %in% samples_lgg]
all(colnames(lgg) %in% samples_lgg) # TRUE

# final dataset for gbm (154 samples)
gbm <- gbm_rawdata$unstranded
samples_gbm <- meta_gbm$barcode
gbm <- gbm[, colnames(gbm) %in% samples_gbm]
all(colnames(gbm) %in% samples_gbm)

# merge the two datasets (both counts and metadata)
glioma <- cbind(lgg, gbm)
meta_glioma <- rbind(meta_lgg, meta_gbm)

# export the final dfs
write.csv(glioma, "glioma_rawcounts.csv", row.names = TRUE)
write.csv(meta_glioma, "glioma_metadata.csv", row.names = FALSE)

# clean and preprocess the data
table(is.na(glioma)) # no NAs
glioma <- TCGAanalyze_Normalization(tabDF = glioma, geneInfo = geneInfoHT, method = "geneLength")
glioma <- TCGAanalyze_Normalization(tabDF = glioma, geneInfo = geneInfoHT, method = "gcContent")
glioma <- TCGAanalyze_Filtering(tabDF = glioma,
                             method = "quantile",
                             qnt.cut = 0.25)

## differential expression analysis 
group.wt <- meta_glioma$barcode[meta_glioma$IDH.status == "WT"]
group.mut <- meta_glioma$barcode[meta_glioma$IDH.status == "Mutant"]

DEA <- TCGAanalyze_DEA(mat1 = glioma[ , group.wt], 
                       mat2 = glioma[ , group.mut],
                       Cond1type = "WT",
                       Cond2type = "Mutant",
                       pipeline = "edgeR")
DEA.Level <- 
  TCGAanalyze_LevelTab(DEA, "WT", "Mutant",
                       glioma[ , group.wt],
                       glioma[ , group.mut])

# annotation
gene_names <- rownames(DEA.Level)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
annot <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = gene_names,
  mart = mart)
gene_with_symbols <- merge(DEA.Level, annot, 
                                by.x = "mRNA", 
                                by.y = "ensembl_gene_id", 
                                all.x = TRUE)

# volcano plot
genes_to_label <- subset(gene_with_symbols, abs(logFC) > 5 & FDR < 0.05)

ggplot(DEA.Level, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = ifelse(FDR < 0.05 & abs(logFC) > 1,
                                ifelse(logFC > 1, "Upregulated", "Downregulated"), "Not significant")), size = 2) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") + 
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "black") + 
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not significant" = "gray")) +
  theme_minimal() +
  labs(x = "Log2 Fold Change", y = "-Log10 FDR", title = "Volcano Plot of WT vs Mut", color = "Gene Regulation") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14), 
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "top") +
  ggrepel::geom_label_repel(data = genes_to_label, aes(label = hgnc_symbol))


# heatmap of DEGs
# top differentially expressed genes 
DEGs <- subset(DEA.Level, abs(logFC) > 3 & FDR < 0.01) # 905 DEGs   
heat.DEGs <- glioma[rownames(DEGs), ]

# color code
ccodes <- ifelse(colnames(heat.DEGs) %in% group.wt, "red", "blue") # wt in red, mut in blue

# heatmap
red_green_palette <- colorRampPalette(c("green", "red"))(10)
heatmap.2(x = as.matrix(heat.DEGs),
          col = red_green_palette,
          Rowv = F, Colv = T,
          scale = 'row',
          sepcolor = 'black',
          trace = "none",
          key = TRUE,
          dendrogram = "col",
          cexRow = 0.9, cexCol = 0.7,
          main = "Heatmap of WT vs Mut",
          na.color = 'black',
          ColSideColors = ccodes,
          margins = c(11,10))
legend("left",                       
       legend = c("WT", "Mut"),     
       col = c("red", "blue"),           
       lty = 1,                          
       lwd = 4,                         
       cex = 0.7,
       xpd = TRUE,
       inset = c(-0.1, 0))
