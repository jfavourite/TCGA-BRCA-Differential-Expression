sessionInfo()
#if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
#}

#BiocManager::install(
#  c(
    "TCGAbiolinks",
    "SummarizedExperiment",
    "DESeq2",
    "edgeR",
    "limma",
    "EnhancedVolcano",
    "pheatmap",
    "apeglm"
  ),
  ask = FALSE,
  update = FALSE
)
#

#Loading the various libraries

library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(apeglm)
library(pheatmap)

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

#Insecting to check the query
query

table(getResults(query)$sample_type)


#Selecting subsets of the whole files for reproducability and control of the wor environment

results <- getResults(query)

tumor_ids <- results$cases[results$sample_type == "Primary Tumor"][1:30]
normal_ids <- results$cases[results$sample_type == "Solid Tissue Normal"][1:30]

selected_barcodes <- c(tumor_ids, normal_ids)

length(selected_barcodes)

#Requery only selected samples(subsets)

query_subset <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode =  selected_barcodes
)

query_subset

#Downloading files 

GDCdownload(query_subset)

#Preparing the data

brca_se <- GDCprepare(query_subset)
brca_se

#Extracting raw counts
counts <- assays(brca_se)$unstranded
dim(counts)

#Building the sample information and defining 'normal' and 'Tumor'

sample_info <- colData(brca_se)

sample_info$group <- ifelse(
  sample_info$sample_type == "Primary Tumor",
  "Tumor",
  "Normal"
)

sample_info$group <- factor(sample_info$group, levels = c("Normal", "Tumor"))

table(sample_info$sample_type, sample_info$group)


#Creating the DESeq2 dataset for statistics

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sample_info,
  design = ~ group
)

#Filtering low expression genes(Noise)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dim(dds)

#Running a DESeq Diferential expression

dds <- DESeq(dds)

#Extracting the results

res <- results(dds, contrast = c("group", "Tumor", "Normal"))
View(res)
summary(res)

#Shrinking logfolds changes because LFCs exaggerate low count genes


resLFC <- lfcShrink(
  dds,
  coef = "group_Tumor_vs_Normal",
  type = "apeglm"
)

head(resLFC[order(resLFC$padj), ])


#Making a volcanoe plot

# 1. Get all names from the original object
all_names <- rowData(brca_se)$gene_name

# 2. Make sure the names are 'named' with their ENSG IDs so R can match them
names(all_names) <- rownames(brca_se)

# 3. Filter those names to match only the genes in your results
matched_labels <- all_names[rownames(resLFC)]

# 4. Now running EnhancedVolcano with the matched labels
EnhancedVolcano(
  resLFC,
  lab = matched_labels, # Use the matched labels here!
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2,
  labSize = 3,
  title = "TCGA-BRCA: Tumor vs Normal",
  subtitle = "RNA-seq Differential Expression",
  caption = "30 Tumor vs 30 Normal (TCGA-BRCA)"
)

 # Making a gggplot duplicate of the enhancedvolcanoe 
res_df <- as.data.frame(resLFC) #changing to dataframe

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 1), alpha = 0.5) +
  scale_color_manual(values = c("black", "red")) + # Red for significant
  theme_minimal() +
  labs(title = "Volcano Plot (ggplot2)",
       x = "Log2 Fold Change",
       y = "-log10 Adjusted P-value") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") 
#save image as pdf and png in the plots section. 

#Selecting the high confidence genes

res_df <- as.data.frame(resLFC)
res_ordered <- res_df[order(res_df$padj), ]


res_ordered$Symbol <- all_names[rownames(res_ordered)]

# View the Top 10 most significant genes
head(res_ordered, 10)

View(res_ordered)
#Save the table as a CSV for Excel
write.csv(res_ordered, "TCGA_BRCA_Top_Genes.csv")


#Creating an heatmap for the top 20 genes
top_genes <- rownames(res_ordered)[1:20]

vsd <- vst(dds, blind=FALSE)
heatmap_matrix <- assay(vsd)[top_genes, ]

rownames(heatmap_matrix) <- res_ordered$Symbol[1:20]

annotation_data <- as.data.frame(colData(dds)[, c("group")])
colnames(annotation_data) <- "Status"
rownames(annotation_data) <- colnames(heatmap_matrix)

pheatmap(heatmap_matrix, 
         annotation_col = annotation_data, 
         show_colnames = FALSE,         
         scale = "row",                 
         clustering_distance_rows = "euclidean", 
         main = "Top 20 Differentially Expressed Genes in TCGA-BRCA",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         filename = "TCGA_BRCA_Heatmap.png", 
         width = 8, height = 10)



