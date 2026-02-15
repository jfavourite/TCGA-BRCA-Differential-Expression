#GSEA
#PREPARING RANKED GENE LIST

# Prepare the initial Dataframe
res_df <- as.data.frame(resLFC)
res_df <- res_df[!is.na(res_df$log2FoldChange), ]

# create a new column 'ENSEMBL' to make merging easier later
res_df$ENSEMBL <- gsub("\\..*", "", rownames(res_df))

# Translate Ensembl to Entrez 
gene_df <- bitr(
  res_df$ENSEMBL,
  fromType = "ENSEMBL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

# Merge Fold Changes with the new Entrez IDs
merged_df <- merge(res_df, gene_df, by = "ENSEMBL")

# Handling many-to-one mapping
# If multiple Ensembl IDs map to 1 Entrez ID, we take the average Log2FoldChange
# This ensures every name in our final vector is UNIQUE.
final_df <- aggregate(log2FoldChange ~ ENTREZID, data = merged_df, FUN = mean)

# Creating the Ranked List (The 'geneList' required)
gene_list <- final_df$log2FoldChange
names(gene_list) <- final_df$ENTREZID

gene_list <- sort(gene_list, decreasing = TRUE)

# adding a tiny bit of noise to break ties
set.seed(42) # This makes the 'random' noise the same every time you run it
gene_list <- gene_list + rnorm(length(gene_list), mean = 0, sd = 1e-12)

# Re-sorting one last time to ensure PERFECT decreasing order
gene_list <- sort(gene_list, decreasing = TRUE) 

# Run GSEA (Biological Process)
gsea_bp <- gseGO(
  geneList      = gene_list,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",         # Biological Process
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  verbose       = FALSE
)

gsea_results_df <- as.data.frame(gsea_bp)
head(gsea_results_df[, c("Description", "NES", "p.adjust")])

gseaplot2(gsea_bp, 
          geneSetID = 1, 
          title = gsea_bp$Description[1], 
          pvalue_table = TRUE)


# Show both Up-regulated (Positive NES) and Down-regulated (Negative NES)
dotplot(gsea_bp, showCategory = 10, split = ".sign") + 
  facet_grid(.~.sign) +
  theme_bw() +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40))

ggsave("GSEA_Dotplot_BRCA_1.png", width = 10, height = 8, dpi = 300)


ridgeplot(gsea_bp, showCategory = 15) + 
  labs(x = "Enrichment Distribution (Log2 Fold Change)") +
  theme_minimal()
ggsave("GSEA_Ridgeplot_BRCA_1.png", width = 12, height = 10, dpi = 300)


