#GO FUNCTIONAL ENRICHMENT

sig_genes <- res_df[
  res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,
]

nrow(sig_genes)

#Converting ENSEMBL → ENTREZ IDs

ensembl_ids <- gsub("\\..*", "", rownames(sig_genes))

entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = ensembl_ids,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

entrez_ids <- na.omit(entrez_ids)
View(entrez_ids)

#Running GO Biological Enrichment function

ego_bp <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

head(ego_bp)

#Visualization

dotplot(ego_bp, showCategory = 15) +
  ggtitle("GO Biological Process Enrichment: TCGA-BRCA Tumor vs Normal")

ggsave(
  "TCGA_BRCA_GO_BP_Dotplot.png",
  width = 8,
  height = 6
)
