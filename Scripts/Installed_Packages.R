sessionInfo()
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

install.packages(c(
  "igraph", "reshape2", "cowplot", "fastmatch", "scatterpie", 
  "ggnewscale", "yulab.utils", "gson", "aplot", "ggfun", 
  "ggtangle", "tidydr", "ape", "ggiraph"
))

BiocManager::install(
  c(
    "TCGAbiolinks",
    "SummarizedExperiment",
    "DESeq2",
    "edgeR",
    "limma",
    "EnhancedVolcano",
    "pheatmap",
    "apeglm",
    "clusterProfiler",
    "org.Hs.eg.db",
    "enrichplot", 
    "DOSE",
    "fgsea",
    "GOSemSim",
    "qvalue",
    "ggtree",
    "treeio",
    "GO.db",
    "tidytree"
  ),
  ask = FALSE,
  update = FALSE
)

#Loading the various libraries

library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(apeglm)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(stringr)
