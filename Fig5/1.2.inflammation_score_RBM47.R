rm(list = ls())
gc()

library(Seurat)
library(escape)
library(ggplot2)
library(clusterProfiler)
library(ggprism)
library(dplyr)
library(GSVA)
library(data.table)

D_path = "/Users/MyFolder/Code/LLN-RBM47/CITE-seq_cytokine"
F_path = "/Users/MyFolder/Code/LLN-RBM47/CITE-seq_cytokine/AML-8_script/V3"


## inflammation geneset
gmt_folder <- "/Users/MyFolder/Code/LLN-RBM47/CITE-seq_cytokine/AML-8_script/V3/gmt"
gmt_files <- list.files(gmt_folder, pattern = "\\.gmt$", full.names = TRUE)

gene_sets_list <- list()

for (f in gmt_files) {
  gmt_df <- read.gmt(f)
  
  # convert to named list
  gs_list <- split(gmt_df$gene, gmt_df$term)
  
  # merge
  gene_sets_list <- c(gene_sets_list, gs_list)
}

str(gene_sets_list)

## scRNA-seq data
sc <- readRDS("/Users/MyFolder/Code/LLN-RBM47/CITE-seq_cytokine/AML-8/AML-8_PeterVG_anno_cytotrace_ADT_V5.rds")

sc


## AddModuleScore
sc_feature <- intersect(row.names(sc), gene_sets_list$GOBP_INFLAMMATORY_RESPONSE)
cd_features <- list(GOBP_INFLAMMATORY_RESPONSE = sc_feature)

GS.hallmark <- list(GOBP_INFLAMMATORY_RESPONSE = sc_feature,
                    HALLMARK_INFLAMMATORY_RESPONSE = gene_sets_list$GOBP_INFLAMMATORY_RESPONSE)

if (T) {
    sc <- AddModuleScore(
    object = sc,
    features = cd_features,
    name = 'GOBP_INFLAMMATORY_RESPONSE'
    )
    head(x = sc[[]])
}

p2 <- FeaturePlot(sc, reduction = "umap", features = "GOBP_INFLAMMATORY_RESPONSE1", 
            order = TRUE, pt.size = 0.4, max.cutoff = NA) +
  scale_colour_gradientn(colors = c("blue", "skyblue","gray", "pink", "magenta", "red")) + 
  ggtitle('GOBP_INFLAMMATORY_RESPONSE') +
      theme_few(base_size = 14)

p3 <- FeaturePlot(sc, reduction = "umap", features = "RBM47", 
                  order = TRUE, pt.size = 0.4, max.cutoff = NA) +
  scale_colour_gradientn(colors = c("blue", "skyblue","gray", "pink", "magenta", "red")) + 
  ggtitle('RBM47') +
  theme_few(base_size = 14)

ggsave(filename = paste0(F_path, "/sc_InflammationScore.pdf"), 
p2 + p3, width = 4.5, height = 7.5)
