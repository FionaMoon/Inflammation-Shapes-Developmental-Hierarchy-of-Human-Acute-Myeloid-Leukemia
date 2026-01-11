## Celltracer2 for stem score

rm(list = ls())
gc()

## ---- load libraries ----
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(CytoTRACE2)
library(ggthemes)

## ---- set paths ----
D_path = "/Users/MyFolder/Code/LLN-RBM47/CITE-seq_cytokine"
F_path = "/Users/MyFolder/Code/LLN-RBM47/CITE-seq_cytokine/ZYG_script/V3"


## ---- load Seurat object ----
seurat_obj <- readRDS("/Users/MyFolder/Code/LLN-RBM47/CITE-seq_cytokine/ZYG/ZYG_PeterVG_anno_cytotrace_ADT_V5.rds")
seurat_obj

## ---- run CytoTRACE2 (human, Seurat) ----
set.seed(14)

exp <- GetAssayData(seurat_obj, layer = "data")
cytotrace2_result <- cytotrace2(
    exp,
    species = "human",           # human species
    # is_seurat = TRUE,            # input is Seurat object
    # slot_type = "counts",        # use raw counts (or "data" for normalized)
    batch_size = 10000,
    smooth_batch_size = 1000,
    parallelize_models = TRUE,
    parallelize_smoothing = TRUE,
    seed = 14,
    ncores = 1
)

head(cytotrace2_result)
if(all(rownames(seurat_obj@meta.data) == rownames(cytotrace2_result))){
  seurat_obj@meta.data <- cbind(seurat_obj@meta.data, cytotrace2_result)
} else {
  stop("\n Not euqal\n")
}
  

p1 <- FeaturePlot(
  seurat_obj,
  features = "CytoTRACE2_Score",
  reduction = "umap",
  split.by = "condition",
  ncol = 4
) & scale_colour_gradientn(colors = c("blue", "skyblue","lightblue","gray", "pink", "magenta", "red")) & 
  theme_few(base_size = 14)

ggsave(file.path(F_path, "CytoTRACE2_Score_UMAP.pdf"), p1, width = 19, height = 4)


# generate prediction and phenotype association plots with plotData function
plots <- plotData(cytotrace2_result = seurat_obj, 
                  is_seurat = T
                  )
p1 <- plots$CytoTRACE2_Potency_UMAP
p2 <- plots$CytoTRACE2_UMAP
p3 <- plots$CytoTRACE2_Relative_UMAP
ggsave(paste0(F_path, "/cytotrace.pdf"), p1+p2+p3, width = 10,height = 3)

