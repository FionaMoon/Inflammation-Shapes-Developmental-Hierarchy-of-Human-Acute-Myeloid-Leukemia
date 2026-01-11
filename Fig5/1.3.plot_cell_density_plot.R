## plot cell density plot
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
F_path = "/Users/MyFolder/Code/LLN-RBM47/CITE-seq_cytokine/ZYG_script/V3"

source("/Users/MyFolder/Code/LLN-RBM47/CITE-seq_cytokine/ZYG_script/V3/FACS_function_Dimplot.R")


obj <- readRDS("/Users/MyFolder/Code/LLN-RBM47/CITE-seq_cytokine/ZYG/ZYG_PeterVG_anno_cytotrace_ADT_V5.rds")

obj$PeterVG_anno_broad <- as.character(obj$PeterVG_anno)

obj$PeterVG_anno_broad[grep(pattern = "B|NK|T", obj$PeterVG_anno_broad)] <- "Other"

Idents(obj) <- "PeterVG_anno_broad"
obj <- subset(obj, idents = "Other", invert = T)

obj$condition <- factor(as.character(obj$condition), levels = c("CT", "IL1B", "TNFA", "Combo"))

Idents(obj) <- "PeterVG_anno"
obj <- subset(obj, downsample = 500)

obj_sub <- SplitObject(obj, split.by = "condition")

plot_list <- lapply(
  names(obj_sub),
  function(cond) {
    galaxy_DimPlot_Seurat(
      seurat_obj = obj_sub[[cond]],
      reduction = "umap"
    ) + ggtitle(cond)
  }
)

names(plot_list) <- names(obj_sub)
plot_list <- plot_list[c("CT", "IL1B", "TNFA", "Combo")]

pdf(paste0(F_path, "/UMAP_by_condition.pdf"), width = 16, height = 3.8)
  wrap_plots(plot_list, ncol = 4)
dev.off()
