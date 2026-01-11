## galaxy-plot for ADT

rm(list = ls())
gc()

## load packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(scCustomize)
library(RColorBrewer)
library(ggthemes)
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 10 * 1e9)

## paths
D_path <- "/Users/MyFolder/Code/LLN-RBM47/CITE-seq_cytokine"
F_path <- "/Users/MyFolder/Code/LLN-RBM47/CITE-seq_cytokine/ZYG_script/V3"

## load object
SeuratObj <- readRDS(
  paste0(D_path, "/ZYG/ZYG_PeterVG_anno_cytotrace_ADT_V5.rds")
)

galaxy_ADT_FeaturePlot <- function(
  seurat_obj,
  feature,
  assay        = "ADT",
  layer        = "data",
  cutoff_mult  = 2,
  cols         = c("lightgrey", "#e89a73ff"),
  pt.size      = 0.05,
  alpha        = 0.8,
  border_width = 1.0
) {

  ## ---------------------------
  ## Safety checks
  ## ---------------------------
  stopifnot(feature %in% rownames(seurat_obj[[assay]]))

  DefaultAssay(seurat_obj) <- assay

  ## ---------------------------
  ## Fetch expression
  ## ---------------------------
  expr <- FetchData(
    object = seurat_obj,
    vars   = feature,
    assay  = assay,
    layer  = layer
  )

  ## ---------------------------
  ## Compute threshold
  ## ---------------------------
  q25 <- quantile(expr[[feature]], 0.25, na.rm = TRUE)
  q75 <- quantile(expr[[feature]], 0.75, na.rm = TRUE)
  med <- median(expr[[feature]], na.rm = TRUE)

  threshold_lower <- 0.25 * q25 + 0.25 * med + 0.5 * q75

  if(threshold_lower == 0){
    threshold_lower = NA
  }

  ## ---------------------------
  ## Plot
  ## ---------------------------
  p <- FeaturePlot(
    object      = seurat_obj,
    features    = feature,
    keep.scale  = "feature",
    slot        = layer,
    cols        = cols,
    min.cutoff  = threshold_lower,
    max.cutoff  = cutoff_mult * threshold_lower,
    order       = TRUE,
    alpha       = alpha,
    pt.size     = pt.size
  ) &
    theme_few() &
    theme(
      panel.border = element_rect(
        colour   = "black",
        fill     = NA,
        linewidth = border_width
      ),
      axis.text  = element_blank(),
      axis.ticks = element_blank()
    )

  ## ---------------------------
  ## Return
  ## ---------------------------
  return(list(
    plot      = p,
    threshold = threshold_lower,
    expr_rng  = range(expr[[feature]], na.rm = TRUE)
  ))
}

res_0 <- galaxy_ADT_FeaturePlot(
  seurat_obj = SeuratObj,
  feature    = "TotalSeq-CD34",
  assay        = "ADT",
  layer        = "data",
  cutoff_mult  = 2,
  cols         = c("lightgrey", "#516be2ff"),
  pt.size      = 0.1,
  alpha        = 0.8,
  border_width = 1.0
)


res_1 <- galaxy_ADT_FeaturePlot(
  seurat_obj = SeuratObj,
  feature    = "MECOM",
  assay        = "RNA",
  layer        = "data",
  cutoff_mult  = 2,
  cols         = c("lightgrey", "#516be2ff"),
  pt.size      = 0.1,
  alpha        = 0.8,
  border_width = 1.0
)


res_2 <- galaxy_ADT_FeaturePlot(
  seurat_obj = SeuratObj,
  feature    = "TotalSeq-CD14",
  assay        = "ADT",
  layer        = "data",
  cutoff_mult  = 2,
  cols         = c("lightgrey", "#e26b51ff"),
  pt.size      = 0.1,
  alpha        = 0.8,
  border_width = 1.0
)


res_3 <- galaxy_ADT_FeaturePlot(
  seurat_obj = SeuratObj,
  feature    = "CD14",
  assay        = "RNA",
  layer        = "data",
  cutoff_mult  = 2,
  cols         = c("lightgrey", "#e26b51ff"),
  pt.size      = 0.1,
  alpha        = 0.8,
  border_width = 1.0
)

res_4 <- galaxy_ADT_FeaturePlot(
  seurat_obj = SeuratObj,
  feature    = "TotalSeq-CD11b",
  assay        = "ADT",
  layer        = "data",
  cutoff_mult  = 2,
  cols         = c("lightgrey", "#e28651ff"),
  pt.size      = 0.1,
  alpha        = 0.8,
  border_width = 1.0
)


res_5 <- galaxy_ADT_FeaturePlot(
  seurat_obj = SeuratObj,
  feature    = "ITGAM",
  assay        = "RNA",
  layer        = "data",
  cutoff_mult  = 2,
  cols         = c("lightgrey", "#d5924aff"),
  pt.size      = 0.1,
  alpha        = 0.8,
  border_width = 1.0
)


res_6 <- galaxy_ADT_FeaturePlot(
  seurat_obj = SeuratObj,
  feature    = "TotalSeq-CD117",
  assay        = "ADT",
  layer        = "data",
  cutoff_mult  = 2,
  cols         = c("lightgrey", "#7461b5ff"),
  pt.size      = 0.1,
  alpha        = 0.8,
  border_width = 1.0
)


res_7 <- galaxy_ADT_FeaturePlot(
  seurat_obj = SeuratObj,
  feature    = "KIT",
  assay        = "RNA",
  layer        = "data",
  cutoff_mult  = 2,
  cols         = c("lightgrey", "#7461b5ff"),
  pt.size      = 0.1,
  alpha        = 0.8,
  border_width = 1.0
)

pdf(paste0(F_path, "/Featureplot.pdf"), width = 3.8, height = 3.2)
  res_0$plot # CD34
  res_1$plot
  res_2$plot # CD14
  res_3$plot
  res_4$plot # CD1b
  res_5$plot
  res_6$plot # CD1b
  res_7$plot
dev.off()
