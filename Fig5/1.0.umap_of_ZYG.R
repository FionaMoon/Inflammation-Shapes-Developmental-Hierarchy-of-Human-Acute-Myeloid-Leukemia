## step1: change color of celltypes

## after clustering_2
rm(list=ls())
gc() ## release memory

## load packages
library(Seurat)
library(dplyr)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(cowplot)
library(progress)
renv::activate(project = "/Users/SCP_env/renv")
library(SCP)
# library(scplotter) # plotthis
options(Seurat.object.assay.version = "v5")
# needs to be set for large dataset analysis
options(future.globals.maxSize = 10 * 1e9)


## add cell hashing tag and merge annotation
D_path = "/Users/MyFolder/Code/LLN-RBM47/CITE-seq_cytokine"
F_path = "/Users/MyFolder/Code/LLN-RBM47/CITE-seq_cytokine/ZYG_script/V3"

SeuratObj <- readRDS("/Users/MyFolder/Code/LLN-RBM47/CITE-seq_cytokine/ZYG/ZYG_PeterVG_anno_cytotrace_ADT_V5.rds")

head(SeuratObj)
table(is.na(SeuratObj$PeterVG_anno))


SeuratObj$condition <- factor(as.character(SeuratObj$condition),
levels = c("CT", "IL1B", "TNFA", "Combo"))



	# Your desired ordering of cell types:
	celltype_order <- c(
  	"LSPC-Quiescent","LSPC-Primed", "LSPC-Cycle", "MEP",
  	"GMP-like", "ProMono-like", "Mono-like", "cDC-like", "pDC-like",
  	"B", "pre/pro−B",
 	"NK", "T cell"
	)

	# Make sure your metadata is a factor with that order:
	unique(SeuratObj$PeterVG_anno) %in% celltype_order
	SeuratObj$PeterVG_anno <- factor(SeuratObj$PeterVG_anno, levels = celltype_order)

	base_cols <- c("#454fd2e8", "#698ff7ff", "#6bbbd8ff","#b5cf23c0",
  	"#ded232ff","#cd57b5ff","#cf465bc4", "#ed977aff", "#e86a40ff",
  	"#84C868", "#00c456",
	"#a588b1", "#965ac1ff")

	names(celltype_order) <- base_cols

	# Then plot with SCP::CellDimPlot using your custom colors
pdf(paste0(F_path, "/annoplot_split_withlym.pdf"), width = 12, height = 10)
	CellDimPlot(
  	SeuratObj,
  	group.by = "PeterVG_anno",
  	reduction = "umap",
  	split.by = "condition",
  	palcolor = base_cols,
  	bg_color = "grey90",
			pt.alpha = 0.8,
	pt.size = 0.2
	)
dev.off()

pdf(paste0(F_path, "/annoplot_all_withlym.pdf"), width = 5, height = 4)
	CellDimPlot(
	  SeuratObj,
	  group.by = "PeterVG_anno",
	  reduction = "umap",
	  palcolor = base_cols,
	  bg_color = "grey90",
			pt.alpha = 0.8,
	pt.size = 0.2
	)
dev.off()

################### make lym as others
SeuratObj$PeterVG_anno_mye <- as.character(SeuratObj$PeterVG_anno)
SeuratObj$PeterVG_anno_mye[grep(pattern = "B|T|NK", SeuratObj$PeterVG_anno_mye)] <- "Lym"

table(SeuratObj$PeterVG_anno_mye)

# Your desired ordering of cell types:
celltype_order <- c(
  	"LSPC-Quiescent", "LSPC-Primed", "LSPC-Cycle", "MEP",
  	"GMP-like", "ProMono-like", "Mono-like", "cDC-like", "pDC-like",
  	"Lym"
	)

	# Make sure your metadata is a factor with that order:
unique(SeuratObj$PeterVG_anno_mye) %in% celltype_order
SeuratObj$PeterVG_anno_mye <- factor(SeuratObj$PeterVG_anno_mye, levels = celltype_order)

base_cols <- c("#454fd2e8", "#698ff7ff", "#6bbbd8ff","#b5cf23c0",
  	"#ded232ff","#cd57b5ff","#cf465bc4", "#ed977aff", "#e86a40ff",
  	"lightgrey") ## change all lym to grey color

names(celltype_order) <- base_cols

	# Then plot with SCP::CellDimPlot using your custom colors
pdf(paste0(F_path, "/annoplot_split.pdf"), width = 12, height = 10)
		CellDimPlot(
  		SeuratObj,
  		group.by = "PeterVG_anno_mye",
  		reduction = "umap",
  		split.by = "condition",
  		palcolor = base_cols,
  		bg_color = "grey90",
			pt.alpha = 0.8,
	pt.size = 0.2
		)
	dev.off()

pdf(paste0(F_path, "/annoplot_all.pdf"), width = 5, height = 4)
	CellDimPlot(
	  SeuratObj,
	  group.by = "PeterVG_anno_mye",
	  reduction = "umap",
	  palcolor = base_cols,
	  bg_color = "grey90",
	pt.alpha = 0.8,
	pt.size = 0.2
	)
dev.off()