## Dotplot of marker

## load packages
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
library(ggthemes)
options(Seurat.object.assay.version = "v5")
# needs to be set for large dataset analysis
options(future.globals.maxSize = 10 * 1e9)

## set path
D_path = "D:/supple_Fig"
F_path = "D:/supple_Fig"


## load data
obj <- readRDS("D:/Sequencing_Analysis/2023_Pei_CD/ssp_scArches_with_MACS.Rds")
Idents(obj) <- "scArches_Cluster"
obj <- subset(obj, idents = "unknown", invert =T)
celltype_order <- c('LSPC-Quiescent',  'LSPC-Primed', 'LSPC-Cycle',
'GMP-like', 'ProMono-like', 'Mono-like', 'cDC-like', 'MEP',
'pre/pro-B', 'B', 'Plasma',
'CD4 T', 'CD8 T', 'NK')
obj$scArches_Cluster <- factor(as.character(obj$scArches_Cluster), 
levels = rev(celltype_order))              

DefaultAssay(object = obj) <- "RNA"


## function of dotplot
T_cell = c('PTPRC','CD3D', 'CD3E','CD3G','CD8A','FOXP3') # T cell 'PTPRC' #免疫细胞
B_plasma = c('TCF3', 'EBF1', "PAX5",'CD19', 'CD79A','MS4A1','MZB1') # B cell
Erythrocyte = c('TFRC', 'HBA1', 'ALAS2','GATA1','KLF1') # erythroid cells
Monocyte = c('CD33', 'LYZ','S100A9', 'S100A8', 'CD14','FCGR3A')# monocyte
dendritic_cells = c('IL3RA',  'CD1C','HLA-DPB1','HLA-DQA1' ,'ITGAX','FCER1A', 'CST3','IRF7') ## DC 
Natural_killer_cell = c('GZMB',"NKG7",'NCAM1','GNLY','NCR1') # NK 
Palate = c('PPBP','PF4','GATA2','ITGA2B','ITGB3') # Mk  # Palate
stem_cell = c('CD34','CD38', 'SPINK2','CRHBP', 'NR4A1', 'RGS18', 'CNRIP1') # HSPC
GMP = c('AZU1', 'MPO', 'ELANE')
## MkP / ErP / MEP
## SPINK2: CMP +
## CRHBP: CMP / MEP +
## NR4A1: MEP
## RGS18: MEP / MkP
## CNRIP1: MEP / ErP

genes_to_check = list(
  HSPC = stem_cell,
        GMP = GMP,
        Mono = Monocyte,
        DC = dendritic_cells,
        Palate_Mk = Palate,
        Ery = Erythrocyte,
        B_plasma = B_plasma,
        T = T_cell,
        NK = Natural_killer_cell
        )

manual_add_marker <- function(x, geneSet = "GS_provided",fig.dir = "./Seurat/fig",data.dir = "./Seurat/data", assays = "SCT", select_ident = 'SCT_snn_res.1.2', samplename, reduction = "umap"){
    if(geneSet == 'GS_provided'){
        genes_to_check = genes_to_check
    } else {
        genes_to_check = get(geneSet)
    }
    # 判断这些基因是否在表达矩阵中
    genes <- unlist(genes_to_check, use.names = FALSE)
    data = x
    Idents(data) = select_ident
    allGenes = row.names(data)
    cat("Genes which are not contain in our sample:")
  cat(genes[!genes %in% allGenes],"\n")
    # 画图探究
    p_umap <- DimPlot(data, reduction = reduction,label = TRUE,repel = TRUE, group.by = select_ident) + labs(title = select_ident)
    p1 <- DotPlot(data, assay = assays, features = genes_to_check, group.by = select_ident, dot.scale = 5, 
            cols = "RdYlBu") + 
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10, face="bold"),
                axis.text.y = element_text(size=12, face="bold"),
                legend.text = element_text(size=12, face="bold"),
                axis.text = element_text(size=9, face="bold"),
                axis.title = element_text(size=10, face="bold"),
                axis.title.x = element_blank(),
                axis.title.y = element_blank())
    SN=paste0(samplename, "_TFsdot_marker.pdf")
    p_name = paste(fig.dir, SN, sep = "/")
    p <- p_umap + p1 + plot_layout(widths = c(1, 2))
    ggsave(p_name, p ,width =22, height = 5)
}


ID = "GSE232559"
manual_add_marker(obj, fig.dir = F_path, data.dir = D_path, 
  assays = "RNA", select_ident = 'scArches_Cluster', samplename = ID,
  reduction = "umap")
table(obj$scArches_Cluster)

## DoHeatmap
# DoHeatmap(subset(obj, downsample = 100), features = features, size = 3)


## DEG of each cluster 





