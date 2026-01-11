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
D_path = "/Users/MyFolder/Code/LLN-RBM47/CITE-seq_cytokine"
F_path = "/Users/MyFolder/Code/LLN-RBM47/CITE-seq_cytokine/ZYG_script/V3"

obj <- readRDS("/Users/MyFolder/Code/LLN-RBM47/CITE-seq_cytokine/ZYG/ZYG_PeterVG_anno_cytotrace_ADT_V5.rds")


DefaultAssay(object = obj) <- "RNA"

TFs_RNA_3 <- list(Prim_based = c("RUNX1", "CEBPA", "MYB","ETV6"),
  Mature_based = c("TNFRSF1B", "IL1R1", "IL1R2", "CEBPB", "SPI1", "MAFB", "CSF1R") ## "ETV6" removed
  )

TFs <- unlist(TFs_RNA_3, use.names = FALSE)

## dotplot

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
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 14, face="bold"),
                axis.text.y = element_text(size=14, face="bold"),
                legend.text = element_text(size=14, face="bold"),
                axis.text = element_text(size=10, face="bold"),
                axis.title = element_text(size=12, face="bold"),
                axis.title.x = element_blank(),
                axis.title.y = element_blank())
    SN=paste0(samplename, "_TFsdot_marker.pdf")
    p_name = paste(fig.dir, SN, sep = "/")
    p <- p_umap + p1 + plot_layout(widths = c(1, 2))
    ggsave(p_name, p ,width = 8, height = 5)
}


ID = "ZYG"
manual_add_marker(obj, geneSet = "TFs_RNA_3", fig.dir = F_path, data.dir = D_path, 
  assays = "RNA", select_ident = 'condition', samplename = ID,
  reduction = "umap")
