##################### only save DEG which is RBP

rm(list = ls())
gc()


library(ggvenn)
library(data.table)
library(dplyr)

## set working dir
setwd("D:/备份盘/PeiLab/工作/刘丽娜老师/Sequencing_Analysis")

# Load RBPs
RBPs <- fread("./LLN_Data/490_RBP_CRISPR_target.csv")$`Gene name`
RBPs <- unique(RBPs)


# List of DEG data frames from different cohorts
DEG_TCGA <- read.csv("./炎症相关/TCGALAML_Inflammation_High_vs_Low_DEG_all.csv")
colnames(DEG_TCGA)[1] <- "gene"
DEG_BEATAML <- read.csv("./BeatAML/BEATAML_Inflammation_High_vs_Low_DEG_all.csv")
colnames(DEG_BEATAML)[1] <- "gene"
DEG_scRNA <- read.csv("./2023_Pei_CD/DEGs_MonoLike_LSPC-Quie.csv")
colnames(DEG_scRNA)[2] <- "log2FoldChange"


# Prepare a named list of upregulated genes
venn_input <- list(
  TCGA = DEG_TCGA %>% filter(abs(log2FoldChange) > 1 & pvalue < 0.05) %>% filter(gene %in% RBPs) %>% pull(gene),
  BEATAML = DEG_BEATAML %>% filter(abs(log2FoldChange) > 1 & pvalue < 0.05) %>% filter(gene %in% RBPs) %>% pull(gene),
  `CITE-seq_SSP` = DEG_scRNA %>% filter(abs(log2FoldChange) > 1 & p_val_adj < 0.05) %>% filter(gene %in% RBPs) %>% pull(gene)
)

co_DEG <- Reduce(intersect, venn_input)
"RBM47" %in% co_DEG
TCGA_RBM47 <- DEG_TCGA %>% filter(abs(log2FoldChange) > 1 & pvalue < 0.05) %>% filter(gene %in% co_DEG) %>% arrange(gene)
BEATAML = DEG_BEATAML %>% filter(abs(log2FoldChange) > 1 & pvalue < 0.05) %>% filter(gene %in% co_DEG) %>% arrange(gene)
`CITE-seq_SSP` = DEG_scRNA %>% filter(abs(log2FoldChange) > 1 & p_val_adj < 0.05) %>% filter(gene %in% co_DEG) %>% arrange(gene)

df <- data.frame(
  Gene = co_DEG[order(co_DEG)],
  TCGA_Log2FC = TCGA_RBM47$log2FoldChange,
  BEATAML_Log2FC = BEATAML$log2FoldChange,
  CITEseq_Log2FC = `CITE-seq_SSP`$log2FoldChange
)

fwrite(df, "./LLN_Data/Venn_co_DEG.csv")


pdf("./LLN_Data/VennPlot_DEG_RBPs.pdf", width = 4, height = 4)
ggvenn(venn_input, 
       fill_color = c("red", "blue", "darkgreen"),
       stroke_size = 0.5, 
       text_size = 4)
dev.off()

