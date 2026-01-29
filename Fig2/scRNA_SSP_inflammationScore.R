# inflammation Score has been calculated by previous like the following
rm(list = ls())
gc()

library(Seurat)
library(ggplot2)
library(ggprism)
library(dplyr)
library(data.table)
library(scCustomize)


setwd("D:/Sequencing_Analysis/2023_Pei_CD")
list.files(pattern = ".Rds")

## scRNA-seq data
sc <- readRDS("./ssp_scArches_with_MACS_Myeloid.Rds")

# scRNA-seq
# DEG LSPC-Quiescent vs Mono-like
DEGs <- FindMarkers(sc_mye, ident.1 = "Mono-like",
                    ident.2 = "LSPC-Quiescent",logfc.threshold = 0,
                    min.pct = 0)
DEGs <- DEGs %>% mutate(gene = rownames(DEGs)) %>%
  arrange(desc(avg_log2FC))
DEGs <- DEGs[-grep(pattern = "^LINC", DEGs$gene),]
DEGs <- DEGs[-grep(pattern = "[.]", DEGs$gene),]
logFC_cutoff <- with(DEGs, mean(abs(avg_log2FC)) + 2*sd(abs(avg_log2FC)))
logFC_cutoff <- round(logFC_cutoff, 2)
DEGs$change <- factor(ifelse(DEGs$p_val_adj < 0.05 & abs(DEGs$avg_log2FC) > logFC_cutoff,
                             ifelse(DEGs$avg_log2FC > logFC_cutoff, 'UP', 'DOWN'), 'STABLE'))
head(DEGs)
fwrite(DEGs, "./DEGs_MonoLike_LSPC-Quie.csv")
DEGs["RBM47",]
