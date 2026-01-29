rm(list = ls())
gc()

## load packages
library(data.table)
library(dplyr)

## set working dir
setwd("F:/xx/Sequencing_Analysis")

# Load RBPs
RBPs <- fread("./LLN_Data/490_RBP_CRISPR_target.csv")$`Gene name`
RBPs <- unique(RBPs)

# List of DEG data frames from different cohorts
DEG_TCGA <- read.csv("./TCGALAML_Inflammation_High_vs_Low_DEG_all.csv")
colnames(DEG_TCGA)[1] <- "gene"
DEG_BEATAML <- read.csv("./BeatAML/BEATAML_Inflammation_High_vs_Low_DEG_all.csv")
colnames(DEG_BEATAML)[1] <- "gene"
DEG_scRNA <- read.csv("./2023_Pei_CD/DEGs_MonoLike_LSPC-Quie.csv")
colnames(DEG_scRNA)[2] <- "log2FoldChange"

DEG_list <- list(DEG_TCGA, DEG_BEATAML, DEG_scRNA)  # replace with your actual data frame names
names(DEG_list) <- c("DEG_TCGA", "DEG_BEATAML", "DEG_scRNA")

# Specify RBP gene(s) to highlight
RBP_labels <- "RBM47"
cohort_titles <- c("TCGA", "BEATAML", "CITE-seq_SSP")
highlight_colors <- c("red", "blue", "darkgreen")  # Different color for each cohort


# Output to PDF
pdf(file = "./LLN_Data/DEG-RBPs_x4.pdf", width = 5, height = 3)  
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1), oma = c(0, 0, 0, 0)) ## inner margins (bottom, left, top, right)

for (i in seq_along(DEG_list)) {
  DEG_RBP <- DEG_list[[i]] %>% 
    filter(gene %in% RBPs) %>%
    arrange(log2FoldChange)
  DEG_RBP$rank <- seq_len(nrow(DEG_RBP))
  fwrite(DEG_RBP , paste0("RBPs_rank_", names(DEG_list)[i], ".csv"))
  
  plot(DEG_RBP$rank, DEG_RBP$log2FoldChange,
       pch = 19,
       xlab = "Rank",
       ylab = "Avg log2FC",
       main = cohort_titles[i],
       col = ifelse(DEG_RBP$gene %in% RBP_labels, 
                    adjustcolor(highlight_colors[i], alpha.f = 0.7), 
                    adjustcolor("gray", alpha.f = 0.7)),
       cex = 0.9)
  
  text(DEG_RBP$rank, DEG_RBP$log2FoldChange,
       labels = ifelse(DEG_RBP$gene %in% RBP_labels, DEG_RBP$gene, ""), 
       pos = 2,
       cex = 1.2,
       offset = 0.3, 
       col = highlight_colors[i])
}

dev.off()
