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

setwd("/Users/liangyue/MyFolder/Code/LLN-RBM47/2023_Pei_CD")

## inflammation geneset
gmtfile = "/Users/liangyue/MyFolder/Code/LLN-RBM47/2023_Pei_CD/GOBP_INFLAMMATORY_RESPONSE.v2024.1.Hs.gmt"
GOBP_inflammation <- read.gmt(gmtfile)
GS.hallmark <- getGeneSets(library = "H")

## scRNA-seq data
sc <- readRDS("/Users/liangyue/MyFolder/Code/LLN-RBM47/2023_Pei_CD/ssp_scArches_with_MACS.Rds")
sc

## 更改2023_CD_scRNA-seq的配色
sc@meta.data$scArches_Cluster <- factor(sc@meta.data$scArches_Cluster,
                   levels = c('LSPC-Quiescent','LSPC-Primed','LSPC-Cycle', 
                              'GMP-like','ProMono-like','Mono-like','cDC-like',
                              "MEP","pre/pro-B","B","Plasma","CD4 T","CD8 T","NK","unknown"))
Idents(sc) <- "scArches_Cluster"

palette_clusters <- c('#e41a1c', '#3e8c3b', '#ff7f00', '#377eb8',
                    "#f781bf","#984ea3","#a65628","#aac9e7","#f6f7a1",
                    "#dbc902","#ffe4ca","#8ad587","#b2df8a","#87a86a","gray")

color_idents <- data.frame(ident = levels(Idents(sc)), color = palette_clusters)
fwrite(color_idents, "./color_idents.csv")

p1 <- DimPlot(sc, reduction = "umap", pt.size = 0.4, alpha = 0.7,  raster =F,
              cols = palette_clusters, label = T,label.size = 3,
              label.box = F, repel =T, group.by = "scArches_Cluster") +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        legend.position="right",
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()) 
ggsave(filename = "./scArches.pdf", p1, width = 5, height = 4)

## AddModuleScore
sc_feature <- intersect(row.names(sc), GOBP_inflammation$gene)
cd_features <- list(GOBP_INFLAMMATORY_RESPONSE = sc_feature)
GS.hallmark <- list(GOBP_INFLAMMATORY_RESPONSE = sc_feature,
                    HALLMARK_INFLAMMATORY_RESPONSE = GS.hallmark$`HALLMARK-INFLAMMATORY-RESPONSE`)

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
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        legend.position="right",
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())
p3 <- FeaturePlot(sc, reduction = "umap", features = "RBM47", 
                  order = TRUE, pt.size = 0.4, max.cutoff = NA) +
  scale_colour_gradientn(colors = c("blue", "skyblue","gray", "pink", "magenta", "red")) + 
  ggtitle('RBM47') +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        legend.position="right",
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())

ggsave(filename = "./sc_InflammationScore.pdf", p2 + p3, width = 5, height = 10)

## only conserve myeloid group
groups <- c("LSPC-Quiescent", "LSPC-Primed", "LSPC-Cycle", "GMP-like",
            "ProMono-like", "cDC-like", "Mono-like")
sc_mye <- subset(sc, idents = groups)

## ssGSEA
## https://www.borch.dev/uploads/screpertoire/articles/running_escape
sc_mye <- runEscape(sc_mye,
                           method = "ssGSEA",
                           gene.sets = GS.hallmark,
                           groups = 5000,
                           min.size = 0,
                           new.assay.name = "escape.ssGSEA")
saveRDS(sc_mye, file = "./ssp_scArches_with_MACS_Myeloid.Rds")

# names(GS.hallmark)[grep(pattern = "INFLAMMATORY", names(GS.hallmark))]
sc_mye$scArches_Cluster <- factor(sc_mye$scArches_Cluster, levels = rev(groups))
p1 <- ridgeEnrichment(sc_mye,
                assay = "escape.ssGSEA",
                group.by = "scArches_Cluster",
                gene.set = "GOBP-INFLAMMATORY-RESPONSE",
                add.rug = TRUE,
                scale = TRUE) + 
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size=15, face = "bold"),
        axis.title.y = element_blank(),
        legend.position="right"
        )
  
p2 <- ridgeEnrichment(sc_mye,
                assay = "escape.ssGSEA",
                group.by = "scArches_Cluster",
                gene.set = "HALLMARK-INFLAMMATORY-RESPONSE",
                add.rug = TRUE,
                scale = TRUE) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size=15, face = "bold"),
        axis.title.y = element_blank(),
        legend.position="right")

ggsave(filename = "./ridgePlot.pdf", p1 + p2, height = 6, width = 10)

## m-LSC vs p-LSC
# Set up the plotting layout: 1 row and 2 columns
pdf(file = "./pLSC_mLSC.pdf", width = 8, height = 4)
  par(mfrow = c(1, 2))
  # Plot the first histogram
  hist(sc_mye$`p-LSC score`, breaks = 100, main = "p-LSC Score", xlab = "p-LSC Score", col = "lightblue")
  # Plot the second histogram
  hist(sc_mye$`m-LSC score`, breaks = 100, main = "m-LSC Score", xlab = "m-LSC Score", col = "lightgreen")
dev.off()
# pLSC 0.1, mLSC 0.1
sc_mye$LSC <- "NO"
sc_mye$LSC[sc_mye$`p-LSC score` > 0.2 & sc_mye$`m-LSC score` < 0] <- "pLSC"
sc_mye$LSC[sc_mye$`m-LSC score` > 0.1 & sc_mye$`p-LSC score` < 0] <- "mLSC"
table(sc_mye$LSC, sc_mye$scArches_Cluster)
Idents(sc_mye) <- "LSC"
DimPlot(sc_mye)

sc_LSC <- subset(sc_mye, idents = c("mLSC", "pLSC"))
DimPlot(sc_LSC) + FeaturePlot(sc_LSC, features = "GOBP_INFLAMMATORY_RESPONSE1")
## ssGSEA
## https://www.borch.dev/uploads/screpertoire/articles/running_escape
sc_LSC <- runEscape(sc_LSC,
                    method = "ssGSEA",
                    gene.sets = GS.hallmark,
                    groups = 1000,
                    min.size = 0,
                    new.assay.name = "escape.ssGSEA")
saveRDS(sc_LSC, file = "./ssp_scArches_with_MACS_LSC,Rds")
p1 <- ridgeEnrichment(sc_LSC,
                      assay = "escape.ssGSEA",
                      group.by = "LSC",
                      gene.set = "GOBP-INFLAMMATORY-RESPONSE",
                      add.rug = TRUE,
                      scale = TRUE) + 
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size=15, face = "bold"),
        axis.title.y = element_blank(),
        legend.position="right"
  )

p2 <- ridgeEnrichment(sc_LSC,
                      assay = "escape.ssGSEA",
                      group.by = "LSC",
                      gene.set = "HALLMARK-INFLAMMATORY-RESPONSE",
                      add.rug = TRUE,
                      scale = TRUE) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size=15, face = "bold"),
        axis.title.y = element_blank(),
        legend.position="right")
library(scCustomize)
VlnPlot_scCustom(seurat_object = sc_LSC, features = "GOBP_INFLAMMATORY_RESPONSE1")
ggsave(filename = "./ridgePlot_LSC.pdf", p1 + p2, height = 6, width = 10)

DEGs_LSC <- FindMarkers(sc_LSC, ident.1 = "mLSC", ident.2 = "pLSC")
DEGs_LSC <- DEGs_LSC %>% arrange(desc(avg_log2FC))
DEGs_LSC["RBM47",]

## DEG LSPC-Quiescent vs Mono-like
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


## load RBPs
RBPs <- fread("../LLN_Data/490_RBP_CRISPR_target.csv")
# RBPs <- readxl::read_xlsx("../LLN_Data/A_cencus_of_human_RBPs.xlsx")
RBPs <- unique(RBPs$`Gene name`)
DEG_RBP <- DEGs[DEGs$gene %in% RBPs,] %>% 
  arrange(avg_log2FC)
DEG_RBP$rank <- 1:nrow(DEG_RBP)

# Specify the RBP gene(s) you want to highlight
 RBP_labels <- "RBM47"

pdf(file = "../LLN_Data/DEG-RBPs.pdf", width = 3.5, height = 4) 
# Plotting the data with transparency
 plot(DEG_RBP$rank, DEG_RBP$avg_log2FC,
      pch = 19, # Solid circle for points
      xlab = "Rank", # Custom x-axis label
      ylab = "Average log2 Fold Change", # Custom y-axis label
      col = ifelse(DEG_RBP$gene %in% RBP_labels, 
                   rgb(1, 0, 0, alpha = 0.5), # Red with alpha 0.5 for RBP
                   rgb(0.5, 0.5, 0.5, alpha = 0.5))) # Grey with alpha 0.5 for other points
# Add text labels with transparency
 text(DEG_RBP$rank, DEG_RBP$avg_log2FC,
      labels = ifelse(DEG_RBP$gene %in% RBP_labels, DEG_RBP$gene, ""), # Label only RBP genes
      pos = 2, # Position the labels at the bottom
      cex = 0.7, # Size of the text
      offset = 0.5, 
      col = "red") # Red with alpha 0.7 for text labels
dev.off()
