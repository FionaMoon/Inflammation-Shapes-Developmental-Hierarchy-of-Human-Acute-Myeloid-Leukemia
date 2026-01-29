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

setwd("/Users/Code/LLN-RBM47/2023_Pei_CD")

## inflammation geneset
gmtfile = "/Users/Code/LLN-RBM47/2023_Pei_CD/GOBP_INFLAMMATORY_RESPONSE.v2024.1.Hs.gmt"
GOBP_inflammation <- read.gmt(gmtfile)
GS.hallmark <- getGeneSets(library = "H")

## scRNA-seq data
sc <- readRDS("/Users/l Code/LLN-RBM47/2023_Pei_CD/ssp_scArches_with_MACS.Rds")
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
  scale_colour_gradientn(colors = c("blue", "skyblue", "gray", "pink", "magenta", "red")) + 
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
