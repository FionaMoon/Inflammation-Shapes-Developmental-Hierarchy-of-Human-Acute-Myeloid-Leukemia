## cellchat for prim/mature subtype interaction

rm(list = ls())
gc()

## ---- libraries ----
library(Seurat)
library(CellChat)
library(dplyr)

setwd("/Users/MyFolder/Code/LLN-RBM47/CITE-seq_cytokine/AML-8")

cellchat.list <- readRDS("./AML-8_chatList.rds")
cellchat.merged <- readRDS("./AML-8_cellchat.merged.rds")

#####==== Visuallization ====#####
cellchat.merged@meta$datasets <- factor(
  cellchat.merged@meta$datasets,
  levels = c("CT", "TNFA", "IL1B", "Combo")
)

levels(cellchat.merged@idents$joint)

pdf("../AML-8_script/V3/circle_group_1.pdf")
weight.max <- getMaxWeight(cellchat.list, attribute = c("idents","count"))
par(mfrow = c(2,2), xpd=TRUE)
for (i in 1:length(cellchat.list)) {
  netVisual_circle(cellchat.list[[i]]@net$count, weight.scale = T, 
   edge.weight.max = weight.max[2],  #label.edge= T,
   sources.use = c("ProMono-like", "Mono-like", "cDC-like", "pDC-like"),
  targets.use = c("LSPC-Primed", "LSPC-Quiescent", "LSPC-Cycle"),
    arrow.width = 1, arrow.size = 0.4,
    title.name = names(cellchat.list)[i])
}
dev.off()


pdf("../AML-8_script/V3/chord_group.pdf", width = 5,height = 5)
par(mfrow = c(1,1), xpd=TRUE)
object.list <- cellchat.list
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = unique(object.list[[i]]@netP$pathways), 
  thresh = 0.01,
    layout = "chord",  
    sources.use = c("ProMono-like", "Mono-like", "cDC-like", "pDC-like"),
  targets.use = c("LSPC-Primed", "LSPC-Quiescent", "LSPC-Cycle"),
    signaling.name = names(object.list)[i], 
    remove.isolate = TRUE)
}
dev.off()

class(cellchat.merged@netP$CT$prob)
