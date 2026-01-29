rm(list = ls())
gc()

setwd("D:/Sequencing_Analysis/BeatAML")

############################
## Load libraries
############################
library(ComplexHeatmap)
library(dplyr)
library(grid)
library(data.table)
library(RColorBrewer)

############################
## Load data
############################
data <- fread("./decon_clinical_order.csv", data.table = FALSE)

############################
## Inflammation quantiles (replace FAB)
############################
inflam <- data$GOBP_INFLAMMATORY_RESPONSE

Inflam_group <- cut(
  inflam,
  breaks = quantile(inflam, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE),
  include.lowest = TRUE,
  labels = c("0–25%", "25–50%", "50–75%", "75–100%")
)

data$Inflam_group <- Inflam_group
Inflam_color <- brewer.pal(9, "PuRd")[c(2,3,5,6)]
# Inflam_color <- c("#FFEDA0", "#FEB24C", "#FD8D3C", "#FC4E2A")
names(Inflam_color) <- levels(Inflam_group)

############################
## Remove unused immune types
############################
data <- data %>% select(!c("NK","MEP_Ery", "CD4 T", "CD8 T", "B"))

############################
## Rename immune columns
############################
colnames(data)[5:9] <- c("LSPC", "DC", "GMP-like", "ProMono-like", "Mono")

############################
## Immune cell colors (unchanged)
############################
immune_colors <- c(
  "DC"           = "#8E7CC3",
  "Mono"         = "#AC71B5",
  "ProMono-like" = "#CCB7E2",
  "GMP-like"     = "#377eb8",
  "LSPC"         = "#EE7576"
)

############################
## Desired order
############################
desired_order <- c("DC", "Mono", "ProMono-like", "GMP-like", "LSPC")

immune_cell_data <- data[, desired_order]
immune_cell_data <- round(immune_cell_data / rowSums(immune_cell_data), 4)

############################
## Stacked barplot + inflammation line
############################
stacked_barplot <- rowAnnotation(
  Immune = anno_barplot(
    immune_cell_data,
    gp = gpar(fill = immune_colors),
    bar_width = 1,
    ylim = c(0, 1)
  ),

  Inflammation = anno_lines(
    inflam,
    add_points = TRUE,
    size = unit(1.2, "mm"),
    pt_gp = gpar(col = Inflam_color[Inflam_group]),
    gp = gpar(col = "black"),
    smooth = TRUE
  ),

  height = unit(16, "cm"),
  width  = unit(4, "cm")
)

############################
## Inflammation group annotation (FAB-style replacement)
############################
Inflam_anno <- rowAnnotation(
  Inflam_group = as.factor(data$Inflam_group),
  col = list(Inflam_group = Inflam_color)
)

############################
## Combine annotations
############################
ha_list <- stacked_barplot + Inflam_anno

############################
## Draw
############################
pdf(file = "./immune_BEATAML_20260128_test.pdf", width = 6, height = 8)
draw(ha_list)
dev.off()
