rm(list = ls())
gc()

setwd("D:/备份盘/PeiLab/工作/刘丽娜老师/Sequencing_Analysis/BeatAML")


# Load necessary libraries
library(ComplexHeatmap)
library(dplyr)
library(grid)
library(data.table)


data <- fread("./decon_clinical_order.csv", data.table = F)

FAB_color = c(
  "M0" = "#1f77b4",  # Blue
  "M1" = "#ff7f0e",  # Orange
  "M2" = "#2ca02c",  # Green
  "M3" = "#d62728",  # Red
  "M4" = "#9467bd",  # Purple
  "M5" = "#8c564b",  # Brown
  "M6" = "#e377c2",  # Pink
  "M7" = "#7f7f7f"   # Gray
)

# Create bar annotations for Inflammation, with FAB color applied
names(data$GOBP_INFLAMMATORY_RESPONSE) <- data$fabBlastMorphology


# Create the named color annotation for FAB
FAB <- rowAnnotation(
  FAB = as.factor(data$fabBlastMorphology),
  col = list(FAB = FAB_color)
)

colnames(data)

data <- data %>% select(!c("NK","MEP_Ery", "CD4 T", "CD8 T", "B"))

colnames(data)[5:9] <- c("LSPC", "DC", "GMP-like", "ProMono-like", "Mono")

# Create stacked barplot annotations using the immune cell data
immune_colors <- c(
  "DC" = "#8E7CC3",          # Gray for DCs
  "Mono" = "#AC71B5",       # Pink for Monocytes
  "ProMono-like" = "#CCB7E2", # Purple for ProMono-like cells
  "GMP-like" = "#377eb8",   # Orange for GMP-like
  "LSPC" = "#EE7576"       # Brown for LSPC
)

# Define the desired order of cell types
desired_order <- c("DC", "Mono", "ProMono-like","GMP-like", "LSPC")


# Reorder the columns of immune_cell_data
immune_cell_data <- data[, desired_order]
immune_cell_data <- round(immune_cell_data/rowSums(immune_cell_data), 4)

# Create stacked barplot annotation with custom colors
stacked_barplot <- rowAnnotation(
 foo = anno_barplot(
    immune_cell_data,
    gp = gpar(fill = immune_colors),  # Apply the custom colors
    # axis_param = list(direction = "reverse"),
    bar_width = 1,
    ylim = c(0, 1)  # Set the y-axis limits for proportion data
  ),
   Inflammation = anno_lines(
    data$GOBP_INFLAMMATORY_RESPONSE,  # Inflammation data for lines
    add_points = TRUE, 
    size = unit(1.2, "mm"),
    pt_gp = gpar(col = FAB_color[as.factor(data$fabBlastMorphology)]),  # Use gpar to set point colors
    gp = gpar(col = "black"),  # Line color (you can customize this)
    smooth = TRUE  # Smooth the lines
  ),
  height = unit(12, "cm"), ## only change height of complex annotation, simple anno won't work
  width = unit(4, "cm")
)

# Combine the annotations using c()
# ha_list <-   stacked_barplot + FAB + ha
ha_list <-   stacked_barplot + FAB

# Draw the combined annotations
pdf(file = "./immune_BEATAML_20251107.pdf", width = 6, height = 8)
  draw(ha_list)
dev.off()
