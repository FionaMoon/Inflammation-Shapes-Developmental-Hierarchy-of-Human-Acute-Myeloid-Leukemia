## =========================================================
## Correlation of cell type fractions and inflammation score
## =========================================================

rm(list = ls())
gc()

setwd("D:/备份盘/PeiLab/工作/刘丽娜老师/Sequencing_Analysis/BeatAML")

## Load packages
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)


## Load data
data <- fread("./decon_clinical_order.csv", data.table = F)

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
celltypes <- c("LSPC","GMP-like",  "ProMono-like", "DC", "Mono")
inflam_var <- "GOBP_INFLAMMATORY_RESPONSE"


## -------------------------------
## Preprocess data
## -------------------------------
plot_df <- data %>%
  pivot_longer(
    cols = all_of(celltypes),
    names_to = "CellType",
    values_to = "Fraction"
  )

cor_df <- plot_df %>%
  group_by(CellType) %>%
  summarise(
    cor = cor(Fraction, .data[[inflam_var]], method = "spearman"),
    p = cor.test(
      Fraction,
      .data[[inflam_var]],
      method = "spearman",
      exact = FALSE
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    logP = -log10(p),
    signif = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

cor_df

ggplot(
  cor_df,
  aes(
    x = CellType,
    y = inflam_var,
    color = cor,
    size = logP
  )
) +
  geom_point() +
  geom_text(aes(label = signif), vjust = -1.2, size = 5) +
  scale_color_gradient2(
    low = "#2C7BB6",
    mid = "white",
    high = "#D7191C",
    midpoint = 0,
    name = "Spearman r"
  ) +
  scale_size(range = c(3, 8), name = "-log10(P)") +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12)
  )


plot_df$CellType <- factor(plot_df$CellType, levels = celltypes)


pdf("./celltype_correlation_BEATAML.pdf", width = 12, height = 4)
ggplot(
  plot_df,
  aes(x = .data[[inflam_var]], y = Fraction)
) +
  geom_point(aes(color = CellType), size = 2, alpha = 0.7) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linetype = "dashed",
    color = "black"
  ) +
  stat_cor(
    method = "spearman",
    size = 4,
    label.x.npc = "left",
    label.y.npc = "top"
  ) +
  scale_color_manual(values = immune_colors, guide = "none") +
  facet_wrap(~ CellType, scales = "free_x", nrow = 1) +
  theme_classic() +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12)
  ) +
  labs(
    x = "GOBP inflammatory response score",
    y = "Cell fraction"
  )
dev.off()
