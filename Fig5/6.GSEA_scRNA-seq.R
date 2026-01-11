## GSEA for scRNA-seq
rm(list = ls())
gc()
packages <- c("GSVA", "Seurat", "dplyr", "tidyr", "ggplot2", "GSEABase", "ggridges", "ComplexHeatmap", "ggthemes")
invisible(lapply(packages, function(pkg) suppressPackageStartupMessages(library(pkg, character.only = TRUE))))
setwd("/Users/MyFolder/Code/LLN-RBM47/CITE-seq_cytokine/AML-8")



seurat_obj <- readRDS("/Users/MyFolder/Code/LLN-RBM47/CITE-seq_cytokine/AML-8/AML-8_PeterVG_anno_cytotrace_ADT_V5.rds")


gmt_folder <- "/Users/MyFolder/Code/LLN-RBM47/CITE-seq_cytokine/AML-8_script/V2/gmt"
# List all files in the folder (assuming only one GMT file, or choose the one you want)
gmt_files <- list.files(gmt_folder, pattern = "\\.gmt$", full.names = TRUE)

# Loop over files to read and combine
gene_sets_list <- list()

for (f in gmt_files) {
  gs <- getGmt(f)
  # Convert GeneSetCollection to named list
  gs_list <- lapply(gs, geneIds)
  names(gs_list) <- names(gs)
  
  # Combine into one list
  gene_sets_list <- c(gene_sets_list, gs_list)
}
str(gene_sets_list)

expr_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")  # log-normalized

ssgsea_scores <- gsvaParam(
  expr = expr_mat,
  geneSets = gene_sets_list
)

gsvaranks <- gsvaRanks(ssgsea_scores)
gsvaranks
es <- gsvaScores(gsvaranks)
es[1:3,1:3]

if(all(colnames(es) == rownames(seurat_obj@meta.data))){
  seurat_obj@meta.data <- cbind(seurat_obj@meta.data, t(es))
}


df <- seurat_obj@meta.data %>% dplyr::select(condition, rownames(es)) %>%
  pivot_longer(
    cols = -condition,
    names_to = "pathway",
    values_to = "score"
  )

# Ridge plot
ggplot(df, aes(x = score, y = pathway, fill = condition)) +
  geom_density_ridges(
    aes(height = ..scaled..),   # scales each ridge to 0–1
    stat = "density",
    alpha = 0.7,
    scale = 1.0
  ) +
  scale_fill_manual(values = color.3) +
  theme_few() +
  labs(x = "ssGSEA score", y = "Pathway") +
  theme(axis.text.y = element_text(size = 14))



#################
# Optional: order pathways by median across all conditions
pathway_order <- df %>%
  group_by(pathway) %>%
  summarize(median_score = median(score, na.rm = TRUE)) %>%
  arrange(desc(median_score)) %>%
  pull(pathway)

df$pathway <- factor(df$pathway, levels = pathway_order)

# Compute medians per pathway and condition
# Assign colors manually to conditions
color.3 <- c("royalblue", "orange1", "springgreen4", "pink")
names(color.3) <- levels(df$condition)  # make sure the names match condition levels

# Compute medians
medians <- df %>%
  group_by(pathway, condition) %>%
  summarize(median_score = median(score, na.rm = TRUE), .groups = "drop")

# Boxplot with median lines and custom colors
p <- ggplot(df, aes(x = condition, y = score, fill = condition)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7, width = 0.7, color = "black") +
  geom_line(data = medians, aes(x = condition, y = median_score, group = pathway),
            color = "black", size = 0.7) +
  geom_point(data = medians, aes(x = condition, y = median_score),
             color = "black", shape = 1, size = 2, stroke = 1) +
  facet_wrap(~pathway, scales = "free_y") +
  scale_fill_manual(values = color.3) +
  theme_few() +
  labs(x = "Condition", y = "ssGSEA score") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10)
  )

ggsave(filename = "../AML-8_script/V3/enrichment.pdf", p, width = 6, height = 3.75)
