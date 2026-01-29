## PCA + Hierarchical clustering
## AML cell lines (gene-restricted)

library(data.table)
library(dplyr)
library(PCAtools)
library(ggthemes)

setwd("D:/备份盘/PeiLab/工作/刘丽娜老师/AML-CellLine")

## -------------------------------
## Load expression data
## -------------------------------
data_logTPM <- fread(
  "Expression_Public_25Q3_subsetted.csv",
  data.table = FALSE
)

data_mye <- data_logTPM %>%
  filter(lineage_1 == "Myeloid")

## -------------------------------
## Cell lines (K562 excluded)
## -------------------------------
cell_select <- c(
  "KASUMI1", "HL60", "KG1",
  "MV411", "MOLM13","K562",
  "U937", "THP1", "OCIAML3"
)

data_select <- data_mye %>%
  filter(cell_line_display_name %in% cell_select) %>%
  select(2, 8:ncol(.))

## -------------------------------
## Load gene set (Cluster 1 genes)
## -------------------------------
cluster1 <- fread(
  "D:/备份盘/PeiLab/工作/刘丽娜老师/Sequencing_Analysis/2023_Pei_CD/cluste_all_DH_up_gene-ordered_mannually.csv",
  data.table = FALSE
)

## -------------------------------
## Expression matrix (genes × cell lines)
## -------------------------------
data_select <- data.frame(t(data_select), check.names = FALSE)
colnames(data_select) <- data_select[1, ]
data_select <- data_select[-1, ]

## enforce numeric
data_select[] <- lapply(data_select, as.numeric)

## remove zero-variance genes
data_select <- data_select[apply(data_select, 1, sd) > 0, ]

## restrict to Cluster 1 gene set
co_gene <- intersect(rownames(data_select), cluster1$cluster_1)
data_select <- data_select[co_gene, ]

## -------------------------------
## PCA using PCAtools
## -------------------------------
p <- pca(
  data_select,
  scale = TRUE,
  center = TRUE,
  removeVar = 0.1
)

## -------------------------------
## Hierarchical clustering (Ward.D2)
## -------------------------------
d <- dist(p$rotated[, 1:9], method = "euclidean")
hc <- hclust(d, method = "ward.D2")

hc_cluster <- factor(cutree(hc, k = 3))

p$metadata <- data.frame(
  CellLine = rownames(p$rotated),
  Cluster = hc_cluster
)

## -------------------------------
## PCA plot (PC1 / PC2 variance shown)
## -------------------------------

pdf("./PCA.pdf", width = 4.5, height = 4)
biplot(
  p,
  colby = "Cluster",
  lab = p$metadata$CellLine,
  pointSize = 4,
  labSize = 4,
  legendPosition = "right",
  # encircle config
  encircle = TRUE,
  encircleFill = T,
  encircleLineSize = 1
) +
  guides(fill = "none") +
  theme_few(base_size = 14) +
  ggtitle("PCA of AML Cell Lines (Hierarchical clustering)")
dev.off()
