gc()

library(ComplexHeatmap)
library(circlize)
library(readxl)

setwd("G:/Liangy-temp/skipper")
list.files(pattern = ".xlsx")

data <- read_xlsx("222_OE HEATMAP.xlsx")
head(data)

data <- as.data.frame(data)
rownames(data) <- data$Gene.name
data <- data[,-1]

## Your selected genes to annotate
selected_genes <- c("CCL8", "IL1B", "IL1A", "CCL2", "CXCL10", "CCL8", "TNF", "IFNGR2")
selected_genes %in% rownames(data)  # check if all selected genes are in data
selected_genes <- intersect(selected_genes, rownames(data))  # ensure exist


group <- factor(c(rep("Vector",3), rep("OE",3)), levels = c("Vector", "OE"))

col_anno <- HeatmapAnnotation(
  Group = group,
  col = list(Group = c(Vector = "light grey", OE = "#ff1616"))
)

# Row annotation with anno_mark: mark selected genes
row_anno <- rowAnnotation(
  Selected = anno_mark(
    at = which(rownames(data) %in% selected_genes),
    labels = selected_genes
  )
)

data <- as.data.frame(t(scale(t(data))))
colnames(data) <- c(paste("Vector", 1:3, sep ="_"), paste("OE", 1:3, sep ="_"))

col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# Plot heatmap
pdf("./heatmap_short.pdf", width = 5, height = 4)
Heatmap(
  data,
  name = "Expression",
  col = col_fun,
  top_annotation = col_anno,
  right_annotation = row_anno,
  show_row_names = FALSE,  # hide all row names to reduce clutter
  cluster_rows = T,
  cluster_columns = F
)
dev.off()
