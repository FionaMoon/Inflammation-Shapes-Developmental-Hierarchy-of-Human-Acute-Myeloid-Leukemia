## Time Series Analysis

rm(list = ls())
gc()

file_path = "D:/xx/Sequencing_Analysis/2023_Pei_CD"
setwd(file_path)

library(Seurat)
library(dplyr)
library(Mfuzz)
library(Biobase)
library(grDevices)

# Extract the pseudo-bulk expression matrix
data <- readRDS("./ssp_scArches_with_MACS_pseudo.Rds")
data
# Get pseudobulk matrix from the 'Integr' assay
pseudo_mat <- GetAssayData(data, assay = "Integr", slot = "counts")

# Check dimensions and format
dim(pseudo_mat)
head(colnames(pseudo_mat))

# Get metadata for sample-celltype mapping
meta <- data@meta.data

head(meta)

# Ensure the rownames match column names in expression matrix
stopifnot(identical(colnames(pseudo_mat), rownames(meta)))

# Make sure pseudo_cluster is set correctly
meta$pseudo_cluster <- sapply(strsplit(split = "-", meta$pseudo_cluster), "[", 1)
meta$pseudo_cluster <- factor(meta$pseudo_cluster, levels = c("LSPC", "GMP", "ProMono", "Mono"))

# Aggregate (average) by pseudo_cluster (i.e., per hematopoietic stage)
clusters <- levels(meta$pseudo_cluster)
avg_expr <- sapply(clusters, function(cluster) {
  samples <- rownames(meta)[meta$pseudo_cluster == cluster]
  rowMeans(pseudo_mat[, samples, drop = FALSE])
})

avg_expr <- avg_expr[rowSums(avg_expr) > 1, ]

# Convert to ExpressionSet
pheno_data <- data.frame(row.names = colnames(avg_expr))
eset <- ExpressionSet(assayData = as.matrix(avg_expr), phenoData = AnnotatedDataFrame(pheno_data))

# Standardize
eset <- standardise(eset)
summary(exprs(eset))

# Estimate fuzzification parameter
m <- mestimate(eset)

# Cluster
cl <- mfuzz(eset, c = 10, m = m)  # Try c = 4~15


# Plot
hematopo_order <- c("LSPC", "GMP", "ProMono", "Mono")

# Create a gradient from grey to blue with 50 colors
custom_colors <- colorRampPalette(c("#f8f8f8","#c8c8c8fe","#a0a0a0fe", "#e4d479","#cd9a35","#de6017"))(100)

pdf("./mfuzz.plot.pdf", width = 18, height = 8)
  par(cex.axis = 1.6, cex.lab = 1.6, cex = 1.6, las = 1, cex.main=2)  # Adjust sizes her
  mfuzz.plot(eset, cl = cl, mfrow = c(2,5),
    time.labels = hematopo_order,
    new.window = FALSE,
    colo = custom_colors)
dev.off()

# Extract increasing expression genes (e.g., cluster 2)
increasing_genes <- names(cl$cluster[cl$cluster == 9]) ## cl$cluster: cluster named by gene name 
increasing_genes <- increasing_genes[-grep(pattern = "^LINC",increasing_genes)]
increasing_genes <- increasing_genes[-grep(pattern = "[.]",increasing_genes)]
"RBM47" %in% increasing_genes
saveRDS(increasing_genes, "./cluste10_DH_up_gene.Rds")


