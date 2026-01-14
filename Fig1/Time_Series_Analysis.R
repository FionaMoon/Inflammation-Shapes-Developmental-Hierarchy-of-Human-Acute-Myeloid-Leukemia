## Time Series Analysis

rm(list = ls())
gc()

file_path = "D:/备份盘/PeiLab/工作/刘丽娜老师/Sequencing_Analysis/2023_Pei_CD"
setwd(file_path)

library(Seurat)
library(dplyr)
library(Mfuzz)
library(Biobase)
library(grDevices)
library(data.table)

# Extract the pseudo-bulk expression matrix
data <- readRDS("./ssp_scArches_with_MACS_pseudo.Rds")
data
# Get pseudobulk matrix from the 'Integr' assay
pseudo_mat <- GetAssayData(data, assay = "Integr", layer = "counts")

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
cl <- mfuzz(eset, c = 10, m = m)  # Try c = 4~6



# pdf("./mfuzz_plot1.pdf", width = 18, height = 8)
#   par(cex.axis = 1.6, cex.lab = 1.6, cex = 1.6, las = 1, cex.main=2)  # Adjust sizes her
#   mfuzz.plot(eset, cl = cl, mfrow = c(2,5),
#     time.labels = hematopo_order,
#     new.window = FALSE,
#     colo = custom_colors)
# dev.off()


my_mfuzz_plot <- function(eset, cl, mfrow=c(2,2),
                          time.labels=NULL,
                          new.window=TRUE,
                          colo,
                          min.mem=0,
                          ylim=NULL){
  if(is.null(ylim)){
    ylim <- range(exprs(eset), na.rm=TRUE)
  }

  if(new.window) dev.new()
  oldpar <- par(mfrow = mfrow)
  on.exit(par(oldpar))

  expr <- exprs(eset)
  cnum <- max(cl$cluster)

  for(i in 1:cnum){

    members <- expr[cl$cluster == i, , drop=FALSE]
    memval  <- cl$membership[cl$cluster == i, i]

    keep <- memval >= min.mem
    members <- members[keep, , drop=FALSE]
    memval  <- memval[keep]

    if(nrow(members) == 0){
      plot(0,0,type="n", main=paste("Cluster",i))
      next
    }

    ## ----- KEY PART -----
    ## order lines by membership so similar colors group
    ord <- order(memval)
    members <- members[ord, , drop=FALSE]
    memval  <- memval[ord]

    ## map membership -> color index
    idx <- round(memval * (length(colo)-1)) + 1
    line_cols <- colo[idx]
    ## ---------------------

    matplot(
      t(members),
      type="l",
      lty=1,
      col=line_cols,
      xaxt="n",
      ylim=ylim,
      main=paste("Cluster", i)
    )

    axis(1, at=1:ncol(members), labels=time.labels)
    axis(2)
  }
}


# Plot
hematopo_order <- c("LSPC", "GMP", "ProMono", "Mono")

# Create a gradient from grey to blue with 50 colors
custom_colors1 <- colorRampPalette(c("#f8f8f8","#c8c8c8fe","#a0a0a0fe", "#e4d479","#cd9a35","#de6017"))(100)
custom_colors2 <- colorRampPalette(c("#f8f8f8","#c8c8c8fe","#a0a0a0fe","#bdd9ff","#5898f7ff","#4354c6ff"))(100)

yr <- round(range(exprs(eset), na.rm = TRUE), digits = 1)


pdf("./mfuzz.plot.fixedaxis1.pdf", width=18, height=8)
par(cex.axis=1.6, cex.lab=1.6, cex=1.6, las=1, cex.main=2)

my_mfuzz_plot(
  eset,
  cl,
  mfrow=c(2,5),
  time.labels=hematopo_order,
  colo=custom_colors1,
  ylim=yr,
  new.window=FALSE
)
dev.off()

pdf("./mfuzz.plot.fixedaxis2.pdf", width=18, height=8)
par(cex.axis=1.6, cex.lab=1.6, cex=1.6, las=1, cex.main=2)

my_mfuzz_plot(
  eset,
  cl,
  mfrow=c(2,5),
  time.labels=hematopo_order,
  colo=custom_colors2,
  ylim=yr,
  new.window=FALSE
)
dev.off()


# Extract increasing expression genes (e.g., cluster 2)
dataList <- list()
for(i in sort(unique(cl$cluster))){
  increasing_genes <- names(cl$cluster[cl$cluster == i])
  increasing_genes <- increasing_genes[-grep(pattern = "^LINC",increasing_genes)]
  increasing_genes <- increasing_genes[-grep(pattern = "[.]",increasing_genes)]
  dataList[[i]] <- increasing_genes
  names(dataList)[i] <- paste0("cluster_", i)
}

max_len <- max(lengths(dataList))
dataList <- lapply(dataList, `length<-`, max_len)
data <- do.call(cbind, dataList)

fwrite(data, "./cluste_all_DH_up_gene.csv")
