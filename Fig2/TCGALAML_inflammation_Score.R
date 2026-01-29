## inflammation score for LAML / BEATAML / scRNA-seq

rm(list = ls())
gc()

pg <- c('ggplot2', 'dplyr', 'patchwork',  'data.table', 'patchwork', 'DESeq2', 
  'edgeR', 'GSVA', 'GSEABase', 'mclust', 'ggridges', 'Seurat', 'ggthemes')
for(i in pg){
    suppressMessages(library(i, character.only = T))
}

setwd("/Volumes/UBUNTU 22_0")

## get inflammation Score of TCGA LAML
TCGA_clinical <- fread("./Sequencing_Analysis/LAML/ssgsea_clnical.csv", data.table =F)
TCGA_clinical <- TCGA_clinical[,-1]
rownames(TCGA_clinical) <- TCGA_clinical$Specimen

TCGA_counts <- fread("./Sequencing_Analysis/LAML/TCGA_LAML_151_RNA-seq_rawCounts.csv", data.table =F)
rownames(TCGA_counts) <- TCGA_counts$V1
TCGA_counts <- TCGA_counts[,-1]
TCGA_counts[1:4,1:4]
raw_counts <- TCGA_counts[apply(TCGA_counts,1,sum) > 3,] # rm low expression genes 
colnames(raw_counts) <- substr(colnames(raw_counts), 1, 12)


samples_to_keep <- intersect(TCGA_clinical$Specimen, colnames(raw_counts))  # Replace with actual column name

# Subset raw counts to match
counts_subset <- raw_counts[, samples_to_keep]
TCGA_clinical <- TCGA_clinical[samples_to_keep,]
dim(counts_subset);dim(TCGA_clinical)


## get group of inflammation high and low
head(TCGA_clinical)
summary(TCGA_clinical$HALLMARK_INFLAMMATORY_RESPONSE)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.1461  0.1790  0.2959  0.3155  0.4498  0.8539 
# Assuming TCGA_clinical is a data.frame and has the column HALLMARK_INFLAMMATORY_RESPONSE
q1 <- quantile(TCGA_clinical$HALLMARK_INFLAMMATORY_RESPONSE, 0.25)
q3 <- quantile(TCGA_clinical$HALLMARK_INFLAMMATORY_RESPONSE, 0.75)

TCGA_clinical$Inflammatory_Status <- with(TCGA_clinical, ifelse(
  HALLMARK_INFLAMMATORY_RESPONSE < q1, "Low",
  ifelse(HALLMARK_INFLAMMATORY_RESPONSE > q3, "High", "Median")
))


## DEG of high and low
inflammation_subset <- TCGA_clinical %>% 
dplyr::filter(Inflammatory_Status %in% c("High", "Low"))
table(inflammation_subset$Inflammatory_Status) ## check if top 25%, bottom 25%
rownames(inflammation_subset) <- inflammation_subset$Specimen
counts_subset <- counts_subset[,rownames(inflammation_subset)]

# Ensure clinical data is in the same order as columns in counts
dds <- DESeqDataSetFromMatrix(countData = counts_subset,
                              colData = inflammation_subset,
                              design = ~ Inflammatory_Status)

# Run DESeq2
dds <- DESeq(dds)

# Extract High vs Low results
res <- results(dds, contrast = c("Inflammatory_Status", "High", "Low"))

head(as.data.frame(res))
res["RBM47",]

# Filter for significant DEGs
res_sig <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]

# View top DEGs
head(res_sig[order(res_sig$padj), ])
write.csv(as.data.frame(res), "./Sequencing_Analysis/TCGALAML_Inflammation_High_vs_Low_DEG_all.csv", row.names = T)
write.csv(as.data.frame(res_sig), "./Sequencing_Analysis/TCGALAML_Inflammation_High_vs_Low_DEG_significant.csv",row.names = T)


