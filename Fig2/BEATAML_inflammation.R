## inflammation score for LAML / BEATAML / scRNA-seq

rm(list = ls())
gc()

pg <- c('ggplot2', 'dplyr', 'patchwork',  'data.table', 'patchwork','DESeq2', 
  'edgeR', 'GSVA', 'GSEABase', 'mclust', 'ggridges', 'Seurat', 'ggthemes')
for(i in pg){
    suppressMessages(library(i, character.only = T))
}

setwd("D:/xx/Sequencing_Analysis/BeatAML")


list.files(pattern = ".txt|csv")

## get clinical data
BEATAML_clinical <- fread("./beataml_wv1to4_clinical.csv", data.table =F)
length(unique(BEATAML_clinical$dbgap_rnaseq_sample))

## get expression data
BEATAML_counts <- fread("./beataml_waves1to4_counts_dbgap.txt", data.table =F)
BEATAML_counts <- BEATAML_counts[!duplicated(BEATAML_counts$display_label),]
rownames(BEATAML_counts) <- BEATAML_counts$display_label
BEATAML_counts <- BEATAML_counts[,-c(1:4)]
raw_counts <- BEATAML_counts[apply(BEATAML_counts,1,sum) > 3,] # rm low expression genes 


## match expression data and clinical data
samples_to_keep <- intersect(BEATAML_clinical$dbgap_rnaseq_sample, colnames(raw_counts))  # 671

# Subset raw counts to match
counts_subset <- raw_counts[, samples_to_keep]
BEATAML_clinical_R <- BEATAML_clinical[BEATAML_clinical$dbgap_rnaseq_sample %in% samples_to_keep,]
dim(counts_subset);dim(BEATAML_clinical_R)
rownames(BEATAML_clinical_R) <- BEATAML_clinical_R$dbgap_rnaseq_sample
BEATAML_clinical_R[1:4,1:4]
BEATAML_clinical_R <- BEATAML_clinical_R[colnames(counts_subset),]

## calculate inflammation score
# Convert to DGEList object
dge <- DGEList(counts = counts_subset)
cpm_matrix <- edgeR::cpm(dge, log = FALSE) # Compute Counts Per Million (CPM)
log2_cpm_matrix <- log2(cpm_matrix + 1) # Apply log2 transformation
## log2(FPKM/TPM/RSEM/CPM + 1) is recommended, cpm is highly recommended
## raw counts or data without log transformation is not recommend

# Read GMT file (HALLMARK_INFLAMMATORY_RESPONSE.gmt)
gmt_file <- "./HALLMARK_INFLAMMATORY_RESPONSE.v2024.1.Hs.gmt"
gene_set <- getGmt(gmt_file)

## Run GSVA
gsva_param <- ssgseaParam(log2_cpm_matrix, gene_set) # Create a GSVA parameter object using ssgseaParam() for ssGSEA

# Run GSVA using the parameter object
gsva_scores <- gsva(gsva_param)
class(gsva_scores)

# Save GSVA results
gsva_scores <- as.data.frame(t(gsva_scores))
head(gsva_scores)

## combine clinical & inflammation score
if(all(rownames(gsva_scores) == rownames(BEATAML_clinical_R))){
  BEATAML_clinical_R$HALLMARK_INFLAMMATORY_RESPONSE <- gsva_scores$HALLMARK_INFLAMMATORY_RESPONSE
} else {
  stop("\n Not equal!!! \n")
}

## save Data
gmt_file <- "./GOBP_INFLAMMATORY_RESPONSE.v2024.1.Hs.gmt"
gene_set <- getGmt(gmt_file)

## Run GSVA
gsva_param <- ssgseaParam(log2_cpm_matrix, gene_set) # Create a GSVA parameter object using ssgseaParam() for ssGSEA

# Run GSVA using the parameter object
gsva_scores <- gsva(gsva_param)
class(gsva_scores)

# Save GSVA results
gsva_scores <- as.data.frame(t(gsva_scores))
head(gsva_scores)

## combine clinical & inflammation score
if(all(rownames(gsva_scores) == rownames(BEATAML_clinical_R))){
  BEATAML_clinical_R$GOBP_INFLAMMATORY_RESPONSE <- gsva_scores$GOBP_INFLAMMATORY_RESPONSE
} else {
  stop("\n Not equal!!! \n")
}


fwrite(BEATAML_clinical_R, "./BEATAML_clinical_RNA_filtered.csv")
write.csv(counts_subset, "./BEATAML_counts_with_clinical.csv", row.names = T)






## get group of inflammation high and low
summary(BEATAML_clinical_R$HALLMARK_INFLAMMATORY_RESPONSE)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.1461  0.1790  0.2959  0.3155  0.4498  0.8539 
# Assuming BEATAML_clinical_R is a data.frame and has the column HALLMARK_INFLAMMATORY_RESPONSE
q1 <- quantile(BEATAML_clinical_R$HALLMARK_INFLAMMATORY_RESPONSE, 0.25)
q3 <- quantile(BEATAML_clinical_R$HALLMARK_INFLAMMATORY_RESPONSE, 0.75)

BEATAML_clinical_R$Inflammatory_Status <- with(BEATAML_clinical_R, ifelse(
  HALLMARK_INFLAMMATORY_RESPONSE < q1, "Low",
  ifelse(HALLMARK_INFLAMMATORY_RESPONSE > q3, "High", "Median")
))


## DEG of high and low
inflammation_subset <- BEATAML_clinical_R %>% 
dplyr::filter(Inflammatory_Status %in% c("High", "Low"))
table(inflammation_subset$Inflammatory_Status) ## check if top 25%, bottom 25%
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
write.csv(as.data.frame(res), "./BEATAML_Inflammation_High_vs_Low_DEG_all.csv", row.names = T)
write.csv(as.data.frame(res_sig), "./BEATAML_Inflammation_High_vs_Low_DEG_significant.csv",row.names = T)


