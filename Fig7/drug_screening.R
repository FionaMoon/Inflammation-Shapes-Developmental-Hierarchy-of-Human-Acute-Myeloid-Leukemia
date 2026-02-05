## Run GSEA
rm(list = ls())
gc()

pg <- c('ggplot2', 'dplyr', 'patchwork',  'data.table', 'patchwork', 
  'edgeR', 'GSVA', 'GSEABase', 'mclust', 'ggridges', 'ggpubr','ggthemes')
for(i in pg){
    suppressMessages(library(i, character.only = T))
}

D_path = "D:/xx/Sequencing_Analysis/BeatAML"

Exp_auc_prot <- fread("D:/GitHub/EAGLE/data/BeatAML/RNA-seq_raw_has_Ven_AUC_ProteinCoding.csv", data.table = F)
rownames(Exp_auc_prot) <- Exp_auc_prot$display_label
Exp_auc_prot[1:3,1:6]
Exp_auc_prot <- Exp_auc_prot[,-c(1:4)]
dim(Exp_auc_prot)
table(rowSums(Exp_auc_prot==0)==ncol(Exp_auc_prot)) # check genes which not express
raw_counts <- Exp_auc_prot[apply(Exp_auc_prot,1,sum) > 3,] # rm low expression genes 
dim(raw_counts) ## check group and expression genes


# Convert to DGEList object
dge <- DGEList(counts = raw_counts)

# Compute Counts Per Million (CPM)
cpm_matrix <- edgeR::cpm(dge, log = FALSE)

# Apply log2 transformation
log2_cpm_matrix <- log2(cpm_matrix + 1)

# Read GMT file (HALLMARK_INFLAMMATORY_RESPONSE.gmt)
gmt_file <- "D:/xx/Sequencing_Analysis/BeatAML/HALLMARK_INFLAMMATORY_RESPONSE.v2024.1.Hs.gmt"
gene_set <- getGmt(gmt_file)

# Run GSVA
# Create a GSVA parameter object using ssgseaParam() for ssGSEA
gsva_param <- ssgseaParam(log2_cpm_matrix, gene_set)

# Run GSVA using the parameter object
gsva_scores <- gsva(gsva_param)
class(gsva_scores)

# Save GSVA results
gsva_scores <- as.data.frame(t(gsva_scores))
head(gsva_scores)
write.csv(gsva_scores, 
"D:/xx/Sequencing_Analysis/BeatAML/GSVA_HALLMARK_INFLAMMATORY_RESPONSE.csv", 
row.names = TRUE)

# Convert GSVA scores to a vector for histogram plotting
gsva_values <- as.vector(gsva_scores$HALLMARK_INFLAMMATORY_RESPONSE)

# Fit a Gaussian Mixture Model (GMM) to identify 2 peaks
gmm <- Mclust(gsva_values, G = 2)  # Specify G = 2 for two components

gsva_scores$Cluster <- as.factor(gmm$classification)

gsva_scores <- gsva_scores %>%
  mutate(Inflammatory_Status = ifelse(Cluster == 1, "Negative", "Positive"))

gsva_scores$project <- "BeatAML"

summary(gsva_values)
# >   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.1291  0.3803  0.4956  0.5160  0.6399  1.1291

# Plot the density ridges with the GMM classification
pdf("D:/xx/Sequencing_Analysis/BeatAML/GMM.pdf", width = 6, height = 4)
  ggplot(gsva_scores, aes(x = HALLMARK_INFLAMMATORY_RESPONSE, y = project, fill = Inflammatory_Status)) +
  geom_density_ridges(aes(point_color = Inflammatory_Status, point_fill = Inflammatory_Status, point_shape = Inflammatory_Status),
    scale = 1, alpha = 0.5, na.rm = TRUE, point_alpha = 1, jittered_points = F) +
  scale_fill_manual(values = c("Negative" = "lightblue", "Positive" = "orange")) +
  scale_point_color_hue(l = 40) +
  scale_x_continuous(
    breaks = seq(0, 1.2, by = 0.2),  # Set the x-axis breaks at intervals of 0.1
    labels = scales::label_number(accuracy = 0.1)  # Format axis labels as numbers
  ) +
  theme_few() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 12),
        legend.background = element_rect(color = NA)) +
  labs(title = "Inflammatory_Score", fill = "Inflammatory_Status")
dev.off()


# ggplot(gsva_scores, aes(x = HALLMARK_INFLAMMATORY_RESPONSE, y = project, fill = Inflammatory_Status)) +
#   geom_density_ridges(
#     aes(point_color = Inflammatory_Status, point_fill = Inflammatory_Status, point_shape = Inflammatory_Status),
#     scale = 1, alpha = 0.5, na.rm = TRUE, point_alpha = 1, jittered_points = TRUE, point_size = 0.6
#   ) +
#   scale_fill_manual(values = c("Negative" = "lightblue", "Positive" = "orange")) +
#   scale_point_color_hue(l = 40) +
#   scale_discrete_manual(aesthetics = "point_shape", values = c(22, 24)) +
#   scale_discrete_manual(aesthetics = "point_color", values = c("Negative" = "lightblue", "Positive" = "orange")) +
#   scale_discrete_manual(aesthetics = "point_fill", values = c("Negative" = "lightblue", "Positive" = "orange")) +
#   scale_x_continuous(
#     breaks = seq(0, 0.8, by = 0.1),  # Set the x-axis breaks at intervals of 0.1
#     labels = scales::label_number(accuracy = 0.1)  # Format the labels to show one decimal place
#   ) +
#   theme_new(border = TRUE) +
#   coord_cartesian(clip = "off") +
#   labs(title = "Inflammatory_Score", fill = "Inflammatory_Status")
  
# Plot histogram with GMM density curve overlay
pdf("D:/xx/Sequencing_Analysis/BeatAML/histogram.pdf")
# Plot histogram with GMM density and show 2 peaks in different colors
  ggplot(gsva_scores, aes(x=HALLMARK_INFLAMMATORY_RESPONSE)) + 
      geom_histogram(aes(y=..density..), colour="grey", fill="white", bins = 100, alpha = 0.6)+
      geom_density(alpha=.2, fill="#FF6666") +
      theme_few() +
      theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size=14),
        axis.title  = element_text(size=16))
dev.off()


## difine neg < 0.4, pos > 0.8
gsva_scores$mannual_define_info <- ifelse(gsva_scores$HALLMARK_INFLAMMATORY_RESPONSE < 0.4, "Negative",
                               ifelse(gsva_scores$HALLMARK_INFLAMMATORY_RESPONSE > 0.8, "Positive", "Intermediate"))
table(gsva_scores$mannual_define_info)
# Intermediate     Negative     Positive
#          225          111           28

gsva_scores <- gsva_scores %>% dplyr::filter(mannual_define_info %in% c("Negative", "Positive"))

Drug_IC_RNA_VEN <- fread("D:/GitHub/EAGLE/data/BeatAML/Drug_auc_Venetoclax.csv", data.table = F)
rownames(Drug_IC_RNA_VEN) <- Drug_IC_RNA_VEN$R_ID
meta <- Drug_IC_RNA_VEN[rownames(gsva_scores),]
table(meta$Venetoclax)

if(all(rownames(gsva_scores) == meta$R_ID)){
  meta <- cbind(meta, gsva_scores)
  head(meta)
}


pdf("D:/xx/Sequencing_Analysis/BeatAML/boxplot.pdf", width = 6,height =4)
  ggviolin(meta, 
  x = "mannual_define_info", 
  y = "auc",
  fill = "mannual_define_info",
  add = c("jitter", "mean_sd"),
  shape = "mannual_define_info",
  palette = c("lightblue", "orange"),
  alpha = 0.7) +
  ylab("Venetoclax AUC") +
  xlab("Inflammatory_Score") +
  theme_few() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 12),
        legend.background = element_rect(color = NA)) +
        stat_compare_means(comparisons = list(c("Negative", "Positive")), 
        label = "p.signif")
dev.off()

