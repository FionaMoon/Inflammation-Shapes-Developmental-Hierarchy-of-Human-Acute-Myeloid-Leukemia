gc()

library(ggplot2)
library(data.table)
library(dplyr)
library(ggrepel)
library(RColorBrewer)

c35 <- c(
  '#FFC312','#12CBC4','#ED4C67','#A3CB38',
  '#F79F1F','#1289A7','#D980FA','#B53471','#cd6133',
  '#EE5A24','#009432','#0652DD','#9980FA','#833471', 
  '#EA2027','#006266','#5758BB','#d1ccc0','#C4E538',
  '#40407a','#706fd3','#34ace0','#FDA7DF',
  '#474787','#aaa69d','#227093','#218c74',
  '#ff5252','#ff793f','#ffda79','#ccae62',
  '#b33939','#84817a','#cc8e35','#33d9b2'
)

setwd("D:/eCLIP/20250904-RBM47_reproducible_enriched_windows")

list.files(pattern = ".csv")

data <- fread("./DOX_RBM47.reproducible_enriched_windows.csv", data.table = F)
colnames(data)
# # table(data$feature_types)
# data$feature_first <- sapply(strsplit(data$feature_types, ":"), "[",1)
# table(data$feature_first)
# table(data$feature_type_top) ## first and top are equal

# log2fc mea>=3, q max<=0.05
data_filter <- data %>% filter(enrichment_l2or_mean >= 3 & q_max < 0.05)
# > dim(data)
# [1] 15426    30
# > dim(data_filter)
# [1] 2436   30

#################################################################################
df <- as.data.frame(table(data_filter$feature_type_top))
colnames(df) <- c("region", "count")

# Add percentage labels
df$percent <- round(100 * df$count / sum(df$count), 1)
df$label <- paste0(df$region, "\n", df$percent, "%")
pie_colors <- c35[1:length(unique(data_filter$feature_types))]

pdf("./filter_pie_window.pdf", width = 6, height = 8)
pie(df$count, labels = df$label, col = pie_colors,
    main = "Features Pie Chart")
legend("right", legend = df$region, fill = pie_colors, cex = 1.2, bty = "n")
dev.off()


#################################################################################
data2 <- fread("D:/20250818_18-08-2025-02-53-23/output/all_reads/DOX_RBM47-18-08-2025-02-52-42.RBM47_IP_1.all_reads_fractions_feature_data.tsv", data.table =F)
sum(data2$input);sum(data2$clip)

# Add percentage labels
data2$input_percent <- round(100 * data2$input, 2)
data2$input_label <- paste0(data2$feature_type_top, "\n", data2$input_percent, "%")
pie_colors <- c35[1:length(unique(data2$feature_type_top))]

pdf("./input_feature.pdf", width = 6, height = 4)
pie(data2$input_percent, labels = data2$input_label, col = pie_colors,
    main = "Input Features")
    # Add legend
legend("right", legend = data2$feature_type_top, fill = pie_colors, cex = 1.2, bty = "n")
dev.off()


# Add percentage labels
data2$clip_percent <- round(100 * data2$clip, 2)
data2$clip_label <- paste0(data2$feature_type_top, "\n", data2$clip_percent, "%")
pie_colors <- c35[1:length(unique(data2$feature_type_top))]

pdf("./clip_feature.pdf", width = 6, height = 4)
pie(data2$clip_percent, labels = data2$clip_label, col = pie_colors,
    main = "clip Features")
    # Add legend
legend("right", legend = data2$feature_type_top, fill = pie_colors, cex = 1.2, bty = "n")
dev.off()


##########################################################
library(tidyr)
library(ggthemes)

data2 <- fread(
  "D:/20250818_18-08-2025-02-53-23/output/all_reads/DOX_RBM47-18-08-2025-02-52-42.RBM47_IP_1.all_reads_fractions_feature_data.tsv", 
  data.table = FALSE
)

# Check totals
sum(data2$input)
sum(data2$clip)

# Reshape wide → long
df_long <- data2 %>%
  pivot_longer(
    cols = c(input, clip),
    names_to = "Sample",
    values_to = "Count"
  ) %>%
  group_by(Sample) %>%
  mutate(Percent = 100 * Count / sum(Count))

# Colors for features
pie_colors <- c35[1:length(unique(data2$feature_type_top))]

# Plot stacked barplot
pdf("./stacked_barplot.pdf", width = 4, height = 6)
ggplot(df_long, aes(x = Sample, y = Percent, fill = feature_type_top)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = pie_colors) +
  labs(x = "", y = "Percentage", fill = "Feature type") +
  theme_few(base_size = 14) +
  coord_polar("y", start=0)

dev.off()
