rm(list = ls())
gc()

####################==== strackbarplot ====####################

source("/Users/MyFolder/Code/LLN-RBM47/CITE-seq_cytokine/Alluvial_stack_barplot.R")
D_path = "/Users/MyFolder/Code/LLN-RBM47/CITE-seq_cytokine"
F_path = "/Users/MyFolder/Code/LLN-RBM47/CITE-seq_cytokine/ZYG_script/V3"


Trt <- readRDS("/Users/MyFolder/Code/LLN-RBM47/CITE-seq_cytokine/ZYG/ZYG_PeterVG_anno_cytotrace_ADT_V5.rds")

## col of celltype and condition can be conserved
Trt_1 <- Trt@meta.data %>% select(PeterVG_anno, condition)

meta <- Trt_1
meta$Alluvial <- as.character(meta$PeterVG_anno)
unique(meta$Alluvial)
meta$Alluvial[grep(pattern = "B|T|NK", meta$Alluvial)] <- "Others"

table(meta$Alluvial)

celltype_order <- c(
  "LSPC-Quiescent", "LSPC-Primed", "LSPC-Cycle", "MEP",
  "GMP-like", "ProMono-like", "Mono-like", "cDC-like", "pDC-like",
  "Other"
)

meta$Alluvial <- factor(meta$Alluvial, levels = celltype_order)

base_cols <- c("#454fd2e8", "#698ff7ff", "#6bbbd8ff","#b5cf23c0",
  	"#ded232ff","#cd57b5ff","#cf465bc4", "#ed977aff", "#e86a40ff",
  	"lightgrey")

names(base_cols) <- celltype_order

p <- sc_Alluvial_pl(meta, 
  x_key = "condition",
  order_x = c("CT", "IL1B","TNFA", "Combo"),
  color_use = base_cols,
stratum_group = "Alluvial")

ggsave(
  paste0(F_path,'/Trt_alluvial.pdf'),
  p ,
  height = 4, width = 5
)
