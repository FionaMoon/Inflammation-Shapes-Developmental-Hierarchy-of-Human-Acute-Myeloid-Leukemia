# Load package
library(SingleCellExperiment)
library(ggplot2)
library(viridis)
library(ggpubr)
library(Rtsne)
library(uwot)
library(ggthemes)
library(data.table)
library(dplyr)
library(patchwork)


# build theme function
galaxyTheme_black = function(base_size = 12, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Specify axis options
      axis.line = element_blank(),  
      axis.text.x = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.text.y = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.ticks = element_line(color = "white", linewidth = 0.2),  
      axis.title.x = element_text(size = base_size, color = "white", margin = margin(0, 10, 0, 0)),  
      axis.title.y = element_text(size = base_size, color = "white", angle = 90, margin = margin(0, 10, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "black"),  
      legend.key = element_rect(color = "white",  fill = "black"),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*0.8, color = "white"),  
      legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),  
      legend.position = "none",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL, 
      # Specify panel options
      panel.background = element_rect(fill = "black", color  =  NA),  
      panel.border = element_rect(fill = NA, color = "white"),  
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank(),  
      panel.spacing = unit(0.5, "lines"),   
      # Specify facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),  
      strip.text.x = element_text(size = base_size*0.8, color = "white"),  
      strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "black", fill = "black"),  
      plot.title = element_text(size = base_size*1.2, color = "white"),  
      plot.margin = unit(rep(1, 4), "lines")
    )
}


galaxy_DimPlot_Seurat <- function(
  seurat_obj,
  reduction = c("umap", "tsne"),
  split.by = NULL,
  facet_ncol = NULL,
  facet_title_size = 14
) {
  reduction <- match.arg(reduction)

  # --- extract embeddings ---
  coords <- as.data.frame(Embeddings(seurat_obj, reduction = reduction))
  colnames(coords)[1:2] <- c("Dim_1", "Dim_2")

  # --- add metadata ---
  meta <- seurat_obj@meta.data
  data <- cbind(coords, meta)

  p <- ggplot(data, aes(x = Dim_1, y = Dim_2)) +
    stat_density_2d(
      aes(fill = after_stat(density)),
      geom = "raster",
      contour = FALSE
    ) +
    geom_point(color = "white", size = 0.02) +
    scale_fill_viridis(option = "magma") +
    galaxyTheme_black() +
    theme_void()

  if (!is.null(split.by) && split.by %in% colnames(meta)) {
    p <- p +
      facet_wrap(as.formula(paste("~", split.by)), ncol = facet_ncol) +
      theme(
        strip.background = element_blank(),
        strip.text = element_text(
          size = facet_title_size,
          face = "bold",
          colour = "white"
        )
      )
  }

  return(p)
}


DimPlot_Seurat <- function(
  seurat_obj,
  reduction = c("umap", "tsne"),
  split.by = NULL,
  facet_ncol = NULL,
  facet_title_size = 12
) {
  reduction <- match.arg(reduction)

  # --- extract embeddings ---
  coords <- as.data.frame(Embeddings(seurat_obj, reduction = reduction))
  colnames(coords)[1:2] <- c("Dim_1", "Dim_2")

  # --- add metadata ---
  meta <- seurat_obj@meta.data
  data <- cbind(coords, meta)

  colorp <- c(
    "#FFFFFF", "#0028F6", "#027DFE", "#08CAFF", "#03FFE8",
    "#04FF99", "#00FF23", "#B4FF06", "#F3FC05", "#FE6700", "#FA0200"
  )

  p <- ggplot(data, aes(x = Dim_1, y = Dim_2)) +
    stat_density_2d(
      aes(fill = after_stat(density)),
      geom = "raster",
      n = 200,
      contour = FALSE
    ) +
    scale_fill_gradientn(colors = colorp, na.value = "#FFFFFF") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_few(base_size = 16) +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white")
    )

  if (!is.null(split.by) && split.by %in% colnames(meta)) {
    p <- p +
      facet_wrap(as.formula(paste("~", split.by)), ncol = facet_ncol) +
      theme(
        strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(
          size = facet_title_size,
          face = "bold"
        )
      )
  }

  return(p)
}
