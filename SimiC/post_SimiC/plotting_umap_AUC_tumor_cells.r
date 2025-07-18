library(reticulate)
reticulate::use_python("/usr/bin/python3")
library(Seurat)
library(cowplot)
library(dplyr)
library(future)
library(viridis)
library(ggplot2)
library(plyr)
library(reshape2)
library(gridExtra)
library(ggridges)
library(stringr)
library(hues)
library(pheatmap)
library(Rtsne)



get_plot_dims <- function(heat_map)
{
  plot_height <- sum(sapply(heat_map$gtable$heights, grid::convertHeight, "in"))
  plot_width  <- sum(sapply(heat_map$gtable$widths, grid::convertWidth, "in"))
  dev.off()
  return(list(height = plot_height, width = plot_width))
}



plot_root_dir <- "/home/workdir/HDAC_tfm/NEW_MAT/K25L/Tumor/outputSimic/figures/"
annotations <- read.csv('/home/workdir/HDAC_tfm/NEW_MAT/K25L/Tumor/inputFiles/all_phenotype_annotation.csv')
annotations$celltype = "Cancer cells"
# data.frame having the annotation per cell with columns: .id = phenotypes, cell_id = cell identifier, cluster_id = cluster of the cell
plasma <- viridis(50, direction = 1, option = "C")

#### NEED THE INPUT OF THE SEURAT DATA WITH THE NAME OF THE PHENOTYPE
seurat_data <- readRDS('/home/workdir/HDAC_tfm/data/tumor_cells_adjusted_13032025.rds')
# The phenotypes studied, be aware that must have the same order to the Simic input 
phenotypes <- c('KAP25L_control','KAP25L_PD-L1','KAP25L_DAC', 'KAP25L_Combination')

# binaryze the phenotypes with base 0
assignment <- as.character(seq(0,length(phenotypes)-1))
names(assignment) <- phenotypes

# Load the weigths from SimiC
weights_file <- '/root/workdir/HDAC_tfm/NEW_MAT/KUV/Tumor/outputSimic/matrices/Output_Tumor_KAP25L_all_L1_0.01_L2_0.001_2_Ws_filtered_BIC.pickle'
SimiC_weights <- py_load_object(filename =weights_file)



auc_dir <- '/home/workdir/HDAC_tfm/NEW_MAT/K25L/Tumor/outputSimic/matrices/Output_Tumor_K25L_all_L1_0.01_L2_0.001_2_filtered_AUCs.pickle'

# Cargar el archivo pickle directamente
SimiC_AUCs <- py_load_object(filename = auc_dir)
AUCs <- list()

for(i in seq_along(phenotypes)){ 
    dims <- dim(py_to_r(SimiC_AUCs[[i]]))
    tmp_df <- py_to_r(SimiC_AUCs[[i]])
    tmp_list <- list()
    for (j in seq_along(tmp_df)){
       tmp_list[[j]] <-   tmp_df[[j]]$tolist()
    }
    cell_indexes <- attr(tmp_df, "pandas.index")
    cellindexes_pylist <- py_to_r(cell_indexes$values)
    cell_ids_list <- list()
    for (j in seq_along(cellindexes_pylist)){
       cell_idx  <- j-1
       cell_ids_list[[j]] <- as.character(cellindexes_pylist[[cell_idx]])
    }
    df_auc_all <- do.call(cbind, tmp_list)
    cell_ids <- do.call(c, cell_ids_list)
    colnames(df_auc_all) <- names(tmp_df)
    rownames(df_auc_all) <- cell_ids
    cell_ids_real <- annotations$Cell[annotations$phenotype == phenotypes[i]]
    df_auc_all <- df_auc_all[rownames(df_auc_all) %in% cell_ids_real, ]
    tmp <- t(df_auc_all)
    tmp <- melt(tmp, variable.name = 'value', stringsAsFactors =F)
    tmp$cluster_id <- phenotypes[i]
    colnames(tmp) <- c('driver', 'cell_id', 'value', "cluster_id")
    AUCs[[i]] <- tmp
}
df_AUC <- as.data.frame(do.call(rbind,AUCs))

#####



colnames(df_AUC)[which(names(df_AUC) == "cluster_id")] <- "phenotype"
df_AUC$cluster_id <- NA 
df_AUC <- df_AUC %>%
  left_join(annotations %>%
              select(Cell), by = c("cell_id" = "Cell")) 

df_AUC$phenotype <- gsub("KAP25L", "KPB25L", df_AUC$phenotype)




MinMax_clust<-NULL
tmp <- df_AUC
MinMax_val <- NULL
clusters <- unique(annotations$celltype)
for(cluster in clusters){
  tmp <- df_AUC  #[df_AUC$cluster_id==cluster,]
  MinMax_val <- NULL
  for (tf in unique(df_AUC$driver)){
      n_breaks = 100
      Wauc_dist <- list()
      for (phenotype in phenotypes){
          Wauc_dist[[phenotype]] <- hist(tmp[tmp$phenotype == phenotype & tmp$driver == tf, 'value'], breaks = c(seq(0,1, 1/n_breaks)), plot=FALSE)$density
      }
      mat <-  do.call('rbind', Wauc_dist)
      mat <- mat[complete.cases(mat),,drop=FALSE]
      minmax_diff <- apply(na.omit(mat), 2, max) - apply(na.omit(mat), 2, min)
      variant <- sum(abs(minmax_diff)) / n_breaks
      variant <- variant/ sum(rowSums(mat)!=0)
      MinMax_val <- append(MinMax_val,variant)
  }
  MinMax_val <- setNames(as.data.frame(MinMax_val), c(cluster))
  rownames(MinMax_val) <-  unique(tmp$driver)
  if (is.null(MinMax_clust)){
    MinMax_clust <- MinMax_val
  }else{
    MinMax_clust<-cbind(MinMax_clust,MinMax_val)
  }
}


for(cluster in clusters){
  tmp <- df_AUC #[df_AUC$cluster_id==cluster,]
  pdf(paste0(plot_root_dir,'Simic_Auc_Cluster_ctrl_all_final_',cluster,'_.pdf'), width = 15, onefile = TRUE)
    plot_counter <- 1
    for (tf in unique(tmp$driver)){
      assign( paste0('p', plot_counter), 
      ggplot(tmp[tmp$driver == tf,], aes(x=value, fill=phenotype)) + 
      geom_density(alpha = 0.6, adjust = 1/8) + theme_classic() + 
      scale_fill_iwanthue() +
      theme(legend.position = 'top')+ geom_rug() + 
      ggtitle(paste0(tf, '   ', MinMax_clust[rownames(MinMax_clust) == tf, ])))
      #  )
      if(plot_counter == 2){
        grid.arrange(p1, p2, ncol=2)
        plot_counter <- 1
      }else{
        plot_counter <- plot_counter +1
      }
    }
    grid.arrange(p1, p2 , ncol=2)
  dev.off()
}

#KPB25L-UV
custom_colors <- c(
  "KPB25L-UV_control" = "#B0B0B0",     # Darker gray
  "KPB25L-UV_DAC" = "#E57373",         # Darker pink
  "KPB25L-UV_PD-L1" = "#648FFF",       # Darker blue
  "KPB25L-UV_Combination" = "#7E57C2"  # Darker purple
)

drivers_of_interest <- c("Zfp827", "Nfkb1")
tmp <- df_AUC[df_AUC$driver %in% drivers_of_interest, ]

# Asegurar que los niveles del factor coinciden con los del vector de colores
df_AUC$phenotype <- factor(df_AUC$phenotype, levels = names(custom_colors))


for (tf in drivers_of_interest){
    p <- ggplot(tmp[tmp$driver == tf, ], aes(x = value, fill = phenotype)) +
    geom_density(alpha = 0.6, adjust = 1/8) +
    theme_classic(base_family = "Times New Roman") +
    scale_fill_manual(values = custom_colors) +
    theme(
        legend.position = "top",
        text = element_text(family = "Times New Roman"),
        axis.title = element_text(family = "Times New Roman", size= 18),
        axis.text = element_text(family = "Times New Roman",size= 18),
        legend.text = element_text(family = "Times New Roman", size= 18),
        legend.title = element_text(family = "Times New Roman", size= 18),
        plot.title = element_text(family = "Times New Roman", hjust = 0, size= 20, face = "bold")
    ) +
    geom_rug() +
    ggtitle(paste0(tf, "   "))
  # Guardar imagen PNG
  ggsave(filename = paste0(plot_root_dir, "Simic_Auc_", tf, "_final_1707.png"),
         plot = p, width = 13, height = 10, dpi= 600)
}


genes_2_plot <- c("Zfp827", "Nfkb1")
phenotypes <- c('KPB25L-UV_control','KPB25L-UV_PD-L1','KPB25L-UV_DAC', 'KPB25L-UV_Combination')

umap_coords <- as.data.frame(seurat_data@reductions$umap@cell.embeddings)
umap_coords$cell_id <- rownames(umap_coords)
clusters <- unique(annotations$celltype)
# Agregar los clusters de resolución 0.095
umap_coords$cluster_id <- seurat_data@meta.data[umap_coords$cell_id, "Adjusted_Clusters"]
#umap_coords <- umap_coords[umap_coords$cluster_id %in% clusters, ]

# Definir directorio donde guardar las imágenes
plot_dir_umap <- "/home/workdir/HDAC_tfm/NEW_MAT/KUV/Tumor/outputSimic/figures/"
levels_order <- c('KPB25L-UV_control','KPB25L-UV_PD-L1','KPB25L-UV_DAC', 'KPB25L-UV_Combination')
for (gene in genes_2_plot) {
  png_filename <- paste0(plot_dir_umap, gene, "_UMAP_ctrl_all.png")
  
  # Filtrar solo células presentes en df_AUC
  plotter_tmp <- umap_coords[umap_coords$cell_id %in% df_AUC[df_AUC$phenotype %in% phenotypes, 'cell_id'], ]
  
  # Unir con la expresión del gen
  plotter_tmp <- merge(
    plotter_tmp,
    df_AUC[df_AUC$driver == gene, c('cell_id', 'value', 'phenotype')],
    by = 'cell_id'
  )

  if (nrow(plotter_tmp) > 0) {
    # Asegurar orden correcto del facet
    plotter_tmp$phenotype <- factor(plotter_tmp$phenotype, levels = levels_order)

    p <- ggplot(plotter_tmp) +
      geom_point(size = 0.8, alpha = 0.4, aes(x = umap_1, y = umap_2, color = value)) +
      theme_classic(base_family = "Times New Roman") +
      scale_color_viridis(option = 'inferno') +
      facet_wrap(~phenotype, ncol = 2, nrow = 2) +
      labs(x = '', y = gene, color = "Activity score") +
      theme(
        text = element_text(family = "Times New Roman", size = 18),
        axis.title = element_text(family = "Times New Roman", size = 18),
        axis.text = element_text(family = "Times New Roman", size = 18),
        legend.text = element_text(family = "Times New Roman", size = 18),
        legend.title = element_text(family = "Times New Roman", size = 18),
        strip.text = element_text(family = "Times New Roman",, size = 18, face = "bold")
      )

    ggsave(filename = png_filename, plot = p, width = 13, height = 8, dpi = 600)
  } else {
    print(paste("No data available for gene:", gene))
  }
}

## KPB25L

custom_colors <- c(
 "KPB25L_control" = "#B0B0B0",        # Light gray
 "KPB25L_DAC" = "#FFB6B6",            # Light pink
 "KPB25L_PD-L1" = "#A8C8FF",          # Light blue
 "KPB25L_Combination" = "#C1A9E0"    # Light purple
)

drivers_of_interest <- c("Bbx", "Foxp1", "Zeb2")
tmp <- df_AUC[df_AUC$driver %in% drivers_of_interest, ]

# Asegurar que los niveles del factor coinciden con los del vector de colores
df_AUC$phenotype <- factor(df_AUC$phenotype, levels = names(custom_colors))

for (tf in drivers_of_interest){
    p <- ggplot(tmp[tmp$driver == tf, ], aes(x = value, fill = phenotype)) +
    geom_density(alpha = 0.8, adjust = 1/8) +
    theme_classic(base_family = "Times New Roman") +
    scale_fill_manual(values = custom_colors) +
    theme(
        legend.position = "top",
        text = element_text(family = "Times New Roman"),
        axis.title = element_text(family = "Times New Roman", size= 18),
        axis.text = element_text(family = "Times New Roman",size= 18),
        legend.text = element_text(family = "Times New Roman", size= 18),
        legend.title = element_text(family = "Times New Roman", size= 18),
        plot.title = element_text(family = "Times New Roman", hjust = 0, size= 20, face = "bold")
    ) +
    geom_rug() +
    ggtitle(paste0(tf, "   "))
  # Guardar imagen PNG
  ggsave(filename = paste0(plot_root_dir, "Simic_Auc_", tf, "_final_1707.png"),
         plot = p, width = 13, height = 10, dpi= 600)
}

for (gene in genes_2_plot) {
  png_filename <- paste0(plot_dir_umap, gene, "_UMAP_ctrl_all_no_sep.png")
  
  # Filtrar solo células presentes en df_AUC
  plotter_tmp <- umap_coords[umap_coords$cell_id %in% df_AUC[df_AUC$phenotype %in% phenotypes, 'cell_id'], ]
  
  # Unir con la expresión del gen
  plotter_tmp <- merge(
    plotter_tmp,
    df_AUC[df_AUC$driver == gene, c('cell_id', 'value', 'phenotype')],
    by = 'cell_id'
  )

  if (nrow(plotter_tmp) > 0) {
    # Asegurar orden correcto del facet
    plotter_tmp$phenotype <- factor(plotter_tmp$phenotype, levels = levels_order)

    p <- ggplot(plotter_tmp) +
      geom_point(size = 0.8, alpha = 0.4, aes(x = umap_1, y = umap_2, color = value)) +
      theme_classic(base_family = "Times New Roman") +
      scale_color_viridis(option = 'inferno') +
      labs(x = '', y = gene, color = "Activity score") +
      theme(
        text = element_text(family = "Times New Roman", size = 18),
        axis.title = element_text(family = "Times New Roman", size = 18),
        axis.text = element_text(family = "Times New Roman", size = 18),
        legend.text = element_text(family = "Times New Roman", size = 18),
        legend.title = element_text(family = "Times New Roman", size = 18),
        strip.text = element_text(family = "Times New Roman",, size = 18, face = "bold")
      )

    ggsave(filename = png_filename, plot = p, width = 13, height = 8, dpi = 600)
  } else {
    print(paste("No data available for gene:", gene))
  }
}


## HEATMAP

MinMax_clust <- MinMax_clust[order(rowSums(MinMax_clust), decreasing = TRUE), , drop = FALSE]

p <- pheatmap(MinMax_clust,color=plasma, fontsize=5, angle_col =45, cellwidth=40, cluster_rows = FALSE, cluster_cols = FALSE)
plot_dims <- get_plot_dims(p)
      
pdf(paste0(plot_root_dir,'Simic_HeatMaps_RegDissScore_ctrl_all_final.pdf'), height = plot_dims$height, width = plot_dims$width )
print(p)
dev.off()


### UMAPs

genes_2_plot <- c("Bbx", "Foxp1", "Zeb2")
phenotypes <- c('KPB25L_control','KPB25L_PD-L1','KPB25L_DAC', 'KPB25L_Combination')

umap_coords <- as.data.frame(seurat_data@reductions$umap@cell.embeddings)
umap_coords$cell_id <- rownames(umap_coords)
clusters <- unique(annotations$celltype)
# Agregar los clusters de resolución 0.095
umap_coords$cluster_id <- seurat_data@meta.data[umap_coords$cell_id, "Adjusted_Clusters"]
#umap_coords <- umap_coords[umap_coords$cluster_id %in% clusters, ]

# Definir directorio donde guardar las imágenes
plot_dir_umap <- "/home/workdir/HDAC_tfm/NEW_MAT/K25L/Tumor/outputSimic/figures/"
levels_order <- c("KPB25L_control", "KPB25L_PD-L1", "KPB25L_DAC", "KPB25L_Combination")
for (gene in genes_2_plot) {
  png_filename <- paste0(plot_dir_umap, gene, "_UMAP_ctrl_all.png")
  
  # Filtrar solo células presentes en df_AUC
  plotter_tmp <- umap_coords[umap_coords$cell_id %in% df_AUC[df_AUC$phenotype %in% phenotypes, 'cell_id'], ]
  
  # Unir con la expresión del gen
  plotter_tmp <- merge(
    plotter_tmp,
    df_AUC[df_AUC$driver == gene, c('cell_id', 'value', 'phenotype')],
    by = 'cell_id'
  )

  if (nrow(plotter_tmp) > 0) {
    # Asegurar orden correcto del facet
    plotter_tmp$phenotype <- factor(plotter_tmp$phenotype, levels = levels_order)

    p <- ggplot(plotter_tmp) +
      geom_point(size = 0.8, alpha = 0.4, aes(x = umap_1, y = umap_2, color = value)) +
      theme_classic(base_family = "Times New Roman") +
      scale_color_viridis(option = 'inferno') +
      facet_wrap(~phenotype, ncol = 2, nrow = 2) +
      labs(x = '', y = gene, color = "Activity score") +
      theme(
        text = element_text(family = "Times New Roman", size = 18),
        axis.title = element_text(family = "Times New Roman", size = 18),
        axis.text = element_text(family = "Times New Roman", size = 18),
        legend.text = element_text(family = "Times New Roman", size = 18),
        legend.title = element_text(family = "Times New Roman", size = 18),
        strip.text = element_text(family = "Times New Roman",, size = 18, face = "bold")
      )

    ggsave(filename = png_filename, plot = p, width = 13, height = 8, dpi = 600)
  } else {
    print(paste("No data available for gene:", gene))
  }
}


for (gene in genes_2_plot) {
  png_filename <- paste0(plot_dir_umap, gene, "_UMAP_ctrl_all_no_sep.png")
  
  # Filtrar solo células presentes en df_AUC
  plotter_tmp <- umap_coords[umap_coords$cell_id %in% df_AUC[df_AUC$phenotype %in% phenotypes, 'cell_id'], ]
  
  # Unir con la expresión del gen
  plotter_tmp <- merge(
    plotter_tmp,
    df_AUC[df_AUC$driver == gene, c('cell_id', 'value', 'phenotype')],
    by = 'cell_id'
  )

  if (nrow(plotter_tmp) > 0) {
    # Asegurar orden correcto del facet
    plotter_tmp$phenotype <- factor(plotter_tmp$phenotype, levels = levels_order)

    p <- ggplot(plotter_tmp) +
      geom_point(size = 0.8, alpha = 0.4, aes(x = umap_1, y = umap_2, color = value)) +
      theme_classic(base_family = "Times New Roman") +
      scale_color_viridis(option = 'inferno') +
      labs(x = '', y = gene, color = "Activity score") +
      theme(
        text = element_text(family = "Times New Roman", size = 18),
        axis.title = element_text(family = "Times New Roman", size = 18),
        axis.text = element_text(family = "Times New Roman", size = 18),
        legend.text = element_text(family = "Times New Roman", size = 18),
        legend.title = element_text(family = "Times New Roman", size = 18),
        strip.text = element_text(family = "Times New Roman",, size = 18, face = "bold")
      )

    ggsave(filename = png_filename, plot = p, width = 13, height = 8, dpi = 600)
  } else {
    print(paste("No data available for gene:", gene))
  }
}
