################################################################################
# Title: Subclustering of tumor cells 
# Author: Paula
# Description:
#   Subclustering with Seurat tools in order to find specific subpopulations
#   in tumor cells
#
################################################################################

library(Seurat)
library(ggplot2)
library(tidyverse)
set.seed(1234)
setwd("/home/workdir/analysis_augie")
immune_combined <- readRDS("./20250120_labeled_immune.combined_edit.RDS")

table(immune_combined@meta.data$sample, immune_combined@meta.data$immuno )


tumor_cells <- subset(immune_combined, subset = immuno %in% c("Cancer cells", "Unannotated"))
tumor_cells <- FindVariableFeatures(tumor_cells, selection.method = "vst", nfeatures = 3000)

VariableFeaturePlot(tumor_cells)

tumor_cells <- ScaleData(tumor_cells, features = VariableFeatures(tumor_cells))
tumor_cells <- RunPCA(tumor_cells, features = VariableFeatures(tumor_cells), npcs = 30)
library(clustree)
tumor_cells <- FindNeighbors(tumor_cells, dims = 1:30)  # Recalcular los vecinos
tumor_cells <- FindClusters(tumor_cells, resolution = res)
tumor_cells <- RunUMAP(tumor_cells, dims = 1:30)


# CLUSTREE ----------------------------------------------------------------
setwd("/home/workdir/")

# Crear un dataframe vacío para almacenar los resultados
cluster_results <- data.frame()

# Realizar clustering en diferentes resoluciones
for (res in seq(0.025, 0.5, by = 0.025)) {
  tumor_cells <- FindNeighbors(tumor_cells, dims = 1:30)  # Recalcular los vecinos
  tumor_cells <- FindClusters(tumor_cells, resolution = res)
  
  # Agregar los resultados a la lista
  cluster_results <- rbind(cluster_results, 
                           data.frame(resolution = res, 
                                      cluster = tumor_cells$seurat_clusters, 
                                      cell = rownames(tumor_cells@meta.data)))
  
  # También podrías visualizar la UMAP para cada resolución, como ya lo estás haciendo
  tumor_cells <- RunUMAP(tumor_cells, dims = 1:30)
  p <- DimPlot(tumor_cells, reduction = "umap", label = TRUE) + ggtitle(paste0("Subclustering cancer cells with resolution: ", res))
  ggsave(paste0("./HDAC_tfm/fig/sub_clust_tumor/subclustering_v2/UMAP_subclustering_tumor_cancer_unann_",res,".png"), plot = p, width = 10, height = 8, dpi=600)
  print("Saved")
}

colnames(cluster_results) <- c("resolution", "cluster", "cell")
cluster_results_wide <- cluster_results %>%
  pivot_wider(names_from = resolution, values_from = cluster, names_prefix = "res_")

# Crear el gráfico de clustree
p <- clustree(cluster_results_wide, prefix = "res_", layout = "tree", node_text_size = 4) +
  theme(
    legend.position = "bottom",
    text = element_text(family = "Times New Roman"),
    plot.title = element_text(family = "Times New Roman"),
    axis.title = element_text(family = "Times New Roman"),
    axis.text = element_text(family = "Times New Roman")
  )

p <- clustree(cluster_results_wide, prefix = "res_", layout = "tree", node_text_size = 6) +
  theme_minimal(base_size = 14)  # sin fuente personalizada

# Guardar el gráfico en un archivo
ggsave("./HDAC_tfm/fig/sub_clust_tumor/subclustering_v2/clustree_res_tumor_cancer_unann_subclustering_plot_2.png", plot = p, width = 18, height = 15, dpi=600)


# SILHOUTTE -----------------------------------------------------------------
library(cluster)
# Calcular matriz de distancias entre células
cell_dists <- dist(tumor_cells@reductions$pca@cell.embeddings, method = "euclidean")

cluster_info <- tumor_cells@meta.data[,grepl(paste0(DefaultAssay(tumor_cells),"_snn_res"),
                                             colnames(tumor_cells@meta.data))] %>%
  dplyr::mutate_all(as.character) %>%
  dplyr::mutate_all(as.numeric)

plot_shilouette <- function(sil.obj, snn, print.summary = TRUE, ...){
  df <- as.data.frame(sil.obj[, 1:3], stringsAsFactors = TRUE)
  # order by cluster and by sil_width
  df <- df[order(df$cluster, df$sil_width), ]
  if(!is.null(rownames(df))){ df$name <- factor(rownames(df), levels = rownames(df))
  }else {df$name <- as.factor(1:nrow(df))}
  df$cluster <- as.factor(df$cluster)
  mapping <- aes(x = sil_width, y = name, fill = cluster)
  
  p <- ggplot(df, mapping) +
    geom_bar(stat="identity") +
    labs(y = "", x = "Silhouette width Si",
         title = paste0("Clusters silhouette plot for SNN resolution of ", snn),
         subtitle =   paste0("\n Average silhouette width: ", 
                             round(mean(df$sil_width), 4)))+
    geom_vline(xintercept = mean(df$sil_width), linetype = "dashed", color = "red" )+
    theme(axis.text.y=element_blank(),axis.ticks=element_blank())
  
  recap <-summary(sil.obj)
  n <- as.matrix(summary(sil.obj)$clus.sizes)
  sil.sum <- data.frame(cluster = names(recap$clus.avg.widths), size =  as.matrix(recap$clus.sizes),
                        ave.sil.width = round(recap$clus.avg.widths,4), stringsAsFactors = TRUE)
  if(print.summary) {print(sil.sum)}
  message( "SNN resolution of ", snn, " average silhouette width: ", round(recap$avg.width,4))
  return(p)
}


pdf("./HDAC_tfm/fig/sub_clust_tumor/subclustering_v2/silhouette_plots.pdf", width = 8, height = 6)

# Generar los plots
lapply(colnames(cluster_info), function(x) {
  if(length(unique(cluster_info[,x])) == 1) {
    return(0)
  }
  si <- silhouette(cluster_info[,x], cell_dists)
  plot_shilouette(si, x)
})

# Cerrar el dispositivo PDF para finalizar el archivo
dev.off()


# RESTO DE PLOTS ----------------------------------------------------------
res= 0.05
p2 <- UMAPPlot(tumor_cells, group.by = "RNA_snn_res.0.075", split.by = "sample", ncol = 4)
ggsave(paste0("./HDAC_tfm/fig/sub_clust_tumor/subclustering_v2/UMAPbytreatments_",res,"_tumor_unann.png"), plot = p2, width = 10, height = 8, dpi=300)


Idents(tumor_cells) <- tumor_cells@meta.data$RNA_snn_res.0.05


tumor_cells <- JoinLayers(tumor_cells)
markers <- FindAllMarkers(tumor_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

res = 0.05
library(dplyr)

top_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))  # Ordenar los genes dentro de cada cluster por avg_log2FC


plot <- DotPlot(tumor_cells, features = unique(top_markers$gene))
ggsave(paste0("./HDAC_tfm/fig/sub_clust_tumor/subclustering_v2/Dotplot_subclustering_tumor_unnan_",res,".png"), plot = plot + coord_flip(), width = 10, height = 8, dpi=600)


p <- DoHeatmap(tumor_cells, features = top_markers$gene) + NoLegend()
ggsave(paste0("./HDAC_tfm/fig/sub_clust_tumor/subclustering_v2/Heatmap_subclustering_tumor_unnan_",res,".png"), plot = p, width = 10, height = 8, dpi=600)


library(UpSetR)

markers_by_cluster <- split(markers$gene, markers$cluster)

# Guardar el gráfico utilizando png o pdf
png(paste0("./HDAC_tfm/fig/sub_clust_tumor/subclustering_v2/venn_upset_tumor_unann", res, ".png"), width = 25, height = 10, units = "in", res = 600)
upset(
  fromList(markers_by_cluster),
  sets = names(markers_by_cluster),
  order.by = "freq",
  main.bar.color = "blue",
  sets.bar.color = "green",
  nsets = length(markers_by_cluster),
  nintersects = NA,
  text.scale = c(2, 2, 2, 2) 

dev.off()  


# SIGNATURES ------------------------------------------------------------------

library(IOBR)
tumor_cells_copy <- tumor_cells
rownames(tumor_cells_copy) <- toupper(rownames(tumor_cells))


gene_sets <- list(
  "WNT_target" = signature_tme$WNT_target,
  "Cell_cycle" = signature_tme$Cell_cycle, 
  "CellCycle_Reg" = signature_tme$CellCycle_Reg,
  "FGFR3_related" = signature_tme$FGFR3_related,
  "Mismatch_Repair" = signature_tme$Mismatch_Repair,
  "Homologous_recombination" = signature_tme$Homologous_recombination,
  "Nucleotide_excision_repair" = signature_tme$Nucleotide_excision_repair,
  "DNA_replication" = signature_tme$DNA_replication,
  "Base_excision_repair" = signature_tme$Base_excision_repair,
  "TMEscoreA_CIR" = signature_tme$TMEscoreA_CIR,
  "TMEscoreB_CIR" = signature_tme$TMEscoreB_CIR,
  "Co_inhibition_T_cell_Rooney_et_al" = signature_tme$Co_inhibition_T_cell_Rooney_et_al,
  "TMEscoreA_plus" = signature_tme$TMEscoreA_plus,
  "TMEscoreB_plus" = signature_tme$TMEscoreB_plus,
  "Interferon_Receptor_Li_et_al" = signature_tme$Interferon_Receptor_Li_et_al,
  "TGFb_Family_Member_Li_et_al" = signature_tme$TGFb_Family_Member_Li_et_al,
  "TGFb_Family_Member_Receptor_Li_et_al" = signature_tme$TGFb_Family_Member_Receptor_Li_et_al
  
)

gene_sets_tumor2 <- list(
  "Nature_metabolism_Hypoxia" = signature_tumor$Nature_metabolism_Hypoxia,
  "Winter_hypoxia_signature" = signature_tumor$Winter_hypoxia_signature,
  "Hu_hypoxia_signature" = signature_tumor$Hu_hypoxia_signature,
  "Molecular_Cancer_m6A" = signature_tumor$Molecular_Cancer_m6A,
  "MT_exosome" = signature_tumor$MT_exosome,
  "SR_exosome" = signature_tumor$SR_exosome,
  "Positive_regulation_of_exosomal_secretion" = signature_tumor$Positive_regulation_of_exosomal_secretion,
  "Negative_regulation_of_exosomal_secretion" = signature_tumor$Negative_regulation_of_exosomal_secretion,
  "Exosomal_secretion" = signature_tumor$Exosomal_secretion,
  "Exosome_assembly" = signature_tumor$Exosome_assembly,
  "Extracellular_vesicle_biogenesis" = signature_tumor$Extracellular_vesicle_biogenesis,
  "MC_Review_Exosome1" = signature_tumor$MC_Review_Exosome1,
  "MC_Review_Exosome2" = signature_tumor$MC_Review_Exosome2,
  "CMLS_Review_Exosome" = signature_tumor$CMLS_Review_Exosome,
  "Ferroptosis" = signature_tumor$Ferroptosis,
  "EV_Cell_2020" = signature_tumor$EV_Cell_2020
)
library(pheatmap)
library(grid)  # Para dibujar el heatmap correctamente


signature_plots <- function(database, signature) {
  cell_cycle_genes <- database
  
  # Añadir el score de la firma al objeto Seurat
  tumor_cells_copy <- AddModuleScore(object = tumor_cells_copy, 
                                     features = list(cell_cycle_genes), 
                                     name = signature)
  print(colnames(tumor_cells_copy@meta.data))
  
  
  # Definir la variable con la columna generada
  signature_col <- paste0(signature, "1")
  print(paste0("La signature es: ", signature))
  print(signature_col)
  
  # df <- FetchData(tumor_cells_copy, vars = c(signature_col, "sample"))
  # df$sample <- as.factor(df$sample)
  # # Verificamos las primeras filas de df para confirmar las columnas
  # print(head(df))
  # print(unique(df$sample))
  # str(df)
  # 
  # 
  # 
  # # Asegurarnos de que la columna 'sample' esté correctamente referenciada
  # stats_per_sample <- 
  #   group_by(df, sample) %>%
  #   summarise(
  #     mean_score = mean(df[[1]], na.rm = TRUE),  # Usamos la primera columna para la media
  #     sd_score = sd(df[[1]], na.rm = TRUE)       # Usamos la primera columna para la desviación estándar
  #   )
  # # Verificamos el resultado
  # print(stats_per_sample)
  # 
  # 
  # color_palette <- colorRampPalette(c("purple", "white", "yellow"))(50)
  # 
  # # Definir un rango fijo para la escala de colores
  # limite_max <- max(abs(stats_per_sample$mean_score), na.rm = TRUE)  
  # breaks_range <- seq(-limite_max, limite_max, length.out = 51)  
  # 
  # sample_levels <- c(
  #   "KAP25L_control", "KAP25L_PD-L1", "KAP25L_DAC", "KAP25L_Combination",
  #   "KAP25L-UV_control", "KAP25L-UV_PD-L1", "KAP25L-UV_DAC", "KAP25L-UV_Combination"
  # )
  # 
  # # Ordenar muestras
  # stats_per_sample$sample <- factor(stats_per_sample$sample, levels = sample_levels)
  # stats_per_sample <- stats_per_sample[order(stats_per_sample$sample), ]
  # 
  # # Convertir la tabla en una matriz para pheatmap
  # heatmap_matrix <- as.matrix(stats_per_sample$mean_score)
  # rownames(heatmap_matrix) <- stats_per_sample$sample
  # colnames(heatmap_matrix) <- "Mean Score"
  # 
  # # Generar el heatmap
  # p0 <- pheatmap(heatmap_matrix, 
  #                cluster_rows = FALSE,  # Mantiene el orden de sample_levels
  #                cluster_cols = FALSE,  
  #                color = color_palette,  
  #                breaks = breaks_range,  
  #                display_numbers = TRUE,  
  #                main = paste0("Heatmap de Expresión Media por Sample de ", signature))
  # 
  # # Guardar el heatmap correctamente
  # png(paste0("/home/workdir/HDAC_tfm/fig/sub_clust_tumor/subclustering_v2/signatures_0.05/Heatmap_tumor_cells_bysample_", signature, ".png"), 
  #     width = 10, height = 8, units = "in", res = 600)
  # grid.draw(p0$gtable)
  # dev.off()
  # 
  # Violin Plot
  p1 <- VlnPlot(tumor_cells_copy, features = signature_col, group.by = "sample")
  
  ggsave(paste0("/home/workdir/HDAC_tfm/fig/sub_clust_tumor/subclustering_v2/signatures_0.05/Violin_plot_tumor_cells_juntas_", signature, ".png"), 
         plot = p1, width = 10, height = 8, dpi = 600)
  
  #Feature Plot en UMAP
  p2 <- FeaturePlot(tumor_cells_copy, 
                    features = signature_col, 
                     reduction = "umap",
                     pt.size = 0.5)
  p1 <- VlnPlot(tumor_cells_copy, features = signature_col, group.by = "RNA_snn_res.0.05")
  
  # # Guardar imágenes
  ggsave(paste0("./HDAC_tfm/fig/sub_clust_tumor/subclustering_v2/signatures_0.05/Violin_plot_", signature, ".png"), 
           plot = p1, width = 10, height = 8, dpi = 600)
  ggsave(paste0("./HDAC_tfm/fig/sub_clust_tumor/subclustering_v2/signatures_0.05/Feature_plot_", signature, ".png"), 
           plot = p2, width = 10, height = 8, dpi = 600)
  # 
  # # Extraer datos para ggplot
  df <- FetchData(tumor_cells_copy, vars = c(signature_col, "RNA_snn_res.0.05", "sample"))
  # # Extraer el prefijo K25L o K25L_UV
  df$category <- ifelse(grepl("UV", df$sample), "KAP25L_UV", "KAP25L")
  # 
  # # Ordenar factor por categorías
  df$category <- factor(df$category, levels = c("KAP25L", "KAP25L_UV"))
  # 
  # # Ordenar los samples dentro de cada categoría
  df$sample <- factor(df$sample, levels = c(
    "KAP25L_control", "KAP25L_PD-L1", "KAP25L_DAC", "KAP25L_Combination",
    "KAP25L-UV_control", "KAP25L-UV_PD-L1", "KAP25L-UV_DAC", "KAP25L-UV_Combination"
  ))
  # 
  # # Violin Plot separado por sample y con filas por categoría
  p3 <- ggplot(df, aes(x = RNA_snn_res.0.05, y = .data[[signature_col]], fill = RNA_snn_res.0.05)) +
    geom_violin() +
    facet_wrap(~ sample, ncol = 4) +  # Filas por K25L/K25L_UV y columnas por sample
    theme_minimal() +
    theme(legend.position = "none")
  ggsave(paste0("./HDAC_tfm/fig/sub_clust_tumor/subclustering_v2/signatures_0.05/Violin_plot_", signature, "_bytreats.png"), 
         plot = p3, width = 10, height = 8, dpi = 600)
  # 


}

# Aplicar la función a todas las firmas en gene_sets
lapply(names(gene_sets_tumor2), function(signature_name) {
  signature_plots(gene_sets_tumor2[[signature_name]], signature_name)
})


heatmap_plots <- function(database, signature) {
  cell_cycle_genes <- database
  # Añadir el score de la firma al objeto Seurat
  tumor_cells_copy <- AddModuleScore(object = tumor_cells_copy, 
                                     features = list(cell_cycle_genes), 
                                     name = signature)
  
  signature_col <- paste0(signature, "1")
  print(paste0("La signature es: ", signature))
  
  df <- FetchData(tumor_cells_copy, vars = c(signature_col, "sample"))
  
  # Verificamos las primeras filas de df para confirmar las columnas
  print(head(df))
  print(unique(df$sample))
  str(df)
  
  # Asegurarnos de que la columna 'sample' esté correctamente referenciada
  # Usamos 'signature_col' directamente para referirnos a la columna correcta
  stats_per_sample <- df %>%
    group_by(sample) %>%
    summarise(
      mean_score = mean(df[[signature_col]], na.rm = TRUE),  # Calculamos la media de la firma
      sd_score = sd(df[[signature_col]], na.rm = TRUE)       # Desviación estándar
    )
  
  print(stats_per_sample)
  
  # Aquí puedes agregar el código para generar el heatmap si lo necesitas
}
