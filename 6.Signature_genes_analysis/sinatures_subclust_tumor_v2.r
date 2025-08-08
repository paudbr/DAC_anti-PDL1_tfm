################################################################################
# Title: Tumor signature analysis
# Author: Paula
# Description:
#   This script studies the expression of different immune-tumor signatures and
#   subtypes signatures in tumor cell subclusters.
#
################################################################################


library(Seurat)
library(ggplot2)
library(tidyverse)
library(IOBR)
library(pheatmap)
library(viridis)

#Load seurat object with only tumor cells already clustered:

tumor_cells_clustered <- readRDS("/home/workdir/HDAC_tfm/data/tumor_cells_adjusted_13032025.rds")


tumor_cells_copy <- JoinLayers(tumor_cells_clustered)
rownames(tumor_cells_copy) <- toupper(rownames(tumor_cells_clustered))


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
  
  ggsave(paste0("/home/workdir/HDAC_tfm/fig/sub_clust_tumor/subclustering_v2/signatures_adjusted/Violin_plot_tumor_cells_juntas_", signature, ".png"), 
         plot = p1, width = 10, height = 8, dpi = 600)
  
  #Feature Plot en UMAP
  p2 <- FeaturePlot(tumor_cells_copy, 
                    features = signature_col, 
                    reduction = "umap",
                    pt.size = 0.5)
  p1 <- VlnPlot(tumor_cells_copy, features = signature_col, group.by = "Adjusted_Clusters")
  
  # # Guardar imágenes
  ggsave(paste0("./HDAC_tfm/fig/sub_clust_tumor/subclustering_v2/signatures_adjusted/Violin_plot_", signature, ".png"), 
         plot = p1, width = 10, height = 8, dpi = 600)
  ggsave(paste0("./HDAC_tfm/fig/sub_clust_tumor/subclustering_v2/signatures_adjusted/Feature_plot_", signature, ".png"), 
         plot = p2, width = 10, height = 8, dpi = 600)
  # 
  # # Extraer datos para ggplot
  df <- FetchData(tumor_cells_copy, vars = c(signature_col, "Adjusted_Clusters", "sample"))
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
  p3 <- ggplot(df, aes(x = Adjusted_Clusters, y = .data[[signature_col]], fill = Adjusted_Clusters)) +
    geom_violin() +
    facet_wrap(~ sample, ncol = 4) +  # Filas por K25L/K25L_UV y columnas por sample
    theme_minimal() +
    theme(legend.position = "none")
  ggsave(paste0("./HDAC_tfm/fig/sub_clust_tumor/subclustering_v2/signatures_adjusted/Violin_plot_", signature, "_bytreats.png"), 
         plot = p3, width = 10, height = 8, dpi = 600)
  # 
  
  
}

# Aplicar la función a todas las firmas en gene_sets
lapply(names(gene_sets), function(signature_name) {
  signature_plots(gene_sets[[signature_name]], signature_name)
  
})

# TNBC SUBTYPES -----------------------------------------------------------

LAR_subtype <- c("DHRS2", "GABRP", "AGR2", "PIP", "FOXA1", "PROM1", "TFF1", "NAT1", 
                 "BCL11A", "ESR1", "FOXC1", "CA12", "TFF3", "SCUBE2", "SFRP1", 
                 "ERBB4", "SIDT1", "PSAT1", "CHI3L1", "AR", "SERPINB5", "SOX10", 
                 "PDE9A", "KCNK5", "IRX1", "MIA", "DSC2", "ST3GAL6", "GZMB", 
                 "LAMP3", "GBP1", "CCL5", "EZH2", "PLAT", "TAP2","AIM2")

MES_subtype <- c("CD36", "OGN", "ABCA8", "CFD", "IGF1", "HBB", "CDH1", "MEOX2", 
                 "GPX3", "SCARA5", "PDK4", "ENPP2", "AGTR1", "LEP", "LPL", "DPT", 
                 "TIMP4", "FHL1", "SRPX", "EDNRB", "DSC2", "FGL2", "PTGER4", "SPON1", 
                 "PBK", "EZH2")

BLIS_subtype <-  c("DHRS2", "GABRP", "AGR2", "PIP", "FOXA1", "PROM1", "TFF1", "NAT1", 
                   "BCL11A", "ESR1", "FOXC1", "CA12", "TFF3", "SCUBE2", "SFRP1", 
                   "ERBB4", "SIDT1", "CD36", "OGN", "ABCA8", "CFD", "MEOX2", "ENPP2", 
                   "AGTR1", "DPT", "SERPINB5", "SOX10", "IRX1", "MIA", "DSC2", "TTYH1", 
                   "COL9A3", "FGL2", "RARRES3", "PDE9A", "BST2", "PTGER4", "KCNK5", 
                   "PSMB9", "HLA-DMA", "EPHB3", "IGSF6", "ST3GAL6", "RHOH", "SGPP1", 
                   "CXCL9", "CXCL11", "GBP5", "GZMB", "CCL5", "SLAMF7", "TAP1", "CD2", 
                   "AIM2")


genes <- c(
  "DHRS2", "GABRP", "AGR2", "PIP", "FOXA1", "PROM1", "TFF1", "NAT1", 
  "BCL11A", "ESR1", "FOXC1", "CA12", "TFF3", "SCUBE2", "SFRP1", "ERBB4", 
  "SIDT1", "PSAT1", "CHI3L1", "AR", "CD36", "OGN", "ABCA8", "CFD", "IGF1", 
  "HBB", "CDH1", "MEOX2", "GPX3", "SCARA5", "PDK4", "ENPP2", "AGTR1", "LEP", 
  "LPL", "DPT", "TIMP4", "FHL1", "SRPX", "EDNRB", "SERPINB5", "SOX10", "IRX1", 
  "MIA", "DSC2", "TTYH1", "COL9A3", "FGL2", "RARRES3", "PDE9A", "BST2", "PTGER4", 
  "KCNK5", "PSMB9", "HLA-DMA", "EPHB3", "IGSF6", "ST3GAL6", "RHOH", "SGPP1", 
  "CXCL9", "CXCL11", "GBP5", "GZMB", "LAMP3", "GBP1", "ADAMDEC1", "CCL5", "SPON1", 
  "PBK", "STAT1", "EZH2", "PLAT", "TAP2", "SLAMF7", "HERC5", "SPOCK1", "TAP1", 
  "CD2", "AIM2"
)

BLIA_subtype <- c("TFF3", "SCUBE2", "ERBB4", "CHI3L1", "AR", "OGN", "IGF1", "PDK4", 
                  "PSMB9", "CXCL11", "GBP5", "GZMB", "LAMP3", "GBP1", "ADAMDEC1", 
                  "CCL5", "SPON1", "PBK", "STAT1", "EZH2", "PLAT", "TAP2", "SLAMF7", 
                  "HERC5", "SPOCK1", "TAP1", "CD2", "AIM2")

all_genes <- c(LAR_subtype, MES_subtype, BLIS_subtype, BLIA_subtype)
genes_objeto <- toupper(rownames(GetAssayData(tumor_cells_copy@assays$integrated, layer = "scale.data" )))


genes_comun <- intersect(genes_objeto,genes)
# Extraer la expresión de los genes de tu objeto Seurat
gene_expression_data <- GetAssayData(tumor_cells_copy@assays$integrated, layer = "scale.data" )





gene_expression_data_ordered <- gene_expression_data[, order(clusters)]

# Crear una tabla de anotaciones para los clusters (mantén la misma orden de las células)
cluster_annotation <- data.frame(Cluster = clusters[order(clusters)])

# Definir la paleta de colores
plasma <- viridis(50, direction = 1, option = "C")

# Crear el heatmap con las células ordenadas por los clusters
pheatmap(gene_expression_data_ordered, 
         annotation_col = cluster_annotation,
         show_rownames = TRUE, 
         show_colnames = FALSE, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         color = plasma)
tumor_cells_copy <- AddModuleScore(object = tumor_cells_copy, 
                               features = LAR_subtype , 
                               name = "LAR_subtype")

tumor_cells_copy <- AddModuleScore(object = tumor_cells_copy, 
                              features = MES_subtype , 
                              name = "MES_subtype")

tumor_cells_copy <- AddModuleScore(object = tumor_cells_copy, 
                                   features = BLIS_subtype , 
                                   name = "BLIS_subtype")

tumor_cells_copy <- AddModuleScore(object = tumor_cells_copy, 
                                   features = BLIA_subtype , 
                                   name = "BLIA_subtype")

df <- FetchData(tumor_cells_copy, vars = c("BLIS_subtype1", "BLIA_subtype1", "Adjusted_Clusters", "sample"))



# Convertir filas en listas y combinarlas en una sola, eliminando NA
df$category <- ifelse(grepl("UV", df$sample), "KAP25L_UV", "KAP25L")

# Ordenar factor por categorías
df$category <- factor(df$category, levels = c("KAP25L", "KAP25L_UV"))

# Ordenar los samples dentro de cada categoría
df$sample <- factor(df$sample, levels = c(
  "KAP25L_control", "KAP25L_PD-L1", "KAP25L_DAC", "KAP25L_Combination",
  "KAP25L-UV_control", "KAP25L-UV_PD-L1", "KAP25L-UV_DAC", "KAP25L-UV_Combination"
))

ggplot(df, aes(x = "", y = BLIS_subtype1, fill = Adjusted_Clusters)) + 
  geom_violin() + 
  stat_summary(fun = median, geom = "point", color = "black", size = 2) +  # Marca la mediana
  facet_wrap(~ sample, ncol = 4) +  # Un gráfico por tipo de célula inmune
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),  # Oculta el eje X vacío
        axis.ticks.x = element_blank()) +
  labs(x = NULL, y = "Expresión",title = "BLIS_subtype1")  # Elimina el eje X vacío y etiqueta eje Y

ggplot(df, aes(x = "", y = BLIA_subtype1, fill = Adjusted_Clusters)) +
  geom_violin() +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  theme_minimal() +
  labs(title = "Distribution of BLIA Subtype Signature Expression Across Different Treatment Groups") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Genes PAM50 -------------------------------------------------------------
genes_PAM <- c(
  "SLC39A6", "MDM2", "MYC", "MLPH", "CXXC5", "BAG1", "FOXA1", "CENPF", 
  "MMP11", "BUB1B", "ERBB2", "PHGDH", "RRM2", "ESR1", "KRT17", "BIRC5", 
  "FOXC1", "CCNB1", "PTTG1", "GPR160", "KRT5", "UBE2C", "SFRP1", "ACTR3B", 
  "MKI67", "TYMS", "COX8", "MAPT", "PGR", "CEP55", "ANLN", "BCL2", "MELK", 
  "CDC6", "MYBL2", "KRT14", "EGFR", "FGFR4", "NAT1", "EXO1", "GRB7", 
  "CDC20", "TMEM45B", "ORC6", "MIA", "UBE2T", "CCNE1", "CDCA1", "KNTC2", 
  "KIF2C"
)



genes_objeto <- toupper(rownames(GetAssayData(tumor_cells_copy@assays$RNA, layer = "scale.data")))


genes_comun <- intersect(genes_objeto,genes_PAM)

clusters <- tumor_cells_copy$Adjusted_Clusters

gene_expression_data <- GetAssayData(tumor_cells_copy@assays$RNA, layer = "scale.data" )

filtered_gene_expression <- gene_expression_data[toupper(rownames(gene_expression_data)) %in% genes_comun, ]

filtered_gene_expression_ordered <- filtered_gene_expression[, order(clusters)]

# Crear una tabla de anotaciones para los clusters (mantén la misma orden de las células)
cluster_annotation <- data.frame(Cluster = clusters[order(clusters)])

# Definir la paleta de colores
plasma <- viridis(50, direction = 1, option = "C")

# Crear el heatmap con las células ordenadas por los clusters
pheatmap(filtered_gene_expression_ordered, 
         annotation_col = cluster_annotation,
         show_rownames = TRUE, 
         show_colnames = FALSE, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         color = plasma)



genes_objeto <- toupper(rownames(GetAssayData(tumor_cells_copy@assays$RNA, layer = "scale.data.1")))


genes_comun <- intersect(genes_objeto,genes_PAM)

clusters <- tumor_cells_copy$Adjusted_Clusters

gene_expression_data <- GetAssayData(tumor_cells_copy@assays$RNA, layer = "scale.data" )

filtered_gene_expression <- gene_expression_data[toupper(rownames(gene_expression_data)) %in% genes_comun, ]

filtered_gene_expression_ordered <- filtered_gene_expression[, order(clusters)]

# Crear una tabla de anotaciones para los clusters (mantén la misma orden de las células)
cluster_annotation <- data.frame(Cluster = clusters[order(clusters)])

# Definir la paleta de colores
plasma <- viridis(50, direction = 1, option = "C")

# Crear el heatmap con las células ordenadas por los clusters
pheatmap(filtered_gene_expression_ordered, 
         annotation_col = cluster_annotation,
         show_rownames = TRUE, 
         show_colnames = FALSE, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         color = plasma)


#######PRUEBA


genes_objeto <- toupper(rownames(GetAssayData(tumor_cells_copy@assays$RNA, layer = "counts")))

genes_comun <- intersect(genes_objeto,genes_PAM)

gene_expression_data <- GetAssayData(tumor_cells_copy@assays$RNA, layer = "counts" )[genes_comun,]

filtered_gene_expression_ordered <- gene_expression_data[, order(clusters)]

pheatmap(filtered_gene_expression_ordered,
         annotation_col = cluster_annotation,
         show_rownames = TRUE, 
         show_colnames = FALSE, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE,
         scale_row = TRUE,
         color = viridis(50),
         breaks = c(seq(0, 2, length.out = 20), seq(2, max(filtered_gene_expression_ordered), length.out = 30)))

pheatmap(filtered_gene_expression_ordered,
         annotation_col = cluster_annotation,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "row",
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = viridis(50),
         breaks = unique(c(seq(0, 2, length.out = 20), seq(2, max(filtered_gene_expression_ordered), length.out = 30))) # Uso de unique() para eliminar duplicados
)

png("/home/workdir/HDAC_tfm/heatmap.png")
DoHeatmap(object = tumor_cells_copy,
          features = head(VariableFeatures(tumor_cells_copy), 50),  # Los primeros 20 genes variables
          group.by = "Adjusted_Clusters",    # Agrupar por identidad de células
          group.bar = TRUE,      # Incluir la barra de grupos
          disp.min = -2.5,       # Valor mínimo para la escala de colores
          disp.max = 2.5,        # Valor máximo para la escala de colores
          slot = "scale.data",   # Usar los datos escalados
          label = TRUE,          # Mostrar los nombres de los genes
          size = 5.5,            # Tamaño de las etiquetas
          angle = 45,            # Rotar las etiquetas a 45 grados
          raster = TRUE,         # Usar imágenes rasterizadas
          draw.lines = TRUE,     # Dibujar líneas separadoras entre filas
          group.bar.height = 0.02 # Alto de la barra de grupo
)
dev.off()

tumor_cells_copy$ident <- tumor_cells_copy$Adjusted_Clusters
markers <- FindAllMarkers(
  tumor_cells_copy,                # Objeto Seurat
  only.pos = TRUE,                  # Solo genes con expresión positiva
  min.pct = 0.25,                   # Mínimo porcentaje de células que deben expresar el gen
  logfc.threshold = 0.25            # Umbral de cambio en la expresión (logFC)
)

top_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))  # Ordenar los genes dentro de cada cluster por avg_log2FC

p <- DoHeatmap(tumor_cells_copy, features = top_markers$gene) + NoLegend()
ggsave("./HDAC_tfm/fig/sub_clust_tumor/subclustering_v2/Heatmap_subclustering_tumor_unnan_1503.png", plot = p, width = 10, height = 8, dpi=600)
