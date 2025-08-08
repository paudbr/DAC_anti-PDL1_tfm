################################################################################
# Title: DGEA and PEA between control vs responder , control vs non responder from both cell lines
# Author: Paula
################################################################################


library(Seurat)
library(ggplot2)
library(tidyverse)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ReactomePA)

#Load seurat object with only tumor cells already clustered:

tumor_cells_clustered <- readRDS("/home/workdir/HDAC_tfm/data/tumor_cells_adjusted_13032025.rds")

tumor_cells_copy <- JoinLayers(tumor_cells_clustered)
cells_responsive_intersect_tfs_uv <- readRDS( "/home/workdir/HDAC_tfm/data/cells_responsive_intersect_tfs_kpb25l_uv.rds")
cells_unresponsive_intersect_tfs_uv <- readRDS("/home/workdir/HDAC_tfm/data/cells_unresponsive_intersect_tfs_kpb25l_uv.rds")

cells_responsive_intersect_tfs <- readRDS("/home/workdir/HDAC_tfm/data/cells_responsive_intersect_tfs_kpb25l.rds")
cells_unresponsive_intersect_tfs<- readRDS("/home/workdir/HDAC_tfm/data/cells_unresponsive_intersect_tfs_kpb25l.rds")


tumor_responsive_unresponsive_kpb25l <- subset(tumor_cells_copy, cells = c(cells_responsive_intersect_tfs,cells_unresponsive_intersect_tfs ))


tumor_responsive_unresponsive_kpb25l_uv <- subset(tumor_cells_copy, cells = c(cells_responsive_intersect_tfs_uv,cells_unresponsive_intersect_tfs_uv ))

tumor_cells_control_kpb25l <- colnames(subset(tumor_cells_copy, subset = sample == "KAP25L_control"))


tumor_cells_control_kuv <- colnames(subset(tumor_cells_copy, subset = sample == "KAP25L-UV_control"))


markers_responsive_vs_control_kpb25l <- FindMarkers(tumor_cells_copy, 
                                    ident.1 = cells_responsive_intersect_tfs, #unresponsive
                                    ident.2 = tumor_cells_control_kpb25l)

markers_unresponsive_vs_control_kpb25l <- FindMarkers(tumor_cells_copy, 
                                                    ident.1 = cells_unresponsive_intersect_tfs, #unresponsive
                                                    ident.2 = tumor_cells_control_kpb25l)

markers_responsive_vs_control_kuv <- FindMarkers(tumor_cells_copy, 
                                                    ident.1 = cells_responsive_intersect_tfs_uv, #unresponsive
                                                    ident.2 = tumor_cells_control_kuv)

markers_unresponsive_vs_control_kuv <- FindMarkers(tumor_cells_copy, 
                                                      ident.1 = cells_unresponsive_intersect_tfs_uv, #unresponsive
                                                      ident.2 = tumor_cells_control_kuv)

filtrar_y_guardar_markers <- function(markers_df, nombre_archivo, 
                                      pval_cutoff = 0.05, logfc_cutoff = 1) {
  # Filtrar por p-valor ajustado y logFC si se indica
  markers_filtrados <- subset(markers_df, p_val_adj < pval_cutoff & abs(avg_log2FC) > logfc_cutoff)
  
  # Guardar en CSV
  write.csv(markers_filtrados, file = nombre_archivo, row.names = TRUE)
  
  # Mensaje de confirmación
  cat("Guardado:", nombre_archivo, "\nGenes significativos:", nrow(markers_filtrados), "\n\n")
}

filtrar_y_guardar_markers(markers_responsive_vs_control_kpb25l, "/home/workdir/HDAC_tfm/data/markers_responsive_vs_control_kpb25l.csv")
filtrar_y_guardar_markers(markers_unresponsive_vs_control_kpb25l, "/home/workdir/HDAC_tfm/data/markers_unresponsive_vs_control_kpb25l.csv")
filtrar_y_guardar_markers(markers_responsive_vs_control_kuv, "/home/workdir/HDAC_tfm/data/markers_responsive_vs_control_kuv.csv")
filtrar_y_guardar_markers(markers_unresponsive_vs_control_kuv, "/home/workdir/HDAC_tfm/data/markers_unresponsive_vs_control_kuv.csv")


leer_y_separar_markers <- function(archivo) {
  # Leer archivo CSV
  markers_df <- read.csv(archivo, row.names = 1)
  
  # Filtrar significativos
  markers_filtrados <- markers_df %>%
    dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 1)
  
  # Separar up y down
  upregulated <- markers_filtrados %>% dplyr::filter(avg_log2FC > 1)
  downregulated <- markers_filtrados %>% dplyr::filter(avg_log2FC < -1)
  
  # Devolver lista con los 3 dataframes
  return(list(
    all_filtered = markers_filtrados,
    upregulated = upregulated,
    downregulated = downregulated
  ))
}

archivos <- list(
  responsive_kpb25l = "/home/workdir/HDAC_tfm/data/markers_responsive_vs_control_kpb25l.csv",
  unresponsive_kpb25l = "/home/workdir/HDAC_tfm/data/markers_unresponsive_vs_control_kpb25l.csv",
  responsive_kuv = "/home/workdir/HDAC_tfm/data/markers_responsive_vs_control_kuv.csv",
  unresponsive_kuv = "/home/workdir/HDAC_tfm/data/markers_unresponsive_vs_control_kuv.csv"
)

# Leer y separar todos
resultados_markers <- lapply(archivos, leer_y_separar_markers)



convert_to_entrez <- function(gene_symbols) {
  bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
}

gene_list <- list(
  KPB25L_Up = convert_to_entrez(rownames(resultados_markers$responsive_kpb25l$upregulated)),
  KPB25L_Down = convert_to_entrez(rownames(resultados_markers$responsive_kpb25l$downregulated)),
  KPB25L_UV_Up = convert_to_entrez(rownames(resultados_markers$responsive_kuv$upregulated)),
  KPB25L_UV_Down = convert_to_entrez(rownames(resultados_markers$responsive_kuv$downregulated))
)


entrez_list <- lapply(gene_list, function(x) x$ENTREZID)

go_cluster <- compareCluster(
  geneCluster = entrez_list,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

# KEGG pathways
kegg_cluster <- compareCluster(
  geneCluster = entrez_list,
  fun = "enrichKEGG",
  organism = "mmu",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
reactome_cluster <- compareCluster(
  geneCluster = entrez_list,
  fun = "enrichPathway",
  organism = "mouse",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)


top_terms <- go_cluster %>%
  as.data.frame() %>%
  arrange(desc(GeneRatio), p.adjust) %>%
  slice_head(n = 20)

top_terms$GeneRatio <- sapply(strsplit(as.character(top_terms$GeneRatio), "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]))

p_go<- dotplot(go_cluster, showCategory = 20, split = "cluster") +
  ggtitle("GO BP PEA by cell line in Responder vs control cancer cells") +
  theme_classic() +
  theme(
    text = element_text(family = "Times New Roman", size = 24, color = "black"),
    axis.text.x = element_text(family = "Times New Roman", size = 24, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(family = "Times New Roman", size = 24, color = "black"),
    legend.text = element_text(family = "Times New Roman", size = 24),
    plot.title = element_text(size = 26, face = "bold", hjust=1)
  )


ggsave(
  filename = "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/Dotplot_go_responder_control_cell_line_compareclusters.png",
  plot = p_go,
  height = 45, width = 15, dpi = 700
)


p_kegg <- dotplot(kegg_cluster, showCategory = 20, split = "cluster") +
  ggtitle("KEGG PEA by cell line in Responder vs control cancer cells") +
  theme_classic() +
  theme(
    text = element_text(family = "Times New Roman", size = 24, color = "black"),
    axis.text.x = element_text(family = "Times New Roman", size = 24, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(family = "Times New Roman", size = 24, color = "black"),
    legend.text = element_text(family = "Times New Roman", size = 24),
    plot.title = element_text(size = 26, face = "bold", hjust=1)
  )
ggsave(
  filename = "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/Dotplot_kegg_responder_control_cell_line_compareclusters.png",
  plot = p_kegg,
  height = 30, width = 13, dpi = 700
)




p_react <- dotplot(reactome_cluster, showCategory = 20, split = "cluster") +
  ggtitle("Reactome PEA by cell line in Responder vs control cancer cells") +
  theme_classic() +
  theme(
    text = element_text(family = "Times New Roman", size = 26, color = "black"),
    axis.text.x = element_text(family = "Times New Roman", size = 26, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(family = "Times New Roman", size = 26, color = "black", hjust = 0),
    legend.text = element_text(family = "Times New Roman", size = 26),
    plot.title = element_text(size = 26, face = "bold", hjust = 1)
  )

ggsave(
  filename = "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/Dotplot_react_responder_control_cell_line_compareclusters.png",
  plot = p_react,
  height = 43, width = 13, dpi = 700)


# NO-RESPONDER ------------------------------------------------------------


gene_list <- list(
  KPB25L_Up = convert_to_entrez(rownames(resultados_markers$unresponsive_kpb25l$upregulated)),
  KPB25L_Down = convert_to_entrez(rownames(resultados_markers$unresponsive_kpb25l$downregulated)),
  KPB25L_UV_Up = convert_to_entrez(rownames(resultados_markers$unresponsive_kuv$upregulated)),
  KPB25L_UV_Down = convert_to_entrez(rownames(resultados_markers$unresponsive_kuv$downregulated))
)


entrez_list <- lapply(gene_list, function(x) x$ENTREZID)

go_cluster <- compareCluster(
  geneCluster = entrez_list,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

# KEGG pathways
kegg_cluster <- compareCluster(
  geneCluster = entrez_list,
  fun = "enrichKEGG",
  organism = "mmu",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
reactome_cluster <- compareCluster(
  geneCluster = entrez_list,
  fun = "enrichPathway",
  organism = "mouse",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)


top_terms <- go_cluster %>%
  as.data.frame() %>%
  arrange(desc(GeneRatio), p.adjust) %>%
  slice_head(n = 20)

top_terms$GeneRatio <- sapply(strsplit(as.character(top_terms$GeneRatio), "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]))

p_go<- dotplot(go_cluster, showCategory = 20, split = "cluster") +
  ggtitle("GO BP PEA by cell line in Non-responder vs control cancer cells") +
  theme_classic() +
  theme(
    text = element_text(family = "Times New Roman", size = 24, color = "black"),
    axis.text.x = element_text(family = "Times New Roman", size = 24, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(family = "Times New Roman", size = 24, color = "black"),
    legend.text = element_text(family = "Times New Roman", size = 24),
    plot.title = element_text(size = 26, face = "bold", hjust=1)
  )


ggsave(
  filename = "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/Dotplot_go_unresponder_control_cell_line_compareclusters.png",
  plot = p_go,
  height = 45, width = 15, dpi = 700
)


p_kegg <- dotplot(kegg_cluster, showCategory = 20, split = "cluster") +
  ggtitle("KEGG PEA by cell line in Non-responder vs control cancer cells") +
  theme_classic() +
  theme(
    text = element_text(family = "Times New Roman", size = 24, color = "black"),
    axis.text.x = element_text(family = "Times New Roman", size = 24, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(family = "Times New Roman", size = 24, color = "black"),
    legend.text = element_text(family = "Times New Roman", size = 24),
    plot.title = element_text(size = 26, face = "bold", hjust=1)
  )
ggsave(
  filename = "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/Dotplot_kegg_unresponder_control_cell_line_compareclusters.png",
  plot = p_kegg,
  height = 30, width = 13, dpi = 700
)




p_react <- dotplot(reactome_cluster, showCategory = 20, split = "cluster") +
  ggtitle("Reactome PEA by cell line in Non-responder vs control cancer cells") +
  theme_classic() +
  theme(
    text = element_text(family = "Times New Roman", size = 26, color = "black"),
    axis.text.x = element_text(family = "Times New Roman", size = 26, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(family = "Times New Roman", size = 26, color = "black", hjust = 0),
    legend.text = element_text(family = "Times New Roman", size = 26),
    plot.title = element_text(size = 26, face = "bold", hjust = 1)
  )

ggsave(
  filename = "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/Dotplot_react_unresponder_control_cell_line_compareclusters.png",
  plot = p_react,
  height = 43, width = 13, dpi = 700)


# P27 MAKER QUIESCENCE ----------------------------------------------------

tumor_responsive_unresponsive_kpb25l <- AddModuleScore(tumor_responsive_unresponsive_kpb25l, features = "Cdkn1b")

ViolinPlot(tumor_responsive_unresponsive_kpb25l@)

tumor_responsive_unresponsive_kpb25l$responsive_unresponsive <- NA  # inicializar columna

# Asignar "responsive"
tumor_responsive_unresponsive_kpb25l$responsive_unresponsive[
  colnames(tumor_responsive_unresponsive_kpb25l) %in% cells_responsive_intersect_tfs
] <- "responder"

# Asignar "unresponsive"
tumor_responsive_unresponsive_kpb25l$responsive_unresponsive[
  colnames(tumor_responsive_unresponsive_kpb25l) %in% cells_unresponsive_intersect_tfs
] <- "non-responder"

VlnPlot(
  tumor_responsive_unresponsive_kpb25l,
  features = "Cdkn1b",
  group.by = "responsive_unresponsive",
)

library(ggpubr)

df <- FetchData(tumor_responsive_unresponsive_kpb25l, vars = c("Cdkn1b", "Cdkn1a", "Cdkn1c", "Cdk4", "Cdk6", "Cdk2","Ccnd1", "Ccnd2", "Ccnd3","Ccne1", "responsive_unresponsive"))
df_long <- reshape2::melt(df, id.vars = "responsive_unresponsive", variable.name = "gene", value.name = "expression")

# Hacer plots separados para cada gen
library(dplyr)
genes <- unique(df_long$gene)

plots <- lapply(genes, function(g) {
  ggviolin(
    df_long %>% filter(gene == g), 
    x = "responsive_unresponsive", 
    y = "expression",
    add = "boxplot",
    fill = "responsive_unresponsive",
    palette = c("responder" = "yellow", "non-responder" = "blue")
  ) +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    ggtitle(g)
})

plots[[10]]

combined_plot <- wrap_plots(plots) + 
  plot_layout(ncol = 2) + 
  plot_annotation(
    title = "Violin Plots of expression of key regulation actors of quiescence entrance in responder vs non-responder tumor cells in KPB25L",
    theme = theme(
      plot.title = element_text(family = "Times New Roman", size = 14, face = "bold", hjust = 0.5)
    )
  )
print(combined_plot)

ggsave( 
  filename = "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/VP_responder_nonresponder_quiescen_k25l.png",
  plot = combined_plot,
  height = 20, width = 15, dpi = 700
)




### KUV 


tumor_responsive_unresponsive_kpb25l_uv$responsive_unresponsive <- NA  # inicializar columna

# Asignar "responsive"
tumor_responsive_unresponsive_kpb25l_uv$responsive_unresponsive[
  colnames(tumor_responsive_unresponsive_kpb25l_uv) %in% cells_responsive_intersect_tfs_uv
] <- "responder"

# Asignar "unresponsive"
tumor_responsive_unresponsive_kpb25l_uv$responsive_unresponsive[
  colnames(tumor_responsive_unresponsive_kpb25l_uv) %in% cells_unresponsive_intersect_tfs_uv
] <- "non-responder"

VlnPlot(
  tumor_responsive_unresponsive_kpb25l_uv,
  features = "Cdkn1b",
  group.by = "responsive_unresponsive",
)

library(ggpubr)

df <- FetchData(tumor_responsive_unresponsive_kpb25l_uv, vars = c("Cdkn1b", "Cdkn1a", "Cdkn1c", "Cdk4", "Cdk6", "Cdk2", "Ccnd1", "Ccnd2", "Ccnd3","Ccne1", "responsive_unresponsive"))
df_long <- reshape2::melt(df, id.vars = "responsive_unresponsive", variable.name = "gene", value.name = "expression")
df_long$responsive_unresponsive <- factor(df_long$responsive_unresponsive, levels = c("responder", "non-responder"))

# Hacer plots separados para cada gen
library(dplyr)
genes <- unique(df_long$gene)

plots <- lapply(genes, function(g) {
  ggviolin(
    df_long %>% filter(gene == g), 
    x = "responsive_unresponsive", 
    y = "expression",
    add = "boxplot",
    fill = "responsive_unresponsive",
    palette = c("responder" = "yellow", "non-responder" = "blue")
  ) +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    ggtitle(g)
})

# Para mostrar el primero
print(plots[[10]])


combined_plot <- wrap_plots(plots) + plot_layout(ncol = 2)  # pones 2 columnas, por ejemplo

print(combined_plot)

combined_plot <- wrap_plots(plots) + 
  plot_layout(ncol = 4) + 
  plot_annotation(
    title = "Violin Plots of expression of key regulation actors of quiescence entrance in responder vs non-responder tumor cells in KPB25L-UV",
    theme = theme(
      plot.title = element_text(family = "Times New Roman", size = 14, face = "bold", hjust = 0.5)
    )
  )
print(combined_plot)

ggsave( 
  filename = "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/VP_responder_nonresponder_quiescen_kuv.png",
        plot = combined_plot,
        height = 20, width = 15, dpi = 700
)
