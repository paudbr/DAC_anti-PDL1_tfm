library(Seurat)
library(ggplot2)
immune_cells <- readRDS("/home/workdir/analysis_augie/20250120_labeled_immune.combined_edit.RDS")
immune_cells_joined <- JoinLayers(immune_cells)

p <- DimPlot(immune_cells_joined, group.by = "immuno", label = FALSE) +
  ggtitle("UMAP of annotated clustering of single cell data") +  # Cambia el título
  theme(
    plot.title = element_text(family = "Times New Roman", face = "bold", size = 16),
    axis.title = element_text(family = "Times New Roman", size = 14),
    axis.text = element_text(family = "Times New Roman", size = 12),
    legend.text = element_text(family = "Times New Roman", size = 12),
    legend.title = element_text(family = "Times New Roman", size = 13)
  )

install.packages("scico")
library(scico)

DimPlot(immune_cells_joined, group.by = "immuno", label = FALSE) +
  ggtitle("UMAP of annotated clustering of single cell data") +
  scico::scale_color_scico_d(palette = "hawaii") +  # Puedes probar: "batlow", "roma", "vik", "hawaii"
  theme(
    plot.title = element_text(family = "Times New Roman", face = "bold", size = 16),
    axis.title = element_text(family = "Times New Roman", size = 14),
    axis.text = element_text(family = "Times New Roman", size = 12),
    legend.text = element_text(family = "Times New Roman", size = 12),
    legend.title = element_text(family = "Times New Roman", size = 13)
  )

library(RColorBrewer)
niveles <- levels(factor(immune_cells_joined$immuno))
colores <- setNames(brewer.pal(n = length(niveles), name = "Set3"), niveles)


p1 <- DimPlot(immune_cells_joined, group.by = "immuno", label = FALSE) +
  ggtitle("UMAP projection of annotated cell clusters") +
  scale_color_manual(values = colores) +
  theme(
    plot.title = element_text(family = "Times New Roman", face = "bold", size = 13),
    axis.title = element_text(family = "Times New Roman", size = 14),
    axis.text = element_text(family = "Times New Roman", size = 12),
    legend.text = element_text(family = "Times New Roman", size = 12),
    legend.title = element_text(family = "Times New Roman", size = 10)
  )
ggsave("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/Initial_clustering_sc.png", plot=p1, height=4, width=7, dpi=600)