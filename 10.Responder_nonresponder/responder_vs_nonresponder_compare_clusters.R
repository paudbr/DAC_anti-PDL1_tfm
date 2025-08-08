################################################################################
# Title: DGEA and PEA between responder and non -responder cells
# Author: Paula
################################################################################


library(Seurat)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(tidyverse)
library(dplyr)
library(patchwork)

# KPB25L ------------------------------------------------------------------

markers_kpb25l <- readRDS("/home/workdir/HDAC_tfm/data/marker_genes_between_unresponsive_responisve_tumor_cells_SimiC_intersection_k25l.RSD")

# Filtrar genes significativamente regulados
markers_kpb25l_filtered <- markers_kpb25l %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 1)

# Separar up y down
up_kpb25l <- markers_kpb25l_filtered %>% filter(avg_log2FC > 1)
down_kpb25l <- markers_kpb25l_filtered %>% filter(avg_log2FC < -1)

# KPB25L-UV --------------------------------------------------------------

markers_kpb25l_uv <- readRDS("/home/workdir/HDAC_tfm/data/marker_genes_between_unresponsive_responisve_tumor_cells_SimiC_intersection_kuv.RSD")

markers_kpb25l_uv_filtered <- markers_kpb25l_uv %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 1)

up_kpb25l_uv <- markers_kpb25l_uv_filtered %>% filter(avg_log2FC > 1)
down_kpb25l_uv <- markers_kpb25l_uv_filtered %>% filter(avg_log2FC < -1)



convert_to_entrez <- function(gene_df) {
  bitr(gene_df$gen, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
}

gene_list <- list(
  KPB25L_Up = convert_to_entrez(up_kpb25l),
  KPB25L_Down = convert_to_entrez(down_kpb25l),
  KPB25L_UV_Up = convert_to_entrez(up_kpb25l_uv),
  KPB25L_UV_Down = convert_to_entrez(down_kpb25l_uv)
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
  ggtitle("GO BP PEA by cell line in Non-responder vs Responder cancer cells") +
  theme_classic() +
  theme(
    text = element_text(family = "Times New Roman", size = 24, color = "black"),
    axis.text.x = element_text(family = "Times New Roman", size = 24, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(family = "Times New Roman", size = 24, color = "black"),
    legend.text = element_text(family = "Times New Roman", size = 24),
    plot.title = element_text(size = 26, face = "bold", hjust=1)
  )

ggsave(
  filename = "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/Dotplot_go_responder_nonresponder_cell_line_compareclusters.png",
  plot = p_go,
  height = 45, width = 15, dpi = 700
)


p_kegg <- dotplot(kegg_cluster, showCategory = 20, split = "cluster") +
  ggtitle("KEGG PEA by cell line in Non-responder vs Responder cancer cells") +
  theme_classic() +
  theme(
    text = element_text(family = "Times New Roman", size = 24, color = "black"),
    axis.text.x = element_text(family = "Times New Roman", size = 24, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(family = "Times New Roman", size = 24, color = "black"),
    legend.text = element_text(family = "Times New Roman", size = 24),
    plot.title = element_text(size = 26, face = "bold", hjust=1)
  )
ggsave(
  filename = "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/Dotplot_kegg_responder_nonresponder_cell_line_compareclusters.png",
  plot = p_kegg,
  height = 30, width = 13, dpi = 700
)

reactome_cluster@compareClusterResult$Description <- gsub(
  "Insulin-like Growth Factor Binding Proteins",
  "",
  reactome_cluster@compareClusterResult$Description
)

# Luego eliminar la más corta
reactome_cluster@compareClusterResult$Description <- gsub(
  "Insulin-like Growth Factor",
  "",
  reactome_cluster@compareClusterResult$Description
)

reactome_cluster@compareClusterResult$Description <- gsub(
  "Insulin-like Growth Factor Binding Proteins",
  "",
  reactome_cluster@compareClusterResult$Description
)

# Luego eliminar la más corta
reactome_cluster@compareClusterResult$Description <- gsub(
  ")",
  "",
  reactome_cluster@compareClusterResult$Description
)


desc <- reactome_cluster@compareClusterResult$Description

desc <- gsub("signaling pathway", "sig. path.", desc, ignore.case = TRUE)
desc <- gsub("regulation of", "reg. of", desc, ignore.case = TRUE)
desc <- gsub("positive", "pos.", desc, ignore.case = TRUE)
desc <- gsub("negative", "neg.", desc, ignore.case = TRUE)
desc <- gsub("Insulin-like growth factor binding proteins", "", desc, ignore.case = TRUE)
desc <- gsub("Insulin-like growth factor", "IGF", desc, ignore.case = TRUE)
desc <- gsub("Holliday Junction Intermediates", "HJI", desc, ignore.case = TRUE)
desc <- gsub("Major pathway of", "", desc, ignore.case = TRUE)
desc <- gsub("signal", "sig", desc, ignore.case = TRUE)

# Limpieza extra: espacios duplicados y trimming
desc <- gsub("\\s{2,}", " ", desc)
desc <- trimws(desc)

reactome_cluster@compareClusterResult$Description <- desc


p_react <- dotplot(reactome_cluster, showCategory = 20, split = "cluster") +
  ggtitle("Reactome PEA by cell line in Non-responder vs Responder cancer cells") +
  theme_classic() +
  theme(
    text = element_text(family = "Times New Roman", size = 26, color = "black"),
    axis.text.x = element_text(family = "Times New Roman", size = 26, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(family = "Times New Roman", size = 26, color = "black", hjust = 0),
    legend.text = element_text(family = "Times New Roman", size = 26),
    plot.title = element_text(size = 26, face = "bold", hjust = 1)
  )

ggsave(
  filename = "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/Dotplot_react_responder_nonresponder_cell_line_compareclusters.png",
  plot = p_react,
  height = 43, width = 13, dpi = 700)



final_plot <- p_go + p_kegg + p_react + 
  plot_layout(nrow = 1, ncol = 3) + 
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      text = element_text(family = "Times New Roman"),
      plot.tag = element_text(family = "Times New Roman", size = 24, face = "bold")
    )
  )


ggsave(
  filename = "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/Dotplot_responder_nonresponder_cell_line_compareclusters.png",
  plot = final_plot,
  height = 40, width = 35, dpi = 700
)
