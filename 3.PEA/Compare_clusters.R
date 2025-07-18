library(dplyr)
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)


deg_files_cancer <- list(
  markers_ctrl_pdl1_k25l  = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kpb25l_ctrl_pdl1.csv",
  markers_ctrl_pdl1_kuv   = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kuv_ctrl_pdl1.csv",
  markers_ctrl_dac_k25l   = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kpb25l_ctrl_dac.csv",
  markers_ctrl_dac_kuv    = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kuv_ctrl_dac.csv",
  markers_ctrl_comb_k25l  = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kpb25l_ctrl_comb.csv",
  markers_ctrl_comb_kuv   = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kuv_ctrl_comb.csv",
  markers_dac_comb_k25l   = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kpb25l_dac_comb.csv",
  markers_dac_comb_kuv    = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kuv_dac_comb.csv"
)

deg_files_unann <- list(
  markers_ctrl_pdl1_k25l  = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kpb25l_ctrl_pdl1_unann.csv",
  markers_ctrl_pdl1_kuv   = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kuv_ctrl_pdl1_unann.csv",
  markers_ctrl_dac_k25l   = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kpb25l_ctrl_dac_unann.csv",
  markers_ctrl_dac_kuv    = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kuv_ctrl_dac_unann.csv",
  markers_ctrl_comb_k25l  = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kpb25l_ctrl_comb_unann.csv",
  markers_ctrl_comb_kuv   = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kuv_ctrl_comb_unann.csv",
  markers_dac_comb_k25l   = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kpb25l_dac_comb_unann.csv",
  markers_dac_comb_kuv    = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kuv_dac_comb_unann.csv"
)

get_positive_genes_from_files <- function(files_list, pattern_line_cell) {
  # Definir los tratamientos que quieres buscar
  tratamientos <- c("ctrl_pdl1", "ctrl_dac", "ctrl_comb", "dac_comb")
  
  # Crear una lista vacía para guardar los genes por tratamiento
  genes_por_tratamiento <- list()
  
  for (tto in tratamientos) {
    genes_all <- c()
    for (f in files_list) {
      if (grepl(pattern_line_cell, f) && grepl(tto, f)) {
        df <- read.csv(f, row.names = 1)
        if ("avg_log2FC" %in% colnames(df)) {
          df_pos <- df[df$avg_log2FC > 1, ]
          genes_all <- c(genes_all, rownames(df_pos))
        }
      }
    }
    genes_por_tratamiento[[tto]] <- unique(genes_all)
  }
  
  return(genes_por_tratamiento)
}


get_negative_genes_from_files <- function(files_list, pattern_line_cell) {
  # Definir los tratamientos que quieres buscar
  tratamientos <- c("ctrl_pdl1", "ctrl_dac", "ctrl_comb", "dac_comb")
  
  # Crear una lista vacía para guardar los genes por tratamiento
  genes_por_tratamiento <- list()
  
  for (tto in tratamientos) {
    genes_all <- c()
    for (f in files_list) {
      if (grepl(pattern_line_cell, f) && grepl(tto, f)) {
        df <- read.csv(f, row.names = 1)
        if ("avg_log2FC" %in% colnames(df)) {
          df_pos <- df[df$avg_log2FC < -1, ]
          genes_all <- c(genes_all, rownames(df_pos))
        }
      }
    }
    genes_por_tratamiento[[tto]] <- unique(genes_all)
  }
  
  return(genes_por_tratamiento)
}

cancer_k25l_genes_pos <- get_positive_genes_from_files(deg_files_cancer, "kpb25l")
unann_k25l_genes_pos  <- get_positive_genes_from_files (deg_files_unann, "kpb25l")

cancer_k25l_genes_neg <- get_negative_genes_from_files(deg_files_cancer, "kpb25l")
unann_k25l_genes_neg  <- get_negative_genes_from_files (deg_files_unann, "kpb25l")

cancer_kuv_genes_neg <- get_negative_genes_from_files(deg_files_cancer, "kuv")
unann_kuv_genes_neg  <- get_negative_genes_from_files(deg_files_unann, "kuv")

cancer_kuv_genes_pos <- get_positive_genes_from_files(deg_files_cancer, "kuv")
unann_kuv_genes_pos  <- get_positive_genes_from_files(deg_files_unann, "kuv")



k25lgene_list_up_tumor_vs_unann <- list(

  PDL1_cancer = cancer_k25l_genes_pos$ctrl_pdl1,
  PDL1_Unann = unann_k25l_genes_pos$ctrl_pdl1,
  DAC_cancer = cancer_k25l_genes_pos$ctrl_dac,
  DAC_Unann = unann_k25l_genes_pos$ctrl_dac,
  Combo_cancer = cancer_k25l_genes_pos$ctrl_comb,
  Combo_Unann = unann_k25l_genes_pos$ctrl_comb
)

kuvgene_list_up_tumor_vs_unann <- list(

  PDL1_cancer = cancer_kuv_genes_pos$ctrl_pdl1,
  PDL1_Unann= unann_kuv_genes_pos$ctrl_pdl1,
  DAC_cancer = cancer_kuv_genes_pos$ctrl_dac,
  DAC_Unann = unann_kuv_genes_pos$ctrl_dac,
  Combo_cancer = cancer_kuv_genes_pos$ctrl_comb,
  Combo_Unann = unann_kuv_genes_pos$ctrl_comb
)


k25lgene_list_down_tumor_vs_unann <- list(

  PDL1_cancer = cancer_k25l_genes_neg$ctrl_pdl1,
  PDL1_Unann = unann_k25l_genes_neg$ctrl_pdl1,
  DAC_cancer = cancer_k25l_genes_neg$ctrl_dac,
  DAC_Unann = unann_k25l_genes_neg$ctrl_dac,
  Combo_cancer = cancer_k25l_genes_neg$ctrl_comb,
  Combo_Unann= unann_k25l_genes_neg$ctrl_comb
)

kuvgene_list_down_tumor_vs_unann <- list(

  PDL1_cancer = cancer_kuv_genes_neg$ctrl_pdl1,
  PDL1_Unann = unann_kuv_genes_neg$ctrl_pdl1,
  DAC_cancer = cancer_kuv_genes_neg$ctrl_dac,
  DAC_Unann = unann_kuv_genes_neg$ctrl_dac,
  Combo_cancer = cancer_kuv_genes_neg$ctrl_comb,
  Combo_Unann= unann_kuv_genes_neg$ctrl_comb
)
  
tratamiento_colors <- c(
  "DAC_cancer"     = "#FF6F61",  # verde oscuro
  "DAC_Unann"      = "#FFD1CC",  # verde claro
  "PDL1_cancer"    = "#4682B4",  # morado oscuro
  "PDL1_Unann"     = "#B0C4DE",  # morado claro
  "Combo_cancer"   = "#7570B3",  # naranja oscuro
  "Combo_Unann"    = "#C2C2E2"   # naranja claro
)

generate_barplot <- function(gene_list, title, filename, pathway = c("KEGG", "Reactome", "GO")) {
  pathway <- match.arg(pathway)
  
  if (pathway == "KEGG") {
    enrichment_result <- compareCluster(
      geneCluster = gene_list,
      fun = "enrichKEGG",
      organism = "mmu",  # código para ratón
      keyType = "kegg",  # importante: puede que necesites conversión
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH"
    )
  } else if (pathway == "Reactome") {
    enrichment_result <- compareCluster(
      geneCluster = gene_list,
      fun = "enrichPathway",
      organism = "mouse",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH"
    )
  } else if (pathway == "GO") {
    enrichment_result <- compareCluster(
      geneCluster = gene_list,
      fun = "enrichGO",
      OrgDb = org.Mm.eg.db,
      keyType = "SYMBOL",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05
    )
    
  }
  
  # Extraer tabla de resultados
  df <- as.data.frame(enrichment_result)
  
  if (nrow(df) == 0) {
    warning("No enriched terms found.")
    return(NULL)
  }
  
  # Asegurar que RichFactor esté disponible
  if (!"RichFactor" %in% colnames(df)) {
    df$RichFactor <- df$Count / as.numeric(sub("/.*", "", df$BgRatio))
  }
  
  # Seleccionar top 10 términos por grupo
  top_terms <- df %>%
    group_by(Cluster) %>%
    slice_max(order_by = RichFactor, n = 10) %>%
    ungroup()
  
  # Para que el eje X esté ordenado por RichFactor
  top_terms$Description <- factor(top_terms$Description, levels = rev(unique(top_terms$Description)))
  
  # Crear gráfico
  p <- ggplot(top_terms, aes(x = Description, y = RichFactor, fill = Cluster)) +
    geom_bar(stat = "identity", position = position_stack()) +
    coord_flip() +
    theme_minimal(base_family = "Times New Roman") +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 1),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      panel.background = element_rect(fill = "white", colour = "white"),
      plot.background = element_rect(fill = "white", colour = "white"),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    ) +
    labs(title = title, x = "Pathway", y = "RichFactor") +
    scale_fill_manual(values = tratamiento_colors)
  
  print(p)
  

  ggsave(filename, plot = p, width = 12, height = 8, dpi = 600)
}





# Generar y guardar dotplot para los genes UP y DOWN comparando Cancer cells vs Unannotated
generate_barplot(k25lgene_list_up_tumor_vs_unann, "GO Enrichment for UP DEGs in KPB25L ( Cancer cells vs Unannotated)", "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/GO_k25l_up_tumor_vs_unann.png", "GO")
generate_barplot(k25lgene_list_down_tumor_vs_unann, "GO Enrichment for DOWN DEGs in KPB25L ( Cancer cells vs Unannotated)", "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/GO_k25l_down_tumor_vs_unann.png", "GO")
generate_barplot(kuvgene_list_up_tumor_vs_unann, "GO Enrichment for UP DEGs in KPB25L-UV ( Cancer cells vs Unannotated)", "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/GO_kuv_up_tumor_vs_unann.png", "GO")
generate_barplot(kuvgene_list_down_tumor_vs_unann, "GO Enrichment for DOWN DEGs in KPB25L-UV ( Cancer cells vs Unannotated)", "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/GO_kuv_down_tumor_vs_unann.png", "GO")



k25lgene_list_up_tumor_vs_unann_entrez <- lapply(k25lgene_list_up_tumor_vs_unann, function(g) {
  bitr(g, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID
})
k25lgene_list_down_tumor_vs_unann_entrez <- lapply(k25lgene_list_down_tumor_vs_unann, function(g) {
  bitr(g, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID
})

kuvgene_list_up_tumor_vs_unann_entrez <- lapply(kuvgene_list_up_tumor_vs_unann, function(g) {
  bitr(g, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID
})
kuvgene_list_down_tumor_vs_unann_entrez <- lapply(kuvgene_list_down_tumor_vs_unann, function(g) {
  bitr(g, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID
})

generate_barplot(k25lgene_list_up_tumor_vs_unann_entrez, "KEGG Enrichment for UP DEGs in KPB25L ( Cancer cells vs Unannotated)", "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/Keg_k25l_up_tumor_vs_unann.png", "KEGG")
generate_barplot(k25lgene_list_down_tumor_vs_unann_entrez, "KEGG Enrichment for DOWN DEGs in in KPB25L ( Cancer cells vs Unannotated)", "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/Keg_k25l_down_tumor_vs_unann.png", "KEGG")
generate_barplot(kuvgene_list_up_tumor_vs_unann_entrez, "KEGG Enrichment for UP DEGs in KPB25L-UV ( Cancer cells vs Unannotated)", "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/Keg_kuv_up_tumor_vs_unann.png", "KEGG")
generate_barplot(kuvgene_list_down_tumor_vs_unann_entrez, "KEGG Enrichment for DOWN DEGs in in KPB25L-UV ( Cancer cells vs Unannotated)", "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/Keg_kuv_down_tumor_vs_unann.png", "KEGG")


generate_barplot(k25lgene_list_up_tumor_vs_unann_entrez, "Reactome Enrichment for UP DEGs in KPB25L ( Cancer cells vs Unannotated)", "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/react_k25l_up_tumor_vs_unann.png", "Reactome")
generate_barplot(k25lgene_list_down_tumor_vs_unann_entrez, "Reactome Enrichment for DOWN DEGs in in KPB25L ( Cancer cells vs Unannotated)", "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/react_k25l_down_tumor_vs_unann.png", "Reactome")
generate_barplot(kuvgene_list_up_tumor_vs_unann_entrez, "Reactome Enrichment for UP DEGs in KPB25L-UV ( Cancer cells vs Unannotated)", "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/react_kuv_up_tumor_vs_unann.png", "Reactome")
generate_barplot(kuvgene_list_down_tumor_vs_unann_entrez, "Reactome Enrichment for DOWN DEGs in in KPB25L-UV ( Cancer cells vs Unannotated)", "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/react_kuv_down_tumor_vs_unann.png", "Reactome")

