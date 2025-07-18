library(dplyr)
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ReactomePA)
library(ggplot2)

deg_files <- list(
  markers_ctrl_pdl1_k25l  = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kpb25l_ctrl_pdl1.csv",
  markers_ctrl_pdl1_kuv   = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kuv_ctrl_pdl1.csv",
  markers_ctrl_dac_k25l   = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kpb25l_ctrl_dac.csv",
  markers_ctrl_dac_kuv    = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kuv_ctrl_dac.csv",
  markers_ctrl_comb_k25l  = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kpb25l_ctrl_comb.csv",
  markers_ctrl_comb_kuv   = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kuv_ctrl_comb.csv",
  markers_dac_comb_k25l   = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kpb25l_dac_comb.csv",
  markers_dac_comb_kuv    = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kuv_dac_comb.csv"
)

# Genes upregulados (log2FC > 1)
deg_upregulated <- lapply(deg_files, function(file) {
  df <- read.csv(file, row.names = 1)
  rownames(df[df$avg_log2FC > 1, ])
})

# Genes downregulados (log2FC < -1)
deg_downregulated <- lapply(deg_files, function(file) {
  df <- read.csv(file, row.names = 1)
  rownames(df[df$avg_log2FC < -1, ])
})


convert_to_entrez <- function(gene_list) {
  bitr(gene_list,
       fromType = "SYMBOL",
       toType = "ENTREZID",
       OrgDb = org.Mm.eg.db) %>% 
    pull(ENTREZID) %>% 
    unique()
}

# Convertir todas las listas
entrez_up <- lapply(deg_upregulated, convert_to_entrez)
entrez_down <- lapply(deg_downregulated, convert_to_entrez)


# UPREGULATED -------------------------------------------------------------


go_results_up <- lapply(entrez_up, function(gene_ids) {
  enrichGO(gene = gene_ids,
           OrgDb = org.Mm.eg.db,
           ont = "BP",
           pAdjustMethod = "BH",
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.2,
           readable = TRUE)
})

kegg_results_up <- lapply(entrez_up, function(gene_ids) {
  enrichKEGG(gene = gene_ids,
             organism = "mmu",
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05)
})

react_results_up <- lapply(entrez_up, function(gene_ids){
  enrichPathway(gene = gene_ids, 
              organism = "mouse", 
              pAdjustMethod = "BH", 
              pvalueCutoff = 0.05)
})


# DOWNREGULATED -----------------------------------------------------------


go_results_down <- lapply(entrez_down, function(gene_ids) {
  enrichGO(gene = gene_ids,
           OrgDb = org.Mm.eg.db,
           ont = "BP",
           pAdjustMethod = "BH",
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.2,
           readable = TRUE)
})

# ORA con KEGG
kegg_results_down <- lapply(entrez_down, function(gene_ids) {
  enrichKEGG(gene = gene_ids,
             organism = "mmu",
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05)
})

react_results_down <- lapply(entrez_down, function(gene_ids){
  enrichPathway(gene = gene_ids, 
                organism = "mouse", 
                pAdjustMethod = "BH", 
                pvalueCutoff = 0.05)
})


# PLOT --------------------------------------------------------------------


plot_dotplots <- function(result_list, result_type, regulation) {
  names(result_list) %>%
    walk(function(name) {
      result <- result_list[[name]]
      if (!is.null(result) && nrow(result) > 0) {
        
        # Detección de tratamientos
        treatments <- c()
        if (grepl("ctrl", name, ignore.case = TRUE)) treatments <- c(treatments, "Control")
        if (grepl("pdl1", name, ignore.case = TRUE)) treatments <- c(treatments, "PD-L1")
        if (grepl("dac", name, ignore.case = TRUE)) treatments <- c(treatments, "DAC")
        if (grepl("comb", name, ignore.case = TRUE)) treatments <- c(treatments, "Combination")
        treatment_label <- paste(treatments, collapse = " vs ")
        
        # Detección de línea celular
        cell_line <- ifelse(grepl("k25l", name, ignore.case = TRUE), "KPB25L",
                            ifelse(grepl("kuv", name, ignore.case = TRUE), "KPB25L-UV", "Unknown"))
        
        # Título final
        title <- paste(cell_line, "-", treatment_label, "-", result_type, "-", regulation)
        filename <- paste0(gsub("-", "_", title), ".png")
        
        # Plot
        p <- dotplot(result, showCategory = 15) +
          ggtitle(title) +
          theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
       
        outdir <- "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_ORA/Cancer_cells/"
        ggsave(
          filename = file.path(outdir, filename),
          plot = p,
          width = 10, height = 8, dpi=600
        )
      }
    })
}


# Upregulated
plot_dotplots(go_results_up, "GO_BP", "UP")
plot_dotplots(kegg_results_up, "KEGG", "UP")
plot_dotplots(react_results_up, "REACTOME", "UP")

# Downregulated
plot_dotplots(go_results_down, "GO_BP", "DOWN")
plot_dotplots(kegg_results_down, "KEGG", "DOWN")
plot_dotplots(react_results_down, "REACTOME", "DOWN")





extract_enrichment_df <- function(result_list, result_type, regulation) {
  all_dfs <- list()
  
  for (name in names(result_list)) {
    res <- result_list[[name]]
    if (!is.null(res) && nrow(res) > 0) {
      df <- as.data.frame(res)
      
      # Extraer tratamientos
      treatments <- c()
      if (grepl("ctrl", name, ignore.case = TRUE)) treatments <- c(treatments, "Control")
      if (grepl("pdl1", name, ignore.case = TRUE)) treatments <- c(treatments, "PD-L1")
      if (grepl("dac", name, ignore.case = TRUE)) treatments <- c(treatments, "DAC")
      if (grepl("comb", name, ignore.case = TRUE)) treatments <- c(treatments, "Combination")
      treatment_label <- paste(treatments, collapse = " vs ")
      
      # Extraer línea celular
      cell_line <- ifelse(grepl("k25l", name, ignore.case = TRUE), "KPB25L",
                          ifelse(grepl("kuv", name, ignore.case = TRUE), "KPB25L-UV", "Unknown"))
      
      df$treatment <- treatment_label
      df$cell_line <- cell_line
      df$regulation <- regulation
      df$source <- result_type
      df$comparison <- name  # para trazabilidad
      
      all_dfs[[name]] <- df
    }
  }
  
  bind_rows(all_dfs)
}

all_enrichments_df <- bind_rows(
  extract_enrichment_df(go_results_up, "GO_BP", "UP"),
  extract_enrichment_df(go_results_down, "GO_BP", "DOWN"),
  extract_enrichment_df(kegg_results_up, "KEGG", "UP"),
  extract_enrichment_df(kegg_results_down, "KEGG", "DOWN"),
  extract_enrichment_df(react_results_up, "Reactome", "UP"),
  extract_enrichment_df(react_results_down, "Reactome", "DOWN")
)

all_enrichments_df$GeneRatio <- sapply(all_enrichments_df$GeneRatio, function(x) {
  if (is.numeric(x)) return(x)
  if (grepl("/", x)) eval(parse(text = x)) else as.numeric(x)
})


library(forcats)
plot_enrichment_summary <- function(df, cell_line, source_filter = "GO_BP", regulation_filter = "UP", top_terms = 15) {
  df_filtered <- df %>%
    filter(cell_line == !!cell_line,
           source == source_filter,
           regulation == regulation_filter) %>%
    group_by(treatment, Description) %>%
    slice_min(p.adjust, n = 1) %>%
    group_by(treatment) %>%
    slice_min(p.adjust, n = top_terms) %>%
    ungroup() %>%
    mutate(Description = fct_reorder(Description, GeneRatio))
  desired_order <- c("Control vs PD-L1", "Control vs DAC", "Control vs Combination", "DAC vs Combination")
  df_filtered$treatment <- factor(df_filtered$treatment, levels = desired_order)
  
  
  # Define título dinámico
  reg_label <- ifelse(regulation_filter == "UP", "Upregulated", "Downregulated")
  
  ggplot(df_filtered, aes(x = GeneRatio, y = Description, 
                          color = p.adjust, size = Count)) +
    geom_point() +
    scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
    scale_size_continuous(name = "Count") +
    facet_wrap(~treatment, scales = "free") +
    labs(
      title = paste(cell_line, "-", reg_label, " DEGs -", source_filter, "Pathway Enrichment in Cancer Cells"),
      x = "Gene Ratio", y = NULL
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(family = "Times New Roman", face = "bold", size = 11),
      axis.text.y = element_text(family = "Times New Roman", size = 13),
      plot.title = element_text(hjust = 0.5, family = "Times New Roman"),
      axis.title = element_text(),
      axis.text.x = element_text(size = 10),
      legend.title = element_text(),
      legend.text = element_text()
    )
}



plot_enrichment_summary(all_enrichments_df, "KPB25L", source_filter = "GO_BP")
plot_enrichment_summary(all_enrichments_df, "KPB25L-UV", source_filter = "GO_BP")
plot_enrichment_summary(all_enrichments_df, "KPB25L", source_filter = "KEGG")
plot_enrichment_summary(all_enrichments_df, "KPB25L-UV", source_filter = "KEGG")

plot_enrichment_summary(all_enrichments_df, "KPB25L", source_filter = "Reactome")
plot_enrichment_summary(all_enrichments_df, "KPB25L-UV", source_filter = "Reactome")



plot_enrichment_summary(all_enrichments_df, "KPB25L", regulation_filter = "DOWN", source_filter = "GO_BP")
plot_enrichment_summary(all_enrichments_df, "KPB25L-UV",  regulation_filter = "DOWN",source_filter = "GO_BP")
plot_enrichment_summary(all_enrichments_df, "KPB25L", regulation_filter = "DOWN", source_filter = "KEGG")
plot_enrichment_summary(all_enrichments_df, "KPB25L-UV", regulation_filter = "DOWN", source_filter = "KEGG")

plot_enrichment_summary(all_enrichments_df, "KPB25L", regulation_filter = "DOWN", source_filter = "Reactome")
plot_enrichment_summary(all_enrichments_df, "KPB25L-UV",  regulation_filter = "DOWN", source_filter = "Reactome")


cell_lines <- c("KPB25L", "KPB25L-UV")
sources <- c("GO_BP", "KEGG", "Reactome")
regulations <- c("UP", "DOWN")

for (cl in cell_lines) {
  for (src in sources) {
    for (reg in regulations) {
      outdir <- "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_ORA/Cancer_cells/"
      p <- plot_enrichment_summary(all_enrichments_df, cell_line = cl, source_filter = src, regulation_filter = reg)
      filename <- paste0(outdir,cl, "_", src, "_", reg, ".png")
      ggsave(filename, plot = p, width = 15, height = 9, dpi= 600)
    }
  }
}


# UNANNOTATED -------------------------------------------------------------


deg_files <- list(
  markers_ctrl_pdl1_k25l  = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kpb25l_ctrl_pdl1_unann.csv",
  markers_ctrl_pdl1_kuv   = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kuv_ctrl_pdl1_unann.csv",
  markers_ctrl_dac_k25l   = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kpb25l_ctrl_dac_unann.csv",
  markers_ctrl_dac_kuv    = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kuv_ctrl_dac_unann.csv",
  markers_ctrl_comb_k25l  = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kpb25l_ctrl_comb_unann.csv",
  markers_ctrl_comb_kuv   = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kuv_ctrl_comb_unann.csv",
  markers_dac_comb_k25l   = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kpb25l_dac_comb_unann.csv",
  markers_dac_comb_kuv    = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kuv_dac_comb_unann.csv"
)

# Genes upregulados (log2FC > 1)
deg_upregulated <- lapply(deg_files, function(file) {
  df <- read.csv(file, row.names = 1)
  rownames(df[df$avg_log2FC > 1, ])
})

# Genes downregulados (log2FC < -1)
deg_downregulated <- lapply(deg_files, function(file) {
  df <- read.csv(file, row.names = 1)
  rownames(df[df$avg_log2FC < -1, ])
})


convert_to_entrez <- function(gene_list) {
  bitr(gene_list,
       fromType = "SYMBOL",
       toType = "ENTREZID",
       OrgDb = org.Mm.eg.db) %>% 
    pull(ENTREZID) %>% 
    unique()
}

# Convertir todas las listas
entrez_up <- lapply(deg_upregulated, convert_to_entrez)
entrez_down <- lapply(deg_downregulated, convert_to_entrez)


# UPREGULATED -------------------------------------------------------------


go_results_up <- lapply(entrez_up, function(gene_ids) {
  enrichGO(gene = gene_ids,
           OrgDb = org.Mm.eg.db,
           ont = "BP",
           pAdjustMethod = "BH",
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.2,
           readable = TRUE)
})

kegg_results_up <- lapply(entrez_up, function(gene_ids) {
  enrichKEGG(gene = gene_ids,
             organism = "mmu",
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05)
})

react_results_up <- lapply(entrez_up, function(gene_ids){
  enrichPathway(gene = gene_ids, 
                organism = "mouse", 
                pAdjustMethod = "BH", 
                pvalueCutoff = 0.05)
})


# DOWNREGULATED -----------------------------------------------------------


go_results_down <- lapply(entrez_down, function(gene_ids) {
  enrichGO(gene = gene_ids,
           OrgDb = org.Mm.eg.db,
           ont = "BP",
           pAdjustMethod = "BH",
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.2,
           readable = TRUE)
})

# ORA con KEGG
kegg_results_down <- lapply(entrez_down, function(gene_ids) {
  enrichKEGG(gene = gene_ids,
             organism = "mmu",
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05)
})

react_results_down <- lapply(entrez_down, function(gene_ids){
  enrichPathway(gene = gene_ids, 
                organism = "mouse", 
                pAdjustMethod = "BH", 
                pvalueCutoff = 0.05)
})


# PLOT --------------------------------------------------------------------


plot_dotplots <- function(result_list, result_type, regulation) {
  names(result_list) %>%
    walk(function(name) {
      result <- result_list[[name]]
      if (!is.null(result) && nrow(result) > 0) {
        
        # Detección de tratamientos
        treatments <- c()
        if (grepl("ctrl", name, ignore.case = TRUE)) treatments <- c(treatments, "Control")
        if (grepl("pdl1", name, ignore.case = TRUE)) treatments <- c(treatments, "PD-L1")
        if (grepl("dac", name, ignore.case = TRUE)) treatments <- c(treatments, "DAC")
        if (grepl("comb", name, ignore.case = TRUE)) treatments <- c(treatments, "Combination")
        treatment_label <- paste(treatments, collapse = " vs ")
        
        # Detección de línea celular
        cell_line <- ifelse(grepl("k25l", name, ignore.case = TRUE), "KPB25L",
                            ifelse(grepl("kuv", name, ignore.case = TRUE), "KPB25L-UV", "Unknown"))
        
        # Título final
        title <- paste(cell_line, "-", treatment_label, "-", result_type, "-", regulation, "Unannotated")
        filename <- paste0(gsub("-", "_", title), "1707.png")
        
        # Plot
        p <- dotplot(result, showCategory = 15) +
          ggtitle(title) +
          theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
        
        outdir <- "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_ORA/Unann/"
        ggsave(
          filename = file.path(outdir, filename),
          plot = p,
          width = 15, height = 10, dpi=600
        )
      }
    })
}


# Upregulated
plot_dotplots(go_results_up, "GO_BP", "UP")
plot_dotplots(kegg_results_up, "KEGG", "UP")
plot_dotplots(react_results_up, "REACTOME", "UP")

# Downregulated
plot_dotplots(go_results_down, "GO_BP", "DOWN")
plot_dotplots(kegg_results_down, "KEGG", "DOWN")
plot_dotplots(react_results_down, "REACTOME", "DOWN")





extract_enrichment_df <- function(result_list, result_type, regulation) {
  all_dfs <- list()
  
  for (name in names(result_list)) {
    res <- result_list[[name]]
    if (!is.null(res) && nrow(res) > 0) {
      df <- as.data.frame(res)
      
      # Extraer tratamientos
      treatments <- c()
      if (grepl("ctrl", name, ignore.case = TRUE)) treatments <- c(treatments, "Control")
      if (grepl("pdl1", name, ignore.case = TRUE)) treatments <- c(treatments, "PD-L1")
      if (grepl("dac", name, ignore.case = TRUE)) treatments <- c(treatments, "DAC")
      if (grepl("comb", name, ignore.case = TRUE)) treatments <- c(treatments, "Combination")
      treatment_label <- paste(treatments, collapse = " vs ")
      
      # Extraer línea celular
      cell_line <- ifelse(grepl("k25l", name, ignore.case = TRUE), "KPB25L",
                          ifelse(grepl("kuv", name, ignore.case = TRUE), "KPB25L-UV", "Unknown"))
      
      df$treatment <- treatment_label
      df$cell_line <- cell_line
      df$regulation <- regulation
      df$source <- result_type
      df$comparison <- name  # para trazabilidad
      
      all_dfs[[name]] <- df
    }
  }
  
  bind_rows(all_dfs)
}

all_enrichments_df <- bind_rows(
  extract_enrichment_df(go_results_up, "GO_BP", "UP"),
  extract_enrichment_df(go_results_down, "GO_BP", "DOWN"),
  extract_enrichment_df(kegg_results_up, "KEGG", "UP"),
  extract_enrichment_df(kegg_results_down, "KEGG", "DOWN"),
  extract_enrichment_df(react_results_up, "Reactome", "UP"),
  extract_enrichment_df(react_results_down, "Reactome", "DOWN")
)

all_enrichments_df$GeneRatio <- sapply(all_enrichments_df$GeneRatio, function(x) {
  if (is.numeric(x)) return(x)
  if (grepl("/", x)) eval(parse(text = x)) else as.numeric(x)
})


library(forcats)
plot_enrichment_summary <- function(df, cell_line, source_filter = "GO_BP", regulation_filter = "UP", top_terms = 15) {
  df_filtered <- df %>%
    filter(cell_line == !!cell_line,
           source == source_filter,
           regulation == regulation_filter) %>%
    group_by(treatment, Description) %>%
    slice_min(p.adjust, n = 1) %>%
    group_by(treatment) %>%
    slice_min(p.adjust, n = top_terms) %>%
    ungroup() %>%
    mutate(Description = fct_reorder(Description, GeneRatio))
  desired_order <- c("Control vs PD-L1", "Control vs DAC", "Control vs Combination", "DAC vs Combination")
  df_filtered$treatment <- factor(df_filtered$treatment, levels = desired_order)
  
  
  # Define título dinámico
  reg_label <- ifelse(regulation_filter == "UP", "Upregulated", "Downregulated")
  
  ggplot(df_filtered, aes(x = GeneRatio, y = Description, 
                          color = p.adjust, size = Count)) +
    geom_point() +
    scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
    scale_size_continuous(name = "Count") +
    facet_wrap(~treatment, scales = "free") +
    labs(
      title = paste(cell_line, "-", reg_label, " DEGs -", source_filter, "Pathway Enrichment in Unannotated"),
      x = "Gene Ratio", y = NULL
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(family = "Times New Roman", face = "bold", size = 11),
      axis.text.y = element_text(family = "Times New Roman", size = 13),
      plot.title = element_text(hjust = 0.5, family = "Times New Roman"),
      axis.title = element_text(),
      axis.text.x = element_text(size = 10),
      legend.title = element_text(),
      legend.text = element_text()
    )
}



plot_enrichment_summary(all_enrichments_df, "KPB25L", source_filter = "GO_BP")
plot_enrichment_summary(all_enrichments_df, "KPB25L-UV", source_filter = "GO_BP")
plot_enrichment_summary(all_enrichments_df, "KPB25L", source_filter = "KEGG")
plot_enrichment_summary(all_enrichments_df, "KPB25L-UV", source_filter = "KEGG")

plot_enrichment_summary(all_enrichments_df, "KPB25L", source_filter = "Reactome")
plot_enrichment_summary(all_enrichments_df, "KPB25L-UV", source_filter = "Reactome")



plot_enrichment_summary(all_enrichments_df, "KPB25L", regulation_filter = "DOWN", source_filter = "GO_BP")
plot_enrichment_summary(all_enrichments_df, "KPB25L-UV",  regulation_filter = "DOWN",source_filter = "GO_BP")
plot_enrichment_summary(all_enrichments_df, "KPB25L", regulation_filter = "DOWN", source_filter = "KEGG")
plot_enrichment_summary(all_enrichments_df, "KPB25L-UV", regulation_filter = "DOWN", source_filter = "KEGG")

plot_enrichment_summary(all_enrichments_df, "KPB25L", regulation_filter = "DOWN", source_filter = "Reactome")
plot_enrichment_summary(all_enrichments_df, "KPB25L-UV",  regulation_filter = "DOWN", source_filter = "Reactome")


cell_lines <- c("KPB25L", "KPB25L-UV")
sources <- c("GO_BP", "KEGG", "Reactome")
regulations <- c("UP", "DOWN")

for (cl in cell_lines) {
  for (src in sources) {
    for (reg in regulations) {
      outdir <- "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_ORA/Unann/"
      p <- plot_enrichment_summary(all_enrichments_df, cell_line = cl, source_filter = src, regulation_filter = reg)
      filename <- paste0(outdir,cl, "_", src, "_", reg, "1707.png")
      ggsave(filename, plot = p, width = 20, height = 10, dpi= 600)
    }
  }
}


# MACROPHAGES -------------------------------------------------------------


deg_files <- list(
  markers_ctrl_pdl1_k25l  = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kpb25l_ctrl_pdl1_macrop.csv",
  markers_ctrl_pdl1_kuv   = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kuv_ctrl_pdl1_macrop.csv",
  markers_ctrl_dac_k25l   = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kpb25l_ctrl_dac_macrop.csv",
  markers_ctrl_dac_kuv    = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kuv_ctrl_dac_macrop.csv",
  markers_ctrl_comb_k25l  = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kpb25l_ctrl_comb_macrop.csv",
  markers_ctrl_comb_kuv   = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kuv_ctrl_comb_macrop.csv",
  markers_dac_comb_k25l   = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kpb25l_dac_comb_macrop.csv",
  markers_dac_comb_kuv    = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kuv_dac_comb_macrop.csv"
)

# Genes upregulados (log2FC > 1)
deg_upregulated <- lapply(deg_files, function(file) {
  df <- read.csv(file, row.names = 1)
  rownames(df[df$avg_log2FC > 1, ])
})

# Genes downregulados (log2FC < -1)
deg_downregulated <- lapply(deg_files, function(file) {
  df <- read.csv(file, row.names = 1)
  rownames(df[df$avg_log2FC < -1, ])
})


convert_to_entrez <- function(gene_list) {
  bitr(gene_list,
       fromType = "SYMBOL",
       toType = "ENTREZID",
       OrgDb = org.Mm.eg.db) %>% 
    pull(ENTREZID) %>% 
    unique()
}

# Convertir todas las listas
entrez_up <- lapply(deg_upregulated, convert_to_entrez)
entrez_down <- lapply(deg_downregulated, convert_to_entrez)


# UPREGULATED -------------------------------------------------------------


go_results_up <- lapply(entrez_up, function(gene_ids) {
  enrichGO(gene = gene_ids,
           OrgDb = org.Mm.eg.db,
           ont = "BP",
           pAdjustMethod = "BH",
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.2,
           readable = TRUE)
})

kegg_results_up <- lapply(entrez_up, function(gene_ids) {
  enrichKEGG(gene = gene_ids,
             organism = "mmu",
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05)
})

react_results_up <- lapply(entrez_up, function(gene_ids){
  enrichPathway(gene = gene_ids, 
                organism = "mouse", 
                pAdjustMethod = "BH", 
                pvalueCutoff = 0.05)
})


# DOWNREGULATED -----------------------------------------------------------


go_results_down <- lapply(entrez_down, function(gene_ids) {
  enrichGO(gene = gene_ids,
           OrgDb = org.Mm.eg.db,
           ont = "BP",
           pAdjustMethod = "BH",
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.2,
           readable = TRUE)
})

# ORA con KEGG
kegg_results_down <- lapply(entrez_down, function(gene_ids) {
  enrichKEGG(gene = gene_ids,
             organism = "mmu",
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05)
})

react_results_down <- lapply(entrez_down, function(gene_ids){
  enrichPathway(gene = gene_ids, 
                organism = "mouse", 
                pAdjustMethod = "BH", 
                pvalueCutoff = 0.05)
})


# PLOT --------------------------------------------------------------------


plot_dotplots <- function(result_list, result_type, regulation) {
  names(result_list) %>%
    walk(function(name) {
      result <- result_list[[name]]
      if (!is.null(result) && nrow(result) > 0) {
        
        # Detección de tratamientos
        treatments <- c()
        if (grepl("ctrl", name, ignore.case = TRUE)) treatments <- c(treatments, "Control")
        if (grepl("pdl1", name, ignore.case = TRUE)) treatments <- c(treatments, "PD-L1")
        if (grepl("dac", name, ignore.case = TRUE)) treatments <- c(treatments, "DAC")
        if (grepl("comb", name, ignore.case = TRUE)) treatments <- c(treatments, "Combination")
        treatment_label <- paste(treatments, collapse = " vs ")
        
        # Detección de línea celular
        cell_line <- ifelse(grepl("k25l", name, ignore.case = TRUE), "KPB25L",
                            ifelse(grepl("kuv", name, ignore.case = TRUE), "KPB25L-UV", "Unknown"))
        
        # Título final
        title <- paste(cell_line, "-", treatment_label, "-", result_type, "-", regulation, "Macrophages")
        filename <- paste0(gsub("-", "_", title), ".png")
        
        # Plot
        p <- dotplot(result, showCategory = 15) +
          ggtitle(title) +
          theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
        
        outdir <- "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_ORA/Macrophages/"
        ggsave(
          filename = file.path(outdir, filename),
          plot = p,
          width = 10, height = 8, dpi=600
        )
      }
    })
}


# Upregulated
plot_dotplots(go_results_up, "GO_BP", "UP")
plot_dotplots(kegg_results_up, "KEGG", "UP")
plot_dotplots(react_results_up, "REACTOME", "UP")

# Downregulated
plot_dotplots(go_results_down, "GO_BP", "DOWN")
plot_dotplots(kegg_results_down, "KEGG", "DOWN")
plot_dotplots(react_results_down, "REACTOME", "DOWN")





extract_enrichment_df <- function(result_list, result_type, regulation) {
  all_dfs <- list()
  
  for (name in names(result_list)) {
    res <- result_list[[name]]
    if (!is.null(res) && nrow(res) > 0) {
      df <- as.data.frame(res)
      
      # Extraer tratamientos
      treatments <- c()
      if (grepl("ctrl", name, ignore.case = TRUE)) treatments <- c(treatments, "Control")
      if (grepl("pdl1", name, ignore.case = TRUE)) treatments <- c(treatments, "PD-L1")
      if (grepl("dac", name, ignore.case = TRUE)) treatments <- c(treatments, "DAC")
      if (grepl("comb", name, ignore.case = TRUE)) treatments <- c(treatments, "Combination")
      treatment_label <- paste(treatments, collapse = " vs ")
      
      # Extraer línea celular
      cell_line <- ifelse(grepl("k25l", name, ignore.case = TRUE), "KPB25L",
                          ifelse(grepl("kuv", name, ignore.case = TRUE), "KPB25L-UV", "Unknown"))
      
      df$treatment <- treatment_label
      df$cell_line <- cell_line
      df$regulation <- regulation
      df$source <- result_type
      df$comparison <- name  # para trazabilidad
      
      all_dfs[[name]] <- df
    }
  }
  
  bind_rows(all_dfs)
}

all_enrichments_df <- bind_rows(
  extract_enrichment_df(go_results_up, "GO_BP", "UP"),
  extract_enrichment_df(go_results_down, "GO_BP", "DOWN"),
  extract_enrichment_df(kegg_results_up, "KEGG", "UP"),
  extract_enrichment_df(kegg_results_down, "KEGG", "DOWN"),
  extract_enrichment_df(react_results_up, "Reactome", "UP"),
  extract_enrichment_df(react_results_down, "Reactome", "DOWN")
)

all_enrichments_df$GeneRatio <- sapply(all_enrichments_df$GeneRatio, function(x) {
  if (is.numeric(x)) return(x)
  if (grepl("/", x)) eval(parse(text = x)) else as.numeric(x)
})


library(forcats)
plot_enrichment_summary <- function(df, cell_line, source_filter = "GO_BP", regulation_filter = "UP", top_terms = 10) {
  df_filtered <- df %>%
    filter(cell_line == !!cell_line,
           source == source_filter,
           regulation == regulation_filter) %>%
    group_by(treatment, Description) %>%
    slice_min(p.adjust, n = 1) %>%
    group_by(treatment) %>%
    slice_min(p.adjust, n = top_terms) %>%
    ungroup() %>%
    mutate(Description = fct_reorder(Description, GeneRatio))
  desired_order <- c("Control vs PD-L1", "Control vs DAC", "Control vs Combination", "DAC vs Combination")
  df_filtered$treatment <- factor(df_filtered$treatment, levels = desired_order)
  
  
  # Define título dinámico
  reg_label <- ifelse(regulation_filter == "UP", "Upregulated", "Downregulated")
  
  ggplot(df_filtered, aes(x = GeneRatio, y = Description, 
                          color = p.adjust, size = Count)) +
    geom_point() +
    scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
    scale_size_continuous(name = "Count") +
    facet_wrap(~treatment, scales = "free") +
    labs(
      title = paste(cell_line, "-", reg_label, " DEGs -", source_filter, "Pathway Enrichment in Macrophages"),
      x = "Gene Ratio", y = NULL
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(family = "Times New Roman", face = "bold", size = 11),
      axis.text.y = element_text(family = "Times New Roman", size = 13),
      plot.title = element_text(hjust = 0.5, family = "Times New Roman"),
      axis.title = element_text(),
      axis.text.x = element_text(size = 10),
      legend.title = element_text(),
      legend.text = element_text()
    )
}



plot_enrichment_summary(all_enrichments_df, "KPB25L", source_filter = "GO_BP")
plot_enrichment_summary(all_enrichments_df, "KPB25L-UV", source_filter = "GO_BP")
plot_enrichment_summary(all_enrichments_df, "KPB25L", source_filter = "KEGG")
plot_enrichment_summary(all_enrichments_df, "KPB25L-UV", source_filter = "KEGG")

plot_enrichment_summary(all_enrichments_df, "KPB25L", source_filter = "Reactome")
plot_enrichment_summary(all_enrichments_df, "KPB25L-UV", source_filter = "Reactome")



plot_enrichment_summary(all_enrichments_df, "KPB25L", regulation_filter = "DOWN", source_filter = "GO_BP")
plot_enrichment_summary(all_enrichments_df, "KPB25L-UV",  regulation_filter = "DOWN",source_filter = "GO_BP")
plot_enrichment_summary(all_enrichments_df, "KPB25L", regulation_filter = "DOWN", source_filter = "KEGG")
plot_enrichment_summary(all_enrichments_df, "KPB25L-UV", regulation_filter = "DOWN", source_filter = "KEGG")

plot_enrichment_summary(all_enrichments_df, "KPB25L", regulation_filter = "DOWN", source_filter = "Reactome")
plot_enrichment_summary(all_enrichments_df, "KPB25L-UV",  regulation_filter = "DOWN", source_filter = "Reactome")


cell_lines <- c("KPB25L", "KPB25L-UV")
sources <- c("GO_BP", "KEGG", "Reactome")
regulations <- c("UP", "DOWN")

for (cl in cell_lines) {
  for (src in sources) {
    for (reg in regulations) {
      outdir <- "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_ORA/Macrophages/"
      p <- plot_enrichment_summary(all_enrichments_df, cell_line = cl, source_filter = src, regulation_filter = reg)
      filename <- paste0(outdir,cl, "_", src, "_", reg, ".png")
      ggsave(filename, plot = p, width = 18, height = 12, dpi= 600)
    }
  }
}



# NK ----------------------------------------------------------------------

deg_files <- list(
  markers_ctrl_pdl1_k25l  = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kpb25l_ctrl_pdl1_nk.csv",
  markers_ctrl_pdl1_kuv   = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kuv_ctrl_pdl1_nk.csv",
  markers_ctrl_dac_k25l   = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kpb25l_ctrl_dac_nk.csv",
  markers_ctrl_dac_kuv    = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kuv_ctrl_dac_nk.csv",
  markers_ctrl_comb_k25l  = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kpb25l_ctrl_comb_nk.csv",
  markers_ctrl_comb_kuv   = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kuv_ctrl_comb_nk.csv",
  markers_dac_comb_k25l   = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kpb25l_dac_comb_nk.csv",
  markers_dac_comb_kuv    = "/home/workdir/HDAC_tfm/data/DEGs_1806/significant/DEGS_kuv_dac_comb_nk.csv"
)

# Genes upregulados (log2FC > 1)
deg_upregulated <- lapply(deg_files, function(file) {
  df <- read.csv(file, row.names = 1)
  rownames(df[df$avg_log2FC > 1, ])
})

# Genes downregulados (log2FC < -1)
deg_downregulated <- lapply(deg_files, function(file) {
  df <- read.csv(file, row.names = 1)
  rownames(df[df$avg_log2FC < -1, ])
})


convert_to_entrez <- function(gene_list) {
  bitr(gene_list,
       fromType = "SYMBOL",
       toType = "ENTREZID",
       OrgDb = org.Mm.eg.db) %>% 
    pull(ENTREZID) %>% 
    unique()
}

# Convertir todas las listas
entrez_up <- lapply(deg_upregulated, convert_to_entrez)
entrez_down <- lapply(deg_downregulated, convert_to_entrez)


# UPREGULATED -------------------------------------------------------------


go_results_up <- lapply(entrez_up, function(gene_ids) {
  enrichGO(gene = gene_ids,
           OrgDb = org.Mm.eg.db,
           ont = "BP",
           pAdjustMethod = "BH",
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.2,
           readable = TRUE)
})

kegg_results_up <- lapply(entrez_up, function(gene_ids) {
  enrichKEGG(gene = gene_ids,
             organism = "mmu",
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05)
})

react_results_up <- lapply(entrez_up, function(gene_ids){
  enrichPathway(gene = gene_ids, 
                organism = "mouse", 
                pAdjustMethod = "BH", 
                pvalueCutoff = 0.05)
})


# DOWNREGULATED -----------------------------------------------------------


go_results_down <- lapply(entrez_down, function(gene_ids) {
  enrichGO(gene = gene_ids,
           OrgDb = org.Mm.eg.db,
           ont = "BP",
           pAdjustMethod = "BH",
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.2,
           readable = TRUE)
})

# ORA con KEGG
kegg_results_down <- lapply(entrez_down, function(gene_ids) {
  enrichKEGG(gene = gene_ids,
             organism = "mmu",
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05)
})

react_results_down <- lapply(entrez_down, function(gene_ids){
  enrichPathway(gene = gene_ids, 
                organism = "mouse", 
                pAdjustMethod = "BH", 
                pvalueCutoff = 0.05)
})


# PLOT --------------------------------------------------------------------


plot_dotplots <- function(result_list, result_type, regulation) {
  names(result_list) %>%
    walk(function(name) {
      result <- result_list[[name]]
      if (!is.null(result) && nrow(result) > 0) {
        
        # Detección de tratamientos
        treatments <- c()
        if (grepl("ctrl", name, ignore.case = TRUE)) treatments <- c(treatments, "Control")
        if (grepl("pdl1", name, ignore.case = TRUE)) treatments <- c(treatments, "PD-L1")
        if (grepl("dac", name, ignore.case = TRUE)) treatments <- c(treatments, "DAC")
        if (grepl("comb", name, ignore.case = TRUE)) treatments <- c(treatments, "Combination")
        treatment_label <- paste(treatments, collapse = " vs ")
        
        # Detección de línea celular
        cell_line <- ifelse(grepl("k25l", name, ignore.case = TRUE), "KPB25L",
                            ifelse(grepl("kuv", name, ignore.case = TRUE), "KPB25L-UV", "Unknown"))
        
        # Título final
        title <- paste(cell_line, "-", treatment_label, "-", result_type, "-", regulation, "Natural killers")
        filename <- paste0(gsub("-", "_", title), ".png")
        
        # Plot
        p <- dotplot(result, showCategory = 15) +
          ggtitle(title) +
          theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
        
        outdir <- "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_ORA/NK/"
        ggsave(
          filename = file.path(outdir, filename),
          plot = p,
          width = 10, height = 8, dpi=600
        )
      }
    })
}


# Upregulated
plot_dotplots(go_results_up, "GO_BP", "UP")
plot_dotplots(kegg_results_up, "KEGG", "UP")
plot_dotplots(react_results_up, "REACTOME", "UP")

# Downregulated
plot_dotplots(go_results_down, "GO_BP", "DOWN")
plot_dotplots(kegg_results_down, "KEGG", "DOWN")
plot_dotplots(react_results_down, "REACTOME", "DOWN")





extract_enrichment_df <- function(result_list, result_type, regulation) {
  all_dfs <- list()
  
  for (name in names(result_list)) {
    res <- result_list[[name]]
    if (!is.null(res) && nrow(res) > 0) {
      df <- as.data.frame(res)
      
      # Extraer tratamientos
      treatments <- c()
      if (grepl("ctrl", name, ignore.case = TRUE)) treatments <- c(treatments, "Control")
      if (grepl("pdl1", name, ignore.case = TRUE)) treatments <- c(treatments, "PD-L1")
      if (grepl("dac", name, ignore.case = TRUE)) treatments <- c(treatments, "DAC")
      if (grepl("comb", name, ignore.case = TRUE)) treatments <- c(treatments, "Combination")
      treatment_label <- paste(treatments, collapse = " vs ")
      
      # Extraer línea celular
      cell_line <- ifelse(grepl("k25l", name, ignore.case = TRUE), "KPB25L",
                          ifelse(grepl("kuv", name, ignore.case = TRUE), "KPB25L-UV", "Unknown"))
      
      df$treatment <- treatment_label
      df$cell_line <- cell_line
      df$regulation <- regulation
      df$source <- result_type
      df$comparison <- name  # para trazabilidad
      
      all_dfs[[name]] <- df
    }
  }
  
  bind_rows(all_dfs)
}

all_enrichments_df <- bind_rows(
  extract_enrichment_df(go_results_up, "GO_BP", "UP"),
  extract_enrichment_df(go_results_down, "GO_BP", "DOWN"),
  extract_enrichment_df(kegg_results_up, "KEGG", "UP"),
  extract_enrichment_df(kegg_results_down, "KEGG", "DOWN"),
  extract_enrichment_df(react_results_up, "Reactome", "UP"),
  extract_enrichment_df(react_results_down, "Reactome", "DOWN")
)

all_enrichments_df$GeneRatio <- sapply(all_enrichments_df$GeneRatio, function(x) {
  if (is.numeric(x)) return(x)
  if (grepl("/", x)) eval(parse(text = x)) else as.numeric(x)
})


library(forcats)
plot_enrichment_summary <- function(df, cell_line, source_filter = "GO_BP", regulation_filter = "UP", top_terms = 15) {
  df_filtered <- df %>%
    filter(cell_line == !!cell_line,
           source == source_filter,
           regulation == regulation_filter) %>%
    group_by(treatment, Description) %>%
    slice_min(p.adjust, n = 1) %>%
    group_by(treatment) %>%
    slice_min(p.adjust, n = top_terms) %>%
    ungroup() %>%
    mutate(Description = fct_reorder(Description, GeneRatio))
  desired_order <- c("Control vs PD-L1", "Control vs DAC", "Control vs Combination", "DAC vs Combination")
  df_filtered$treatment <- factor(df_filtered$treatment, levels = desired_order)
  
  
  # Define título dinámico
  reg_label <- ifelse(regulation_filter == "UP", "Upregulated", "Downregulated")
  
  ggplot(df_filtered, aes(x = GeneRatio, y = Description, 
                          color = p.adjust, size = Count)) +
    geom_point() +
    scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
    scale_size_continuous(name = "Count") +
    facet_wrap(~treatment, scales = "free") +
    labs(
      title = paste(cell_line, "-", reg_label, " DEGs -", source_filter, "Pathway Enrichment in Natural killers"),
      x = "Gene Ratio", y = NULL
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(family = "Times New Roman", face = "bold", size = 11),
      axis.text.y = element_text(family = "Times New Roman", size = 13),
      plot.title = element_text(hjust = 0.5, family = "Times New Roman"),
      axis.title = element_text(),
      axis.text.x = element_text(size = 10),
      legend.title = element_text(),
      legend.text = element_text()
    )
}



plot_enrichment_summary(all_enrichments_df, "KPB25L", source_filter = "GO_BP")
plot_enrichment_summary(all_enrichments_df, "KPB25L-UV", source_filter = "GO_BP")
plot_enrichment_summary(all_enrichments_df, "KPB25L", source_filter = "KEGG")
plot_enrichment_summary(all_enrichments_df, "KPB25L-UV", source_filter = "KEGG")

plot_enrichment_summary(all_enrichments_df, "KPB25L", source_filter = "Reactome")
plot_enrichment_summary(all_enrichments_df, "KPB25L-UV", source_filter = "Reactome")



plot_enrichment_summary(all_enrichments_df, "KPB25L", regulation_filter = "DOWN", source_filter = "GO_BP")
plot_enrichment_summary(all_enrichments_df, "KPB25L-UV",  regulation_filter = "DOWN",source_filter = "GO_BP")
plot_enrichment_summary(all_enrichments_df, "KPB25L", regulation_filter = "DOWN", source_filter = "KEGG")
plot_enrichment_summary(all_enrichments_df, "KPB25L-UV", regulation_filter = "DOWN", source_filter = "KEGG")

plot_enrichment_summary(all_enrichments_df, "KPB25L", regulation_filter = "DOWN", source_filter = "Reactome")
plot_enrichment_summary(all_enrichments_df, "KPB25L-UV",  regulation_filter = "DOWN", source_filter = "Reactome")


cell_lines <- c("KPB25L", "KPB25L-UV")
sources <- c("GO_BP", "KEGG", "Reactome")
regulations <- c("UP", "DOWN")

for (cl in cell_lines) {
  for (src in sources) {
    for (reg in regulations) {
      outdir <- "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_ORA/NK/"
      p <- plot_enrichment_summary(all_enrichments_df, cell_line = cl, source_filter = src, regulation_filter = reg)
      filename <- paste0(outdir,cl, "_", src, "_", reg, ".png")
      ggsave(filename, plot = p, width = 15, height = 9, dpi= 600)
    }
  }
}



# TUMOR_UNANN -------------------------------------------------------------


deg_files <- list(
  markers_ctrl_pdl1_k25l  = "/home/workdir/HDAC_tfm/data/DEGs_1406/significant/DEGS_kpb25l_ctrl_pdl1.csv",
  markers_ctrl_pdl1_kuv   = "/home/workdir/HDAC_tfm/data/DEGs_1406/significant/DEGS_kuv_ctrl_pdl1.csv",
  markers_ctrl_dac_k25l   = "/home/workdir/HDAC_tfm/data/DEGs_1406/significant/DEGS_kpb25l_ctrl_dac.csv",
  markers_ctrl_dac_kuv    = "/home/workdir/HDAC_tfm/data/DEGs_1406/significant/DEGS_kuv_ctrl_dac.csv",
  markers_ctrl_comb_k25l  = "/home/workdir/HDAC_tfm/data/DEGs_1406/significant/DEGS_kpb25l_ctrl_comb.csv",
  markers_ctrl_comb_kuv   = "/home/workdir/HDAC_tfm/data/DEGs_1406/significant/DEGS_kuv_ctrl_comb.csv",
  markers_dac_comb_k25l   = "/home/workdir/HDAC_tfm/data/DEGs_1406/significant/DEGS_kpb25l_dac_comb.csv",
  markers_dac_comb_kuv    = "/home/workdir/HDAC_tfm/data/DEGs_1406/significant/DEGS_kuv_dac_comb.csv"
)

# Genes upregulados (log2FC > 1)
deg_upregulated <- lapply(deg_files, function(file) {
  df <- read.csv(file, row.names = 1)
  rownames(df[df$avg_log2FC > 1, ])
})

# Genes downregulados (log2FC < -1)
deg_downregulated <- lapply(deg_files, function(file) {
  df <- read.csv(file, row.names = 1)
  rownames(df[df$avg_log2FC < -1, ])
})


convert_to_entrez <- function(gene_list) {
  bitr(gene_list,
       fromType = "SYMBOL",
       toType = "ENTREZID",
       OrgDb = org.Mm.eg.db) %>% 
    pull(ENTREZID) %>% 
    unique()
}

# Convertir todas las listas
entrez_up <- lapply(deg_upregulated, convert_to_entrez)
entrez_down <- lapply(deg_downregulated, convert_to_entrez)


# UPREGULATED -------------------------------------------------------------


go_results_up <- lapply(entrez_up, function(gene_ids) {
  enrichGO(gene = gene_ids,
           OrgDb = org.Mm.eg.db,
           ont = "BP",
           pAdjustMethod = "BH",
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.2,
           readable = TRUE)
})

kegg_results_up <- lapply(entrez_up, function(gene_ids) {
  enrichKEGG(gene = gene_ids,
             organism = "mmu",
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05)
})

react_results_up <- lapply(entrez_up, function(gene_ids){
  enrichPathway(gene = gene_ids, 
                organism = "mouse", 
                pAdjustMethod = "BH", 
                pvalueCutoff = 0.05)
})


# DOWNREGULATED -----------------------------------------------------------


go_results_down <- lapply(entrez_down, function(gene_ids) {
  enrichGO(gene = gene_ids,
           OrgDb = org.Mm.eg.db,
           ont = "BP",
           pAdjustMethod = "BH",
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.2,
           readable = TRUE)
})

# ORA con KEGG
kegg_results_down <- lapply(entrez_down, function(gene_ids) {
  enrichKEGG(gene = gene_ids,
             organism = "mmu",
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05)
})

react_results_down <- lapply(entrez_down, function(gene_ids){
  enrichPathway(gene = gene_ids, 
                organism = "mouse", 
                pAdjustMethod = "BH", 
                pvalueCutoff = 0.05)
})


# PLOT --------------------------------------------------------------------


plot_dotplots <- function(result_list, result_type, regulation) {
  names(result_list) %>%
    walk(function(name) {
      result <- result_list[[name]]
      if (!is.null(result) && nrow(result) > 0) {
        
        # Detección de tratamientos
        treatments <- c()
        if (grepl("ctrl", name, ignore.case = TRUE)) treatments <- c(treatments, "Control")
        if (grepl("pdl1", name, ignore.case = TRUE)) treatments <- c(treatments, "PD-L1")
        if (grepl("dac", name, ignore.case = TRUE)) treatments <- c(treatments, "DAC")
        if (grepl("comb", name, ignore.case = TRUE)) treatments <- c(treatments, "Combination")
        treatment_label <- paste(treatments, collapse = " vs ")
        
        # Detección de línea celular
        cell_line <- ifelse(grepl("k25l", name, ignore.case = TRUE), "KPB25L",
                            ifelse(grepl("kuv", name, ignore.case = TRUE), "KPB25L-UV", "Unknown"))
        
        # Título final
        title <- paste(cell_line, "-", treatment_label, "-", result_type, "-", regulation, " Cancer + Unannotated")
        filename <- paste0(gsub("-", "_", title), ".png")
        
        # Plot
        p <- dotplot(result, showCategory = 15) +
          ggtitle(title) +
          theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
        
        outdir <- "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_ORA/tumor_unann/"
        ggsave(
          filename = file.path(outdir, filename),
          plot = p,
          width = 10, height = 8, dpi=600
        )
      }
    })
}


# Upregulated
plot_dotplots(go_results_up, "GO_BP", "UP")
plot_dotplots(kegg_results_up, "KEGG", "UP")
plot_dotplots(react_results_up, "REACTOME", "UP")

# Downregulated
plot_dotplots(go_results_down, "GO_BP", "DOWN")
plot_dotplots(kegg_results_down, "KEGG", "DOWN")
plot_dotplots(react_results_down, "REACTOME", "DOWN")





extract_enrichment_df <- function(result_list, result_type, regulation) {
  all_dfs <- list()
  
  for (name in names(result_list)) {
    res <- result_list[[name]]
    if (!is.null(res) && nrow(res) > 0) {
      df <- as.data.frame(res)
      
      # Extraer tratamientos
      treatments <- c()
      if (grepl("ctrl", name, ignore.case = TRUE)) treatments <- c(treatments, "Control")
      if (grepl("pdl1", name, ignore.case = TRUE)) treatments <- c(treatments, "PD-L1")
      if (grepl("dac", name, ignore.case = TRUE)) treatments <- c(treatments, "DAC")
      if (grepl("comb", name, ignore.case = TRUE)) treatments <- c(treatments, "Combination")
      treatment_label <- paste(treatments, collapse = " vs ")
      
      # Extraer línea celular
      cell_line <- ifelse(grepl("k25l", name, ignore.case = TRUE), "KPB25L",
                          ifelse(grepl("kuv", name, ignore.case = TRUE), "KPB25L-UV", "Unknown"))
      
      df$treatment <- treatment_label
      df$cell_line <- cell_line
      df$regulation <- regulation
      df$source <- result_type
      df$comparison <- name  # para trazabilidad
      
      all_dfs[[name]] <- df
    }
  }
  
  bind_rows(all_dfs)
}

all_enrichments_df <- bind_rows(
  extract_enrichment_df(go_results_up, "GO_BP", "UP"),
  extract_enrichment_df(go_results_down, "GO_BP", "DOWN"),
  extract_enrichment_df(kegg_results_up, "KEGG", "UP"),
  extract_enrichment_df(kegg_results_down, "KEGG", "DOWN"),
  extract_enrichment_df(react_results_up, "Reactome", "UP"),
  extract_enrichment_df(react_results_down, "Reactome", "DOWN")
)

all_enrichments_df$GeneRatio <- sapply(all_enrichments_df$GeneRatio, function(x) {
  if (is.numeric(x)) return(x)
  if (grepl("/", x)) eval(parse(text = x)) else as.numeric(x)
})


library(forcats)
plot_enrichment_summary <- function(df, cell_line, source_filter = "GO_BP", regulation_filter = "UP", top_terms = 15) {
  df_filtered <- df %>%
    filter(cell_line == !!cell_line,
           source == source_filter,
           regulation == regulation_filter) %>%
    group_by(treatment, Description) %>%
    slice_min(p.adjust, n = 1) %>%
    group_by(treatment) %>%
    slice_min(p.adjust, n = top_terms) %>%
    ungroup() %>%
    mutate(Description = fct_reorder(Description, GeneRatio))
  desired_order <- c("Control vs PD-L1", "Control vs DAC", "Control vs Combination", "DAC vs Combination")
  df_filtered$treatment <- factor(df_filtered$treatment, levels = desired_order)
  
  
  # Define título dinámico
  reg_label <- ifelse(regulation_filter == "UP", "Upregulated", "Downregulated")
  
  ggplot(df_filtered, aes(x = GeneRatio, y = Description, 
                          color = p.adjust, size = Count)) +
    geom_point() +
    scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
    scale_size_continuous(name = "Count") +
    facet_wrap(~treatment, scales = "free") +
    labs(
      title = paste(cell_line, "-", reg_label, " DEGs -", source_filter, "Pathway Enrichment in Cancer + Unnanotated"),
      x = "Gene Ratio", y = NULL
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(family = "Times New Roman", face = "bold", size = 11),
      axis.text.y = element_text(family = "Times New Roman", size = 13),
      plot.title = element_text(hjust = 0.5, family = "Times New Roman"),
      axis.title = element_text(),
      axis.text.x = element_text(size = 10),
      legend.title = element_text(),
      legend.text = element_text()
    )
}



plot_enrichment_summary(all_enrichments_df, "KPB25L", source_filter = "GO_BP")
plot_enrichment_summary(all_enrichments_df, "KPB25L-UV", source_filter = "GO_BP")
plot_enrichment_summary(all_enrichments_df, "KPB25L", source_filter = "KEGG")
plot_enrichment_summary(all_enrichments_df, "KPB25L-UV", source_filter = "KEGG")

plot_enrichment_summary(all_enrichments_df, "KPB25L", source_filter = "Reactome")
plot_enrichment_summary(all_enrichments_df, "KPB25L-UV", source_filter = "Reactome")



plot_enrichment_summary(all_enrichments_df, "KPB25L", regulation_filter = "DOWN", source_filter = "GO_BP")
plot_enrichment_summary(all_enrichments_df, "KPB25L-UV",  regulation_filter = "DOWN",source_filter = "GO_BP")
plot_enrichment_summary(all_enrichments_df, "KPB25L", regulation_filter = "DOWN", source_filter = "KEGG")
plot_enrichment_summary(all_enrichments_df, "KPB25L-UV", regulation_filter = "DOWN", source_filter = "KEGG")

plot_enrichment_summary(all_enrichments_df, "KPB25L", regulation_filter = "DOWN", source_filter = "Reactome")
plot_enrichment_summary(all_enrichments_df, "KPB25L-UV",  regulation_filter = "DOWN", source_filter = "Reactome")


cell_lines <- c("KPB25L", "KPB25L-UV")
sources <- c("GO_BP", "KEGG", "Reactome")
regulations <- c("UP", "DOWN")

for (cl in cell_lines) {
  for (src in sources) {
    for (reg in regulations) {
      outdir <- "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_ORA/tumor_unann/"
      p <- plot_enrichment_summary(all_enrichments_df, cell_line = cl, source_filter = src, regulation_filter = reg)
      filename <- paste0(outdir,cl, "_", src, "_", reg, ".png")
      ggsave(filename, plot = p, width = 15, height = 9, dpi= 600)
    }
  }
}

