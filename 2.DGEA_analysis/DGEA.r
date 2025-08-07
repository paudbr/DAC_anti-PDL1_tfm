###
#author :  Paula de Blas Rioja
# script for obtaining DEGs of each cell type annotated + Violin Plots


# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)

# -------------------------------------------------------------------------
# Load the integrated Seurat object containing labeled immune and tumor cells
rnaseq <- readRDS("/home/workdir/analysis_augie/20250120_labeled_immune.combined_edit.RDS")

# Subset tumor cells only (cluster ID = "1")
tumor_cells <- subset(rnaseq, idents = c("1"))

# Fix sample name if needed (standardize label for downstream grouping)
tumor_cells@meta.data$sample <- gsub("KAP25L", "KPB25L", tumor_cells@meta.data$sample)

# Join split assays/layers (necessary for FindMarkers if layers exist)
tumor_cells_joined <- JoinLayers(tumor_cells)

# -------------------------------------------------------------------------
# Find DEGs between PD-L1 and control treatment in KPB25L (tumor cells)
markers_ctrl_pdl1_k25l <- FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L_PD-L1",
  ident.2 = "KPB25L_control",
  group.by = "sample",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0.1
)
write.csv(markers_ctrl_pdl1_k25l, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_ctrl_pdl1.csv", row.names = TRUE)

# Find DEGs between DAC and control in KPB25L
markers_ctrl_dac_k25l <- FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L_DAC",
  ident.2 = "KPB25L_control",
  group.by = "sample",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0.1
)
write.csv(markers_ctrl_dac_k25l, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_ctrl_dac.csv", row.names = TRUE)

# Find DEGs between Combination and control in KPB25L
markers_ctrl_comb_k25l <- FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L_Combination",
  ident.2 = "KPB25L_control",
  group.by = "sample",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0.1
)
write.csv(markers_ctrl_comb_k25l, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_ctrl_comb.csv", row.names = TRUE)

# Find DEGs between Combination and DAC in KPB25L
markers_dac_comb_k25l <- FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L_Combination",
  ident.2 = "KPB25L_DAC",
  group.by = "sample",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0.1
)
write.csv(markers_dac_comb_k25l, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_dac_comb.csv", row.names = TRUE)

# -------------------------------------------------------------------------
# Repeat the DEG analysis for KPB25L-UV (UV-preconditioned) tumor cells

# PD-L1 vs Control
markers_ctrl_pdl1_kuv <- FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L-UV_PD-L1",
  ident.2 = "KPB25L-UV_control",
  group.by = "sample",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0.1
)

# DAC vs Control
markers_ctrl_dac_kuv <- FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L-UV_DAC",
  ident.2 = "KPB25L-UV_control",
  group.by = "sample",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0.1
)

# Combination vs Control
markers_ctrl_comb_kuv <- FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L-UV_Combination",
  ident.2 = "KPB25L-UV_control",
  group.by = "sample",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0.1
)

# Combination vs DAC
markers_dac_comb_kuv <- FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L-UV_Combination",
  ident.2 = "KPB25L-UV_DAC",
  group.by = "sample",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0.1
)

# Save all DEG results
write.csv(markers_ctrl_pdl1_kuv, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_ctrl_pdl1.csv", row.names = TRUE)
write.csv(markers_ctrl_dac_kuv, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_ctrl_dac.csv", row.names = TRUE)
write.csv(markers_ctrl_comb_kuv, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_ctrl_comb.csv", row.names = TRUE)
write.csv(markers_dac_comb_kuv, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_dac_comb.csv", row.names = TRUE)

# -------------------------------------------------------------------------
# Load all DEG files into a named list
deg_files <- list(
  markers_ctrl_pdl1_k25l  = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_ctrl_pdl1.csv",
  markers_ctrl_pdl1_kuv   = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_ctrl_pdl1.csv",
  markers_ctrl_dac_k25l   = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_ctrl_dac.csv",
  markers_ctrl_dac_kuv    = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_ctrl_dac.csv",
  markers_ctrl_comb_k25l  = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_ctrl_comb.csv",
  markers_ctrl_comb_kuv   = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_ctrl_comb.csv",
  markers_dac_comb_k25l   = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_dac_comb.csv",
  markers_dac_comb_kuv    = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_dac_comb.csv"
)

# Read and store all DEG dataframes
deg_filtered <- lapply(deg_files, function(file) {
  df <- read.csv(file, row.names = 1)
})
list2env(deg_filtered, envir = .GlobalEnv)

# -------------------------------------------------------------------------
# Process DEG tables to prepare for volcano plot visualization
deg_forplot <- lapply(deg_filtered, function(df) {
  df <- df %>%
    tibble::rownames_to_column(var = "gene") %>%
    mutate(
      color = case_when(
        abs(avg_log2FC) > 1 & p_val_adj < 0.05 & avg_log2FC > 0 ~ "Up",
        abs(avg_log2FC) > 1 & p_val_adj < 0.05 & avg_log2FC < 0 ~ "Down",
        TRUE ~ "grey"
      )
    )
  
  # Replace zero adjusted p-values to avoid Inf in log scale
  df$p_val_adj[df$p_val_adj == 0] <- 1e-300

  # Identify top 10 upregulated and top 10 downregulated genes for labeling
  top_up <- df %>% filter(color == "Up") %>%
    arrange(desc(abs(avg_log2FC) * -log10(p_val_adj))) %>%
    slice_head(n = 10)
  
  top_down <- df %>% filter(color == "Down") %>%
    arrange(desc(abs(avg_log2FC) * -log10(p_val_adj))) %>%
    slice_head(n = 10)
  
  top_genes <- c(top_up$gene, top_down$gene)
  df$extreme <- ifelse(df$gene %in% top_genes, df$gene, NA)
  
  return(df)
})

# Assign to global environment with new names
names(deg_forplot) <- paste0(names(deg_filtered), "_forplot")
list2env(deg_forplot, envir = .GlobalEnv)

# -------------------------------------------------------------------------
# Create volcano plots and save each as PNG
volcano_plots <- lapply(names(deg_forplot), function(name) {
  df <- deg_forplot[[name]]
  
  # Detect experimental conditions from name
  treatments <- c()
  if (grepl("ctrl", name, ignore.case = TRUE)) treatments <- c(treatments, "Control")
  if (grepl("pdl1", name, ignore.case = TRUE)) treatments <- c(treatments, "PD-L1")
  if (grepl("dac", name, ignore.case = TRUE)) treatments <- c(treatments, "DAC")
  if (grepl("comb", name, ignore.case = TRUE)) treatments <- c(treatments, "Combination")
  treatment_label <- paste(treatments, collapse = " vs ")

  # Detect cell line
  cell_line <- ifelse(grepl("k25l", name, ignore.case = TRUE), "KPB25L",
                      ifelse(grepl("kuv", name, ignore.case = TRUE), "KPB25L-UV", "Unknown"))

  # Define dynamic y-axis limit
  max_log10p <- min(350, ceiling(max(-log10(df$p_val_adj), na.rm = TRUE)) + 5)

  # Create volcano plot using ggplot2
  p <- ggplot(df, aes(x = avg_log2FC, y = -log10(p_val_adj), color = color)) +
    geom_point(alpha = 0.75, size = 1.8) +
    scale_color_manual(
      values = c("Up" = "#D7301F", "Down" = "#1F78B4", "grey" = "lightgrey"),
      labels = c("Up" = "Upregulated", "Down" = "Downregulated", "grey" = "Non-significant")
    ) +
    geom_text_repel(
      aes(label = extreme),
      na.rm = TRUE,
      size = 3.8,
      box.padding = 0.4,
      point.padding = 0.3,
      segment.color = "grey60",
      color = "black",
      max.overlaps = 20,
      force = 3,
      max.iter = 6000,
      min.segment.length = 0
    ) +
    coord_cartesian(ylim = c(0, max_log10p)) +
    theme_minimal(base_family = "Times New Roman") +
    labs(
      title = paste0("Volcano Plot - ", cell_line, " (", treatment_label, ") Cancer cells"),
      x = expression(log[2]~Fold~Change),
      y = expression(-log[10]~Adjusted~P~value),
      color = "Expression"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
    )

  # Save the plot to file (sanitize label for filename)
  safe_treatment_label <- gsub("[^A-Za-z0-9_]", "_", treatment_label)
  filepath <- paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/VP_", cell_line, "_", safe_treatment_label, "cancer.png")
  ggsave(filepath, plot = p, height = 6, width = 8, dpi = 600)

  return(p)
})


# UNANNOTATED -------------------------------------------------------------

tumor_cells <- subset(rnaseq, idents = c("0"))


tumor_cells@meta.data$sample <- gsub("KAP25L", "KPB25L", tumor_cells@meta.data$sample)

tumor_cells_joined <- JoinLayers(tumor_cells)

markers_ctrl_pdl1_k25l<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L_PD-L1",      
  ident.2 = "KPB25L_control",        
  group.by = "sample",      
  test.use = "wilcox",
  logfc.threshold = 0,      
  min.pct = 0.1               
) 

write.csv(markers_ctrl_pdl1_k25l, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_ctrl_pdl1_unann.csv", row.names = TRUE)


markers_ctrl_dac_k25l<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L_DAC",      # condición 1
  ident.2 = "KPB25L_control",         # condición 2 (comparación de referencia)
  group.by = "sample",      # columna que indica condiciones
  test.use = "wilcox",
  logfc.threshold = 0,      # mínimo log2 Fold Change (opcional)
  min.pct = 0.1                # mínimo % de células que expresan el gen (opcional)
) 

write.csv(markers_ctrl_dac_k25l, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_ctrl_dac_unann.csv", row.names = TRUE)



markers_ctrl_comb_k25l<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L_Combination",      # condición 1
  ident.2 = "KPB25L_control",         # condición 2 (comparación de referencia)
  group.by = "sample",      # columna que indica condiciones
  test.use = "wilcox",
  logfc.threshold = 0,      # mínimo log2 Fold Change (opcional)
  min.pct = 0.1                # mínimo % de células que expresan el gen (opcional)
) 

write.csv(markers_ctrl_comb_k25l, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_ctrl_comb_unann.csv", row.names = TRUE)



markers_dac_comb_k25l<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L_Combination",      # condición 1
  ident.2 = "KPB25L_DAC",         # condición 2 (comparación de referencia)
  group.by = "sample",      # columna que indica condiciones
  test.use = "wilcox",
  logfc.threshold = 0,      # mínimo log2 Fold Change (opcional)
  min.pct = 0.1                # mínimo % de células que expresan el gen (opcional)
) 

write.csv(markers_dac_comb_k25l, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_dac_comb_unann.csv", row.names = TRUE)


#KAP25L-UV

markers_ctrl_pdl1_kuv<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L-UV_PD-L1",      
  ident.2 = "KPB25L-UV_control",        
  group.by = "sample",      
  test.use = "wilcox",
  logfc.threshold = 0,      
  min.pct = 0.1               
) 

markers_ctrl_dac_kuv<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L-UV_DAC",      # condición 1
  ident.2 = "KPB25L-UV_control",         # condición 2 (comparación de referencia)
  group.by = "sample",      # columna que indica condiciones
  test.use = "wilcox",
  logfc.threshold = 0,      # mínimo log2 Fold Change (opcional)
  min.pct = 0.1                # mínimo % de células que expresan el gen (opcional)
) 



markers_ctrl_comb_kuv<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L-UV_Combination",      # condición 1
  ident.2 = "KPB25L-UV_control",         # condición 2 (comparación de referencia)
  group.by = "sample",      # columna que indica condiciones
  test.use = "wilcox",
  logfc.threshold = 0,      # mínimo log2 Fold Change (opcional)
  min.pct = 0.1                # mínimo % de células que expresan el gen (opcional)
) 


markers_dac_comb_kuv<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L-UV_Combination",      # condición 1
  ident.2 = "KPB25L-UV_DAC",         # condición 2 (comparación de referencia)
  group.by = "sample",      # columna que indica condiciones
  test.use = "wilcox",
  logfc.threshold = 0,      # mínimo log2 Fold Change (opcional)
  min.pct = 0.1                # mínimo % de células que expresan el gen (opcional)
) 
write.csv(markers_ctrl_pdl1_kuv, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_ctrl_pdl1_unann.csv", row.names = TRUE)

write.csv(markers_ctrl_dac_kuv, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_ctrl_dac_unann.csv", row.names = TRUE)

write.csv(markers_ctrl_comb_kuv, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_ctrl_comb_unann.csv", row.names = TRUE)

write.csv(markers_dac_comb_kuv, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_dac_comb_unann.csv", row.names = TRUE)

deg_files <- list(
  markers_ctrl_pdl1_k25l  = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_ctrl_pdl1_unann.csv",
  markers_ctrl_pdl1_kuv   = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_ctrl_pdl1_unann.csv",
  markers_ctrl_dac_k25l   = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_ctrl_dac_unann.csv",
  markers_ctrl_dac_kuv    = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_ctrl_dac_unann.csv",
  markers_ctrl_comb_k25l  = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_ctrl_comb_unann.csv",
  markers_ctrl_comb_kuv   = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_ctrl_comb_unann.csv",
  markers_dac_comb_k25l   = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_dac_comb_unann.csv",
  markers_dac_comb_kuv    = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_dac_comb_unann.csv"
)

deg_filtered <- lapply(deg_files, function(file) {
  df <- read.csv(file, row.names = 1)
})

list2env(deg_filtered, envir = .GlobalEnv)

deg_forplot <- lapply(deg_filtered, function(df) {
  df <- df %>%
    tibble::rownames_to_column(var = "gene") %>%
    mutate(
      color = case_when(
        abs(avg_log2FC) > 1 & p_val_adj < 0.05 & avg_log2FC > 0 ~ "Up",
        abs(avg_log2FC) > 1 & p_val_adj < 0.05 & avg_log2FC < 0 ~ "Down",
        TRUE ~ "grey"
      )
    )
  
  # Reemplazar p_val_adj = 0 por un valor pequeño para evitar Inf en log10
  df$p_val_adj[df$p_val_adj == 0] <- 1e-300
  
  # Seleccionar top 10 Up y top 10 Down significativos
  top_up <- df %>%
    filter(color == "Up") %>%
    arrange(desc(abs(avg_log2FC) * -log10(p_val_adj))) %>%
    slice_head(n = 10)
  
  top_down <- df %>%
    filter(color == "Down") %>%
    arrange(desc(abs(avg_log2FC) * -log10(p_val_adj))) %>%
    slice_head(n = 10)
  
  top_genes <- c(top_up$gene, top_down$gene)
  
  # Etiquetar solo los top
  df$extreme <- ifelse(df$gene %in% top_genes, df$gene, NA)
  
  return(df)
})


# Asignar al entorno global con nuevos nombres
names(deg_forplot) <- paste0(names(deg_filtered), "_forplot")
list2env(deg_forplot, envir = .GlobalEnv)


deg_forplot_1 <- deg_forplot[[1]]
# Asumiendo que tienes deg_forplot ya creada con los dataframes preparados


volcano_plots <- lapply(names(deg_forplot), function(name) {
  df <- deg_forplot[[name]]
  
  # Detect all treatments mentioned in the name
  treatments <- c()
  if (grepl("ctrl", name, ignore.case = TRUE)) treatments <- c(treatments, "Control")
  if (grepl("pdl1", name, ignore.case = TRUE)) treatments <- c(treatments, "PD-L1")
  if (grepl("dac", name, ignore.case = TRUE)) treatments <- c(treatments, "DAC")
  if (grepl("comb", name, ignore.case = TRUE)) treatments <- c(treatments, "Combination")
  
  treatment_label <- paste(treatments, collapse = " vs ")
  
  # Detect cell line from name
  cell_line <- ifelse(grepl("k25l", name, ignore.case = TRUE), "KPB25L",
                      ifelse(grepl("kuv", name, ignore.case = TRUE), "KPB25L-UV", "Unknown"))
  
  # Cálculo del máximo para ylim dinámico
  max_log10p <- min(305, ceiling(max(-log10(df$p_val_adj), na.rm = TRUE)) + 5)
  
  ggplot(df, aes(x = avg_log2FC, y = -log10(p_val_adj), color = color)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_manual(
      values = c("Up" = "red", "Down" = "blue", "grey" = "grey"),
      labels = c("Up" = "Upregulated", "Down" = "Downregulated", "grey" = "Non-significant")
    ) +
    theme_minimal() +
    labs(
      title = paste0("Volcano Plot - ", cell_line, " (", treatment_label, ") Unannotated "),
      x = "log2 Fold Change",
      y = "-log10 Adjusted P-value",
      color = "Expression"
    ) +
    geom_text_repel(
      aes(label = extreme),
      na.rm = TRUE,
      size = 3.5,
      box.padding = 0.3,
      point.padding = 0.3,
      segment.color = "grey50",
      color = "black",
      max.overlaps = 15,
      force = 2,
      max.iter = 5000,
      min.segment.length = 0
    ) +
    coord_cartesian(ylim = c(0, max_log10p)) +
    theme(
      text = element_text(family = "Times New Roman"),
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(),
      axis.text = element_text(),
      legend.title = element_text(),
      legend.text = element_text()
    )
})

print(volcano_plots[[1]])



library(ggplot2)
library(ggrepel)
volcano_plots <- lapply(names(deg_forplot), function(name) {
  df <- deg_forplot[[name]]
  
  # Detect treatments
  treatments <- c()
  if (grepl("ctrl", name, ignore.case = TRUE)) treatments <- c(treatments, "Control")
  if (grepl("pdl1", name, ignore.case = TRUE)) treatments <- c(treatments, "PD-L1")
  if (grepl("dac", name, ignore.case = TRUE)) treatments <- c(treatments, "DAC")
  if (grepl("comb", name, ignore.case = TRUE)) treatments <- c(treatments, "Combination")
  treatment_label <- paste(treatments, collapse = " vs ")
  
  # Detect cell line
  cell_line <- ifelse(grepl("k25l", name, ignore.case = TRUE), "KPB25L",
                      ifelse(grepl("kuv", name, ignore.case = TRUE), "KPB25L-UV", "Unknown"))
  
  # Calcular límite superior de Y dinámicamente (máximo -log10(p_val_adj))
  max_log10p <- min(350, ceiling(max(-log10(df$p_val_adj), na.rm = TRUE)) + 5)
  
  p <- ggplot(df, aes(x = avg_log2FC, y = -log10(p_val_adj), color = color)) +
    geom_point(alpha = 0.75, size = 1.8) +
    scale_color_manual(
      values = c("Up" = "#D7301F", "Down" = "#1F78B4", "grey" = "lightgrey"),
      labels = c("Up" = "Upregulated", "Down" = "Downregulated", "grey" = "Non-significant")
    ) +
    geom_text_repel(
      aes(label = extreme),
      na.rm = TRUE,
      size = 3.8,
      box.padding = 0.4,
      point.padding = 0.3,
      segment.color = "grey60",
      color = "black",
      max.overlaps = 20,
      force = 3,
      max.iter = 6000,
      min.segment.length = 0
    ) +
    coord_cartesian(ylim = c(0, max_log10p)) +
    theme_minimal(base_family = "Times New Roman") +
    labs(
      title = paste0("Volcano Plot - ", cell_line, " (", treatment_label, ") Unannotated"),
      x = expression(log[2]~Fold~Change),
      y = expression(-log[10]~Adjusted~P~value),
      color = "Expression"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
    )
  
  # Guardar el gráfico, sanitizando el nombre de archivo
  safe_treatment_label <- gsub("[^A-Za-z0-9_]", "_", treatment_label)
  filepath <- paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/VP_", cell_line, "_", safe_treatment_label, "unann.png")
  ggsave(filepath, plot = p, height = 6, width = 8, dpi = 600)
  
  return(p)
})


# MACROPHAGES---------------------------------------------------------



tumor_cells <- subset(rnaseq, idents = c("2"))


tumor_cells@meta.data$sample <- gsub("KAP25L", "KPB25L", tumor_cells@meta.data$sample)

tumor_cells_joined <- JoinLayers(tumor_cells)

markers_ctrl_pdl1_k25l<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L_PD-L1",      
  ident.2 = "KPB25L_control",        
  group.by = "sample",      
  test.use = "wilcox",
  logfc.threshold = 0,      
  min.pct = 0.1               
) 

write.csv(markers_ctrl_pdl1_k25l, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_ctrl_pdl1_macrop.csv", row.names = TRUE)


markers_ctrl_dac_k25l<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L_DAC",      # condición 1
  ident.2 = "KPB25L_control",         # condición 2 (comparación de referencia)
  group.by = "sample",      # columna que indica condiciones
  test.use = "wilcox",
  logfc.threshold = 0,      # mínimo log2 Fold Change (opcional)
  min.pct = 0.1                # mínimo % de células que expresan el gen (opcional)
) 

write.csv(markers_ctrl_dac_k25l, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_ctrl_dac_macrop.csv", row.names = TRUE)



markers_ctrl_comb_k25l<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L_Combination",      # condición 1
  ident.2 = "KPB25L_control",         # condición 2 (comparación de referencia)
  group.by = "sample",      # columna que indica condiciones
  test.use = "wilcox",
  logfc.threshold = 0,      # mínimo log2 Fold Change (opcional)
  min.pct = 0.1                # mínimo % de células que expresan el gen (opcional)
) 

write.csv(markers_ctrl_comb_k25l, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_ctrl_comb_macrop.csv", row.names = TRUE)



markers_dac_comb_k25l<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L_Combination",      # condición 1
  ident.2 = "KPB25L_DAC",         # condición 2 (comparación de referencia)
  group.by = "sample",      # columna que indica condiciones
  test.use = "wilcox",
  logfc.threshold = 0,      # mínimo log2 Fold Change (opcional)
  min.pct = 0.1                # mínimo % de células que expresan el gen (opcional)
) 

write.csv(markers_dac_comb_k25l, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_dac_comb_macrop.csv", row.names = TRUE)


#KAP25L-UV

markers_ctrl_pdl1_kuv<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L-UV_PD-L1",      
  ident.2 = "KPB25L-UV_control",        
  group.by = "sample",      
  test.use = "wilcox",
  logfc.threshold = 0,      
  min.pct = 0.1               
) 

markers_ctrl_dac_kuv<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L-UV_DAC",      # condición 1
  ident.2 = "KPB25L-UV_control",         # condición 2 (comparación de referencia)
  group.by = "sample",      # columna que indica condiciones
  test.use = "wilcox",
  logfc.threshold = 0,      # mínimo log2 Fold Change (opcional)
  min.pct = 0.1                # mínimo % de células que expresan el gen (opcional)
) 



markers_ctrl_comb_kuv<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L-UV_Combination",      # condición 1
  ident.2 = "KPB25L-UV_control",         # condición 2 (comparación de referencia)
  group.by = "sample",      # columna que indica condiciones
  test.use = "wilcox",
  logfc.threshold = 0,      # mínimo log2 Fold Change (opcional)
  min.pct = 0.1                # mínimo % de células que expresan el gen (opcional)
) 


markers_dac_comb_kuv<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L-UV_Combination",      # condición 1
  ident.2 = "KPB25L-UV_DAC",         # condición 2 (comparación de referencia)
  group.by = "sample",      # columna que indica condiciones
  test.use = "wilcox",
  logfc.threshold = 0,      # mínimo log2 Fold Change (opcional)
  min.pct = 0.1                # mínimo % de células que expresan el gen (opcional)
) 
write.csv(markers_ctrl_pdl1_kuv, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_ctrl_pdl1_macrop.csv", row.names = TRUE)

write.csv(markers_ctrl_dac_kuv, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_ctrl_dac_macrop.csv", row.names = TRUE)

write.csv(markers_ctrl_comb_kuv, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_ctrl_comb_macrop.csv", row.names = TRUE)

write.csv(markers_dac_comb_kuv, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_dac_comb_macrop.csv", row.names = TRUE)

deg_files <- list(
  markers_ctrl_pdl1_k25l  = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_ctrl_pdl1_macrop.csv",
  markers_ctrl_pdl1_kuv   = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_ctrl_pdl1_macrop.csv",
  markers_ctrl_dac_k25l   = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_ctrl_dac_macrop.csv",
  markers_ctrl_dac_kuv    = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_ctrl_dac_macrop.csv",
  markers_ctrl_comb_k25l  = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_ctrl_comb_macrop.csv",
  markers_ctrl_comb_kuv   = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_ctrl_comb_macrop.csv",
  markers_dac_comb_k25l   = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_dac_comb_macrop.csv",
  markers_dac_comb_kuv    = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_dac_comb_macrop.csv"
)

deg_filtered <- lapply(deg_files, function(file) {
  df <- read.csv(file, row.names = 1)
})

list2env(deg_filtered, envir = .GlobalEnv)

deg_forplot <- lapply(deg_filtered, function(df) {
  df <- df %>%
    tibble::rownames_to_column(var = "gene") %>%
    mutate(
      color = case_when(
        abs(avg_log2FC) > 1 & p_val_adj < 0.05 & avg_log2FC > 0 ~ "Up",
        abs(avg_log2FC) > 1 & p_val_adj < 0.05 & avg_log2FC < 0 ~ "Down",
        TRUE ~ "grey"
      )
    )
  
  # Reemplazar p_val_adj = 0 por un valor pequeño para evitar Inf en log10
  df$p_val_adj[df$p_val_adj == 0] <- 1e-300
  
  # Seleccionar top 10 Up y top 10 Down significativos
  top_up <- df %>%
    filter(color == "Up") %>%
    arrange(desc(abs(avg_log2FC) * -log10(p_val_adj))) %>%
    slice_head(n = 10)
  
  top_down <- df %>%
    filter(color == "Down") %>%
    arrange(desc(abs(avg_log2FC) * -log10(p_val_adj))) %>%
    slice_head(n = 10)
  
  top_genes <- c(top_up$gene, top_down$gene)
  
  # Etiquetar solo los top
  df$extreme <- ifelse(df$gene %in% top_genes, df$gene, NA)
  
  return(df)
})


# Asignar al entorno global con nuevos nombres
names(deg_forplot) <- paste0(names(deg_filtered), "_forplot")
list2env(deg_forplot, envir = .GlobalEnv)


deg_forplot_1 <- deg_forplot[[1]]
# Asumiendo que tienes deg_forplot ya creada con los dataframes preparados


volcano_plots <- lapply(names(deg_forplot), function(name) {
  df <- deg_forplot[[name]]
  
  # Detect all treatments mentioned in the name
  treatments <- c()
  if (grepl("ctrl", name, ignore.case = TRUE)) treatments <- c(treatments, "Control")
  if (grepl("pdl1", name, ignore.case = TRUE)) treatments <- c(treatments, "PD-L1")
  if (grepl("dac", name, ignore.case = TRUE)) treatments <- c(treatments, "DAC")
  if (grepl("comb", name, ignore.case = TRUE)) treatments <- c(treatments, "Combination")
  
  treatment_label <- paste(treatments, collapse = " vs ")
  
  # Detect cell line from name
  cell_line <- ifelse(grepl("k25l", name, ignore.case = TRUE), "KPB25L",
                      ifelse(grepl("kuv", name, ignore.case = TRUE), "KPB25L-UV", "Unknown"))
  
  # Cálculo del máximo para ylim dinámico
  max_log10p <- min(305, ceiling(max(-log10(df$p_val_adj), na.rm = TRUE)) + 5)
  
  ggplot(df, aes(x = avg_log2FC, y = -log10(p_val_adj), color = color)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_manual(
      values = c("Up" = "red", "Down" = "blue", "grey" = "grey"),
      labels = c("Up" = "Upregulated", "Down" = "Downregulated", "grey" = "Non-significant")
    ) +
    theme_minimal() +
    labs(
      title = paste0("Volcano Plot - ", cell_line, " (", treatment_label, ") Macrophages "),
      x = "log2 Fold Change",
      y = "-log10 Adjusted P-value",
      color = "Expression"
    ) +
    geom_text_repel(
      aes(label = extreme),
      na.rm = TRUE,
      size = 3.5,
      box.padding = 0.3,
      point.padding = 0.3,
      segment.color = "grey50",
      color = "black",
      max.overlaps = 15,
      force = 2,
      max.iter = 5000,
      min.segment.length = 0
    ) +
    coord_cartesian(ylim = c(0, max_log10p)) +
    theme(
      text = element_text(family = "Times New Roman"),
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(),
      axis.text = element_text(),
      legend.title = element_text(),
      legend.text = element_text()
    )
})

print(volcano_plots[[1]])



library(ggplot2)
library(ggrepel)
volcano_plots <- lapply(names(deg_forplot), function(name) {
  df <- deg_forplot[[name]]
  
  # Detect treatments
  treatments <- c()
  if (grepl("ctrl", name, ignore.case = TRUE)) treatments <- c(treatments, "Control")
  if (grepl("pdl1", name, ignore.case = TRUE)) treatments <- c(treatments, "PD-L1")
  if (grepl("dac", name, ignore.case = TRUE)) treatments <- c(treatments, "DAC")
  if (grepl("comb", name, ignore.case = TRUE)) treatments <- c(treatments, "Combination")
  treatment_label <- paste(treatments, collapse = " vs ")
  
  # Detect cell line
  cell_line <- ifelse(grepl("k25l", name, ignore.case = TRUE), "KPB25L",
                      ifelse(grepl("kuv", name, ignore.case = TRUE), "KPB25L-UV", "Unknown"))
  
  # Calcular límite superior de Y dinámicamente (máximo -log10(p_val_adj))
  max_log10p <- min(350, ceiling(max(-log10(df$p_val_adj), na.rm = TRUE)) + 5)
  
  p <- ggplot(df, aes(x = avg_log2FC, y = -log10(p_val_adj), color = color)) +
    geom_point(alpha = 0.75, size = 1.8) +
    scale_color_manual(
      values = c("Up" = "#D7301F", "Down" = "#1F78B4", "grey" = "lightgrey"),
      labels = c("Up" = "Upregulated", "Down" = "Downregulated", "grey" = "Non-significant")
    ) +
    geom_text_repel(
      aes(label = extreme),
      na.rm = TRUE,
      size = 3.8,
      box.padding = 0.4,
      point.padding = 0.3,
      segment.color = "grey60",
      color = "black",
      max.overlaps = 20,
      force = 3,
      max.iter = 6000,
      min.segment.length = 0
    ) +
    coord_cartesian(ylim = c(0, max_log10p)) +
    theme_minimal(base_family = "Times New Roman") +
    labs(
      title = paste0("Volcano Plot - ", cell_line, " (", treatment_label, ") Macrophages"),
      x = expression(log[2]~Fold~Change),
      y = expression(-log[10]~Adjusted~P~value),
      color = "Expression"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
    )
  
  # Guardar el gráfico, sanitizando el nombre de archivo
  safe_treatment_label <- gsub("[^A-Za-z0-9_]", "_", treatment_label)
  filepath <- paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/VP_", cell_line, "_", safe_treatment_label, "macrop.png")
  ggsave(filepath, plot = p, height = 6, width = 8, dpi = 600)
  
  return(p)
})


# NK ----------------------------------------------------------------------

tumor_cells <- subset(rnaseq, idents = c("3"))


tumor_cells@meta.data$sample <- gsub("KAP25L", "KPB25L", tumor_cells@meta.data$sample)

tumor_cells_joined <- JoinLayers(tumor_cells)

markers_ctrl_pdl1_k25l<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L_PD-L1",      
  ident.2 = "KPB25L_control",        
  group.by = "sample",      
  test.use = "wilcox",
  logfc.threshold = 0,      
  min.pct = 0.1               
) 

write.csv(markers_ctrl_pdl1_k25l, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_ctrl_pdl1_nk.csv", row.names = TRUE)


markers_ctrl_dac_k25l<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L_DAC",      # condición 1
  ident.2 = "KPB25L_control",         # condición 2 (comparación de referencia)
  group.by = "sample",      # columna que indica condiciones
  test.use = "wilcox",
  logfc.threshold = 0,      # mínimo log2 Fold Change (opcional)
  min.pct = 0.1                # mínimo % de células que expresan el gen (opcional)
) 

write.csv(markers_ctrl_dac_k25l, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_ctrl_dac_nk.csv", row.names = TRUE)



markers_ctrl_comb_k25l<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L_Combination",      # condición 1
  ident.2 = "KPB25L_control",         # condición 2 (comparación de referencia)
  group.by = "sample",      # columna que indica condiciones
  test.use = "wilcox",
  logfc.threshold = 0,      # mínimo log2 Fold Change (opcional)
  min.pct = 0.1                # mínimo % de células que expresan el gen (opcional)
) 

write.csv(markers_ctrl_comb_k25l, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_ctrl_comb_nk.csv", row.names = TRUE)



markers_dac_comb_k25l<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L_Combination",      # condición 1
  ident.2 = "KPB25L_DAC",         # condición 2 (comparación de referencia)
  group.by = "sample",      # columna que indica condiciones
  test.use = "wilcox",
  logfc.threshold = 0,      # mínimo log2 Fold Change (opcional)
  min.pct = 0.1                # mínimo % de células que expresan el gen (opcional)
) 

write.csv(markers_dac_comb_k25l, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_dac_comb_nk.csv", row.names = TRUE)


#KAP25L-UV

markers_ctrl_pdl1_kuv<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L-UV_PD-L1",      
  ident.2 = "KPB25L-UV_control",        
  group.by = "sample",      
  test.use = "wilcox",
  logfc.threshold = 0,      
  min.pct = 0.1               
) 

markers_ctrl_dac_kuv<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L-UV_DAC",      # condición 1
  ident.2 = "KPB25L-UV_control",         # condición 2 (comparación de referencia)
  group.by = "sample",      # columna que indica condiciones
  test.use = "wilcox",
  logfc.threshold = 0,      # mínimo log2 Fold Change (opcional)
  min.pct = 0.1                # mínimo % de células que expresan el gen (opcional)
) 



markers_ctrl_comb_kuv<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L-UV_Combination",      # condición 1
  ident.2 = "KPB25L-UV_control",         # condición 2 (comparación de referencia)
  group.by = "sample",      # columna que indica condiciones
  test.use = "wilcox",
  logfc.threshold = 0,      # mínimo log2 Fold Change (opcional)
  min.pct = 0.1                # mínimo % de células que expresan el gen (opcional)
) 


markers_dac_comb_kuv<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L-UV_Combination",      # condición 1
  ident.2 = "KPB25L-UV_DAC",         # condición 2 (comparación de referencia)
  group.by = "sample",      # columna que indica condiciones
  test.use = "wilcox",
  logfc.threshold = 0,      # mínimo log2 Fold Change (opcional)
  min.pct = 0.1                # mínimo % de células que expresan el gen (opcional)
) 
write.csv(markers_ctrl_pdl1_kuv, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_ctrl_pdl1_nk.csv", row.names = TRUE)

write.csv(markers_ctrl_dac_kuv, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_ctrl_dac_nk.csv", row.names = TRUE)

write.csv(markers_ctrl_comb_kuv, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_ctrl_comb_nk.csv", row.names = TRUE)

write.csv(markers_dac_comb_kuv, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_dac_comb_nk.csv", row.names = TRUE)

deg_files <- list(
  markers_ctrl_pdl1_k25l  = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_ctrl_pdl1_nk.csv",
  markers_ctrl_pdl1_kuv   = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_ctrl_pdl1_nk.csv",
  markers_ctrl_dac_k25l   = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_ctrl_dac_nk.csv",
  markers_ctrl_dac_kuv    = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_ctrl_dac_nk.csv",
  markers_ctrl_comb_k25l  = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_ctrl_comb_nk.csv",
  markers_ctrl_comb_kuv   = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_ctrl_comb_nk.csv",
  markers_dac_comb_k25l   = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_dac_comb_nk.csv",
  markers_dac_comb_kuv    = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_dac_comb_nk.csv"
)

deg_filtered <- lapply(deg_files, function(file) {
  df <- read.csv(file, row.names = 1)
})

list2env(deg_filtered, envir = .GlobalEnv)

deg_forplot <- lapply(deg_filtered, function(df) {
  df <- df %>%
    tibble::rownames_to_column(var = "gene") %>%
    mutate(
      color = case_when(
        abs(avg_log2FC) > 1 & p_val_adj < 0.05 & avg_log2FC > 0 ~ "Up",
        abs(avg_log2FC) > 1 & p_val_adj < 0.05 & avg_log2FC < 0 ~ "Down",
        TRUE ~ "grey"
      )
    )
  
  # Reemplazar p_val_adj = 0 por un valor pequeño para evitar Inf en log10
  df$p_val_adj[df$p_val_adj == 0] <- 1e-300
  
  # Seleccionar top 10 Up y top 10 Down significativos
  top_up <- df %>%
    filter(color == "Up") %>%
    arrange(desc(abs(avg_log2FC) * -log10(p_val_adj))) %>%
    slice_head(n = 10)
  
  top_down <- df %>%
    filter(color == "Down") %>%
    arrange(desc(abs(avg_log2FC) * -log10(p_val_adj))) %>%
    slice_head(n = 10)
  
  top_genes <- c(top_up$gene, top_down$gene)
  
  # Etiquetar solo los top
  df$extreme <- ifelse(df$gene %in% top_genes, df$gene, NA)
  
  return(df)
})


# Asignar al entorno global con nuevos nombres
names(deg_forplot) <- paste0(names(deg_filtered), "_forplot")
list2env(deg_forplot, envir = .GlobalEnv)


deg_forplot_1 <- deg_forplot[[4]]
# Asumiendo que tienes deg_forplot ya creada con los dataframes preparados


volcano_plots <- lapply(names(deg_forplot), function(name) {
  df <- deg_forplot[[name]]
  
  # Detect all treatments mentioned in the name
  treatments <- c()
  if (grepl("ctrl", name, ignore.case = TRUE)) treatments <- c(treatments, "Control")
  if (grepl("pdl1", name, ignore.case = TRUE)) treatments <- c(treatments, "PD-L1")
  if (grepl("dac", name, ignore.case = TRUE)) treatments <- c(treatments, "DAC")
  if (grepl("comb", name, ignore.case = TRUE)) treatments <- c(treatments, "Combination")
  
  treatment_label <- paste(treatments, collapse = " vs ")
  
  # Detect cell line from name
  cell_line <- ifelse(grepl("k25l", name, ignore.case = TRUE), "KPB25L",
                      ifelse(grepl("kuv", name, ignore.case = TRUE), "KPB25L-UV", "Unknown"))
  
  # Cálculo del máximo para ylim dinámico
  max_log10p <- min(305, ceiling(max(-log10(df$p_val_adj), na.rm = TRUE)) + 5)
  
  ggplot(df, aes(x = avg_log2FC, y = -log10(p_val_adj), color = color)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_manual(
      values = c("Up" = "red", "Down" = "blue", "grey" = "grey"),
      labels = c("Up" = "Upregulated", "Down" = "Downregulated", "grey" = "Non-significant")
    ) +
    theme_minimal() +
    labs(
      title = paste0("Volcano Plot - ", cell_line, " (", treatment_label, ") Macrophages "),
      x = "log2 Fold Change",
      y = "-log10 Adjusted P-value",
      color = "Expression"
    ) +
    geom_text_repel(
      aes(label = extreme),
      na.rm = TRUE,
      size = 3.5,
      box.padding = 0.3,
      point.padding = 0.3,
      segment.color = "grey50",
      color = "black",
      max.overlaps = 15,
      force = 2,
      max.iter = 5000,
      min.segment.length = 0
    ) +
    coord_cartesian(ylim = c(0, max_log10p)) +
    theme(
      text = element_text(family = "Times New Roman"),
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(),
      axis.text = element_text(),
      legend.title = element_text(),
      legend.text = element_text()
    )
})




library(ggplot2)
library(ggrepel)
volcano_plots <- lapply(names(deg_forplot), function(name) {
  df <- deg_forplot[[name]]
  
  # Detect treatments
  treatments <- c()
  if (grepl("ctrl", name, ignore.case = TRUE)) treatments <- c(treatments, "Control")
  if (grepl("pdl1", name, ignore.case = TRUE)) treatments <- c(treatments, "PD-L1")
  if (grepl("dac", name, ignore.case = TRUE)) treatments <- c(treatments, "DAC")
  if (grepl("comb", name, ignore.case = TRUE)) treatments <- c(treatments, "Combination")
  treatment_label <- paste(treatments, collapse = " vs ")
  
  # Detect cell line
  cell_line <- ifelse(grepl("k25l", name, ignore.case = TRUE), "KPB25L",
                      ifelse(grepl("kuv", name, ignore.case = TRUE), "KPB25L-UV", "Unknown"))
  
  # Calcular límite superior de Y dinámicamente (máximo -log10(p_val_adj))
  max_log10p <- min(350, ceiling(max(-log10(df$p_val_adj), na.rm = TRUE)) + 5)
  
  p <- ggplot(df, aes(x = avg_log2FC, y = -log10(p_val_adj), color = color)) +
    geom_point(alpha = 0.75, size = 1.8) +
    scale_color_manual(
      values = c("Up" = "#D7301F", "Down" = "#1F78B4", "grey" = "lightgrey"),
      labels = c("Up" = "Upregulated", "Down" = "Downregulated", "grey" = "Non-significant")
    ) +
    geom_text_repel(
      aes(label = extreme),
      na.rm = TRUE,
      size = 3.8,
      box.padding = 0.4,
      point.padding = 0.3,
      segment.color = "grey60",
      color = "black",
      max.overlaps = 20,
      force = 3,
      max.iter = 6000,
      min.segment.length = 0
    ) +
    coord_cartesian(ylim = c(0, max_log10p)) +
    theme_minimal(base_family = "Times New Roman") +
    labs(
      title = paste0("Volcano Plot - ", cell_line, " (", treatment_label, ") NK cells"),
      x = expression(log[2]~Fold~Change),
      y = expression(-log[10]~Adjusted~P~value),
      color = "Expression"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
    )
  
  # Guardar el gráfico, sanitizando el nombre de archivo
  safe_treatment_label <- gsub("[^A-Za-z0-9_]", "_", treatment_label)
  filepath <- paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/VP_", cell_line, "_", safe_treatment_label, "nk.png")
  ggsave(filepath, plot = p, height = 6, width = 8, dpi = 600)
  
  return(p)
})


# Tcells ------------------------------------------------------------------

tumor_cells <- subset(rnaseq, idents = c("9"))


tumor_cells@meta.data$sample <- gsub("KAP25L", "KPB25L", tumor_cells@meta.data$sample)

tumor_cells_joined <- JoinLayers(tumor_cells)

markers_ctrl_pdl1_k25l<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L_PD-L1",      
  ident.2 = "KPB25L_control",        
  group.by = "sample",      
  test.use = "wilcox",
  logfc.threshold = 0,      
  min.pct = 0.1               
) 

write.csv(markers_ctrl_pdl1_k25l, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_ctrl_pdl1_t.csv", row.names = TRUE)


markers_ctrl_dac_k25l<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L_DAC",      # condición 1
  ident.2 = "KPB25L_control",         # condición 2 (comparación de referencia)
  group.by = "sample",      # columna que indica condiciones
  test.use = "wilcox",
  logfc.threshold = 0,      # mínimo log2 Fold Change (opcional)
  min.pct = 0.1                # mínimo % de células que expresan el gen (opcional)
) 

write.csv(markers_ctrl_dac_k25l, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_ctrl_dac_t.csv", row.names = TRUE)



markers_ctrl_comb_k25l<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L_Combination",      # condición 1
  ident.2 = "KPB25L_control",         # condición 2 (comparación de referencia)
  group.by = "sample",      # columna que indica condiciones
  test.use = "wilcox",
  logfc.threshold = 0,      # mínimo log2 Fold Change (opcional)
  min.pct = 0.1                # mínimo % de células que expresan el gen (opcional)
) 

write.csv(markers_ctrl_comb_k25l, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_ctrl_comb_t.csv", row.names = TRUE)



markers_dac_comb_k25l<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L_Combination",      # condición 1
  ident.2 = "KPB25L_DAC",         # condición 2 (comparación de referencia)
  group.by = "sample",      # columna que indica condiciones
  test.use = "wilcox",
  logfc.threshold = 0,      # mínimo log2 Fold Change (opcional)
  min.pct = 0.1                # mínimo % de células que expresan el gen (opcional)
) 

write.csv(markers_dac_comb_k25l, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_dac_comb_t.csv", row.names = TRUE)


#KAP25L-UV

markers_ctrl_pdl1_kuv<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L-UV_PD-L1",      
  ident.2 = "KPB25L-UV_control",        
  group.by = "sample",      
  test.use = "wilcox",
  logfc.threshold = 0,      
  min.pct = 0.1               
) 

markers_ctrl_dac_kuv<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L-UV_DAC",      # condición 1
  ident.2 = "KPB25L-UV_control",         # condición 2 (comparación de referencia)
  group.by = "sample",      # columna que indica condiciones
  test.use = "wilcox",
  logfc.threshold = 0,      # mínimo log2 Fold Change (opcional)
  min.pct = 0.1                # mínimo % de células que expresan el gen (opcional)
) 



markers_ctrl_comb_kuv<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L-UV_Combination",      # condición 1
  ident.2 = "KPB25L-UV_control",         # condición 2 (comparación de referencia)
  group.by = "sample",      # columna que indica condiciones
  test.use = "wilcox",
  logfc.threshold = 0,      # mínimo log2 Fold Change (opcional)
  min.pct = 0.1                # mínimo % de células que expresan el gen (opcional)
) 


markers_dac_comb_kuv<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L-UV_Combination",      # condición 1
  ident.2 = "KPB25L-UV_DAC",         # condición 2 (comparación de referencia)
  group.by = "sample",      # columna que indica condiciones
  test.use = "wilcox",
  logfc.threshold = 0,      # mínimo log2 Fold Change (opcional)
  min.pct = 0.1                # mínimo % de células que expresan el gen (opcional)
) 
write.csv(markers_ctrl_pdl1_kuv, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_ctrl_pdl1_t.csv", row.names = TRUE)

write.csv(markers_ctrl_dac_kuv, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_ctrl_dac_t.csv", row.names = TRUE)

write.csv(markers_ctrl_comb_kuv, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_ctrl_comb_t.csv", row.names = TRUE)

write.csv(markers_dac_comb_kuv, file = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_dac_comb_t.csv", row.names = TRUE)

deg_files <- list(
  markers_ctrl_pdl1_k25l  = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_ctrl_pdl1_t.csv",
  markers_ctrl_pdl1_kuv   = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_ctrl_pdl1_t.csv",
  markers_ctrl_dac_k25l   = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_ctrl_dac_t.csv",
  markers_ctrl_dac_kuv    = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_ctrl_dac_t.csv",
  markers_ctrl_comb_k25l  = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_ctrl_comb_t.csv",
  markers_ctrl_comb_kuv   = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_ctrl_comb_t.csv",
  markers_dac_comb_k25l   = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_dac_comb_t.csv",
  markers_dac_comb_kuv    = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_dac_comb_t.csv"
)

deg_filtered <- lapply(deg_files, function(file) {
  df <- read.csv(file, row.names = 1)
})

list2env(deg_filtered, envir = .GlobalEnv)

deg_forplot <- lapply(deg_filtered, function(df) {
  df <- df %>%
    tibble::rownames_to_column(var = "gene") %>%
    mutate(
      color = case_when(
        abs(avg_log2FC) > 1 & p_val_adj < 0.05 & avg_log2FC > 0 ~ "Up",
        abs(avg_log2FC) > 1 & p_val_adj < 0.05 & avg_log2FC < 0 ~ "Down",
        TRUE ~ "grey"
      )
    )
  
  # Reemplazar p_val_adj = 0 por un valor pequeño para evitar Inf en log10
  df$p_val_adj[df$p_val_adj == 0] <- 1e-300
  
  # Seleccionar top 10 Up y top 10 Down significativos
  top_up <- df %>%
    filter(color == "Up") %>%
    arrange(desc(abs(avg_log2FC) * -log10(p_val_adj))) %>%
    slice_head(n = 10)
  
  top_down <- df %>%
    filter(color == "Down") %>%
    arrange(desc(abs(avg_log2FC) * -log10(p_val_adj))) %>%
    slice_head(n = 10)
  
  top_genes <- c(top_up$gene, top_down$gene)
  
  # Etiquetar solo los top
  df$extreme <- ifelse(df$gene %in% top_genes, df$gene, NA)
  
  return(df)
})


# Asignar al entorno global con nuevos nombres
names(deg_forplot) <- paste0(names(deg_filtered), "_forplot")
list2env(deg_forplot, envir = .GlobalEnv)


deg_forplot_1 <- deg_forplot[[4]]
# Asumiendo que tienes deg_forplot ya creada con los dataframes preparados


volcano_plots <- lapply(names(deg_forplot), function(name) {
  df <- deg_forplot[[name]]
  
  # Detect all treatments mentioned in the name
  treatments <- c()
  if (grepl("ctrl", name, ignore.case = TRUE)) treatments <- c(treatments, "Control")
  if (grepl("pdl1", name, ignore.case = TRUE)) treatments <- c(treatments, "PD-L1")
  if (grepl("dac", name, ignore.case = TRUE)) treatments <- c(treatments, "DAC")
  if (grepl("comb", name, ignore.case = TRUE)) treatments <- c(treatments, "Combination")
  
  treatment_label <- paste(treatments, collapse = " vs ")
  
  # Detect cell line from name
  cell_line <- ifelse(grepl("k25l", name, ignore.case = TRUE), "KPB25L",
                      ifelse(grepl("kuv", name, ignore.case = TRUE), "KPB25L-UV", "Unknown"))
  
  # Cálculo del máximo para ylim dinámico
  max_log10p <- min(305, ceiling(max(-log10(df$p_val_adj), na.rm = TRUE)) + 5)
  
  ggplot(df, aes(x = avg_log2FC, y = -log10(p_val_adj), color = color)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_manual(
      values = c("Up" = "red", "Down" = "blue", "grey" = "grey"),
      labels = c("Up" = "Upregulated", "Down" = "Downregulated", "grey" = "Non-significant")
    ) +
    theme_minimal() +
    labs(
      title = paste0("Volcano Plot - ", cell_line, " (", treatment_label, ") Macrophages "),
      x = "log2 Fold Change",
      y = "-log10 Adjusted P-value",
      color = "Expression"
    ) +
    geom_text_repel(
      aes(label = extreme),
      na.rm = TRUE,
      size = 3.5,
      box.padding = 0.3,
      point.padding = 0.3,
      segment.color = "grey50",
      color = "black",
      max.overlaps = 15,
      force = 2,
      max.iter = 5000,
      min.segment.length = 0
    ) +
    coord_cartesian(ylim = c(0, max_log10p)) +
    theme(
      text = element_text(family = "Times New Roman"),
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(),
      axis.text = element_text(),
      legend.title = element_text(),
      legend.text = element_text()
    )
})




library(ggplot2)
library(ggrepel)
volcano_plots <- lapply(names(deg_forplot), function(name) {
  df <- deg_forplot[[name]]
  
  # Detect treatments
  treatments <- c()
  if (grepl("ctrl", name, ignore.case = TRUE)) treatments <- c(treatments, "Control")
  if (grepl("pdl1", name, ignore.case = TRUE)) treatments <- c(treatments, "PD-L1")
  if (grepl("dac", name, ignore.case = TRUE)) treatments <- c(treatments, "DAC")
  if (grepl("comb", name, ignore.case = TRUE)) treatments <- c(treatments, "Combination")
  treatment_label <- paste(treatments, collapse = " vs ")
  
  # Detect cell line
  cell_line <- ifelse(grepl("k25l", name, ignore.case = TRUE), "KPB25L",
                      ifelse(grepl("kuv", name, ignore.case = TRUE), "KPB25L-UV", "Unknown"))
  
  # Calcular límite superior de Y dinámicamente (máximo -log10(p_val_adj))
  max_log10p <- min(350, ceiling(max(-log10(df$p_val_adj), na.rm = TRUE)) + 5)
  
  p <- ggplot(df, aes(x = avg_log2FC, y = -log10(p_val_adj), color = color)) +
    geom_point(alpha = 0.75, size = 1.8) +
    scale_color_manual(
      values = c("Up" = "#D7301F", "Down" = "#1F78B4", "grey" = "lightgrey"),
      labels = c("Up" = "Upregulated", "Down" = "Downregulated", "grey" = "Non-significant")
    ) +
    geom_text_repel(
      aes(label = extreme),
      na.rm = TRUE,
      size = 3.8,
      box.padding = 0.4,
      point.padding = 0.3,
      segment.color = "grey60",
      color = "black",
      max.overlaps = 20,
      force = 3,
      max.iter = 6000,
      min.segment.length = 0
    ) +
    coord_cartesian(ylim = c(0, max_log10p)) +
    theme_minimal(base_family = "Times New Roman") +
    labs(
      title = paste0("Volcano Plot - ", cell_line, " (", treatment_label, ") Naive CD8+ T cells"),
      x = expression(log[2]~Fold~Change),
      y = expression(-log[10]~Adjusted~P~value),
      color = "Expression"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
    )
  
  # Guardar el gráfico, sanitizando el nombre de archivo
  safe_treatment_label <- gsub("[^A-Za-z0-9_]", "_", treatment_label)
  filepath <- paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/VP_", cell_line, "_", safe_treatment_label, "t.png")
  ggsave(filepath, plot = p, height = 6, width = 8, dpi = 600)
  
  return(p)
})



# TUMOR CELLS AND UNANNOTATED ---------------------------------------------

tumor_cells <- readRDS("/home/workdir/HDAC_tfm/data/tumor_cells_adjusted_13032025.rds")
tumor_cells@meta.data$sample <- gsub("KAP25L", "KPB25L", tumor_cells@meta.data$sample)

tumor_cells_joined <- JoinLayers(tumor_cells)

markers_ctrl_pdl1_k25l<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L_PD-L1",      
  ident.2 = "KPB25L_control",        
  group.by = "sample",      
  test.use = "wilcox",
  logfc.threshold = 0,      
  min.pct = 0.1               
) 

write.csv(markers_ctrl_pdl1_k25l, file = "/home/workdir/HDAC_tfm/data/DEGs_1406/DEGS_kpb25l_ctrl_pdl1.csv", row.names = TRUE)


markers_ctrl_dac_k25l<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L_DAC",      # condición 1
  ident.2 = "KPB25L_control",         # condición 2 (comparación de referencia)
  group.by = "sample",      # columna que indica condiciones
  test.use = "wilcox",
  logfc.threshold = 0,      # mínimo log2 Fold Change (opcional)
  min.pct = 0.1                # mínimo % de células que expresan el gen (opcional)
) 

write.csv(markers_ctrl_dac_k25l, file = "/home/workdir/HDAC_tfm/data/DEGs_1406/DEGS_kpb25l_ctrl_dac.csv", row.names = TRUE)



markers_ctrl_comb_k25l<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L_Combination",      # condición 1
  ident.2 = "KPB25L_control",         # condición 2 (comparación de referencia)
  group.by = "sample",      # columna que indica condiciones
  test.use = "wilcox",
  logfc.threshold = 0,      # mínimo log2 Fold Change (opcional)
  min.pct = 0.1                # mínimo % de células que expresan el gen (opcional)
) 

write.csv(markers_ctrl_comb_k25l, file = "/home/workdir/HDAC_tfm/data/DEGs_1406/DEGS_kpb25l_ctrl_comb.csv", row.names = TRUE)



markers_dac_comb_k25l<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L_Combination",      # condición 1
  ident.2 = "KPB25L_DAC",         # condición 2 (comparación de referencia)
  group.by = "sample",      # columna que indica condiciones
  test.use = "wilcox",
  logfc.threshold = 0,      # mínimo log2 Fold Change (opcional)
  min.pct = 0.1                # mínimo % de células que expresan el gen (opcional)
) 

write.csv(markers_dac_comb_k25l, file = "/home/workdir/HDAC_tfm/data/DEGs_1406/DEGS_kpb25l_dac_comb.csv", row.names = TRUE)


#KAP25L-UV

markers_ctrl_pdl1_kuv<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L-UV_PD-L1",      
  ident.2 = "KPB25L-UV_control",        
  group.by = "sample",      
  test.use = "wilcox",
  logfc.threshold = 0,      
  min.pct = 0.1               
) 

markers_ctrl_dac_kuv<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L-UV_DAC",      # condición 1
  ident.2 = "KPB25L-UV_control",         # condición 2 (comparación de referencia)
  group.by = "sample",      # columna que indica condiciones
  test.use = "wilcox",
  logfc.threshold = 0,      # mínimo log2 Fold Change (opcional)
  min.pct = 0.1                # mínimo % de células que expresan el gen (opcional)
) 



markers_ctrl_comb_kuv<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L-UV_Combination",      # condición 1
  ident.2 = "KPB25L-UV_control",         # condición 2 (comparación de referencia)
  group.by = "sample",      # columna que indica condiciones
  test.use = "wilcox",
  logfc.threshold = 0,      # mínimo log2 Fold Change (opcional)
  min.pct = 0.1                # mínimo % de células que expresan el gen (opcional)
) 


markers_dac_comb_kuv<-FindMarkers(
  object = tumor_cells_joined,
  ident.1 = "KPB25L-UV_Combination",      # condición 1
  ident.2 = "KPB25L-UV_DAC",         # condición 2 (comparación de referencia)
  group.by = "sample",      # columna que indica condiciones
  test.use = "wilcox",
  logfc.threshold = 0,      # mínimo log2 Fold Change (opcional)
  min.pct = 0.1                # mínimo % de células que expresan el gen (opcional)
) 
write.csv(markers_ctrl_pdl1_kuv, file = "/home/workdir/HDAC_tfm/data/DEGs_1406/DEGS_kuv_ctrl_pdl1.csv", row.names = TRUE)

write.csv(markers_ctrl_dac_kuv, file = "/home/workdir/HDAC_tfm/data/DEGs_1406/DEGS_kuv_ctrl_dac.csv", row.names = TRUE)

write.csv(markers_ctrl_comb_kuv, file = "/home/workdir/HDAC_tfm/data/DEGs_1406/DEGS_kuv_ctrl_comb.csv", row.names = TRUE)

write.csv(markers_dac_comb_kuv, file = "/home/workdir/HDAC_tfm/data/DEGs_1406/DEGS_kuv_dac_comb.csv", row.names = TRUE)

deg_files <- list(
  markers_ctrl_pdl1_k25l  = "/home/workdir/HDAC_tfm/data/DEGs_1406/DEGS_kpb25l_ctrl_pdl1.csv",
  markers_ctrl_pdl1_kuv   = "/home/workdir/HDAC_tfm/data/DEGs_1406/DEGS_kuv_ctrl_pdl1.csv",
  markers_ctrl_dac_k25l   = "/home/workdir/HDAC_tfm/data/DEGs_1406/DEGS_kpb25l_ctrl_dac.csv",
  markers_ctrl_dac_kuv    = "/home/workdir/HDAC_tfm/data/DEGs_1406/DEGS_kuv_ctrl_dac.csv",
  markers_ctrl_comb_k25l  = "/home/workdir/HDAC_tfm/data/DEGs_1406/DEGS_kpb25l_ctrl_comb.csv",
  markers_ctrl_comb_kuv   = "/home/workdir/HDAC_tfm/data/DEGs_1406/DEGS_kuv_ctrl_comb.csv",
  markers_dac_comb_k25l   = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kpb25l_dac_comb.csv",
  markers_dac_comb_kuv    = "/home/workdir/HDAC_tfm/data/DEGs_1806/DEGS_kuv_dac_comb.csv"
)

deg_filtered <- lapply(deg_files, function(file) {
  df <- read.csv(file, row.names = 1)
})
list2env(deg_filtered, envir = .GlobalEnv)

deg_forplot <- lapply(deg_filtered, function(df) {
  df <- df %>%
    tibble::rownames_to_column(var = "gene") %>%
    mutate(
      color = case_when(
        abs(avg_log2FC) > 1 & p_val_adj < 0.05 & avg_log2FC > 0 ~ "Up",
        abs(avg_log2FC) > 1 & p_val_adj < 0.05 & avg_log2FC < 0 ~ "Down",
        TRUE ~ "grey"
      )
    )
  
  # Reemplazar p_val_adj = 0 por un valor pequeño para evitar Inf en log10
  df$p_val_adj[df$p_val_adj == 0] <- 1e-300
  
  # Seleccionar top 10 Up y top 10 Down significativos
  top_up <- df %>%
    filter(color == "Up") %>%
    arrange(desc(abs(avg_log2FC) * -log10(p_val_adj))) %>%
    slice_head(n = 10)
  
  top_down <- df %>%
    filter(color == "Down") %>%
    arrange(desc(abs(avg_log2FC) * -log10(p_val_adj))) %>%
    slice_head(n = 10)
  
  top_genes <- c(top_up$gene, top_down$gene)
  
  # Etiquetar solo los top
  df$extreme <- ifelse(df$gene %in% top_genes, df$gene, NA)
  
  return(df)
})


# Asignar al entorno global con nuevos nombres
names(deg_forplot) <- paste0(names(deg_filtered), "_forplot")
list2env(deg_forplot, envir = .GlobalEnv)


deg_forplot_1 <- deg_forplot[[1]]
# Asumiendo que tienes deg_forplot ya creada con los dataframes preparados

volcano_plots <- lapply(names(deg_forplot), function(name) {
  df <- deg_forplot[[name]]
  
  # Detect all treatments mentioned in the name
  treatments <- c()
  if (grepl("ctrl", name, ignore.case = TRUE)) treatments <- c(treatments, "Control")
  if (grepl("pdl1", name, ignore.case = TRUE)) treatments <- c(treatments, "PD-L1")
  if (grepl("dac", name, ignore.case = TRUE)) treatments <- c(treatments, "DAC")
  if (grepl("comb", name, ignore.case = TRUE)) treatments <- c(treatments, "Combination")
  
  treatment_label <- paste(treatments, collapse = " vs ")
  
  # Detect cell line from name
  cell_line <- ifelse(grepl("k25l", name, ignore.case = TRUE), "KPB25L",
                      ifelse(grepl("kuv", name, ignore.case = TRUE), "KPB25L-UV", "Unknown"))
  
  ggplot(df, aes(x = avg_log2FC, y = -log10(p_val_adj), color = color)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "grey" = "grey"),
                       labels = c("Up" = "Upregulated", "Down" = "Downregulated", "grey" = "Non-significant")) +
    theme_minimal() +
    labs(
      title = paste0("Volcano Plot - ", cell_line, " (", treatment_label, ")"),
      x = "log2 Fold Change",
      y = "-log10 Adjusted P-value",
      color = "Expression"
    ) +
    geom_text_repel(
      aes(label = extreme),
      na.rm = TRUE,
      size = 2.5,
      box.padding = 0.3,
      point.padding = 0.3,
      segment.color = "grey50",
      color = "black",
      max.overlaps = 15,
      force = 2,
      max.iter = 5000,
      min.segment.length = 0
    ) +
    coord_cartesian(ylim = c(0, 350)) +
    theme(
      text = element_text(family = "Times New Roman"),
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(),
      axis.text = element_text(),
      legend.title = element_text(),
      legend.text = element_text()
    )
})

print(volcano_plots[[5]])



library(ggplot2)
library(ggrepel)

volcano_plots <- lapply(names(deg_forplot), function(name) {
  df <- deg_forplot[[name]]
  
  # Detect treatments
  treatments <- c()
  if (grepl("ctrl", name, ignore.case = TRUE)) treatments <- c(treatments, "Control")
  if (grepl("pdl1", name, ignore.case = TRUE)) treatments <- c(treatments, "PD-L1")
  if (grepl("dac", name, ignore.case = TRUE)) treatments <- c(treatments, "DAC")
  if (grepl("comb", name, ignore.case = TRUE)) treatments <- c(treatments, "Combination")
  treatment_label <- paste(treatments, collapse = " vs ")
  
  # Detect cell line
  cell_line <- ifelse(grepl("k25l", name, ignore.case = TRUE), "KPB25L",
                      ifelse(grepl("kuv", name, ignore.case = TRUE), "KPB25L-UV", "Unknown"))
  
  p<- ggplot(df, aes(x = avg_log2FC, y = -log10(p_val_adj), color = color)) +
    geom_point(alpha = 0.75, size = 1.8) +
    scale_color_manual(
      values = c("Up" = "#D7301F", "Down" = "#1F78B4", "grey" = "lightgrey"),
      labels = c("Up" = "Upregulated", "Down" = "Downregulated", "grey" = "Non-significant")
    ) +
    geom_text_repel(
      aes(label = extreme),
      na.rm = TRUE,
      size = 3.8,
      box.padding = 0.4,
      point.padding = 0.3,
      segment.color = "grey60",
      color = "black",
      max.overlaps = 20,
      force = 3,
      max.iter = 6000,
      min.segment.length = 0
    ) +
    coord_cartesian(ylim = c(0, 350)) +
    theme_minimal(base_family = "Times New Roman") +
    labs(
      title = paste0("Volcano Plot - ", cell_line, " (", treatment_label, ") Cancer cells"),
      x = expression(log[2]~Fold~Change),
      y = expression(-log[10]~Adjusted~P~value),
      color = "Expression"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
    )
  
  ggsave(paste("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/VP_",cell_line,"_", treatment_label,"tumor_and_unann.png"), plot= p, height = 6, width = 8, dpi=600)
})



# SAVING_SIGNIFICANTLY_EXPRESSED ------------------------------------------

input_dir <- "/home/workdir/HDAC_tfm/data/DEGs_1806"
output_dir <- file.path(input_dir, "significant")

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

filter_and_save <- function(file) {

  df <- read.csv(file)
  
  df_filtered <- subset(df, p_val_adj < 0.05 & abs(avg_log2FC) > 1)
  
  file_name <- basename(file)
  
  write.csv(df_filtered, file = file.path(output_dir, file_name), row.names = FALSE)
}

lapply(files, filter_and_save)

input_dir <- "/home/workdir/HDAC_tfm/data/DEGs_1406"

files_joined <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
lapply(files_joined, filter_and_save)
