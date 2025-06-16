library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
# DEGS CANCER CELLS -------------------------------------------------------

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
  markers_ctrl_comb_kuv   = "/home/workdir/HDAC_tfm/data/DEGs_1406/DEGS_kuv_ctrl_comb.csv"
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
      ),
      extreme = ifelse(abs(avg_log2FC) > 2 & -log10(p_val_adj) > 150, gene, NA)
    )
  return(df)
})

# Asignar al entorno global con nuevos nombres
names(deg_forplot) <- paste0(names(deg_filtered), "_forplot")
list2env(deg_forplot, envir = .GlobalEnv)



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
      color = "black"  
    ) +
    coord_cartesian(ylim = c(0, 500)) +
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