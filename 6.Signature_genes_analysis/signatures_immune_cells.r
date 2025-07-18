library(IOBR)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(readxl)
library(dplyr)
library(stringr)
library(purrr)
library(car)
library(ggsignif)

immune_combined <- readRDS("/home/workdir/analysis_augie/20250120_labeled_immune.combined_edit.RDS")

immune_cells <- subset(immune_combined, immuno %in% c("Cancer cells", "Unannotated") == FALSE)

rownames(immune_cells) <- toupper(rownames(immune_cells))

immune_cells$sample <- gsub("KAP25L", "KPB25L", immune_cells$sample)

gene_sets <- list(
  "ICB_resistance_Peng_et_al" = signature_tme$ICB_resistance_Peng_et_al,
  "Cell_cycle" = signature_tme$Cell_cycle, 
  "CellCycle_Reg" = signature_tme$CellCycle_Reg, 
  "TGFb_Family_Member_Li" = signature_tme$TGFb_Family_Member_Li_et_al,
  "TGFb_Family_Member_Receptor" = signature_tme$TGFb_Family_Member_Receptor_Li_et_al,
  "TIP_Killing_of_cancer_cells_1" = signature_tme$TIP_Killing_of_cancer_cells_1,
  "TIP_Killing_of_cancer_cells_2" = signature_tme$TIP_Killing_of_cancer_cells_2,
  "TIP_Infiltration_of_immune_cells_into_tumors_1" = signature_tme$TIP_Infiltration_of_immune_cells_into_tumors_1,
  "TIP_Infiltration_of_immune_cells_into_tumors_2" = signature_tme$TIP_Infiltration_of_immune_cells_into_tumors_2,
  "TIP_Recognition_of_cancer_cells_by_T_cells_1" = signature_tme$TIP_Recognition_of_cancer_cells_by_T_cells_1,
  "TIP_Recognition_of_cancer_cells_by_T_cells_2" = signature_tme$TIP_Recognition_of_cancer_cells_by_T_cells_2,
  "TIP_Priming_and_activation_2" = signature_tme$TIP_Priming_and_activation_2,
  "TIP_Trafficking_of_immune_cells_to_tumors" = signature_tme$TIP_Trafficking_of_immune_cells_to_tumors,
  "TIP_Priming_and_activation_1" = signature_tme$TIP_Priming_and_activation_1,
  "TIP_Cancer_antigen_presentation_1" = signature_tme$TIP_Cancer_antigen_presentation_1,
  "TIP_Cancer_antigen_presentation_2" = signature_tme$TIP_Cancer_antigen_presentation_2,
  "Immune_Checkpoint" = signature_tme$Immune_Checkpoint,
  "TLS_Nature" = signature_tme$TLS_Nature,
  "TMEscoreA_plus" = signature_tme$TMEscoreA_plus,
  "TMEscoreB_plus" = signature_tme$TMEscoreB_plus,
  "TMEscoreA_CIR" = signature_tme$TMEscoreA_CIR,
  "TMEscoreB_CIR" = signature_tme$TMEscoreB_CIR,
  "TIP_Release_of_cancer_cell_antigens" = signature_tme$TIP_Release_of_cancer_cell_antigens,
  "T_cell_accumulation" = signature_tme$T_cell_accumulation_Peng_et_al,
  "T_cell_exhaustion" = signature_tme$T_cell_exhaustion_Peng_et_al,
  "T_cell_regulatory" = signature_tme$T_cell_regulatory_Peng_et_al,
  "Interleukins_Li" = signature_tme$Interleukins_Li_et_al,
  "Interleukins_Receptor_Li" = signature_tme$Interferon_Receptor_Li_et_al,
  "Interferons_Li" = signature_tme$Interferons_Li_et_al,
  "Interferon_Receptor_Li" = signature_tme$Interferon_Receptor_Li_et_al,
  "Cytokines_Li" = signature_tme$Cytokines_Li_et_al,
  "Cytokine_Receptors_Li" = signature_tme$Cytokine_Receptors_Li_et_al,
  "Chemokines_Li" = signature_tme$Chemokines_Li_et_al,
  "Chemokine_Receptors_Li" = signature_tme$Chemokine_Receptors_Li_et_al,
  "TNF_Family_Members_Li" = signature_tme$TNF_Family_Members_Li_et_al,
  "TNF_Family_Members_Receptors_Li" = signature_tme$TNF_Family_Members_Receptors_Li_et_al,
  "ICB_resistance_Peng_et_al" = signature_tme$ICB_resistance_Peng_et_al,
  "MeTIL_Jeschke_et_al" = signature_tme$MeTIL_Jeschke_et_al,
  "T_cell_probe_Jeschke" = signature_tme$T_cell_probe_Jeschke_et_al,
  "Co_stimulation_APC_Rooney" = signature_tme$Co_stimulation_APC_Rooney_et_al,
  "Co_stimulation_T_cell_Rooney" = signature_tme$Co_stimulation_T_cell_Rooney_et_al,
  "Co_inhibition_APC_Rooney" = signature_tme$Co_inhibition_APC_Rooney_et_al,
  "Co_inhibition_T_cell_Rooney" = signature_tme$Co_inhibition_T_cell_Rooney_et_al,
  "Type_I_IFN_Reponse_Rooney" = signature_tme$Type_I_IFN_Reponse_Rooney_et_al,
  "Type_II_IFN_Reponse_Rooney" = signature_tme$Type_II_IFN_Reponse_Rooney_et_al,
  "Cytolytic_Activity_Rooney" = signature_tme$Cytolytic_Activity_Rooney_et_al
)

immune_cells <- JoinLayers(immune_cells)
rownames(immune_cells) <- toupper(rownames(immune_cells))

signature_plots <- function(database, signature){
  cell_cycle_genes <- database
  
  # Añadir el score de la firma al objeto Seurat
  immune_cells <- AddModuleScore(object = immune_cells, 
                                 features = list(cell_cycle_genes), 
                                 name = signature)
  colnames(immune_cells@meta.data)
  
  # Definir la variable con la columna generada
  signature_col <- paste0(signature, "1")
  
  # Violin Plot normal con etiquetas rotadas
  p1 <- VlnPlot(immune_cells, features = signature_col, group.by = "immuno") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  # Feature Plot en UMAP
  p2 <- FeaturePlot(immune_cells, 
                    features = signature_col, 
                    reduction = "umap",
                    pt.size = 0.5)
  
  # Guardar imágenes
  ggsave(paste0("/home/workdir/HDAC_tfm/fig/nodownsampled/signatures_immune/Violin_plot_", signature, ".png"), 
         plot = p1, width = 10, height = 8, dpi = 600)
  ggsave(paste0("/home/workdir/HDAC_tfm/fig/nodownsampled/signatures_immune/Feature_plot_", signature, ".png"), 
         plot = p2, width = 10, height = 8, dpi = 600)
  
  # Extraer datos para ggplot
  df <- FetchData(immune_cells, vars = c(signature_col, "immuno", "sample"))
  # Extraer el prefijo K25L o K25L_UV
  df$category <- ifelse(grepl("UV", df$sample), "KPB25L_UV", "KPB25L")
  
  # Ordenar factor por categorías
  df$category <- factor(df$category, levels = c("KPB25L", "KPB25L_UV"))
  
  # Ordenar los samples dentro de cada categoría
  df$sample <- factor(df$sample, levels = c(
    "KPB25L_control", "KPB25L_PD-L1", "KPB25L_DAC", "KPB25L_Combination",
    "KPB25L-UV_control", "KPB25L-UV_PD-L1", "KPB25L-UV_DAC", "KPB25L-UV_Combination"
  ))
  
  
  df_violin_macrof <- subset(df, immuno == "Macrophages")
  
  df_violin_macrof$sample <- factor(df_violin_macrof$sample, levels = c(
    "KPB25L_control", "KPB25L_PD-L1", "KPB25L_DAC", "KPB25L_Combination",
    "KPB25L-UV_control", "KPB25L-UV_PD-L1", "KPB25L-UV_DAC", "KPB25L-UV_Combination"
  ))
  # Violin Plot separado por sample y con filas por categoría, con etiquetas rotadas
  p3 <- ggplot(df_violin_macrof, aes(x = "", y = .data[[signature_col]], fill = immuno)) + 
    geom_violin() + 
    stat_summary(fun = median, geom = "point", color = "black", size = 2) +  # Marca la mediana
    facet_wrap(~ sample, ncol = 8) +  # Un gráfico por tipo de célula inmune
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_blank(),  # Oculta el eje X vacío
          axis.ticks.x = element_blank()) +
    labs(x = NULL, y = "Expresión",title = signature_col)  # Elimina el eje X vacío y etiqueta eje Y
  
  print(p3)
  
  
  ggsave(paste0("/home/workdir/HDAC_tfm/fig/nodownsampled/signatures_immune/Violin_plot_", signature, "_byimmune.png"), 
         plot = p3, width = 20, height = 12, dpi = 600)
}

lapply(names(gene_sets), function(signature_name) {
  signature_plots(gene_sets[[signature_name]], signature_name)
})



# M1 SIGNATURE ------------------------------------------------------------


immune_signatures_ir <- read_excel("/home/workdir/HDAC_tfm/data/immune_signatures.xlsx")


macropages_m1_sig <- immune_signatures_ir %>%
  filter(str_starts(Celltype_Source_ID, "Macrophages M1"))

macropages_m1_sig <- macropages_m1_sig[, 3:ncol(macropages_m1_sig)]


macropages_m1_sig_full_list <- c()
row_lists <- apply(macropages_m1_sig, 1, as.list)

row_list_unlist <- unlist(row_lists, use.names = FALSE)

macrophages_m1_sig_full_list <- row_list_unlist[!is.na(unlist(row_lists, use.names = FALSE))]


immune_cells <- AddModuleScore(object = immune_cells, 
                               features = list(macrophages_m1_sig_full_list), 
                               name = "Macrophagues_M")

df <- FetchData(immune_cells, vars = c("Macrophagues_M1", "immuno", "sample"))

# Convertir filas en listas y combinarlas en una sola, eliminando NA
df$category <- ifelse(grepl("UV", df$sample), "KPB25L_UV", "KPB25L")

# Ordenar factor por categorías
df$category <- factor(df$category, levels = c("KPB25L", "KPB25L_UV"))

# Ordenar los samples dentro de cada categoría
df$sample <- factor(df$sample, levels = c(
  "KPB25L_control", "KPB25L_PD-L1", "KPB25L_DAC", "KPB25L_Combination",
  "KPB25L-UV_control", "KPB25L-UV_PD-L1", "KPB25L-UV_DAC", "KPB25L-UV_Combination"
))

df_violin_macrof <- subset(df, immuno == "Macrophages")

df_violin_macrof$sample <- factor(df_violin_macrof$sample, levels = c(
  "KPB25L_control", "KPB25L_PD-L1", "KPB25L_DAC", "KPB25L_Combination",
  "KPB25L-UV_control", "KPB25L-UV_PD-L1", "KPB25L-UV_DAC", "KPB25L-UV_Combination"
))
p3 <- ggplot(df_violin_macrof, aes(x = "", y = .data[["Macrophagues_M1"]], fill = immuno)) + 
  geom_violin() + 
  stat_summary(fun = median, geom = "point", color = "black", size = 2) +  # Marca la mediana
  facet_wrap(~ sample, ncol = 8) +  # Un gráfico por tipo de célula inmune
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),  # Oculta el eje X vacío
        axis.ticks.x = element_blank()) +
  labs(x = NULL, y = "Expresión",title = "Macrophagues_M1")  # Elimina el eje X vacío y etiqueta eje Y

print(p3)

ggsave(paste0("/home/workdir/HDAC_tfm/fig/nodownsampled/signatures_immune/Violin_plot_", "Macrophagues_M1", "_bytreat.png"), 
       plot = p3, width = 20, height = 12, dpi = 600)

shapiro.test(df_violin_macrof$Macrophagues_M1[df_violin_macrof$sample == "KPB25L-UV_Combination"])
shapiro.test(df_violin_macrof$Macrophagues_M1[df_violin_macrof$sample == "KPB25L-UV_PD-L1"])
shapiro.test(df_violin_macrof$Macrophagues_M1[df_violin_macrof$sample == "KPB25L-UV_DAC"])
shapiro.test(df_violin_macrof$Macrophagues_M1[df_violin_macrof$sample == "KPB25L-UV_control"])

shapiro.test(df_violin_macrof$Macrophagues_M1[df_violin_macrof$sample == "KPB25L_Combination"])
shapiro.test(df_violin_macrof$Macrophagues_M1[df_violin_macrof$sample == "KPB25L_PD-L1"])
shapiro.test(df_violin_macrof$Macrophagues_M1[df_violin_macrof$sample == "KPB25L_DAC"])
shapiro.test(df_violin_macrof$Macrophagues_M1[df_violin_macrof$sample == "KPB25L_control"])


leveneTest(Macrophagues_M1 ~ sample, data = df_violin_macrof)

kruskal.test(Macrophagues_M1 ~ sample, data = df_violin_macrof)
resultado_test <- pairwise.wilcox.test(df_violin_macrof$Macrophagues_M1, 
                                       df_violin_macrof$sample, 
                                       p.adjust.method = "bonferroni")
p_values <- resultado_test$p.value

df_significancias <- as.data.frame(as.table(p_values))
df_wider <- df_significancias %>%
  pivot_wider(names_from = Var2, values_from = Freq)

write.csv(df_wider, "/home/workdir/HDAC_tfm/fig/nodownsampled/signatures_immune/resultados_wilcoxon_M1.csv", row.names = FALSE)

p_val <- 2e-16  # Valor muy significativo

# Define la etiqueta de significancia
significance_label <- ifelse(p_val < 0.001, "***", 
                             ifelse(p_val < 0.01, "**", 
                                    ifelse(p_val < 0.05, "*", "ns")))


p3 <- ggplot(df_violin_macrof, aes(x = sample, y = Macrophagues_M1, fill = sample)) +
  geom_violin() +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  theme_minimal() +
  labs(title = "Distribution of M1 Gene Signature Expression Across Different Treatment Groups") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(df_violin_macrof, aes(x = sample, y = Macrophagues_M1, fill = sample)) +
  geom_violin() +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  scale_fill_manual(
    values = c(
      # KPB25L (tonos pastel)
      "KPB25L_control" = "#E0E0E0",        # Light gray
      "KPB25L_DAC" = "#FFB6B6",            # Light pink
      "KPB25L_PD-L1" = "#A8C8FF",          # Light blue
      "KPB25L_Combination" = "#C1A9E0",    # Light purple
      
      # KPB25L-UV (versiones más oscuras)
      "KPB25L-UV_control" = "#B0B0B0",     # Darker gray
      "KPB25L-UV_DAC" = "#E57373",         # Darker pink
      "KPB25L-UV_PD-L1" = "#648FFF",       # Darker blue
      "KPB25L-UV_Combination" = "#7E57C2"  # Darker purple
    )
  ) +
  theme_minimal() +
  labs(title = "Distribution of M1 Gene Signature Expression Across Different Treatment Groups") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0("/home/workdir/HDAC_tfm/fig/nodownsampled/signatures_immune/Violin_plot_", "Macrophagues_M1", "_bytreat.png"), 
       plot = p3, width = 20, height = 12, dpi = 600)

p3 <- ggplot(df_violin_macrof, aes(x = sample, y = Macrophagues_M1, fill = sample)) +
  geom_violin() +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  theme_minimal() +
  labs(title = "Distribution of M1 Gene Signature Expression Across Different Treatment Groups",
       x = "Sample", y = "M1 Macrophage Expression") + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Tamaño del texto del eje X
    axis.text.y = element_text(size = 10 ),  # Tamaño del texto del eje Y
    axis.title.x = element_text(size = 10, face = "bold"),  # Tamaño del título del eje X
    axis.title.y = element_text(size = 10, face = "bold"),  # Tamaño del título del eje Y
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  # Tamaño del título principal
    panel.grid = element_blank(),    
    panel.background = element_blank(),
    plot.background = element_blank()
  ) +
  geom_signif(comparisons = list(c("KPB25L_control", "KPB25L_DAC"), 
                                 c("KPB25L_control", "KPB25L_Combination"),
                                 c("KPB25L-UV_control", "KPB25L-UV_DAC"),
                                 c("KPB25L-UV_control", "KPB25L-UV_Combination")),
              annotations = c("***", "***", "***", "***"),
              y_position = c(max(df_violin_macrof$Macrophagues_M1) * 1, 
                             max(df_violin_macrof$Macrophagues_M1) * 1.15,
                             max(df_violin_macrof$Macrophagues_M1) * 1, 
                             max(df_violin_macrof$Macrophagues_M1) * 1.15), 
              tip_length = 0.02, 
              textsize = 5.5)  # Tamaño del texto de los asteriscos

ggsave(paste0("/home/workdir/HDAC_tfm/fig/nodownsampled/signatures_immune/Violin_plot_", "Macrophagues_M1", "_bytreat_withstats.png"), 
       plot = p3, width = 18, height = 14, dpi = 800)


# FINAL FIGURE ------------------------------------------------------------


p_final <- ggplot(df_violin_macrof, aes(x = sample, y = Macrophagues_M1, fill = sample)) +
  geom_violin() +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_signif(
    comparisons = list(
      c("KPB25L_control", "KPB25L_DAC"), 
      c("KPB25L_control", "KPB25L_Combination"),
      c("KPB25L-UV_control", "KPB25L-UV_DAC"),
      c("KPB25L-UV_control", "KPB25L-UV_Combination")
    ),
    annotations = c("***", "***", "***", "***"),
    y_position = c(
      max(df_violin_macrof$Macrophagues_M1) * 1.00, 
      max(df_violin_macrof$Macrophagues_M1) * 1.15,
      max(df_violin_macrof$Macrophagues_M1) * 1.00, 
      max(df_violin_macrof$Macrophagues_M1) * 1.15
    ),
    tip_length = 0.02,
    textsize = 5.5
  ) +
  scale_fill_manual(
    name = "Cell line and treatment",
    values = c(
      "KPB25L_control" = "#E0E0E0",
      "KPB25L_DAC" = "#FFB6B6",
      "KPB25L_PD-L1" = "#A8C8FF",
      "KPB25L_Combination" = "#C1A9E0",
      "KPB25L-UV_control" = "#B0B0B0",
      "KPB25L-UV_DAC" = "#E57373",
      "KPB25L-UV_PD-L1" = "#648FFF",
      "KPB25L-UV_Combination" = "#7E57C2"
    )
  ) +
  theme_minimal(base_family = "Times New Roman") +
  labs(
    title = "Distribution of M1 Gene Signature Expression Across Different Treatment Groups",
    x = "Sample",
    y = "M1 Signature Expression"
  ) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    
    # Tamaño de textos personalizados
    plot.title = element_text(size = 17, hjust = 0.5),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 13, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12)
  ) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 5, color = NA)))

ggsave(paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/Signature_Macrophagues_M1_bytreat_withstats.png"), 
       plot = p_final, width = 10, height = 8, dpi = 800)

