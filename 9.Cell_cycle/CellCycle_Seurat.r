################################################################################
# Title: Cell cycle analysis
# Author: Paula
# Description:
#   Using Seurat CellCycleScoring() function to know in which cell cycle phase is
#   each cell
#
################################################################################

library(Seurat)
library(OmnipathR)
library(magrittr)
library(biomaRt)
library(gprofiler2)

tumor_cells <- readRDS("/home/workdir/HDAC_tfm/data/tumor_cells_adjusted_13032025.rds")

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

tumor_cells_copy <- JoinLayers(tumor_cells)

conversion_s_genes <- gorth(
  query = s.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)

conversion_g2m_genes <- gorth(
  query = g2m.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)


gm2_mouse <- conversion_g2m_genes$ortholog_name
s_mouse <- conversion_s_genes$ortholog_name


tumor_cells_copy <- CellCycleScoring(
  object = tumor_cells_copy,
  s.features = s_mouse,
  g2m.features = gm2_mouse,
  set.ident = TRUE
)

saveRDS(tumor_cells_copy, "/home/workdir/HDAC_tfm/data/tumor_cells_adjusted_13032025_with_cellcycle.rds")

tumor_cells_copy <- readRDS("/home/workdir/HDAC_tfm/data/tumor_cells_adjusted_13032025_with_cellcycle.rds")

df <- tumor_cells_copy@meta.data

# Separar la columna 'sample' en dos columnas nuevas 'cell_line' y 'treat'
df <- df %>%
  separate(col = sample, into = c("cell_line", "treat"), sep = "_", remove = FALSE)

# Actualizamos el meta.data con las nuevas columnas
tumor_cells_copy <- AddMetaData(tumor_cells_copy, metadata = df[, c("cell_line", "treat")])

df <- tumor_cells_copy@meta.data



df$cell_line <- as.character(df$cell_line)  # Asegurar que es carácter

df$cell_line[df$cell_line == "KAP25L"] <- "KPB25L"
df$cell_line[df$cell_line == "KAP25L-UV"] <- "KPB25L-UV"


tumor_cells_copy <- AddMetaData(tumor_cells_copy, metadata = df["cell_line"])
df$treat <- factor(df$treat, levels = c("control", "PD-L1", "DAC", "Combination"))

p <- ggplot(df, aes(x = Phase, fill = Phase)) +
  geom_bar(position = "fill") +   # Proporciones dentro de cada grupo
  ylab("Proporción de células") +
  xlab("Fase del ciclo celular") +
  ggtitle("Distribución del ciclo celular por línea celular y tratamiento") +
  scale_fill_brewer(palette = "Set2") +
  facet_grid(cell_line ~ treat) +  # filas = cell_line, columnas = treat
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p)


p3<- ggplot(df, aes(x = S.Score, fill = Phase)) +
  geom_histogram(binwidth = 0.05, position = "identity", alpha = 0.6) +
  facet_grid(cell_line ~ treat) +
  xlab("S.Score") +
  ylab("Number of cells") +
  ggtitle("Histogram of S.Score for each Cell Cycle Phase in tumor cells by treatment and cell-line") +
  theme_bw()
ggsave(filename = "/home/workdir/HDAC_tfm/fig/CellCycle/Histogram_Sscore_cellcycle_phases.png",
       plot = p3,
       height = 6, width = 12, dpi = 600)

df_bar <- df %>%
  group_by(cell_line, treat, Phase) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cell_line, treat) %>%
  mutate(proportion = n / sum(n))

p3 <- ggplot(df_bar, aes(x = treat, y = proportion, fill = Phase)) +
  geom_bar(stat = "identity") +
  facet_wrap(~cell_line) +
  scale_fill_manual(
    values = c("G1" = "#ff8000", "S" = "#ffd5b1", "G2M" = "#4daf4a"),
    name = "Cell Cycle Phase"  # Aquí defines el título de la leyenda
  ) +
  ylab("Proportion of cells") +
  xlab("Treatment") +
  ggtitle("Cell cycle phase proportions in cancer cells by treatment and cell line") +
  theme_classic2() +
  theme(
    text = element_text(family = "Times New Roman", size = 18),
    plot.title = element_text(family = "Times New Roman", size = 20, face = "bold")
  )

ggsave(
  filename = "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/Barplot_cellcycle_phases_bytreat.png",
  plot = p3,
  height = 6, width = 9, dpi = 600
)


df_summary <- df %>%
  group_by(cell_line, treat, Phase) %>%
  summarise(
    percent_expressing = mean(S.Score > 0) * 100,
    avg_S_score = mean(S.Score[S.Score > 0], na.rm = TRUE),
    .groups = "drop"
  )

df <- df %>%
  group_by(cell_line) %>%
  mutate(S.Score.z = scale(S.Score)) %>%
  ungroup()
df_summary <- df %>%
  filter(S.Score > 0) %>%
  group_by(cell_line, treat, Phase) %>%
  summarise(
    percent_expressing = n() / nrow(df[df$cell_line == unique(cell_line) & df$treat == unique(treat) & df$Phase == unique(Phase), ]) * 100,
    mean_z_S_score = mean(S.Score.z, na.rm = TRUE),
    .groups = "drop"
  )

p3<- ggplot(df_summary, aes(x = treat, y = Phase)) +
  geom_point(aes(size = percent_expressing, color = mean_z_S_score)) +
  facet_wrap(~cell_line) +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  labs(
    title = "DotPlot of S.Score expression in tumor cells by treatment and cell line",
    x = "Treatment", y = "Cell Cycle Phase",
    size = "% cells expressing S.Score > 0",
    color = "Average Z-S.Score"
  ) +
  theme_bw()

ggsave(filename="/home/workdir/HDAC_tfm/fig/CellCycle/Dotplot_S_score.png",plot=p3, height =4, width = 8, dpi=600)


df_summary <- df %>%
  group_by(cell_line, treat, Phase) %>%
  summarise(
    percent_expressing = mean(G2M.Score > 0) * 100,
    avg_S_score = mean(S.Score[G2M.Score > 0], na.rm = TRUE),
    .groups = "drop"
  )

df <- df %>%
  group_by(cell_line) %>%
  mutate(S.Score.z = scale(G2M.Score)) %>%
  ungroup()
df_summary <- df %>%
  filter(G2M.Score > 0) %>%
  group_by(cell_line, treat, Phase) %>%
  summarise(
    percent_expressing = n() / nrow(df[df$cell_line == unique(cell_line) & df$treat == unique(treat) & df$Phase == unique(Phase), ]) * 100,
    mean_z_G2M_score = mean(S.Score.z, na.rm = TRUE),
    .groups = "drop"
  )

p3<- ggplot(df_summary, aes(x = treat, y = Phase)) +
  geom_point(aes(size = percent_expressing, color = mean_z_G2M_score)) +
  facet_wrap(~cell_line) +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  labs(
    title = "DotPlot of G2M.Score expression in tumor cells by treatment and cell line",
    x = "Treatment", y = "Cell Cycle Phase",
    size = "% cells expressing G2M.Score > 0",
    color = "Average Z-G2M.Score"
  ) +
  theme_bw()

ggsave(filename="/home/workdir/HDAC_tfm/fig/CellCycle/Dotplot_G2M_score.png",plot=p3, height =4, width = 8, dpi=600)


p3 <- ggplot(df, aes(x = G2M.Score, fill = Phase)) +
  geom_histogram(binwidth = 0.05, position = "identity", alpha = 0.6) +
  facet_grid(cell_line ~ treat) +
  xlab("S.Score") +
  ylab("Number of cells") +
  ggtitle("Histogram of G2M.Score for each Cell Cycle Phase in tumor cells by treatment and cell-line") +
  theme_bw()
ggsave(filename = "/home/workdir/HDAC_tfm/fig/CellCycle/Histogram_G2Mscore_cellcycle_phases.png",
       plot = p3,
       height = 6, width = 12, dpi = 600)


DimPlot(tumor_cells_copy@meta.data$RNA_snn_res.0.05, group.by = "Phase", pt.size = 0.5) + 
  scale_color_manual(values = c("G1" = "blue", "S" = "green", "G2M" = "red")) +
  ggtitle("UMAP de células tumorales por fase del ciclo celular") +
  theme_minimal() +
  theme(legend.position = "right")


p3 <- UMAPPlot(tumor_cells_copy, group.by = "RNA_snn_res.0.05", split.by = "cell_line", ncol = 4) +
  theme_bw() +
  ggtitle("Cancer cells subclusters by treatment and cell line")

ggsave(filename = "/home/workdir/HDAC_tfm/fig/CellCycle/UMAP_Cancer_cells_subclusters_by_treatment.png",
       plot = p3,
       height = 6, width = 12, dpi = 600)


DimPlot(tumor_cells_copy, 
        group.by = "Phase", 
        split.by = "cell_line", 
        ncol = 4, 
        pt.size = 0.5) +
  scale_color_manual(values = c("G1" = "blue", "S" = "green", "G2M" = "red")) +
  ggtitle("UMAP por fase del ciclo celular separado por sample") +
  theme_minimal()



library(patchwork)

# Plot ajustado de clusters
p1 <- DimPlot(tumor_cells_copy, group.by = "RNA_snn_res.0.05") + ggtitle("Adjusted clusters")

tumor_cells_copy$group <- paste(tumor_cells_copy$cell_line, tumor_cells_copy$treat, sep = "_")

tumor_cells_copy$cell_line <- factor(tumor_cells_copy$cell_line, levels = c("KPB25L", "KPB25L-UV"))
tumor_cells_copy$treat <- factor(tumor_cells_copy$treat, levels = c("control", "PD-L1", "DAC", "Combination"))

# Crear la columna group combinando cell_line y treat
tumor_cells_copy$group <- factor(
  paste(tumor_cells_copy$cell_line, tumor_cells_copy$treat, sep = "_"),
  levels = c(
    "KPB25L_control", "KPB25L_PD-L1", "KPB25L_DAC", "KPB25L_Combination",
    "KPB25L-UV_control", "KPB25L-UV_PD-L1", "KPB25L-UV_DAC", "KPB25L-UV_Combination"
  )
)
p2 <- DimPlot(tumor_cells_copy, group.by = "Phase", split.by = "cell_line", ncol = 4, pt.size = 0.5) + 
  scale_color_manual(values = c(
    "G1" = "#ff8000",
    "S" = "#ffd5b1",
    "G2M" = "#4daf4a"
  )) +
  ggtitle("Cell cycle phases by treatment and cell line") + 
  theme_bw()
print(p1)

ggsave(filename = "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/DimPlot_Cell_cycle_phases_by_cell_line.png",
       plot = p2,
       height = 6, width = 12, dpi = 600)

# Plot fases dividido por sample
p2 <- DimPlot(tumor_cells_copy, group.by = "Phase", split.by = "cell_line~treat", ncol = 4) + 
  scale_color_manual(values = c("G1" = "pink", "S" = "green", "G2M" = "red")) +
  ggtitle("Fases del ciclo celular") + theme_bw()

# Combinar con etiquetas arriba
p1 / p2 + plot_annotation(tag_levels = "A")


p1 <- DimPlot(tumor_cells_copy, group.by = "Phase", split.by = "cell_line", ncol = 2, pt.size = 0.5) + 
  scale_color_manual(values = c(
    "G1" = "#ff8000",
    "S" = "#ffd5b1",
    "G2M" = "#4daf4a"
  )) +
  ggtitle("Cell cycle phases in cancer cells by cell line") + 
  theme_classic() +
  theme(
    text = element_text(family = "Times New Roman", size = 18),
    plot.title = element_text(size = 20, face = "bold")
  )
  
# Asegurar niveles de factores
tumor_cells_copy$cell_line <- factor(tumor_cells_copy$cell_line, levels = c("KPB25L", "KPB25L-UV"))
tumor_cells_copy$treat <- factor(tumor_cells_copy$treat, levels = c("control", "PD-L1", "DAC", "Combination"))

# Crear columna group con niveles ordenados
tumor_cells_copy$group <- factor(
  paste(tumor_cells_copy$cell_line, tumor_cells_copy$treat, sep = "_"),
  levels = c(
    "KPB25L_control", "KPB25L_PD-L1", "KPB25L_DAC", "KPB25L_Combination",
    "KPB25L-UV_control", "KPB25L-UV_PD-L1", "KPB25L-UV_DAC", "KPB25L-UV_Combination"
  )
)

# Gráfico 2: Fases del ciclo celular por línea
p2 <- DimPlot(tumor_cells_copy, group.by = "Phase", split.by = "group", ncol = 4, pt.size = 0.5) + 
  scale_color_manual(values = c("G1" = "#ff8000", "S" = "#ffd5b1", "G2M" = "#4daf4a")) +
  ggtitle("Cell cycle phases in cancer cells by treatment and cell line") +
  theme_classic() +
  theme(
    text = element_text(family = "Times New Roman", size = 18),
    plot.title = element_text(size = 20, face = "bold")
  )

final_plot <- p3 + p1 / p2 + 
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      text = element_text(family = "Times New Roman"),
      plot.tag = element_text(family = "Times New Roman", size = 24, face = "bold")
    )
  )
# Guardar el gráfico combinado
ggsave(
  filename = "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/DimPlot_combined_cell_cycle_times_new_roman_2.png",
  plot = final_plot,
  height = 12, width = 25, dpi = 700
)



# RESPONSIVE_UNRESPONSIVE -------------------------------------------------

cells_responsive_intersect_tfs_uv <- readRDS( "/home/workdir/HDAC_tfm/data/cells_responsive_intersect_tfs_kpb25l_uv.rds")
cells_unresponsive_intersect_tfs_uv <- readRDS("/home/workdir/HDAC_tfm/data/cells_unresponsive_intersect_tfs_kpb25l_uv.rds")

cells_responsive_intersect_tfs <- readRDS("/home/workdir/HDAC_tfm/data/cells_responsive_intersect_tfs_kpb25l.rds")
cells_unresponsive_intersect_tfs<- readRDS("/home/workdir/HDAC_tfm/data/cells_unresponsive_intersect_tfs_kpb25l.rds")


tumor_responsive_unresponsive_kpb25l <- subset(tumor_cells_copy, cells = c(cells_responsive_intersect_tfs,cells_unresponsive_intersect_tfs ))


tumor_responsive_unresponsive_kpb25l_uv <- subset(tumor_cells_copy, cells = c(cells_responsive_intersect_tfs_uv,cells_unresponsive_intersect_tfs_uv ))


labels_kpb25l <- ifelse(Cells(tumor_responsive_unresponsive_kpb25l) %in% cells_responsive_intersect_tfs,
                        "responsive", "unresponsive")

tumor_responsive_unresponsive_kpb25l$TF_response <- labels_kpb25l

labels_kpb25l_uv <- ifelse(Cells(tumor_responsive_unresponsive_kpb25l_uv) %in% cells_responsive_intersect_tfs_uv,
                           "responsive", "unresponsive")

tumor_responsive_unresponsive_kpb25l_uv$TF_response <- labels_kpb25l_uv


#KPB25L

tumor_responsive_unresponsive_kpb25l <- CellCycleScoring(
  object = tumor_responsive_unresponsive_kpb25l,
  s.features = s_mouse,
  g2m.features = gm2_mouse,
  set.ident = TRUE
)

df_unresponsive_responsive <- tumor_responsive_unresponsive_kpb25l@meta.data

p1 <- ggplot(df_unresponsive_responsive, aes(x = S.Score, fill = Phase)) +
  geom_histogram(binwidth = 0.05, position = "identity", alpha = 0.6) +
  facet_grid(~TF_response) +
  xlab("S.Score") +
  ylab("Number of cells") +
  ggtitle("Histogram of S.Score for each Cell Cycle Phase in responsive vs unresponsive cells in KPB25L") +
  theme_bw()

ggsave(filename = "/home/workdir/HDAC_tfm/fig/CellCycle/Histogram_Sscore_cellcycle_phases_unresponsive_responsive_kpb25l.png",
       plot = p1,
       height = 6, width = 12, dpi = 600)


df_bar_tf <- df_unresponsive_responsive %>%
  group_by(TF_response, Phase) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(TF_response) %>%
  mutate(proportion = n / sum(n))

# Plot
p_tf <- ggplot(df_bar_tf, aes(x = TF_response, y = proportion, fill = Phase)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("G1" = "#f781bf", "S" = "#377eb8", "G2M" = "#4daf4a")) +
  ylab("Proportion of cells") +
  xlab("TF_response group") +
  ggtitle("Cell Cycle Phase Proportions in Responsive vs Unresponsive Tumor Cells (KPB25L)") +
  theme_bw() +
  theme(legend.title = element_blank())

ggsave(filename = "/home/workdir/HDAC_tfm/fig/CellCycle/Barplot_unresponsive_responsive_kpb25l.png",
       plot = p_tf,
       height = 6, width = 12, dpi = 600)



p1 <- ggplot(df_unresponsive_responsive, aes(x = G2M.Score, fill = Phase)) +
  geom_histogram(binwidth = 0.05, position = "identity", alpha = 0.6) +
  facet_grid(~TF_response) +
  xlab("S.Score") +
  ylab("Number of cells") +
  ggtitle("Histogram of G2M.Score for each Cell Cycle Phase in responsive vs unresponsive cells in KPB25L") +
  theme_bw()

ggsave(filename = "/home/workdir/HDAC_tfm/fig/CellCycle/Histogram_G2Mscore_cellcycle_phases_unresponsive_responsive_kpb25l.png",
       plot = p1,
       height = 6, width = 12, dpi = 600)


umap_all <- Embeddings(tumor_cells_copy, reduction = "umap")
meta_all <- tumor_cells_copy@meta.data
meta_all$UMAP_1 <- umap_all[, 1]
meta_all$UMAP_2 <- umap_all[, 2]

# Definir TF_response
meta_all$TF_response <- NA
meta_all$TF_response[rownames(meta_all) %in% cells_responsive_intersect_tfs] <- "responsive"
meta_all$TF_response[rownames(meta_all) %in% cells_unresponsive_intersect_tfs] <- "unresponsive"
meta_all$TF_response[is.na(meta_all$TF_response)] <- "background"

# Fijar niveles y formas
meta_all$TF_response <- factor(meta_all$TF_response, levels = c("background", "responsive", "unresponsive"))

# Colores por fase del ciclo
phase_colors <- c("G1" = "#ff8000", "S" = "#ffd5b1", "G2M" = "#4daf4a")  # verde

# Formas por TF_response
shape_values <- c("responsive" = 16,    # círculo
                  "unresponsive" = 17)  # triángulo

# Plot con fondo gris, color por Phase, shape por TF_response
library(ggplot2)

p3 <- ggplot(meta_all, aes(x = UMAP_1, y = UMAP_2)) +
  # Fondo: todas las células en gris claro, sin distinción
  geom_point(data = subset(meta_all, TF_response == "background"), 
             color = "lightgrey", size = 0.3, alpha = 0.5) +
  
  # Encima: responsive/unresponsive con color por Phase y shape diferente
  geom_point(data = subset(meta_all, TF_response != "background"),
             aes(color = Phase, shape = TF_response), size = 1.2, alpha = 0.9) +
  
  scale_color_manual(values = phase_colors) +
  scale_shape_manual(values = shape_values) +
  
  ggtitle("UMAP: Cancer cells responsive vs unresponsive for KPB25L") +
  theme_bw() +
  theme(legend.title = element_blank())

ggsave(filename = "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/UMAP_byphase_unresponsive_responsive_kpb25l.png",
       plot = p3,
       height = 6, width = 12, dpi = 600)


#KPB25L-UV


tumor_responsive_unresponsive_kpb25l <- CellCycleScoring(
  object = tumor_responsive_unresponsive_kpb25l_uv,
  s.features = s_mouse,
  g2m.features = gm2_mouse,
  set.ident = TRUE
)

df_unresponsive_responsive <- tumor_responsive_unresponsive_kpb25l_uv@meta.data

p1 <- ggplot(df_unresponsive_responsive, aes(x = S.Score, fill = Phase)) +
  geom_histogram(binwidth = 0.05, position = "identity", alpha = 0.6) +
  facet_grid(~TF_response) +
  xlab("S.Score") +
  ylab("Number of cells") +
  ggtitle("Histogram of S.Score for each Cell Cycle Phase in responsive vs unresponsive cells in KPB25L-UV") +
  theme_bw()

ggsave(filename = "/home/workdir/HDAC_tfm/fig/CellCycle/Histogram_Sscore_cellcycle_phases_unresponsive_responsive_kpb25_uv.png",
       plot = p1,
       height = 6, width = 12, dpi = 600)

df_bar_tf <- df_unresponsive_responsive %>%
  group_by(TF_response, Phase) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(TF_response) %>%
  mutate(proportion = n / sum(n))

# Plot
p_tf <- ggplot(df_bar_tf, aes(x = TF_response, y = proportion, fill = Phase)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("G1" = "#f781bf", "S" = "#377eb8", "G2M" = "#4daf4a")) +
  ylab("Proportion of cells") +
  xlab("TF_response group") +
  ggtitle("Cell Cycle Phase Proportions in Responsive vs Unresponsive Tumor Cells (KPB25L-UV)") +
  theme_bw() +
  theme(legend.title = element_blank())

ggsave(filename = "/home/workdir/HDAC_tfm/fig/CellCycle/Barplot_unresponsive_responsive_kpb25l_uv.png",
       plot = p_tf,
       height = 6, width = 12, dpi = 600)




p1 <- ggplot(df_unresponsive_responsive, aes(x = G2M.Score, fill = Phase)) +
  geom_histogram(binwidth = 0.05, position = "identity", alpha = 0.6) +
  facet_grid(~TF_response) +
  xlab("S.Score") +
  ylab("Number of cells") +
  ggtitle("Histogram of G2M.Score for each Cell Cycle Phase in responsive vs unresponsive cells in KPB25L-UV") +
  theme_bw()

ggsave(filename = "/home/workdir/HDAC_tfm/fig/CellCycle/Histogram_G2Mscore_cellcycle_phases_unresponsive_responsive_kpb25_uv.png",
       plot = p1,
       height = 6, width = 12, dpi = 600)


umap_all <- Embeddings(tumor_cells_copy, reduction = "umap")
meta_all <- tumor_cells_copy@meta.data
meta_all$UMAP_1 <- umap_all[, 1]
meta_all$UMAP_2 <- umap_all[, 2]

# Definir TF_response
meta_all$TF_response <- NA
meta_all$TF_response[rownames(meta_all) %in% cells_responsive_intersect_tfs_uv] <- "responsive"
meta_all$TF_response[rownames(meta_all) %in% cells_unresponsive_intersect_tfs_uv] <- "unresponsive"
meta_all$TF_response[is.na(meta_all$TF_response)] <- "background"

# Fijar niveles y formas
meta_all$TF_response <- factor(meta_all$TF_response, levels = c("background", "responsive", "unresponsive"))

# Colores por fase del ciclo
phase_colors <- c("G1" = "#f781bf",   # rosa
                  "S" = "#377eb8",    # azul
                  "G2M" = "#4daf4a")  # verde

# Formas por TF_response
shape_values <- c("responsive" = 16,    # círculo
                  "unresponsive" = 17)  # triángulo

# Plot con fondo gris, color por Phase, shape por TF_response
library(ggplot2)

p3 <- ggplot(meta_all, aes(x = UMAP_1, y = UMAP_2)) +
  # Fondo: todas las células en gris claro, sin distinción
  geom_point(data = subset(meta_all, TF_response == "background"), 
             color = "lightgrey", size = 0.3, alpha = 0.5) +
  
  # Encima: responsive/unresponsive con color por Phase y shape diferente
  geom_point(data = subset(meta_all, TF_response != "background"),
             aes(color = Phase, shape = TF_response), size = 1.2, alpha = 0.9) +
  
  scale_color_manual(values = phase_colors) +
  scale_shape_manual(values = shape_values) +
  
  ggtitle("UMAP: Cancer cells responsive vs unresponsive for KPB25L-UV") +
  theme_bw() +
  theme(legend.title = element_blank())

ggsave(filename = "/home/workdir/HDAC_tfm/fig/CellCycle/UMAP_byphase_unresponsive_responsive_kpb25_uv.png",
       plot = p3,
       height = 6, width = 12, dpi = 600)


