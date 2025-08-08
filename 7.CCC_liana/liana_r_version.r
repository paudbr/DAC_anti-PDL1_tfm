################################################################################
# Title: Liana ( R version)
# Author: Paula
# Description:
#   This script studies the Ligand Receptor interactions with Liana. 
#
################################################################################


library(Seurat)
library(SingleCellExperiment)
library(liana)
library(tidyverse)
library(OmnipathR)
library(magrittr)
library(biomaRt)
library(gprofiler2)

immune_cells <- readRDS("/home/workdir/analysis_augie/20250120_labeled_immune.combined_edit.RDS")
immune_cells[["RNA"]] <- JoinLayers(immune_cells[["RNA"]])
sce_object <- as.SingleCellExperiment(immune_cells)



# ORTHOLOGUE CONVERSION ---------------------------------------------------

mouse_genes <- rownames(sce_object@assays@data$logcounts)

conversion <- gorth(
  query = mouse_genes,
  source_organism = "mmusculus",
  target_organism = "hsapiens"
)
logcounts_mat <- assay(sce_object, "logcounts")
counts_mat <- assay(sce_object, "counts")

# Paso 2: Crear el vector de renombramiento
gene_map <- setNames(conversion$ortholog_name, conversion$input)

# Paso 3: Renombrar solo si hay mapeo
new_rownames <- gene_map[rownames(logcounts_mat)]

# Paso 4: Quitar genes sin mapeo
valid_rows <- !is.na(new_rownames)

# Filtrar matrices y nombres
logcounts_mat <- logcounts_mat[valid_rows, ]
counts_mat <- counts_mat[valid_rows, ]
new_rownames <- new_rownames[valid_rows]

# Paso 5: Identificar duplicados en new_rownames
non_duplicated <- !duplicated(new_rownames)

# Filtrar filas duplicadas en ambas matrices y renombrar
logcounts_mat <- logcounts_mat[non_duplicated, ]
counts_mat <- counts_mat[non_duplicated, ]
rownames(logcounts_mat) <- new_rownames[non_duplicated]
rownames(counts_mat) <- new_rownames[non_duplicated]

sce_object_human <- SingleCellExperiment(
  assays = list(logcounts = logcounts_mat),
  colData = colData(sce_object)
)
assay(sce_object_human, "counts") <- counts_mat

sce_object_human$immuno[sce_object_human$immuno == "Unannotated"] <- "Cancer cells"

library(zellkonverter)

# Guardar el objeto sce como .h5ad
writeH5AD(
  sce_object_human,
  file = "/home/workdir/HDAC_tfm/data/sce_object_human.h5ad",
  compression = "gzip"  # opcional: comprime el archivo
  

)

write.csv(as.matrix(assay(sce_object_human, "logcounts")), file = "/home/workdir/HDAC_tfm/data/anndata/logcounts.csv")
write.csv(as.matrix(assay(sce_object_human, "counts")), file = "/home/workdir/HDAC_tfm/data/anndata/counts.csv")

write.csv(as.data.frame(colData(sce_object_human)), "/home/workdir/HDAC_tfm/data/anndata/obs.csv")
write.csv(as.data.frame(rowData(sce_object_human)), "/home/workdir/HDAC_tfm/data/anndata/var.csv")


complex_test <- liana_wrap(
  sce_object_human,
  method = c('natmi', 'sca', 'logfc'),
  resource = 'CellPhoneDB',
  species = 'human',
  expr_prop = 0.1,
  idents_col = "immuno"
)

saveRDS(complex_test,"/home/workdir/HDAC_tfm/data/liana.RSD")

#### POR TRATAMIENTOS
all_resources <- show_resources()
all_methods <- show_methods()
tratamientos <- unique(sce_object_human$sample)

# Lista para guardar resultados por tratamiento
resultados_por_tratamiento <- list()

# Ejecutar liana_wrap por cada tratamiento
for (tt in tratamientos) {
  # Filtrar solo las células de ese tratamiento
  sce_sub <- sce_object_human[, sce_object_human$sample == tt]
  
  # Ejecutar LIANA solo en ese subconjunto
  resultados_por_tratamiento[[tt]] <- liana_wrap(
    sce_sub,
    method = all_methods,
    resource ="Consensus",
    species = 'human',
    expr_prop = 0.1,
    idents_col = "immuno"
  )
}

for (tt in names(resultados_por_tratamiento)) {
    resultados_por_tratamiento[[tt]] <- liana_aggregate(
      resultados_por_tratamiento[[tt]],
      resource = "CellPhoneDB"
    )
}

dplyr::glimpse(resultados_por_tratamiento)


resultados_combinados <- bind_rows(
  lapply(names(resultados_por_tratamiento), function(tt) {
    df <- resultados_por_tratamiento[[tt]]
    df$tratamiento <- tt   # Añade columna con el nombre del tratamiento
    return(df)
  })
)

resultados_separados <- resultados_combinados %>%
  separate(tratamiento, into = c("linea_celular", "tratamiento"), sep = "_", extra = "merge")
### CANCER CELLS AS SOURCE
resultados_kap25l_souce_cancer <- resultados_separados  %>%
  filter(linea_celular == "KAP25L") %>% filter(source == "Cancer cells")
top15_por_tratamiento <- resultados_kap25l_souce_cancer %>%
  group_by(tratamiento) %>%
  arrange(aggregate_rank) %>%
  slice_head(n = 15) %>%
  ungroup()

# Paso 2: Crear columna para la interacción en eje Y (ligand.complex - receptor.complex)
top15_por_tratamiento <- top15_por_tratamiento %>%
  mutate(interaccion_y = paste(ligand.complex, receptor.complex, sep = " - "))

# Paso 3: Reordenar factores para eje Y según aggregate_rank
top15_por_tratamiento <- top15_por_tratamiento %>%
  group_by(tratamiento) %>%
  mutate(interaccion_y = fct_reorder(interaccion_y, aggregate_rank)) %>%
  ungroup()

# Paso 4: Para eje X dejamos el target tal cual (podrías ordenar si quieres)
top15_por_tratamiento <- top15_por_tratamiento %>%
  mutate(target = factor(target))  # opcional: fct_inorder o fct_reorder si quieres otro orden

# Paso 5: Plot con ggplot + facetas por tratamiento, con subtítulo fijo
p <- ggplot(top15_por_tratamiento, aes(x = target, y = interaccion_y,
                                       color = sca.LRscore, size = natmi.edge_specificity)) +
  geom_point(alpha = 0.8) +
  scale_color_viridis_c(option = "plasma", name = "Expression (sca.LRscore)") +
  scale_size(range = c(1, 6), name = "Specificity (edge_specificity)") +
  facet_wrap(~ tratamiento, scales = "free") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 75, hjust = 1, size = 7),
    axis.text.y = element_text(size = 9),
    strip.text = element_text(face = "bold", size = 10)
  ) +
  labs(x = "Target",
       y = "Interacción (Ligand complex - Receptor complex)",
       title = "Top 15 Interacciones en KAP25L por tratamiento",
       subtitle = "Source: Cancer cells",
       color = "Expresión",
       size = "Especificidad")

# Mostrar plot
print(p)

ggsave("/home/workdir/HDAC_tfm/fig/liana/Dotplots_results_KAP25L.png", plot = p, height=15, width=12, dpi=600)


### CANCER CELLS AS TARGET
resultados_kap25l_target_cancer <- resultados_separados  %>%
  filter(linea_celular == "KAP25L") %>% filter(target == "Cancer cells")
top15_por_tratamiento <- resultados_kap25l_target_cancer %>%
  group_by(tratamiento) %>%
  arrange(aggregate_rank) %>%
  slice_head(n = 15) %>%
  ungroup()

# Paso 2: Crear columna para la interacción en eje Y (ligand.complex - receptor.complex)
top15_por_tratamiento <- top15_por_tratamiento %>%
  mutate(interaccion_y = paste(ligand.complex, receptor.complex, sep = " - "))

# Paso 3: Reordenar factores para eje Y según aggregate_rank
top15_por_tratamiento <- top15_por_tratamiento %>%
  group_by(tratamiento) %>%
  mutate(interaccion_y = fct_reorder(interaccion_y, aggregate_rank)) %>%
  ungroup()

# Paso 4: Para eje X dejamos el target tal cual (podrías ordenar si quieres)
top15_por_tratamiento <- top15_por_tratamiento %>%
  mutate(source = factor(source))  # opcional: fct_inorder o fct_reorder si quieres otro orden

# Paso 5: Plot con ggplot + facetas por tratamiento, con subtítulo fijo
p <- ggplot(top15_por_tratamiento, aes(x = target, y = interaccion_y,
                                       color = sca.LRscore, size = natmi.edge_specificity)) +
  geom_point(alpha = 0.8) +
  scale_color_viridis_c(option = "plasma", name = "Expression (sca.LRscore)") +
  scale_size(range = c(1, 6), name = "Specificity (edge_specificity)") +
  facet_wrap(~ tratamiento, scales = "free") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 75, hjust = 1, size = 7),
    axis.text.y = element_text(size = 9),
    strip.text = element_text(face = "bold", size = 10)
  ) +
  labs(x = "Target",
       y = "Interacción (Ligand complex - Receptor complex)",
       title = "Top 15 Interacciones en KAP25L por tratamiento",
       subtitle = "Target: Cancer cells",
       color = "Expresión",
       size = "Especificidad")

# Mostrar plot
print(p)

ggsave("/home/workdir/HDAC_tfm/fig/liana/Dotplots_results_KAP25L_cancer_as_target.png", plot = p, height=15, width=12, dpi=600)


### HACIENDO TODAS LAS COMBINACIONES:

# Función auxiliar para capitalizar
source("/home/workdir/HDAC_tfm/scripts/liana-py/generar_plot.r")
ucfirst <- function(string) {
  paste0(toupper(substr(string, 1, 1)), substr(string, 2, nchar(string)))
}

# Ejecutar para todas las combinaciones posibles
tipos_celulares <- unique(resultados_separados$source)

for (tipo in tipos_celulares) {
  # source = tipo, target = resto
  generar_dotplot(resultados_separados, linea_celular = "KAP25L-UV", cell_type = tipo, modo = "source", output_dir = "/home/workdir/HDAC_tfm/fig/liana/automat/")
  
  # target = tipo, source = resto
  generar_dotplot(resultados_separados, linea_celular = "KAP25L-UV", cell_type = tipo, modo = "target", output_dir = "/home/workdir/HDAC_tfm/fig/liana/automat/")
}