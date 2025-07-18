library(ggplot2)
library(tidyverse)
library(dplyr)
#install.packages("ggVennDiagram")
library(ggVennDiagram)
library(patchwork)

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

get_genes_from_files <- function(files_list, pattern_line_cell) {
  genes_all <- c()
  for (f in files_list) {
    if (grepl(pattern_line_cell, f)) {
      df <- read.csv(f, row.names = 1)  # lee con rownames
      genes_all <- c(genes_all, rownames(df))
    }
  }
  unique(genes_all)
}

# Para KPB25L (cáncer)
cancer_k25l_genes <- get_genes_from_files(deg_files_cancer, "kpb25l")
unann_k25l_genes  <- get_genes_from_files(deg_files_unann, "kpb25l")

# Para KPB25L-UV (unannotated)
cancer_kuv_genes <- get_genes_from_files(deg_files_cancer, "kuv")
unann_kuv_genes  <- get_genes_from_files(deg_files_unann, "kuv")


# Preparar listas de genes
gene_sets_k25l <- list(
  Cancer = cancer_k25l_genes,
  Unann = unann_k25l_genes
)

gene_sets_kuv <- list(
  Cancer = cancer_kuv_genes,
  Unann = unann_kuv_genes
)

base_theme <- theme(
  text = element_text(family = "Times New Roman"),
  plot.title = element_text(size = 16, face = "bold", hjust = 0.2),
  legend.title = element_text(size = 12),
  legend.text = element_text(),
  plot.margin = margin(8, 5, 5, 5),
  legend.position = "right"
)

set_venn_text_family <- function(plot_obj, font_family = "Times New Roman") {
  # Identificar las capas de texto dentro del plot (los números de intersección)
  plot_obj$layers <- lapply(plot_obj$layers, function(layer) {
    if ("GeomText" %in% class(layer$geom)) {
      layer$aes_params$family <- font_family
      layer$geom_params$family <- font_family
    }
    layer
  })
  return(plot_obj)
}

p1 <- ggVennDiagram(gene_sets_k25l, label_alpha = 0) +
  scale_fill_gradient(low = "skyblue", high = "#D8BFD8") +
  labs(title = "DEGs shared: Cancer vs Unann (KPB25L)") +
  base_theme +
  coord_fixed(clip = "off")

p2 <- ggVennDiagram(gene_sets_kuv, label_alpha = 0) +
  scale_fill_gradient(low = "skyblue", high = "#D8BFD8") +
  labs(title = "DEGs shared: Cancer vs Unann (KPB25L-UV)") +
  base_theme +
  coord_fixed(clip = "off")

p1 <- set_venn_text_family(p1, "Times New Roman")
p2 <- set_venn_text_family(p2, "Times New Roman")

combined_plot <- p1 + p2 + plot_layout(ncol = 2, widths = 3)
print(combined_plot)

ggsave("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/VennDiag_cancer_unan.png", plot=combined_plot, height = 6, width = 12, dpi=600)
