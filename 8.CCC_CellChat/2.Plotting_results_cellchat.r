library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(ggraph)
library(igraph)
celltype_colors <- c(
  "Cancer cells"                 = "#FFB3BA",  # pastel coral rosado
  "Endothelial"                 = "#AFCBFF",  # azul pastel claro
  "Macrophages"                 = "#B5EAD7",  # verde menta suave
  "Naive CD8+ T cells"          = "#CBAACB",  # lila suave
  "Natural killer  cells"       = "#FFDAC1",  # melocotón claro
  "Plasmacytoid Dendritic cells"= "#E2F0CB",  # verde lima pastel
  "Pre-B cells"                 = "#D5AAFF",  # lavanda
  "Pro-B cells"                 = "#FFFFD1"   # amarillo suave
)


# KAP25L ------------------------------------------------------------------


# #CONTROL ----------------------------------------------------------------



cellchat_k25l_control<- readRDS("/home/workdir/HDAC_tfm/data/cellchat/cellchat_SS_30_05_KAP25L_control.rds")
cellchat <- cellchat_k25l_control
unique_ligands_sig_control <- unique(unlist(cellchat@LR$LRsig$ligand))
ptm = Sys.time()
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Primer gráfico: número de interacciones
png("/home/workdir/HDAC_tfm/fig/CellChat/KAP25L/control/circle_interactions_count_control.png", 
    width = 2000, height = 2000, res = 300)

netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  color.use = celltype_colors,
  title.name = "Number of interactions"
)

dev.off()

# Segundo gráfico: fuerza de interacción
png("/home/workdir/HDAC_tfm/fig/CellChat/KAP25L/control/circle_interactions_weight_control.png", 
    width = 2000, height = 2000, res = 300)

netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  color.use = celltype_colors,
  title.name = "Interaction weights/strength"
)

dev.off()


pdf("/home/workdir/HDAC_tfm/fig/CellChat/KAP25L/circle_per_signal_control.pdf", 
    width = 16, height = 12)  # Puedes ajustar tamaño según lo necesites

mat <- cellchat@net$weight

par(mfrow = c(4, 2), xpd = TRUE, mar = c(1, 1, 2, 1))  # Márgenes para que quepa el título

for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(
    mat2,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    edge.weight.max = max(mat),
    title.name = rownames(mat)[i]
  )
}

dev.off()


# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0("/home/workdir/HDAC_tfm/fig/CellChat/KAP25L/",pathways.show.all[i], "_L-R_contribution_control.pdf"), plot=gg, width = 10, height = 12 , dpi = 600)
}

cells_source <- levels(cellchat@idents)



for (i in 1:8) {
gg <- netVisual_bubble(cellchat, sources.use = i, targets.use = c(1:8), remove.isolate = FALSE)
ggsave(filename=paste0("/home/workdir/HDAC_tfm/fig/CellChat/KAP25L/control",cells_source[i], "source_L-R_prob_control.png"), plot=gg, width = 10, height = 12 , dpi = 600)

}


netVisual_chord_gene(cellchat, sources.use = 1, targets.use = c(1:8), lab.cex = 0.5,legend.pos.y = 30, color.use = celltype_colors)

netVisual_chord_gene(
  cellchat,
  sources.use = 1,
  targets.use = c(1:8),
  lab.cex = 0.1,   # más pequeño aún
  legend.pos.y = 30,
  color.use = celltype_colors
  
)

pdf("/home/workdir/HDAC_tfm/fig/CellChat/KAP25L/control/chord_gene_source_Cancer_cells.pdf", width = 10, height = 10)

netVisual_chord_gene(
  cellchat,
  sources.use = 1,
  targets.use = c(1:8),
  color.use = celltype_colors,
  lab.cex = 0.7,
  legend.pos.y = 30
)

dev.off()

for (i in 1:7) {
  print(cells_source[i])
  pdf(paste0("/home/workdir/HDAC_tfm/fig/CellChat/KAP25L/control/chord_gene_source_",cells_source[i],"_control_2.pdf"), width = 10, height = 10)
  
  netVisual_chord_gene(
    cellchat,
    sources.use = i,
    targets.use = c(1:8),
    lab.cex = 0.7,
    color.use = celltype_colors,
    legend.pos.y = 30
  )
  
  dev.off()

}



for (i in 1:7) {
  print(i)
  pdf(paste0("/home/workdir/HDAC_tfm/fig/CellChat/KAP25L/control/chord_gene_target_",cells_source[i],"_control_2.pdf"), width = 10, height = 10)
  
  netVisual_chord_gene(
    cellchat,
    sources.use = c(1:8),
    targets.use = i,
    lab.cex = 0.7,
    color.use = celltype_colors,
    legend.pos.y = 30
  )
  
  dev.off()
  
}



# #PD-L1 ----------------------------------------------------------------

cellchat_k25l_pdl1<- readRDS("/home/workdir/HDAC_tfm/data/cellchat/cellchat_SS_30_05_KAP25L_PD-L1.rds")
cellchat <- cellchat_k25l_pdl1
unique_ligands_sig_pdl1 <- unique(unlist(cellchat@LR$LRsig$ligand))
ptm = Sys.time()
groupSize <- as.numeric(table(cellchat_k25l_pdl1@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Primer gráfico: número de interacciones
png("/home/workdir/HDAC_tfm/fig/CellChat/KAP25L/pdl1/circle_interactions_count_pdl1.png", 
    width = 2000, height = 2000, res = 300)

netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  color.use = celltype_colors,
  title.name = "Number of interactions"
)

dev.off()

# Segundo gráfico: fuerza de interacción
png("/home/workdir/HDAC_tfm/fig/CellChat/KAP25L/pdl1/circle_interactions_weight_pdl1.png", 
    width = 2000, height = 2000, res = 300)

netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  color.use = celltype_colors,
  title.name = "Interaction weights/strength"
)

dev.off()


mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}


pdf("/home/workdir/HDAC_tfm/fig/CellChat/KAP25L/pdl1/circle_per_signal_pdl1.pdf", 
    width = 16, height = 12)  # Puedes ajustar tamaño según lo necesites

mat <- cellchat@net$weight

par(mfrow = c(4, 2), xpd = TRUE, mar = c(1, 1, 2, 1))  # Márgenes para que quepa el título

for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(
    mat2,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    edge.weight.max = max(mat),
    title.name = rownames(mat)[i]
  )
}

dev.off()


# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0("/home/workdir/HDAC_tfm/fig/CellChat/KAP25L/pdl1/",pathways.show.all[i], "_L-R_contribution_pdl1.pdf"), plot=gg, width = 10, height = 12 , dpi = 600)
}

cells_source <- levels(cellchat@idents)

for (i in 1:8) {
  gg <- netVisual_bubble(cellchat, sources.use = i, targets.use = c(1:8), remove.isolate = FALSE)
  ggsave(filename=paste0("/home/workdir/HDAC_tfm/fig/CellChat/KAP25L/pdl1/",cells_source[i], "source_L-R_prob_pdl1.png"), plot=gg, width = 10, height = 12 , dpi = 600)
  
}



pdf("/home/workdir/HDAC_tfm/fig/CellChat/KAP25L/chord_gene_source_Cancer_cells.pdf", width = 10, height = 10)

netVisual_chord_gene(
  cellchat,
  sources.use = 1,
  targets.use = c(1:8),
  lab.cex = 0.7,
  legend.pos.y = 30
)

dev.off()

for (i in 1:5) {
  print(i)
  pdf(paste0("/home/workdir/HDAC_tfm/fig/CellChat/KAP25L/pdl1/chord_gene_source_",cells_source[i],"_pdl1.pdf"), width = 10, height = 10)
  
  netVisual_chord_gene(
    cellchat,
    sources.use = i,
    targets.use = c(1:6),
    lab.cex = 0.7,
    legend.pos.y = 30,
    color.use = celltype_colors
  )
  
  dev.off()
  
}

for (i in 1:5) {
  print(i)
  pdf(paste0("/home/workdir/HDAC_tfm/fig/CellChat/KAP25L/pdl1/chord_gene_target_",cells_source[i],"_pdl1.pdf"), width = 10, height = 10)
  
  netVisual_chord_gene(
    cellchat,
    sources.use = c(1:6),
    targets.use = i,
    lab.cex = 0.7,
    color.use = celltype_colors,
    legend.pos.y = 30
  )
  
  dev.off()
  
}


# DAC ---------------------------------------------------------------------


cellchat_k25l_dac<- readRDS("/home/workdir/HDAC_tfm/data/cellchat/cellchat_SS_30_05_KAP25L_DAC.rds")
cellchat <- cellchat_k25l_dac
unique_ligands_sig_dac <- unique(unlist(cellchat@LR$LRsig$ligand))
ptm = Sys.time()
groupSize <- as.numeric(table(cellchat_k25l_dac@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Primer gráfico: número de interacciones
png("/home/workdir/HDAC_tfm/fig/CellChat/KAP25L/dac/circle_interactions_count_dac.png", 
    width = 2000, height = 2000, res = 300)

netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  color.use = celltype_colors,
  title.name = "Number of interactions"
)

dev.off()

# Segundo gráfico: fuerza de interacción
png("/home/workdir/HDAC_tfm/fig/CellChat/KAP25L/dac/circle_interactions_weight_dac.png", 
    width = 2000, height = 2000, res = 300)

netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  color.use = celltype_colors,
  title.name = "Interaction weights/strength"
)

dev.off()


mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}


pdf("/home/workdir/HDAC_tfm/fig/CellChat/KAP25L/dac/circle_per_signal_dac.pdf", 
    width = 16, height = 12)  # Puedes ajustar tamaño según lo necesites

mat <- cellchat@net$weight

par(mfrow = c(4, 2), xpd = TRUE, mar = c(1, 1, 2, 1))  # Márgenes para que quepa el título

for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(
    mat2,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    edge.weight.max = max(mat),
    title.name = rownames(mat)[i]
  )
}

dev.off()


# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0("/home/workdir/HDAC_tfm/fig/CellChat/KAP25L/dac/",pathways.show.all[i], "_L-R_contribution_dac.pdf"), plot=gg, width = 10, height = 12 , dpi = 600)
}

cells_source <- levels(cellchat@idents)

for (i in 1:8) {
  gg <- netVisual_bubble(cellchat, sources.use = i, targets.use = c(1:8), remove.isolate = FALSE)
  ggsave(filename=paste0("/home/workdir/HDAC_tfm/fig/CellChat/KAP25L/dac/",cells_source[i], "source_L-R_prob_dac.png"), plot=gg, width = 10, height = 12 , dpi = 600)
  
}


netVisual_chord_gene(cellchat, sources.use = 1, targets.use = c(1:8), lab.cex = 0.5,legend.pos.y = 30)

netVisual_chord_gene(
  cellchat,
  sources.use = 1,
  targets.use = c(1:2),
  lab.cex = 0.1,   # más pequeño aún
  legend.pos.y = 30
)

pdf("/home/workdir/HDAC_tfm/fig/CellChat/KAP25L/chord_gene_source_Cancer_cells.pdf", width = 10, height = 10)

netVisual_chord_gene(
  cellchat,
  sources.use = 1,
  targets.use = c(1:8),
  color.use = celltype_colors,
  lab.cex = 0.7,
  legend.pos.y = 30
)

dev.off()

for (i in 1:7) {
  print(i)
  pdf(paste0("/home/workdir/HDAC_tfm/fig/CellChat/KAP25L/dac/chord_gene_source_",cells_source[i],"_dac.pdf"), width = 10, height = 10)
  
  netVisual_chord_gene(
    cellchat,
    sources.use = i,
    targets.use = c(1:8),
    lab.cex = 0.7,
    color.use = celltype_colors,
    legend.pos.y = 30
  )
  
  dev.off()
  
}

for (i in 1:7) {
  print(i)
  pdf(paste0("/home/workdir/HDAC_tfm/fig/CellChat/KAP25L/dac/chord_gene_target_",cells_source[i],"_dac.pdf"), width = 10, height = 10)
  
  netVisual_chord_gene(
    cellchat,
    sources.use = c(1:8),
    targets.use = i,
    lab.cex = 0.7,
    color.use = celltype_colors,
    legend.pos.y = 30
  )
  
  dev.off()
  
}



# comb --------------------------------------------------------------------




cellchat_k25l_comb<- readRDS("/home/workdir/HDAC_tfm/data/cellchat/cellchat_SS_30_05_KAP25L_Combination_min_cells_5.rds")
cellchat <- cellchat_k25l_comb

unique_ligands_sig_comb <- unique(unlist(cellchat@LR$LRsig$ligand))
ptm = Sys.time()
groupSize <- as.numeric(table(cellchat_k25l_comb@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Primer gráfico: número de interacciones
png("/home/workdir/HDAC_tfm/fig/CellChat/KAP25L/comb/circle_interactions_count_comb.png", 
    width = 2000, height = 2000, res = 300)

netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  color.use = celltype_colors,
  title.name = "Number of interactions"
)

dev.off()

# Segundo gráfico: fuerza de interacción
png("/home/workdir/HDAC_tfm/fig/CellChat/KAP25L/comb/circle_interactions_weight_comb.png", 
    width = 2000, height = 2000, res = 300)

netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  color.use = celltype_colors,
  title.name = "Interaction weights/strength"
)

dev.off()


mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}


pdf("/home/workdir/HDAC_tfm/fig/CellChat/KAP25L/comb/circle_per_signal_comb.pdf", 
    width = 16, height = 12)  # Puedes ajustar tamaño según lo necesites

mat <- cellchat@net$weight

par(mfrow = c(4, 2), xpd = TRUE, mar = c(1, 1, 2, 1))  # Márgenes para que quepa el título

for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(
    mat2,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    edge.weight.max = max(mat),
    title.name = rownames(mat)[i]
  )
}

dev.off()


# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0("/home/workdir/HDAC_tfm/fig/CellChat/KAP25L/comb/",pathways.show.all[i], "_L-R_contribution_comb.pdf"), plot=gg, width = 10, height = 12 , dpi = 600)
}

cells_source <- levels(cellchat@idents)

for (i in 1:8) {
  gg <- netVisual_bubble(cellchat, sources.use = i, targets.use = c(1:8), remove.isolate = FALSE)
  ggsave(filename=paste0("/home/workdir/HDAC_tfm/fig/CellChat/KAP25L/comb/",cells_source[i], "source_L-R_prob_comb.png"), plot=gg, width = 10, height = 12 , dpi = 600)
  
}


netVisual_chord_gene(cellchat, sources.use = 1, targets.use = c(1:8), lab.cex = 0.5,legend.pos.y = 30)

netVisual_chord_gene(
  cellchat,
  sources.use = 1,
  targets.use = c(1:2),
  lab.cex = 0.1,   # más pequeño aún
  legend.pos.y = 30
)

pdf("/home/workdir/HDAC_tfm/fig/CellChat/KAP25L/chord_gene_source_Cancer_cells.pdf", width = 10, height = 10)

netVisual_chord_gene(
  cellchat,
  sources.use = 1,
  targets.use = c(1:8),
  lab.cex = 0.7,
  legend.pos.y = 30
)

dev.off()

for (i in 1:7) {
  print(i)
  pdf(paste0("/home/workdir/HDAC_tfm/fig/CellChat/KAP25L/comb/chord_gene_source_",cells_source[i],"_comb.pdf"), width = 10, height = 10)
  
  netVisual_chord_gene(
    cellchat,
    sources.use = i,
    targets.use = c(1:8),
    lab.cex = 0.7,
    color.use = celltype_colors,
    legend.pos.y = 30
  )
  
  dev.off()
  
}

for (i in 1:7) {
  print(i)
  pdf(paste0("/home/workdir/HDAC_tfm/fig/CellChat/KAP25L/comb/chord_gene_target_",cells_source[i],"_comb.pdf"), width = 10, height = 10)
  
  netVisual_chord_gene(
    cellchat,
    sources.use = c(1:8),
    targets.use = i,
    lab.cex = 0.7,
    color.use = celltype_colors,
    legend.pos.y = 30
  )
  
  dev.off()
  
}

################


df_net_control <- subsetCommunication(cellchat_k25l_control)
df_net_combo <- subsetCommunication(cellchat_k25l_comb)
df_net_dac <- subsetCommunication(cellchat_k25l_dac)
df_net_pdl1 <- subsetCommunication(cellchat_k25l_pdl1)
diff_combo_control <- setdiff(df_net_combo$interaction_name_2, df_net_control$interaction_name_2)
diff_dac_control <- setdiff(df_net_dac$interaction_name_2, df_net_control$interaction_name_2)
diff_pdl1_control <- setdiff(df_net_pdl1$interaction_name_2, df_net_control$interaction_name_2)

diff_combo <- setdiff(diff_combo_control,diff_dac_control)
diff_dac <-  setdiff(diff_dac_control,diff_combo_control)


df_net_combo <- subsetCommunication(cellchat_k25l_comb)

# Suponiendo que diff_combo contiene nombres en interaction_name_2
df_filtered <- df_net_combo[df_net_combo$interaction_name_2 %in% diff_combo, ]

# Extraemos las vías asociadas a esas interacciones
pathways_to_plot <- unique(df_filtered$pathway_name)

# Ahora sí puedes llamar a la función con esas vías
netVisual_bubble(
  cellchat_k25l_comb,
  signaling = pathways_to_plot,
  remove.isolate = TRUE,
  title.name = "Unique pathways (Combo only)"
)






df_net_dac <- subsetCommunication(cellchat_k25l_dac)

# Suponiendo que diff_combo contiene nombres en interaction_name_2
df_filtered <- df_net_dac[df_net_dac$interaction_name_2 %in% diff_dac, ]

# Extraemos las vías asociadas a esas interacciones
pathways_to_plot <- unique(df_filtered$pathway_name)

# Ahora sí puedes llamar a la función con esas vías
netVisual_bubble(
  cellchat_k25l_dac,
  signaling = pathways_to_plot,
  remove.isolate = TRUE,
  title.name = "Unique pathways (DAC only)"
)





# # KAP25L-UV -------------------------------------------------------------

# Ver todas las interacciones del DB
head(CellChatDB.mouse$interaction)

# Buscar si existe la interacción PD-1 / PD-L1
CellChatDB.mouse$interaction %>%
  dplyr::filter(ligand == "Pdcd1" & receptor == "Cd274")


# #CONTROL ----------------------------------------------------------------

cellchat_kuv_control<- readRDS("/home/workdir/HDAC_tfm/data/cellchat/cellchat_SS_30_05_KAP25L-UV_control.rds")
cellchat <- cellchat_kuv_control
ptm = Sys.time()
groupSize <- as.numeric(table(cellchat_kuv_control@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Primer gráfico: número de interacciones
png("/home/workdir/HDAC_tfm/fig/CellChat/KUV/control/circle_interactions_count_control.png", 
    width = 2000, height = 2000, res = 300)

netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  color.use = celltype_colors,
  title.name = "Number of interactions"
)

dev.off()

# Segundo gráfico: fuerza de interacción
png("/home/workdir/HDAC_tfm/fig/CellChat/KUV/control/circle_interactions_weight_control.png", 
    width = 2000, height = 2000, res = 300)

netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  color.use = celltype_colors,
  title.name = "Interaction weights/strength"
)

dev.off()


mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}


pdf("/home/workdir/HDAC_tfm/fig/CellChat/KUV/control/circle_per_signal_control.pdf", 
    width = 16, height = 12)  # Puedes ajustar tamaño según lo necesites

mat <- cellchat@net$weight

par(mfrow = c(4, 2), xpd = TRUE, mar = c(1, 1, 2, 1))  # Márgenes para que quepa el título

for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(
    mat2,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    edge.weight.max = max(mat),
    title.name = rownames(mat)[i]
  )
}

dev.off()


# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0("/home/workdir/HDAC_tfm/fig/CellChat/KUV/control/",pathways.show.all[i], "_L-R_contribution_control.pdf"), plot=gg, width = 10, height = 12 , dpi = 600)
}

cells_source <- levels(cellchat@idents)

for (i in 1:8) {
  gg <- netVisual_bubble(cellchat, sources.use = i, targets.use = c(1:8), remove.isolate = FALSE)
  ggsave(filename=paste0("/home/workdir/HDAC_tfm/fig/CellChat/KUV/control/",cells_source[i], "source_L-R_prob_control.png"), plot=gg, width = 10, height = 12 , dpi = 600)
  
}


netVisual_chord_gene(cellchat, sources.use = 1, targets.use = c(1:8), lab.cex = 0.5,legend.pos.y = 30)

netVisual_chord_gene(
  cellchat,
  sources.use = 1,
  targets.use = c(1:2),
  lab.cex = 0.1,   # más pequeño aún
  legend.pos.y = 30
)

pdf("/home/workdir/HDAC_tfm/fig/CellChat/KUV/control/chord_gene_source_Cancer_cells.pdf", width = 10, height = 10)

netVisual_chord_gene(
  cellchat,
  sources.use = 1,
  targets.use = c(1:8),
  lab.cex = 0.7,
  legend.pos.y = 30
)

dev.off()

for (i in 1:7) {
  print(i)
  pdf(paste0("/home/workdir/HDAC_tfm/fig/CellChat/KUV/control/chord_gene_source_",cells_source[i],"_control.pdf"), width = 10, height = 10)
  
  netVisual_chord_gene(
    cellchat,
    sources.use = i,
    targets.use = c(1:8),
    lab.cex = 0.7,
    color.use = celltype_colors,
    legend.pos.y = 30
  )
  
  dev.off()
  
}

for (i in 1:7) {
  print(i)
  pdf(paste0("/home/workdir/HDAC_tfm/fig/CellChat/KUV/control/chord_gene_target_",cells_source[i],"_control.pdf"), width = 10, height = 10)
  
  netVisual_chord_gene(
    cellchat,
    sources.use = c(1:8),
    targets.use = i,
    lab.cex = 0.7,
    color.use = celltype_colors,
    legend.pos.y = 30
  )
  
  dev.off()
  
}



# #PD-L1 ----------------------------------------------------------------

cellchat_kuv_pdl1<- readRDS("/home/workdir/HDAC_tfm/data/cellchat/cellchat_SS_30_05_KAP25L-UV_PD-L1.rds")
cellchat <- cellchat_kuv_pdl1
ptm = Sys.time()
groupSize <- as.numeric(table(cellchat_kuv_pdl1@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Primer gráfico: número de interacciones
png("/home/workdir/HDAC_tfm/fig/CellChat/KUV/pdl1/circle_interactions_count_pdl1.png", 
    width = 2000, height = 2000, res = 300)

netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  color.use = celltype_colors,
  title.name = "Number of interactions"
)

dev.off()

# Segundo gráfico: fuerza de interacción
png("/home/workdir/HDAC_tfm/fig/CellChat/KUV/pdl1/circle_interactions_weight_pdl1.png", 
    width = 2000, height = 2000, res = 300)

netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  color.use = celltype_colors,
  title.name = "Interaction weights/strength"
)

dev.off()


mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}


pdf("/home/workdir/HDAC_tfm/fig/CellChat/KUV/pdl1/circle_per_signal_pdl1.pdf", 
    width = 16, height = 12)  # Puedes ajustar tamaño según lo necesites

mat <- cellchat@net$weight

par(mfrow = c(4, 2), xpd = TRUE, mar = c(1, 1, 2, 1))  # Márgenes para que quepa el título

for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(
    mat2,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    edge.weight.max = max(mat),
    title.name = rownames(mat)[i]
  )
}

dev.off()


# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0("/home/workdir/HDAC_tfm/fig/CellChat/KUV/pdl1/",pathways.show.all[i], "_L-R_contribution_pdl1.pdf"), plot=gg, width = 10, height = 12 , dpi = 600)
}

cells_source <- levels(cellchat@idents)

for (i in 1:8) {
  gg <- netVisual_bubble(cellchat, sources.use = i, targets.use = c(1:8), remove.isolate = FALSE)
  ggsave(filename=paste0("/home/workdir/HDAC_tfm/fig/CellChat/KUV/pdl1/",cells_source[i], "source_L-R_prob_pdl1.png"), plot=gg, width = 10, height = 12 , dpi = 600)
  
}


netVisual_chord_gene(cellchat, sources.use = 1, targets.use = c(1:8), lab.cex = 0.5,legend.pos.y = 30)

netVisual_chord_gene(
  cellchat,
  sources.use = 1,
  targets.use = c(1:2),
  lab.cex = 0.1,   # más pequeño aún
  legend.pos.y = 30
)

pdf("/home/workdir/HDAC_tfm/fig/CellChat/KUV/chord_gene_source_Cancer_cells.pdf", width = 10, height = 10)

netVisual_chord_gene(
  cellchat,
  sources.use = 1,
  targets.use = c(1:8),
  lab.cex = 0.7,
  legend.pos.y = 30
)

dev.off()

for (i in 1:3) {
  print(i)
  pdf(paste0("/home/workdir/HDAC_tfm/fig/CellChat/KUV/pdl1/chord_gene_source_",cells_source[i],"_pdl1.pdf"), width = 10, height = 10)
  
  netVisual_chord_gene(
    cellchat,
    sources.use = i,
    targets.use = c(1:8),
    lab.cex = 0.7,
    color.use = celltype_colors,
    legend.pos.y = 30
  )
  
  dev.off()
  
}

for (i in 1:3) {
  print(i)
  pdf(paste0("/home/workdir/HDAC_tfm/fig/CellChat/KUV/pdl1/chord_gene_target_",cells_source[i],"_pdl1.pdf"), width = 10, height = 10)
  
  netVisual_chord_gene(
    cellchat,
    sources.use = c(1:8),
    targets.use = i,
    lab.cex = 0.7,
    color.use = celltype_colors,
    legend.pos.y = 30
  )
  
  dev.off()
  
}


# DAC ---------------------------------------------------------------------


cellchat_kuv_dac<- readRDS("/home/workdir/HDAC_tfm/data/cellchat/cellchat_SS_30_05_KAP25L-UV_DAC.rds")
cellchat <- cellchat_kuv_dac
ptm = Sys.time()
groupSize <- as.numeric(table(cellchat_kuv_dac@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Primer gráfico: número de interacciones
png("/home/workdir/HDAC_tfm/fig/CellChat/KUV/dac/circle_interactions_count_dac.png", 
    width = 2000, height = 2000, res = 300)

netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  color.use = celltype_colors,
  title.name = "Number of interactions"
)

dev.off()

# Segundo gráfico: fuerza de interacción
png("/home/workdir/HDAC_tfm/fig/CellChat/KUV/dac/circle_interactions_weight_dac.png", 
    width = 2000, height = 2000, res = 300)

netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  color.use = celltype_colors,
  title.name = "Interaction weights/strength"
)

dev.off()


mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}


pdf("/home/workdir/HDAC_tfm/fig/CellChat/KUV/dac/circle_per_signal_dac.pdf", 
    width = 16, height = 12)  # Puedes ajustar tamaño según lo necesites

mat <- cellchat@net$weight

par(mfrow = c(4, 2), xpd = TRUE, mar = c(1, 1, 2, 1))  # Márgenes para que quepa el título

for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(
    mat2,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    edge.weight.max = max(mat),
    title.name = rownames(mat)[i]
  )
}

dev.off()


# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0("/home/workdir/HDAC_tfm/fig/CellChat/KUV/dac/",pathways.show.all[i], "_L-R_contribution_dac.pdf"), plot=gg, width = 10, height = 12 , dpi = 600)
}

cells_source <- levels(cellchat@idents)

for (i in 1:8) {
  gg <- netVisual_bubble(cellchat, sources.use = i, targets.use = c(1:8), remove.isolate = FALSE)
  ggsave(filename=paste0("/home/workdir/HDAC_tfm/fig/CellChat/KUV/dac/",cells_source[i], "source_L-R_prob_dac.png"), plot=gg, width = 10, height = 12 , dpi = 600)
  
}


netVisual_chord_gene(cellchat, sources.use = 1, targets.use = c(1:8), lab.cex = 0.5,legend.pos.y = 30)

netVisual_chord_gene(
  cellchat,
  sources.use = 1,
  targets.use = c(1:2),
  lab.cex = 0.1,   # más pequeño aún
  legend.pos.y = 30
)

pdf("/home/workdir/HDAC_tfm/fig/CellChat/KUV/chord_gene_source_Cancer_cells.pdf", width = 10, height = 10)

netVisual_chord_gene(
  cellchat,
  sources.use = 1,
  targets.use = c(1:8),
  lab.cex = 0.7,
  legend.pos.y = 30
)

dev.off()

for (i in 1:3) {
  print(i)
  pdf(paste0("/home/workdir/HDAC_tfm/fig/CellChat/KUV/dac/chord_gene_source_",cells_source[i],"_dac.pdf"), width = 10, height = 10)
  
  netVisual_chord_gene(
    cellchat,
    sources.use = i,
    targets.use = c(1:8),
    lab.cex = 0.7,
    color.use = celltype_colors,
    legend.pos.y = 30
  )
  
  dev.off()
  
}

for (i in 1:7) {
  print(i)
  pdf(paste0("/home/workdir/HDAC_tfm/fig/CellChat/KUV/dac/chord_gene_target_",cells_source[i],"_dac.pdf"), width = 10, height = 10)
  
  netVisual_chord_gene(
    cellchat,
    sources.use = c(1:8),
    targets.use = i,
    lab.cex = 0.7,
    color.use = celltype_colors,
    legend.pos.y = 30
  )
  
  dev.off()
  
}



# comb --------------------------------------------------------------------



cellchat_kuv_comb<- readRDS("/home/workdir/HDAC_tfm/data/cellchat/cellchat_SS_30_05_KAP25L-UV_Combination.rds")
cellchat <- cellchat_kuv_comb
ptm = Sys.time()
groupSize <- as.numeric(table(cellchat_kuv_comb@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Primer gráfico: número de interacciones
png("/home/workdir/HDAC_tfm/fig/CellChat/KUV/comb/circle_interactions_count_comb.png", 
    width = 2000, height = 2000, res = 300)

netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  color.use = celltype_colors,
  title.name = "Number of interactions"
)

dev.off()

# Segundo gráfico: fuerza de interacción
png("/home/workdir/HDAC_tfm/fig/CellChat/KUV/comb/circle_interactions_weight_comb.png", 
    width = 2000, height = 2000, res = 300)

netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  color.use = celltype_colors,
  title.name = "Interaction weights/strength"
)

dev.off()


mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}


pdf("/home/workdir/HDAC_tfm/fig/CellChat/KUV/comb/circle_per_signal_comb.pdf", 
    width = 16, height = 12)  # Puedes ajustar tamaño según lo necesites

mat <- cellchat@net$weight

par(mfrow = c(4, 2), xpd = TRUE, mar = c(1, 1, 2, 1))  # Márgenes para que quepa el título

for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(
    mat2,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    edge.weight.max = max(mat),
    title.name = rownames(mat)[i]
  )
}

dev.off()


# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0("/home/workdir/HDAC_tfm/fig/CellChat/KUV/comb/",pathways.show.all[i], "_L-R_contribution_comb.pdf"), plot=gg, width = 10, height = 12 , dpi = 600)
}

cells_source <- levels(cellchat@idents)

for (i in 1:8) {
  gg <- netVisual_bubble(cellchat, sources.use = i, targets.use = c(1:8), remove.isolate = FALSE)
  ggsave(filename=paste0("/home/workdir/HDAC_tfm/fig/CellChat/KUV/comb/",cells_source[i], "source_L-R_prob_comb.png"), plot=gg, width = 10, height = 12 , dpi = 600)
  
}


netVisual_chord_gene(cellchat, sources.use = 1, targets.use = c(1:8), lab.cex = 0.5,legend.pos.y = 30)

netVisual_chord_gene(
  cellchat,
  sources.use = 1,
  targets.use = c(1:2),
  lab.cex = 0.1,   # más pequeño aún
  legend.pos.y = 30
)

pdf("/home/workdir/HDAC_tfm/fig/CellChat/KUV/chord_gene_source_Cancer_cells.pdf", width = 10, height = 10)

netVisual_chord_gene(
  cellchat,
  sources.use = 1,
  targets.use = c(1:8),
  lab.cex = 0.7,
  legend.pos.y = 30
)

dev.off()

for (i in 1:3) {
  print(cells_source[i])
  pdf(paste0("/home/workdir/HDAC_tfm/fig/CellChat/KUV/comb/chord_gene_source_",cells_source[i],"_comb.pdf"), width = 10, height = 10)
  
  netVisual_chord_gene(
    cellchat,
    sources.use = i,
    targets.use = c(1:8),
    lab.cex = 0.7,
    color.use = celltype_colors,
    legend.pos.y = 30
  )
  
  dev.off()
  
}

for (i in 1:7) {
  print(i)
  pdf(paste0("/home/workdir/HDAC_tfm/fig/CellChat/KUV/comb/chord_gene_target_",cells_source[i],"_comb.pdf"), width = 10, height = 10)
  
  netVisual_chord_gene(
    cellchat,
    sources.use = c(1:8),
    targets.use = i,
    lab.cex = 0.7,
    color.use = celltype_colors,
    legend.pos.y = 30
  )
  
  dev.off()
  
}

