library(Seurat)
library(ggplot2)
library(htmlwidgets)
library(webshot2)
tumor_cells <- readRDS("/home/workdir/HDAC_tfm/data/tumor_cells_adjusted_13032025_with_cellcycle.rds")

cells_responsive_intersect_tfs_uv <- readRDS( "/home/workdir/HDAC_tfm/data/cells_responsive_intersect_tfs_kpb25l_uv.rds")
cells_unresponsive_intersect_tfs_uv <- readRDS("/home/workdir/HDAC_tfm/data/cells_unresponsive_intersect_tfs_kpb25l_uv.rds")

cells_responsive_intersect_tfs <- readRDS("/home/workdir/HDAC_tfm/data/cells_responsive_intersect_tfs_kpb25l.rds")
cells_unresponsive_intersect_tfs<- readRDS("/home/workdir/HDAC_tfm/data/cells_unresponsive_intersect_tfs_kpb25l.rds")


tumor_responsive_unresponsive_kpb25l <- subset(tumor_cells, cells = c(cells_responsive_intersect_tfs,cells_unresponsive_intersect_tfs ))


tumor_responsive_unresponsive_kpb25l_uv <- subset(tumor_cells, cells = c(cells_responsive_intersect_tfs_uv,cells_unresponsive_intersect_tfs_uv ))


labels_kpb25l <- ifelse(Cells(tumor_responsive_unresponsive_kpb25l) %in% cells_responsive_intersect_tfs,
                        "Responder", "Non-responder")

tumor_responsive_unresponsive_kpb25l$TF_response <- labels_kpb25l

labels_kpb25l_uv <- ifelse(Cells(tumor_responsive_unresponsive_kpb25l_uv) %in% cells_responsive_intersect_tfs_uv,
                           "Responder", "Non-responder")

tumor_responsive_unresponsive_kpb25l_uv$TF_response <- labels_kpb25l_uv

tumor_responsive_unresponsive_kpb25l$Seurat_Phase_TF_res <- paste(
  tumor_responsive_unresponsive_kpb25l$Phase,
  tumor_responsive_unresponsive_kpb25l$TF_response,
  sep = "_"
)

tumor_responsive_unresponsive_kpb25l_uv$Seurat_Phase_TF_res <- paste(
  tumor_responsive_unresponsive_kpb25l_uv$Phase,
  tumor_responsive_unresponsive_kpb25l_uv$TF_response,
  sep = "_"
)

df <- data.frame(table(tumor_responsive_unresponsive_kpb25l$TF_response, tumor_responsive_unresponsive_kpb25l$Phase))
df_uv <- data.frame(table(tumor_responsive_unresponsive_kpb25l_uv$TF_response, tumor_responsive_unresponsive_kpb25l_uv$Phase))

colnames(df) <- c("source", "target", "value")
nodes <- data.frame(name = unique(c(df$source, df$target)))

df$source_id <- match(df$source, nodes$name) - 1
df$target_id <- match(df$target, nodes$name) - 1

my_color <- 'd3.scaleOrdinal() .domain(["Non-responder", "Responder","G1", "G2M", "S"]) .range([ "#5A7CA5", "#B7D1F1", "#ff8000", "green", "#ffd5b1" ])'

sankey <- sankeyNetwork(Links = df,
                        Nodes = nodes,
                        Source = "source_id",
                        Target = "target_id",
                        Value = "value",
                        NodeID = "name",
                        fontSize = 30,
                        colourScale= my_color,
                        nodeWidth = 30)

saveWidget(sankey, file = "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/sankey_temp.html", selfcontained = TRUE)

colnames(df_uv) <- c("source", "target", "value")
nodes <- data.frame(name = unique(c(df_uv$source, df_uv$target)))

df_uv$source_id <- match(df$source, nodes$name) - 1
df_uv$target_id <- match(df$target, nodes$name) - 1

my_color <- 'd3.scaleOrdinal() .domain(["Non-responder", "Responder","G1", "G2M", "S"]) .range([ "#aa8caf", "#e2d8e4", "#ff8000", "green", "#ffd5b1" ])'
sankey <- sankeyNetwork(Links = df_uv,
                        Nodes = nodes,
                        Source = "source_id",
                        Target = "target_id",
                        Value = "value",
                        NodeID = "name",
                        fontSize = 30,
                        colourScale=my_color,
                        nodeWidth = 30)

saveWidget(sankey, file = "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/sankey_temp_uv.html", selfcontained = TRUE)