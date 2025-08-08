################################################################################
# Title: Pathway Enrichment analysis (PEA)
# Author: Paula
# Description:
#   This script performs PEA of the significant DEGs expressed in different cell types,
#   using GO BP, KEGG nad Reactome. 
#
################################################################################

library(Seurat)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(tidyverse)
library(dplyr)
library(pheatmap)

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


# Extraemos genes comunes como antes
genes_k25l_pdl1 <- rownames(deg_filtered$markers_ctrl_pdl1_k25l)
genes_k25l_dac  <- rownames(deg_filtered$markers_ctrl_dac_k25l)
genes_k25l_comb <- rownames(deg_filtered$markers_ctrl_comb_k25l)

common_k25l <- Reduce(intersect, list(genes_k25l_pdl1, genes_k25l_dac, genes_k25l_comb))

# Filtrar dataframes para quedarnos solo con genes comunes
df_pdl1 <- deg_filtered$markers_ctrl_pdl1_k25l[common_k25l, ] %>%
  rownames_to_column("gene") %>%
  rename_with(~ paste0(., "_pdl1"), -gene)

df_dac <- deg_filtered$markers_ctrl_dac_k25l[common_k25l, ] %>%
  rownames_to_column("gene") %>%
  rename_with(~ paste0(., "_dac"), -gene)

df_comb <- deg_filtered$markers_ctrl_comb_k25l[common_k25l, ] %>%
  rownames_to_column("gene") %>%
  rename_with(~ paste0(., "_comb"), -gene)

# Unir todo por gene
df_common_k25l <- df_pdl1 %>%
  full_join(df_dac, by = "gene") %>%
  full_join(df_comb, by = "gene")

write_csv(df_common_k25l, "/home/workdir/HDAC_tfm/data/DEGs_1406/k25l_common_DEGs_pdl1_dac_comb.csv")
genes_kuv_pdl1 <- rownames(deg_filtered$markers_ctrl_pdl1_kuv)
genes_kuv_dac  <- rownames(deg_filtered$markers_ctrl_dac_kuv)
genes_kuv_comb <- rownames(deg_filtered$markers_ctrl_comb_kuv)

common_kuv <- Reduce(intersect, list(genes_kuv_pdl1, genes_kuv_dac, genes_kuv_comb))

df_pdl1_kuv <- deg_filtered$markers_ctrl_pdl1_kuv[common_kuv, ] %>%
  rownames_to_column("gene") %>%
  rename_with(~ paste0(., "_pdl1"), -gene)

df_dac_kuv <- deg_filtered$markers_ctrl_dac_kuv[common_kuv, ] %>%
  rownames_to_column("gene") %>%
  rename_with(~ paste0(., "_dac"), -gene)

df_comb_kuv <- deg_filtered$markers_ctrl_comb_kuv[common_kuv, ] %>%
  rownames_to_column("gene") %>%
  rename_with(~ paste0(., "_comb"), -gene)

df_common_kuv <- df_pdl1_kuv %>%
  full_join(df_dac_kuv, by = "gene") %>%
  full_join(df_comb_kuv, by = "gene")

write_csv(df_common_kuv, "/home/workdir/HDAC_tfm/data/DEGs_1406/kuv_common_DEGs_pdl1_dac_comb.csv")

df_common_kuv <- read_csv("/home/workdir/HDAC_tfm/data/DEGs_1406/kuv_common_DEGs_pdl1_dac_comb.csv")
# KPB25L ------------------------------------------------------------------


ego <- enrichGO(gene = df_common_k25l$gene,
                OrgDb = org.Mm.eg.db,
                keyType = "SYMBOL",
                ont = "BP",              # BP: biological process
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

gene_entrez <- bitr(df_common_k25l$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)

ekegg <- enrichKEGG(gene = gene_entrez$ENTREZID,
                    organism = 'mmu', 
                    pvalueCutoff = 0.05)
ereact <- enrichPathway(gene = gene_entrez$ENTREZID,
                        organism = "mouse", 
                        pvalueCutoff = 0.05,
                        readable = TRUE)

go_annot <- ego@result %>%
  dplyr::filter(qvalue < 0.05) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::select(geneID, Description) %>%
  dplyr::group_by(geneID) %>%
  dplyr::summarise(GO_BP = paste(unique(Description), collapse = "; "), .groups = "drop") %>%
  dplyr::rename(gene = geneID)

kegg_annot <- ekegg@result %>%
  dplyr::filter(p.adjust < 0.05) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::select(geneID, Description) %>%
  dplyr::group_by(geneID) %>%
  dplyr::summarise(KEGG = paste(unique(Description), collapse = "; ")) %>%
  dplyr::left_join(gene_entrez, by = c("geneID" = "ENTREZID")) %>%
  dplyr::rename(gene = SYMBOL) %>%
  dplyr::select(gene, KEGG)


reactome_annot <- ereact@result %>%
  dplyr::filter(p.adjust < 0.05) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::select(geneID, Description) %>%
  dplyr::group_by(geneID) %>%
  dplyr::summarise(Reactome = paste(unique(Description), collapse = "; "), .groups = "drop") %>%
  dplyr::rename(gene = geneID)  # geneID ya es símbolo



degs_annotated <- df_common_k25l %>%
  left_join(go_annot, by = "gene") %>%
  left_join(kegg_annot, by = "gene") %>%
  left_join(reactome_annot, by = "gene")

write_csv(degs_annotated , "/home/workdir/HDAC_tfm/data/DEGs_1406/k25l_common_DEGs_pdl1_dac_comb_enriched.csv")

degs_annotated <- read_csv("/home/workdir/HDAC_tfm/data/DEGs_1406/k25l_common_DEGs_pdl1_dac_comb_enriched.csv")
degs_annotated_ordered_up <- degs_annotated %>%
  filter(p_val_adj_comb < 0.05) %>%
  arrange(desc(avg_log2FC_comb))

degs_annotated_ordered_down <- degs_annotated %>%
  filter(p_val_adj_comb < 0.05) %>%
  arrange(avg_log2FC_comb)




## UP -------------------------------------------------------------------


genes_ordered <- degs_annotated_ordered_up$gene 

heatmap_matrix <- degs_annotated %>%
  filter(gene %in% genes_ordered) %>%
  arrange(match(gene, genes_ordered)) %>%  # mantener el orden
  dplyr::select(gene, avg_log2FC_pdl1, avg_log2FC_dac, avg_log2FC_comb) %>%
  as.matrix()

rownames(heatmap_matrix) <- heatmap_matrix[,1]
heatmap_matrix <- heatmap_matrix[,2:4]

row_annot <- degs_annotated %>%
  filter(gene %in% genes_ordered) %>%
  arrange(match(gene, genes_ordered)) %>%
  mutate(
    Reactome_first = ifelse(is.na(Reactome), "", str_split(Reactome, ";") %>% map_chr(~str_trim(.[1]))),
    KEGG_first    = ifelse(is.na(KEGG),    "", str_split(KEGG, ";") %>% map_chr(~str_trim(.[1]))),
    GO_BP_first   = ifelse(is.na(GO_BP),   "", str_split(GO_BP, ";") %>% map_chr(~str_trim(.[1])))
  ) %>%
  mutate(
    Pathways = paste0(
      ifelse(Reactome_first == "", "", paste0("R: ", Reactome_first, "; ")),
      ifelse(KEGG_first    == "", "", paste0("K: ", KEGG_first, "; ")),
      ifelse(GO_BP_first   == "", "", paste0("G: ", GO_BP_first))
    )
  ) %>%
  dplyr::select(gene, Pathways)
rownames(row_annot) <- row_annot[,1]
row_annot<- row_annot %>% dplyr::select(-c(gene))
top_n <- 50
heatmap_matrix <- heatmap_matrix[1:top_n, ]
heatmap_matrix <- as.data.frame(heatmap_matrix)  # si no lo es ya
heatmap_matrix[] <- lapply(heatmap_matrix, as.numeric)  # convierte todas las columnas a numérico
heatmap_matrix <- as.matrix(heatmap_matrix)  # convierte a matriz numérica
row_annot <- row_annot[1:50, , drop=FALSE]
row_annot$Pathways[is.na(row_annot$Pathways) | row_annot$Pathways == ""] <- "None"
png("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/heatmap_common_degs_and_pathways_k25l_up.png", 
    width = 1400, height = 1500, res = 150, family = "Times")
pheatmap(
  mat = heatmap_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  annotation_row = row_annot,
  fontsize = 16,
  border_color = NA,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  main = "log2FC Heatmap - Common DEGs across treatments vs control in KPB25L cell line"
)
dev.off()



## DOWN --------------------------------------------------------------------



genes_ordered <- degs_annotated_ordered_down$gene 

heatmap_matrix <- degs_annotated %>%
  filter(gene %in% genes_ordered) %>%
  arrange(match(gene, genes_ordered)) %>%  # mantener el orden
  dplyr::select(gene, avg_log2FC_pdl1, avg_log2FC_dac, avg_log2FC_comb) %>%
  as.matrix()

rownames(heatmap_matrix) <- heatmap_matrix[,1]
heatmap_matrix <- heatmap_matrix[,2:4]

row_annot <- degs_annotated %>%
  filter(gene %in% genes_ordered) %>%
  arrange(match(gene, genes_ordered)) %>%
  mutate(
    Reactome_first = ifelse(is.na(Reactome), "", str_split(Reactome, ";") %>% map_chr(~str_trim(.[1]))),
    KEGG_first    = ifelse(is.na(KEGG),    "", str_split(KEGG, ";") %>% map_chr(~str_trim(.[1]))),
    GO_BP_first   = ifelse(is.na(GO_BP),   "", str_split(GO_BP, ";") %>% map_chr(~str_trim(.[1])))
  ) %>%
  mutate(
    Pathways = paste0(
      ifelse(Reactome_first == "", "", paste0("R: ", Reactome_first, "; ")),
      ifelse(KEGG_first    == "", "", paste0("K: ", KEGG_first, "; ")),
      ifelse(GO_BP_first   == "", "", paste0("G: ", GO_BP_first))
    )
  ) %>%
  dplyr::select(gene, Pathways)
rownames(row_annot) <- row_annot$gene
row_annot<- row_annot %>% dplyr::select(-c(gene))
top_n <- 50
heatmap_matrix <- heatmap_matrix[1:top_n, ]
heatmap_matrix <- as.data.frame(heatmap_matrix)  # si no lo es ya
heatmap_matrix[] <- lapply(heatmap_matrix, as.numeric)  # convierte todas las columnas a numérico
heatmap_matrix <- as.matrix(heatmap_matrix)  # convierte a matriz numérica
row_annot <- row_annot[1:50, , drop=FALSE]
row_annot$Pathways[is.na(row_annot$Pathways) | row_annot$Pathways == ""] <- "None"
png("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/heatmap_common_degs_and_pathways_k25l_down.png", 
    width = 1400, height = 1500, res = 150, family = "Times")
pheatmap(
  mat = heatmap_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  annotation_row = row_annot,
  fontsize = 16,
  border_color = NA,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  main = "log2FC Heatmap - Common DEGs across treatments vs control in KPB25L cell line"
)
dev.off()



# KPB25L-UV ---------------------------------------------------------------


ego <- enrichGO(gene = df_common_kuv$gene,
                OrgDb = org.Mm.eg.db,
                keyType = "SYMBOL",
                ont = "BP",              # BP: biological process
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

gene_entrez <- bitr(df_common_kuv$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)

ekegg <- enrichKEGG(gene = gene_entrez$ENTREZID,
                    organism = 'mmu', 
                    pvalueCutoff = 0.05)
ereact <- enrichPathway(gene = gene_entrez$ENTREZID,
                        organism = "mouse", 
                        pvalueCutoff = 0.05,
                        readable = TRUE)

go_annot <- ego@result %>%
  dplyr::filter(qvalue < 0.05) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::select(geneID, Description) %>%
  dplyr::group_by(geneID) %>%
  dplyr::summarise(GO_BP = paste(unique(Description), collapse = "; "), .groups = "drop") %>%
  dplyr::rename(gene = geneID)

kegg_annot <- ekegg@result %>%
  dplyr::filter(p.adjust < 0.05) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::select(geneID, Description) %>%
  dplyr::group_by(geneID) %>%
  dplyr::summarise(KEGG = paste(unique(Description), collapse = "; ")) %>%
  dplyr::left_join(gene_entrez, by = c("geneID" = "ENTREZID")) %>%
  dplyr::rename(gene = SYMBOL) %>%
  dplyr::select(gene, KEGG)


reactome_annot <- ereact@result %>%
  dplyr::filter(p.adjust < 0.05) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::select(geneID, Description) %>%
  dplyr::group_by(geneID) %>%
  dplyr::summarise(Reactome = paste(unique(Description), collapse = "; "), .groups = "drop") %>%
  dplyr::rename(gene = geneID)  # geneID ya es símbolo



degs_annotated <- df_common_kuv %>%
  left_join(go_annot, by = "gene") %>%
  left_join(kegg_annot, by = "gene") %>%
  left_join(reactome_annot, by = "gene")

write_csv(degs_annotated , "/home/workdir/HDAC_tfm/data/DEGs_1406/kuv_common_DEGs_pdl1_dac_comb_enriched.csv")

degs_annotated <- read_csv("/home/workdir/HDAC_tfm/data/DEGs_1406/kuv_common_DEGs_pdl1_dac_comb_enriched.csv")
degs_annotated_ordered_up <- degs_annotated %>%
  filter(p_val_adj_comb < 0.05) %>%
  arrange(desc(avg_log2FC_comb))

degs_annotated_ordered_down <- degs_annotated %>%
  filter(p_val_adj_comb < 0.05) %>%
  arrange(avg_log2FC_comb)





## UP -------------------------------------------------------------------


genes_ordered <- degs_annotated_ordered_up$gene 

heatmap_matrix <- degs_annotated %>%
  filter(gene %in% genes_ordered) %>%
  arrange(match(gene, genes_ordered)) %>%  # mantener el orden
  dplyr::select(gene, avg_log2FC_pdl1, avg_log2FC_dac, avg_log2FC_comb) %>%
  as.matrix()

rownames(heatmap_matrix) <- heatmap_matrix[,1]
heatmap_matrix <- heatmap_matrix[,2:4]

row_annot <- degs_annotated %>%
  filter(gene %in% genes_ordered) %>%
  arrange(match(gene, genes_ordered)) %>%
  mutate(
    Reactome_first = ifelse(is.na(Reactome), "", str_split(Reactome, ";") %>% map_chr(~str_trim(.[1]))),
    KEGG_first    = ifelse(is.na(KEGG),    "", str_split(KEGG, ";") %>% map_chr(~str_trim(.[1]))),
    GO_BP_first   = ifelse(is.na(GO_BP),   "", str_split(GO_BP, ";") %>% map_chr(~str_trim(.[1])))
  ) %>%
  mutate(
    Pathways = paste0(
      ifelse(Reactome_first == "", "", paste0("R: ", Reactome_first, "; ")),
      ifelse(KEGG_first    == "", "", paste0("K: ", KEGG_first, "; ")),
      ifelse(GO_BP_first   == "", "", paste0("G: ", GO_BP_first))
    )
  ) %>%
  dplyr::select(gene, Pathways)
row_annot<- row_annot %>% dplyr::select(-c(gene))
top_n <- 50
heatmap_matrix <- heatmap_matrix[1:top_n, ]
heatmap_matrix <- as.data.frame(heatmap_matrix)  # si no lo es ya
heatmap_matrix[] <- lapply(heatmap_matrix, as.numeric)  # convierte todas las columnas a numérico
heatmap_matrix <- as.matrix(heatmap_matrix)  # convierte a matriz numérica
row_annot <- row_annot[1:50, , drop=FALSE]
row_annot$Pathways[is.na(row_annot$Pathways) | row_annot$Pathways == ""] <- "None"
png("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/heatmap_common_degs_and_pathways_kuv_up.png", 
    width = 1400, height = 1500, res = 150, family = "Times")
pheatmap(
  mat = heatmap_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  annotation_row = row_annot,
  fontsize = 16,
  border_color = NA,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  main = "log2FC Heatmap - Common DEGs across treatments vs control in KPB25L-UV cell line"
)
dev.off()



## DOWN --------------------------------------------------------------------



genes_ordered <- degs_annotated_ordered_down$gene 

heatmap_matrix <- degs_annotated %>%
  filter(gene %in% genes_ordered) %>%
  arrange(match(gene, genes_ordered)) %>%  # mantener el orden
  dplyr::select(gene, avg_log2FC_pdl1, avg_log2FC_dac, avg_log2FC_comb) %>%
  as.matrix()

rownames(heatmap_matrix) <- heatmap_matrix[,1]
heatmap_matrix <- heatmap_matrix[,2:4]

row_annot <- degs_annotated %>%
  filter(gene %in% genes_ordered) %>%
  arrange(match(gene, genes_ordered)) %>%
  mutate(
    Reactome_first = ifelse(is.na(Reactome), "", str_split(Reactome, ";") %>% map_chr(~str_trim(.[1]))),
    KEGG_first    = ifelse(is.na(KEGG),    "", str_split(KEGG, ";") %>% map_chr(~str_trim(.[1]))),
    GO_BP_first   = ifelse(is.na(GO_BP),   "", str_split(GO_BP, ";") %>% map_chr(~str_trim(.[1])))
  ) %>%
  mutate(
    Pathways = paste0(
      ifelse(Reactome_first == "", "", paste0("R: ", Reactome_first, "; ")),
      ifelse(KEGG_first    == "", "", paste0("K: ", KEGG_first, "; ")),
      ifelse(GO_BP_first   == "", "", paste0("G: ", GO_BP_first))
    )
  ) %>%
  dplyr::select(gene, Pathways)
rownames(row_annot) <- row_annot$gene
row_annot<- row_annot %>% dplyr::select(-c(gene))
top_n <- 50
heatmap_matrix <- heatmap_matrix[1:top_n, ]
heatmap_matrix <- as.data.frame(heatmap_matrix)  # si no lo es ya
heatmap_matrix[] <- lapply(heatmap_matrix, as.numeric)  # convierte todas las columnas a numérico
heatmap_matrix <- as.matrix(heatmap_matrix)  # convierte a matriz numérica
row_annot <- row_annot[1:50, , drop=FALSE]
row_annot$Pathways[is.na(row_annot$Pathways) | row_annot$Pathways == ""] <- "None"
png("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/heatmap_common_degs_and_pathways_kuv_down.png", 
    width = 1500, height = 1500, res = 150, family = "Times")
pheatmap(
  mat = heatmap_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  annotation_row = row_annot,
  fontsize = 16,
  border_color = NA,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  main = "log2FC Heatmap - Common DEGs across treatments vs control in KPB25L-UV cell line"
)
dev.off()



# UNIQUE DEGS -------------------------------------------------------------

genes_pdl1 <- rownames(deg_filtered$markers_ctrl_pdl1_k25l)
genes_dac  <- rownames(deg_filtered$markers_ctrl_dac_k25l)
genes_comb <- rownames(deg_filtered$markers_ctrl_comb_k25l)

# Genes exclusivos de cada tratamiento
only_pdl1 <- setdiff(genes_pdl1, union(genes_dac, genes_comb))
only_dac  <- setdiff(genes_dac, union(genes_pdl1, genes_comb))
only_comb <- setdiff(genes_comb, union(genes_pdl1, genes_dac))

df_pdl1 <- deg_filtered$markers_ctrl_pdl1_k25l[only_pdl1, ] %>%
  rownames_to_column("gene") %>%
  rename_with(~ paste0(., "_pdl1"), -gene)

df_dac <- deg_filtered$markers_ctrl_dac_k25l[only_dac, ] %>%
  rownames_to_column("gene") %>%
  rename_with(~ paste0(., "_dac"), -gene)

df_comb <- deg_filtered$markers_ctrl_comb_k25l[only_comb, ] %>%
  rownames_to_column("gene") %>%
  rename_with(~ paste0(., "_comb"), -gene)

genes_pdl1_entrez <- bitr(df_pdl1$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")$ENTREZID
genes_dac_entrez  <- bitr(df_dac$gene,  fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")$ENTREZID
genes_comb_entrez <- bitr(df_comb$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")$ENTREZID

comp_GO <- compareCluster(
  geneCluster = list(
    PDL1 = genes_pdl1_entrez,
    DAC = genes_dac_entrez,
    COMB = genes_comb_entrez
  ),
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  ont = "CC",  # Cambia a "MF" o "CC" si quieres explorar otras ontologías
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  readable = TRUE  # Convierte de Entrez a símbolo para mejor legibilidad
)

p<- dotplot(comp_GO, showCategory = 30, split = "Cluster") + ggtitle("GO CC enrichment")
print(p)

ggsave("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/GO_dotplot.pdf", plot = last_plot(), width = 10, height = 6)

comp_reactome <- compareCluster(
  geneCluster = list(
    PDL1 = genes_pdl1_entrez,
    DAC = genes_dac_entrez,
    COMB = genes_comb_entrez
  ),
  fun = "enrichPathway",
  organism = "mouse", 
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  readable = TRUE
)

p_reactome <- dotplot(comp_reactome, showCategory = 30, split = "Cluster") +
  ggtitle("Reactome Pathway Enrichment")

print(p_reactome)

ggsave(
  "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/Reactome_dotplot.pdf",
  plot = p_reactome,
  width = 10,
  height = 6
)

## NO HAY KEGG

#KUV
genes_pdl1 <- rownames(deg_filtered$markers_ctrl_pdl1_kuv)
genes_dac  <- rownames(deg_filtered$markers_ctrl_dac_kuv)
genes_comb <- rownames(deg_filtered$markers_ctrl_comb_kuv)

# Genes exclusivos de cada tratamiento
only_pdl1 <- setdiff(genes_pdl1, union(genes_dac, genes_comb))
only_dac  <- setdiff(genes_dac, union(genes_pdl1, genes_comb))
only_comb <- setdiff(genes_comb, union(genes_pdl1, genes_dac))

df_pdl1 <- deg_filtered$markers_ctrl_pdl1_k25l[only_pdl1, ] %>%
  rownames_to_column("gene") %>%
  rename_with(~ paste0(., "_pdl1"), -gene)

df_dac <- deg_filtered$markers_ctrl_dac_k25l[only_dac, ] %>%
  rownames_to_column("gene") %>%
  rename_with(~ paste0(., "_dac"), -gene)

df_comb <- deg_filtered$markers_ctrl_comb_k25l[only_comb, ] %>%
  rownames_to_column("gene") %>%
  rename_with(~ paste0(., "_comb"), -gene)

genes_pdl1_entrez <- bitr(df_pdl1$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")$ENTREZID
genes_dac_entrez  <- bitr(df_dac$gene,  fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")$ENTREZID
genes_comb_entrez <- bitr(df_comb$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")$ENTREZID

comp_GO <- compareCluster(
  geneCluster = list(
    PDL1 = genes_pdl1_entrez,
    DAC = genes_dac_entrez,
    COMB = genes_comb_entrez
  ),
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  ont = "BP",  # Cambia a "MF" o "CC" si quieres explorar otras ontologías
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  readable = TRUE  # Convierte de Entrez a símbolo para mejor legibilidad
)

p<- dotplot(comp_GO, showCategory = 20, split = "Cluster") + ggtitle("GO BP enrichment")
print(p)

ggsave("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/GO_BP_dotplot_kuv.pdf", plot = p, width = 20, height = 23)

comp_reactome <- compareCluster(
  geneCluster = list(
    PDL1 = genes_pdl1_entrez,
    DAC = genes_dac_entrez,
    COMB = genes_comb_entrez
  ),
  fun = "enrichPathway",
  organism = "mouse", 
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  readable = TRUE
)

p_reactome <- dotplot(comp_reactome, showCategory = 30, split = "Cluster") +
  ggtitle("Reactome Pathway Enrichment")

print(p_reactome)

ggsave(
  "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/Reactome_dotplot_kuv.pdf",
  plot = p_reactome,
  width = 12,
  height = 15
)




comp_reactome <- compareCluster(
  geneCluster = list(
    PDL1 = genes_pdl1_entrez,
    DAC = genes_dac_entrez,
    COMB = genes_comb_entrez
  ),
  fun = "enrichKEGG",
  organism = "mouse", 
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

p_reactome <- dotplot(comp_reactome, showCategory = 30, split = "Cluster") +
  ggtitle("KEGG Pathway Enrichment")

print(p_reactome)

ggsave(
  "/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/KEGG_dotplot_kuv.pdf",
  plot = p_reactome,
  width = 12,
  height = 15
)