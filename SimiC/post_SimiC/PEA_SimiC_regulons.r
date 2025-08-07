library(Seurat)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(tidyverse)
library(dplyr)
library(patchwork)


# PEA REGULON TFs from KPB25L -------------------------------------------


## Bbx ---------------------------------------------------------------------
line <- "KPB25L"
TF <- "Bbx"
cell_type <- "Cancer cells"

regulon <- c(
  "Phf21a", "Cul4b", "Enah", "ENSMUSG00000120668", "Mbd5", "Asph", "Senp7", "Map2k6",
  "Sipa1l1", "Tnrc6b", "Fbxl20", "Sh3rf1", "Vps13a", "Fer", "Sgms1", "Slc23a2",
  "Tgfbr3", "Galk2", "Tnc", "Fgf7", "Pde7a", "Has2", "Gsk3b", "Lyst", "Atr",
  "Smyd3", "Igf2bp2", "Lrrc8c", "Rbms3", "Il1rap", "Nmt2", "Pcgf5", "Ralgapa2",
  "Tanc2", "Gab1", "Spidr", "Hs3st5", "Dzip3", "Dcbld2", "Ephb2", "Rbfox2", "Nf1",
  "Cmip", "Efna5", "Fto", "Dennd1a", "Dym", "Camk2d", "Btbd11", "Cdh13", "Lrig1",
  "Utrn", "Sclt1", "Klhl24", "Csgalnact1", "Dst", "Sdc2", "Map3k5", "Unc5b", "Ddr2",
  "Rbms1", "Sema5a", "Lrrk2", "Trpm6", "Gm50397", "Luc7l2", "PTPRG", "Ghr", "Rnf217",
  "Atosa", "Stag3", "Suz12", "Nova1", "Tanc1", "Cd44", "Tpm1", "Col4a3", "Taf15",
  "Spp1", "Tmc5", "Hdac8", "St3gal5", "Adgrl2", "Ube2h", "Dtd1", "Ide", "Crybg1",
  "Fendrr", "Lpin2", "Atxn7l1", "Nectin3", "Ncapd3", "Lhfpl2", "Nav1", "Dmxl1",
  "Igf1r", "Npepps", "Mtbp", "Clasp2", "Kank1"
)

genes_entrez <- bitr(regulon, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
enrich_go_result <- enrichGO(gene = genes_entrez$ENTREZID, 
                             OrgDb = org.Mm.eg.db, 
                             ont = "BP",  
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.05)

top_terms <- enrich_go_result %>%
  as.data.frame() %>%
  arrange(desc(GeneRatio), p.adjust) %>%
  slice_head(n = 15)

top_terms$GeneRatio <- sapply(strsplit(as.character(top_terms$GeneRatio), "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]))

p_go <- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                      color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
  scale_size_continuous(name = "Count") +
  labs(
    title = paste("GO BP enrichment for", TF, "'s regulon in", line, cell_type),
    x = "Gene Ratio", y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(family = "Times New Roman", face = "bold", size = 18),
    axis.text.y = element_text(family = "Times New Roman", size = 16, color="black"),
    axis.text.x = element_text(family = "Times New Roman", size = 16),
    axis.title = element_text(family = "Times New Roman", size = 16),
    plot.title = element_text(hjust = 1, family = "Times New Roman", size = 16,face = "bold"),
    legend.title = element_text(family = "Times New Roman", size = 16),
    legend.text = element_text(family = "Times New Roman", size = 14)
  )

ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","GO_enrich_",TF,".png"), plot = p_go, width = 10, height = 7, dpi = 600)



enrich_KEGG <- enrichKEGG(gene = genes_entrez$ENTREZID, 
                                       organism = "mouse", 
                                       pAdjustMethod = "BH", 
                                       pvalueCutoff = 0.05)

top_terms <- enrich_KEGG %>%
  as.data.frame() %>%
  arrange(desc(GeneRatio), p.adjust) %>%
  slice_head(n = 15)


top_terms$GeneRatio <- sapply(strsplit(as.character(top_terms$GeneRatio), "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]))

p_k<- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                      color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
  scale_size_continuous(name = "Count") +
  labs(
    title = paste("KEGG enrichment for", TF, "'s regulon in", line, cell_type),
    x = "Gene Ratio", y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(family = "Times New Roman", face = "bold", size = 18),
    axis.text.y = element_text(family = "Times New Roman", size = 16, color="black"),
    axis.text.x = element_text(family = "Times New Roman", size = 16),
    axis.title = element_text(family = "Times New Roman", size = 16),
    plot.title = element_text(hjust = 1, family = "Times New Roman", size = 16,face = "bold"),
    legend.title = element_text(family = "Times New Roman", size = 16),
    legend.text = element_text(family = "Times New Roman", size = 14))
    
ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","KEGG_enrich_",TF,".png"), plot = p_k, width = 8, height = 7, dpi = 600)
    
    


enrich_pathway_result <- enrichPathway(gene = genes_entrez$ENTREZID, 
                                       organism = "mouse", 
                                       pAdjustMethod = "BH", 
                                       pvalueCutoff = 0.2)



top_terms <- enrich_pathway_result %>%
  as.data.frame() %>%
  arrange(desc(GeneRatio), p.adjust) %>%
  slice_head(n = 15)


top_terms$GeneRatio <- sapply(strsplit(as.character(top_terms$GeneRatio), "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]))


p_r<- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                          color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
  scale_size_continuous(name = "Count") +
  labs(
    title = paste("Reactome enrichment for", TF, "'s regulon in", line, cell_type),
    x = "Gene Ratio", y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(family = "Times New Roman", face = "bold", size = 18),
    axis.text.y = element_text(family = "Times New Roman", size = 16, color="black"),
    axis.text.x = element_text(family = "Times New Roman", size = 16),
    axis.title = element_text(family = "Times New Roman", size = 16),
    plot.title = element_text(hjust = 1, family = "Times New Roman", size = 16,face = "bold"),
    legend.title = element_text(family = "Times New Roman", size = 16),
    legend.text = element_text(family = "Times New Roman", size = 14))

ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","React_enrich_",TF,".png"), plot = p, width = 8, height = 7, dpi = 600)

# NO HAY 


combined_plot <- p_go + p_k + plot_layout(ncol = 2)

ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","Combined_enrich_",TF,".png"), plot = combined_plot, width = 18, height = 7, dpi = 600)



## Foxp1 ---------------------------------------------------------------------
line <- "KPB25L"
TF <- "Foxp1"
cell_type <- "Cancer cells"

regulon <-  c(
  "Fndc3b", "Ext1", "Osmr", "Diaph2", "Ctif", "Phldb2", "Pde4d", "Nav2", "Efnb2", "Cpne8",
  "App", "Kcnk1", "Fto", "Lrp1", "Nedd4", "Cald1", "Itgb5", "Dennd1a", "Lpp", "Aig1",
  "Nrp2", "Rbms3", "Ptprk", "Trio", "Aebp1", "Wwox", "B4galnt3", "Smyd3", "Tns3", "Asap1",
  "Col5a2", "Gmds", "Slc4a4", "Maml3", "5430437J10Rik", "Fgfr1", "Cblb", "Sema4d", "Atp8a1",
  "Peak1", "Inpp4b", "Gm49890", "Exoc4", "Igf1r", "Arhgap42", "Lhfpl2", "Fstl1", "Cacna1c",
  "Prickle1", "Casp4", "Large1", "Mical2", "2610307P16Rik", "Vps13b", "Uxs1", "Col3a1",
  "Map1b", "Vgll4", "Map4k3", "Col6a1", "Usp10", "Samd4", "Cmip", "Morc4", "Cdkn2a",
  "Hs3st5", "Antxr2", "Egfr", "Dock1", "Parva", "Mgat4a", "Cdh1", "4930467D21Rik",
  "Pdgfra", "Notch2", "Nisch", "Pdgfrb", "Nhsl1", "Patj", "Bmper", "Denn2b", "Magi1",
  "Senp7", "Kmt2e", "Ptges", "Ddx60", "Fbxl17", "Arhgap22", "Mical3", "Stk3", "Akap13",
  "Ubash3b", "Rgs22", "Fgf7", "Myo1d", "Cntn5", "Gm26652", "Rin2", "Airn", "Gm14636"
)
genes_entrez <- bitr(regulon, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
enrich_go_result <- enrichGO(gene = genes_entrez$ENTREZID, 
                             OrgDb = org.Mm.eg.db, 
                             ont = "BP",  
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.05)

top_terms <- enrich_go_result %>%
  as.data.frame() %>%
  arrange(desc(GeneRatio), p.adjust) %>%
  slice_head(n = 15)


top_terms$GeneRatio <- sapply(strsplit(as.character(top_terms$GeneRatio), "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]))

p_go <- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                           color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
  scale_size_continuous(name = "Count") +
  labs(
    title = paste("GO BP enrichment for", TF, "'s regulon in", line, cell_type),
    x = "Gene Ratio", y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(family = "Times New Roman", face = "bold", size = 18),
    axis.text.y = element_text(family = "Times New Roman", size = 16, color="black"),
    axis.text.x = element_text(family = "Times New Roman", size = 16),
    axis.title = element_text(family = "Times New Roman", size = 16),
    plot.title = element_text(hjust = 1, family = "Times New Roman", size = 16,face = "bold"),
    legend.title = element_text(family = "Times New Roman", size = 16),
    legend.text = element_text(family = "Times New Roman", size = 14)
  )

ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","GO_enrich_",TF,".png"), plot = p, width = 10, height = 7, dpi = 600)



enrich_KEGG <- enrichKEGG(gene = genes_entrez$ENTREZID, 
                          organism = "mouse", 
                          pAdjustMethod = "BH", 
                          pvalueCutoff = 0.05)

top_terms <- enrich_KEGG %>%
  as.data.frame() %>%
  arrange(desc(GeneRatio), p.adjust) %>%
  slice_head(n = 15)

top_terms$GeneRatio <- sapply(strsplit(as.character(top_terms$GeneRatio), "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]))

p_k<- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                          color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
  scale_size_continuous(name = "Count") +
  labs(
    title = paste("KEGG enrichment for", TF, "'s regulon in", line, cell_type),
    x = "Gene Ratio", y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(family = "Times New Roman", face = "bold", size = 18),
    axis.text.y = element_text(family = "Times New Roman", size = 16, color = "black"),
    axis.text.x = element_text(family = "Times New Roman", size = 16),
    axis.title = element_text(family = "Times New Roman", size = 16),
    plot.title = element_text(hjust = 1, family = "Times New Roman", size = 16, face = "bold"),
    legend.title = element_text(family = "Times New Roman", size = 16),
    legend.text = element_text(family = "Times New Roman", size = 14)
  )
ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","KEGG_enrich_",TF,".png"), plot = p, width = 8, height = 7, dpi = 600)




enrich_pathway_result <- enrichPathway(gene = genes_entrez$ENTREZID, 
                                       organism = "mouse", 
                                       pAdjustMethod = "BH", 
                                       pvalueCutoff = 0.05)




top_terms <- enrich_pathway_result %>%
  as.data.frame() %>%
  arrange(desc(GeneRatio), p.adjust) %>%
  slice_head(n = 15)

top_terms$GeneRatio <- sapply(strsplit(as.character(top_terms$GeneRatio), "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]))

p_r<- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                          color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
  scale_size_continuous(name = "Count") +
  labs(
    title = paste("Reactome enrichment for", TF, "'s regulon in", line, cell_type),
    x = "Gene Ratio", y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(family = "Times New Roman", face = "bold", size = 18),
    axis.text.y = element_text(family = "Times New Roman", size = 16, color="black"),
    axis.text.x = element_text(family = "Times New Roman", size = 16),
    axis.title = element_text(family = "Times New Roman", size = 16),
    plot.title = element_text(hjust = 1, family = "Times New Roman", size = 16,face = "bold"),
    legend.title = element_text(family = "Times New Roman", size = 16),
    legend.text = element_text(family = "Times New Roman", size = 14))

ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","React_enrich_",TF,".png"), plot = p, width = 8, height = 7, dpi = 600)


combined_plot <- p_go + p_k + p_r + plot_layout(ncol = 3)

ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","Combined_enrich_",TF,".png"), plot = combined_plot, width = 30, height = 7, dpi = 600)



## Zeb2 -------------------------------------------------------------------

line <- "KPB25L"
TF <- "Zeb2"
cell_type <- "Cancer cells"

regulon <-   c(
  "Spp1", "Slit2", "Casp12", "Psd3", "Arhgap44", "Mgat4a", "Lhfpl2", "Tgfbr3",
  "Susd6", "Kifc3", "Usp46", "Prickle1", "Nipal2", "Ccnd3", "Slco3a1", "Calcb",
  "Mme", "Fgd4", "Skap2", "Zcchc24", "Sema5a", "Mdfic", "Vcan", "Kank1",
  "Shroom3", "Irak2", "Mindy3", "Luc7l2", "Kif26b", "Nedd4l", "Vav3", "Arhgap6",
  "Peak1", "Dpysl3", "Mical2", "Ank", "Ccn5", "A830082K12Rik", "Casp4", "Ldlrad3",
  "Fgfr1", "Wls", "Pik3cb", "Fstl1", "Cemip2", "Ust", "Abcc4", "Itgb5", "Ube2h",
  "Aig1", "Col6a1", "Pcdh7", "Map4k3", "Diaph2", "Cgn", "Acvr1", "Shtn1", "Hgf",
  "Ctps", "Mir99ahg", "Syt17", "Col5a2", "Senp7", "Map7", "Adgrl2", "Hipk2",
  "Epb41l2", "Myo18a", "Ide", "Dennd1b", "Nhsl1", "Lgr6", "Rin2", "Antxr2",
  "Il34", "Lrig1", "Ctnnb1", "Rbms1", "Kansl1l", "Herc3", "Dock1", "Greb1l",
  "2310002F09Rik", "Nova1", "Plekha7", "Mrps28", "Fndc3b", "Morrbid", "Col3a1",
  "Patj", "Krt7", "Rbfox1", "Dtd1", "Unc13b", "Mrps6", "2900026A02Rik", "Herc6",
  "Pag1", "Sh3pxd2b", "Acan"
)
genes_entrez <- bitr(regulon, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
enrich_go_result <- enrichGO(gene = genes_entrez$ENTREZID, 
                             OrgDb = org.Mm.eg.db, 
                             ont = "BP",  
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.05)

top_terms <- enrich_go_result %>%
  as.data.frame() %>%
  arrange(desc(GeneRatio), p.adjust) %>%
  slice_head(n = 15)


top_terms$GeneRatio <- sapply(strsplit(as.character(top_terms$GeneRatio), "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]))

p_go <- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                              color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
  scale_size_continuous(name = "Count") +
  labs(
    title = paste("GO BP enrichment for", TF, "'s regulon in", line, cell_type),
    x = "Gene Ratio", y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(family = "Times New Roman", face = "bold", size = 18),
    axis.text.y = element_text(family = "Times New Roman", size = 16, color="black"),
    axis.text.x = element_text(family = "Times New Roman", size = 16),
    axis.title = element_text(family = "Times New Roman", size = 16),
    plot.title = element_text(hjust = 1, family = "Times New Roman", size = 16,face = "bold"),
    legend.title = element_text(family = "Times New Roman", size = 16),
    legend.text = element_text(family = "Times New Roman", size = 14)
  )

ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","GO_enrich_",TF,".png"), plot = p_go, width = 10, height = 7, dpi = 600)



enrich_KEGG <- enrichKEGG(gene = genes_entrez$ENTREZID, 
                          organism = "mouse", 
                          pAdjustMethod = "BH", 
                          pvalueCutoff = 0.05)

top_terms <- enrich_KEGG %>%
  as.data.frame() %>%
  arrange(desc(GeneRatio), p.adjust) %>%
  slice_head(n = 15)

top_terms$GeneRatio <- sapply(strsplit(as.character(top_terms$GeneRatio), "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]))

p_k<- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                            color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
  scale_size_continuous(name = "Count") +
  labs(
    title = paste("KEGG enrichment for", TF, "'s regulon in", line, cell_type),
    x = "Gene Ratio", y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(family = "Times New Roman", face = "bold", size = 18),
    axis.text.y = element_text(family = "Times New Roman", size = 16, color = "black"),
    axis.text.x = element_text(family = "Times New Roman", size = 16),
    axis.title = element_text(family = "Times New Roman", size = 16),
    plot.title = element_text(hjust = 1, family = "Times New Roman", size = 16, face = "bold"),
    legend.title = element_text(family = "Times New Roman", size = 16),
    legend.text = element_text(family = "Times New Roman", size = 14)
  )
ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","KEGG_enrich_",TF,".png"), plot = p_k, width = 8, height = 7, dpi = 600)




enrich_pathway_result <- enrichPathway(gene = genes_entrez$ENTREZID, 
                                       organism = "mouse", 
                                       pAdjustMethod = "BH", 
                                       pvalueCutoff = 0.05)




top_terms <- enrich_pathway_result %>%
  as.data.frame() %>%
  arrange(desc(GeneRatio), p.adjust) %>%
  slice_head(n = 15)

top_terms$GeneRatio <- sapply(strsplit(as.character(top_terms$GeneRatio), "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]))

p_r<- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                            color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
  scale_size_continuous(name = "Count") +
  labs(
    title = paste("Reactome enrichment for", TF, "'s regulon in", line, cell_type),
    x = "Gene Ratio", y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(family = "Times New Roman", face = "bold", size = 18),
    axis.text.y = element_text(family = "Times New Roman", size = 16, color="black"),
    axis.text.x = element_text(family = "Times New Roman", size = 16),
    axis.title = element_text(family = "Times New Roman", size = 16),
    plot.title = element_text(hjust = 1, family = "Times New Roman", size = 16,face = "bold"),
    legend.title = element_text(family = "Times New Roman", size = 16),
    legend.text = element_text(family = "Times New Roman", size = 14))

ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","React_enrich_",TF,".png"), plot = p_r, width = 8, height = 7, dpi = 600)


combined_plot <- p_go + p_k + p_r + plot_layout(ncol = 3)

ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","Combined_enrich_",TF,".png"), plot = combined_plot, width = 25, height = 7, dpi = 600)




# KPB25L-UV ---------------------------------------------------------------


## Zfp827 ---------------------------------------------------------------------

line <- "KPB25L-UV"
TF <- "Zfp827"
cell_type <- "Cancer cells"

regulon <-  c(
  "Gdpd5", "Ntn1", "Pacs1", "Slc4a4", "Kifc3", "Pxdn", "Ppp1r12b", "Fendrr", "Synj2",
  "Gpr39", "Pde7a", "Fam107b", "Usp46", "Atp9a", "Slc10a7", "Fam168a", "Gm49797",
  "Nisch", "Kirrel", "Rbm47", "Il18r1", "Pak1", "Fer", "Lypd6b", "Cdcp1", "Vwa8",
  "Sh3d19", "Rin2", "Supt3", "Alcam", "Bicd1", "Fbxo32", "Sorcs2", "Fhl2", "Pced1b",
  "Masp1", "Rtp4", "Itpr3", "Prkd1", "Skap2", "Tra2b", "Atp11a", "Tspan5", "Znrf1",
  "Osbpl3", "Dock4", "Denn2b", "Enah", "Pkp1", "Arhgap31", "Atp13a3", "Wls",
  "Phf21a", "Hspd1", "Hipk2", "Rngtt", "Ipo7", "Itgb1", "Pappa", "Rap1gds1",
  "Hpcal1", "Nrip1", "Tenm1", "Pstpip2", "Nop58", "Slco3a1", "Rps6ka2", "Dtx4",
  "Ptprk", "Nol4", "Prickle1", "Gspt1", "Col18a1", "Lrba", "Gm2682", "Dock7",
  "Patj", "Rps6ka3", "Cadm1", "Rsu1", "Pak3", "Bcl2", "Chd9", "Igf2bp2", "Lrrk1",
  "Stau2", "Gm26652", "Ccbe1", "Fubp1", "Rims2", "Syt17", "Tenm4", "Arhgap24",
  "Ephb2", "Exoc6b", "Maml3", "Pak2", "Pdlim5", "Tbck", "St3gal"
)
genes_entrez <- bitr(regulon, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
enrich_go_result <- enrichGO(gene = genes_entrez$ENTREZID, 
                             OrgDb = org.Mm.eg.db, 
                             ont = "BP",  
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.05)

top_terms <- enrich_go_result %>%
  as.data.frame() %>%
  arrange(desc(GeneRatio), p.adjust) %>%
  slice_head(n = 15)


top_terms$GeneRatio <- sapply(strsplit(as.character(top_terms$GeneRatio), "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]))

p_go <- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                              color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
  scale_size_continuous(name = "Count") +
  labs(
    title = paste("GO BP enrichment for", TF, "'s regulon in", line, cell_type),
    x = "Gene Ratio", y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(family = "Times New Roman", face = "bold", size = 18),
    axis.text.y = element_text(family = "Times New Roman", size = 16, color="black"),
    axis.text.x = element_text(family = "Times New Roman", size = 16),
    axis.title = element_text(family = "Times New Roman", size = 16),
    plot.title = element_text(hjust = 1, family = "Times New Roman", size = 16,face = "bold"),
    legend.title = element_text(family = "Times New Roman", size = 16),
    legend.text = element_text(family = "Times New Roman", size = 14)
  )

ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","GO_enrich_",TF,".png"), plot = p_go, width = 10, height = 7, dpi = 600)



enrich_KEGG <- enrichKEGG(gene = genes_entrez$ENTREZID, 
                          organism = "mouse", 
                          pAdjustMethod = "BH", 
                          pvalueCutoff = 0.05)

top_terms <- enrich_KEGG %>%
  as.data.frame() %>%
  arrange(desc(GeneRatio), p.adjust) %>%
  slice_head(n = 15)

top_terms$GeneRatio <- sapply(strsplit(as.character(top_terms$GeneRatio), "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]))

p_k<- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                            color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
  scale_size_continuous(name = "Count") +
  labs(
    title = paste("KEGG enrichment for", TF, "'s regulon in", line, cell_type),
    x = "Gene Ratio", y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(family = "Times New Roman", face = "bold", size = 18),
    axis.text.y = element_text(family = "Times New Roman", size = 16, color = "black"),
    axis.text.x = element_text(family = "Times New Roman", size = 16),
    axis.title = element_text(family = "Times New Roman", size = 16),
    plot.title = element_text(hjust = 1, family = "Times New Roman", size = 16, face = "bold"),
    legend.title = element_text(family = "Times New Roman", size = 16),
    legend.text = element_text(family = "Times New Roman", size = 14)
  )
ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","KEGG_enrich_",TF,".png"), plot = p_k, width = 8, height = 7, dpi = 600)




enrich_pathway_result <- enrichPathway(gene = genes_entrez$ENTREZID, 
                                       organism = "mouse", 
                                       pAdjustMethod = "BH", 
                                       pvalueCutoff = 0.05)




top_terms <- enrich_pathway_result %>%
  as.data.frame() %>%
  arrange(desc(GeneRatio), p.adjust) %>%
  slice_head(n = 15)

top_terms$GeneRatio <- sapply(strsplit(as.character(top_terms$GeneRatio), "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]))

p_r<- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                            color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
  scale_size_continuous(name = "Count") +
  labs(
    title = paste("Reactome enrichment for", TF, "'s regulon in", line, cell_type),
    x = "Gene Ratio", y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(family = "Times New Roman", face = "bold", size = 18),
    axis.text.y = element_text(family = "Times New Roman", size = 16, color="black"),
    axis.text.x = element_text(family = "Times New Roman", size = 16),
    axis.title = element_text(family = "Times New Roman", size = 16),
    plot.title = element_text(hjust = 1, family = "Times New Roman", size = 16,face = "bold"),
    legend.title = element_text(family = "Times New Roman", size = 16),
    legend.text = element_text(family = "Times New Roman", size = 14))

ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","React_enrich_",TF,".png"), plot = p_r, width = 8, height = 7, dpi = 600)


combined_plot <- p_go + p_k + p_r + plot_layout(ncol = 3)

ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","Combined_enrich_",TF,".png"), plot = combined_plot, width = 25, height = 7, dpi = 600)


## Nfkb1 -------------------------------------------------------------------
line <- "KPB25L-UV"
TF <- "Nfkb1"
cell_type <- "Cancer cells"

regulon <-  c(
  "Ash1l", "Rap1gds1", "Mindy3", "Sec24b", "Usp6nl", "Gm14636", "Lypd6b", "Rabgap1",
  "Stard13", "Gsn", "Gapvd1", "Tiam2", "AU020206", "Magi3", "Ptbp2", "Dpysl3",
  "Ehmt1", "Fnbp1l", "Sh3d19", "Cd44", "Fubp1", "Dock5", "Parp8", "Styk1", "Trim24",
  "Antxr2", "Myof", "Cacna2d1", "Myo1e", "Map4k4", "Sh3kbp1", "Efl1", "Arhgef3",
  "Wwc1", "Pdlim5", "Traf3", "Tnfaip8", "Anxa1", "Hpcal1", "Cobll1", "Osbpl3",
  "Myo10", "Cdk6", "Mlkl", "Notch2", "Frrs1", "Fnbp1", "Strbp", "Ksr1", "Phip",
  "Ralgps1", "Cnot1", "Cdh13", "Ralgps2", "Wnk1", "Eepd1", "Asap2", "Inpp4b",
  "Cpe", "Sh3pxd2a", "Atp2b1", "Camk2d", "Tmem131", "Ppp3ca", "Igsf11", "Akap13",
  "Arhgap22", "Scn5a", "Lrp1", "Lamc2", "Klhl24", "Pik3cb", "Rasa2", "Acvr1",
  "Fndc3b", "Dennd1b", "Abtb2", "Nmt2", "Trim33", "Btaf1", "Gm12648", "Kif26b",
  "Epg5", "Eya1", "Hace1", "C9orf72", "Cd9", "Fn1", "Map3k5", "Samd4", "Hs2st1",
  "Nlk", "Ldlrad4", "Erbin", "Abi1", "Pstpip2", "Itga2", "Gm12610", "Tgfbr3",
  "Epb41l4"
)
genes_entrez <- bitr(regulon, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
enrich_go_result <- enrichGO(gene = genes_entrez$ENTREZID, 
                             OrgDb = org.Mm.eg.db, 
                             ont = "BP",  
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.05)

top_terms <- enrich_go_result %>%
  as.data.frame() %>%
  arrange(desc(GeneRatio), p.adjust) %>%
  slice_head(n = 15)


top_terms$GeneRatio <- sapply(strsplit(as.character(top_terms$GeneRatio), "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]))

p_go <- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                              color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
  scale_size_continuous(name = "Count") +
  labs(
    title = paste("GO BP enrichment for", TF, "'s regulon in", line, cell_type),
    x = "Gene Ratio", y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(family = "Times New Roman", face = "bold", size = 18),
    axis.text.y = element_text(family = "Times New Roman", size = 16, color="black"),
    axis.text.x = element_text(family = "Times New Roman", size = 16),
    axis.title = element_text(family = "Times New Roman", size = 16),
    plot.title = element_text(hjust = 1, family = "Times New Roman", size = 16,face = "bold"),
    legend.title = element_text(family = "Times New Roman", size = 16),
    legend.text = element_text(family = "Times New Roman", size = 14)
  )

ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","GO_enrich_",TF,".png"), plot = p_go, width = 10, height = 7, dpi = 600)



enrich_KEGG <- enrichKEGG(gene = genes_entrez$ENTREZID, 
                          organism = "mouse", 
                          pAdjustMethod = "BH", 
                          pvalueCutoff = 0.05)

top_terms <- enrich_KEGG %>%
  as.data.frame() %>%
  arrange(desc(GeneRatio), p.adjust) %>%
  slice_head(n = 15)

top_terms$GeneRatio <- sapply(strsplit(as.character(top_terms$GeneRatio), "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]))

p_k<- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                            color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
  scale_size_continuous(name = "Count") +
  labs(
    title = paste("KEGG enrichment for", TF, "'s regulon in", line, cell_type),
    x = "Gene Ratio", y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(family = "Times New Roman", face = "bold", size = 18),
    axis.text.y = element_text(family = "Times New Roman", size = 16, color = "black"),
    axis.text.x = element_text(family = "Times New Roman", size = 16),
    axis.title = element_text(family = "Times New Roman", size = 16),
    plot.title = element_text(hjust = 1, family = "Times New Roman", size = 16, face = "bold"),
    legend.title = element_text(family = "Times New Roman", size = 16),
    legend.text = element_text(family = "Times New Roman", size = 14)
  )
ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","KEGG_enrich_",TF,".png"), plot = p_k, width = 8, height = 7, dpi = 600)




enrich_pathway_result <- enrichPathway(gene = genes_entrez$ENTREZID, 
                                       organism = "mouse", 
                                       pAdjustMethod = "BH", 
                                       pvalueCutoff = 0.1)




top_terms <- enrich_pathway_result %>%
  as.data.frame() %>%
  arrange(desc(GeneRatio), p.adjust) %>%
  slice_head(n = 15)

top_terms$GeneRatio <- sapply(strsplit(as.character(top_terms$GeneRatio), "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]))

p_r<- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                            color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
  scale_size_continuous(name = "Count") +
  labs(
    title = paste("Reactome enrichment for", TF, "'s regulon in", line, cell_type),
    x = "Gene Ratio", y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(family = "Times New Roman", face = "bold", size = 18),
    axis.text.y = element_text(family = "Times New Roman", size = 16, color="black"),
    axis.text.x = element_text(family = "Times New Roman", size = 16),
    axis.title = element_text(family = "Times New Roman", size = 16),
    plot.title = element_text(hjust = 1, family = "Times New Roman", size = 16,face = "bold"),
    legend.title = element_text(family = "Times New Roman", size = 16),
    legend.text = element_text(family = "Times New Roman", size = 14))

ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","React_enrich_",TF,".png"), plot = p_r, width = 8, height = 7, dpi = 600)


combined_plot <- p_go + p_k + p_r + plot_layout(ncol = 3)

ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","Combined_enrich_",TF,".png"), plot = combined_plot, width = 25, height = 7, dpi = 600)



# NK ----------------------------------------------------------------------


## Rfx3 --------------------------------------------------------------------


line <- "NK"
TF <- "Rfx3"
cell_type <- "NK cells"

regulon <-  c(
  "Supt3", "Babam2", "Raly", "Mbd5", "Tbc1d5", "Crppa", "Ric8b", "Ptprd", "Ahnak", "Cbfa2t2",
  "Cul1", "Mib1", "Casp12", "Dleu2", "Cog5", "Myh9", "Robo1", "Nav3", "Lrp1", "Neo1",
  "Stim2", "Dst", "Acvr1", "Ssh2", "Col1a1", "Col12a1", "Kansl1l", "Slit3", "Gm4258", "Fap",
  "Neat1", "Srrm2", "Stard13", "Malat1", "Tnc", "Brinp3", "Abi3bp", "Col1a2", "Col3a1", "Kif26b",
  "Slit2", "Spon1", "Kcnma1", "Sema3a", "Pcsk5", "Tenm3", "Col8a1", "Cdh11", "Adam12", "Pde7b",
  "Clstn2", "Itgbl1", "Adgrl3", "Plcl1", "Tgfbr3", "Adamts6", "Palld", "Cacna1c", "Cald1", "Ncam1",
  "Fhod3", "Meg3", "Serpine1", "Prkg1", "Nrg1", "Nxn", "Thbs2", "Pdzrn3", "Hs6st2", "Adamts12",
  "Postn", "Rian", "Col7a1", "Col5a2", "Epha3", "Ror1", "Lama2", "Fn1", "Magi2", "Prickle1",
  "Auts2", "Piezo2", "Lrmda", "Mir99ahg", "Ifi205", "Slc8a1", "Adamts16", "Svep1", "Tln2", "Col6a3",
  "Svil", "Tpm1", "Col5a3", "Parp8", "Thbs1", "Pde10a", "Plcb1", "Nox4", "Igfbp7", "Col5a1"
)
genes_entrez <- bitr(regulon, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
enrich_go_result <- enrichGO(gene = genes_entrez$ENTREZID, 
                             OrgDb = org.Mm.eg.db, 
                             ont = "BP",  
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.05)

top_terms <- enrich_go_result %>%
  as.data.frame() %>%
  arrange(desc(GeneRatio), p.adjust) %>%
  slice_head(n = 15)


top_terms$GeneRatio <- sapply(strsplit(as.character(top_terms$GeneRatio), "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]))

p_go <- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                              color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
  scale_size_continuous(name = "Count") +
  labs(
    title = paste("GO BP enrichment for", TF, "'s regulon in", cell_type),
    x = "Gene Ratio", y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(family = "Times New Roman", face = "bold", size = 18),
    axis.text.y = element_text(family = "Times New Roman", size = 16, color="black"),
    axis.text.x = element_text(family = "Times New Roman", size = 16),
    axis.title = element_text(family = "Times New Roman", size = 16),
    plot.title = element_text(hjust = 1, family = "Times New Roman", size = 16,face = "bold"),
    legend.title = element_text(family = "Times New Roman", size = 16),
    legend.text = element_text(family = "Times New Roman", size = 14)
  )

ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","GO_enrich_",TF,".png"), plot = p_go, width = 10, height = 7, dpi = 600)



enrich_KEGG <- enrichKEGG(gene = genes_entrez$ENTREZID, 
                          organism = "mouse", 
                          pAdjustMethod = "BH", 
                          pvalueCutoff = 0.05)

top_terms <- enrich_KEGG %>%
  as.data.frame() %>%
  arrange(desc(GeneRatio), p.adjust) %>%
  slice_head(n = 15)

top_terms$GeneRatio <- sapply(strsplit(as.character(top_terms$GeneRatio), "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]))

p_k<- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                            color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
  scale_size_continuous(name = "Count") +
  labs(
    title = paste("KEGG enrichment for", TF, "'s regulon in", cell_type),
    x = "Gene Ratio", y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(family = "Times New Roman", face = "bold", size = 18),
    axis.text.y = element_text(family = "Times New Roman", size = 16, color = "black"),
    axis.text.x = element_text(family = "Times New Roman", size = 16),
    axis.title = element_text(family = "Times New Roman", size = 16),
    plot.title = element_text(hjust = 1, family = "Times New Roman", size = 16, face = "bold"),
    legend.title = element_text(family = "Times New Roman", size = 16),
    legend.text = element_text(family = "Times New Roman", size = 14)
  )
ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","KEGG_enrich_",TF,".png"), plot = p_k, width = 8, height = 7, dpi = 600)




enrich_pathway_result <- enrichPathway(gene = genes_entrez$ENTREZID, 
                                       organism = "mouse", 
                                       pAdjustMethod = "BH", 
                                       pvalueCutoff = 0.05)




top_terms <- enrich_pathway_result %>%
  as.data.frame() %>%
  arrange(desc(GeneRatio), p.adjust) %>%
  slice_head(n = 15)

top_terms$GeneRatio <- sapply(strsplit(as.character(top_terms$GeneRatio), "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]))

p_r<- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                            color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
  scale_size_continuous(name = "Count") +
  labs(
    title = paste("Reactome enrichment for", TF, "'s regulon in",cell_type),
    x = "Gene Ratio", y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(family = "Times New Roman", face = "bold", size = 18),
    axis.text.y = element_text(family = "Times New Roman", size = 16, color="black"),
    axis.text.x = element_text(family = "Times New Roman", size = 16),
    axis.title = element_text(family = "Times New Roman", size = 16),
    plot.title = element_text(hjust = 1, family = "Times New Roman", size = 16,face = "bold"),
    legend.title = element_text(family = "Times New Roman", size = 16),
    legend.text = element_text(family = "Times New Roman", size = 14))

ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","React_enrich_",TF,".png"), plot = p_r, width = 8, height = 7, dpi = 600)


combined_plot <- p_go + p_k + p_r + plot_layout(ncol = 3)

ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","Combined_enrich_",TF,".png"), plot = combined_plot, width = 30, height = 7, dpi = 600)



## Hoxc6 -------------------------------------------------------------------

line <- "NK"
TF <- "Hoxc6"
cell_type <- "NK cells"

regulon <- c(
  "Azin2", "Col1a1", "Fgfr2", "Col24a1", "Ttc28", "Cdh11", "Tnc", "Ncam1", "Prkd1", "Adamts9",
  "Neo1", "Itga11", "Col12a1", "Gabrb3", "Nuak1", "Tanc2", "Spats2l", "Malat1", "Brinp3", "Abi3bp",
  "Col1a2", "Col3a1", "Kif26b", "Slit2", "Spon1", "Kcnma1", "Sema3a", "Nav3", "Pcsk5", "Tenm3",
  "Col8a1", "Adam12", "Pde7b", "Clstn2", "Itgbl1", "Adgrl3", "Plcl1", "Tgfbr3", "Adamts6", "Palld",
  "Cacna1c", "Cald1", "Fhod3", "Meg3", "Serpine1", "Prkg1", "Slit3", "Nrg1", "Nxn", "Thbs2",
  "Pdzrn3", "Hs6st2", "Adamts12", "Postn", "Rian", "Col7a1", "Col5a2", "Epha3", "Ror1", "Lama2",
  "Fn1", "Magi2", "Prickle1", "Auts2", "Piezo2", "Lrmda", "Mir99ahg", "Ifi205", "Slc8a1", "Adamts16",
  "Svep1", "Tln2", "Col6a3", "Svil", "Tpm1", "Col5a3", "Parp8", "Thbs1", "Robo1", "Pde10a",
  "Plcb1", "Nox4", "Myh9", "Igfbp7", "Col5a1", "Igf1", "Dlc1", "Flnb", "Celf2", "Chl1",
  "Antxr1", "Angpt1", "Lpp", "Lhfp", "Fndc1", "Adam19", "Neat1", "Fbn2", "9530026P05Rik", "Ccbe1"
)


genes_entrez <- bitr(regulon, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
enrich_go_result <- enrichGO(gene = genes_entrez$ENTREZID, 
                             OrgDb = org.Mm.eg.db, 
                             ont = "BP",  
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.05)

top_terms <- enrich_go_result %>%
  as.data.frame() %>%
  arrange(desc(GeneRatio), p.adjust) %>%
  slice_head(n = 15)


top_terms$GeneRatio <- sapply(strsplit(as.character(top_terms$GeneRatio), "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]))

p_go <- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                              color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
  scale_size_continuous(name = "Count") +
  labs(
    title = paste("GO BP enrichment for", TF, "'s regulon in", cell_type),
    x = "Gene Ratio", y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(family = "Times New Roman", face = "bold", size = 18),
    axis.text.y = element_text(family = "Times New Roman", size = 16, color="black"),
    axis.text.x = element_text(family = "Times New Roman", size = 16),
    axis.title = element_text(family = "Times New Roman", size = 16),
    plot.title = element_text(hjust = 1, family = "Times New Roman", size = 16,face = "bold"),
    legend.title = element_text(family = "Times New Roman", size = 16),
    legend.text = element_text(family = "Times New Roman", size = 14)
  )

ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","GO_enrich_",TF,".png"), plot = p_go, width = 10, height = 7, dpi = 600)



enrich_KEGG <- enrichKEGG(gene = genes_entrez$ENTREZID, 
                          organism = "mouse", 
                          pAdjustMethod = "BH", 
                          pvalueCutoff = 0.05)

top_terms <- enrich_KEGG %>%
  as.data.frame() %>%
  arrange(desc(GeneRatio), p.adjust) %>%
  slice_head(n = 15)

top_terms$GeneRatio <- sapply(strsplit(as.character(top_terms$GeneRatio), "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]))

p_k<- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                            color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
  scale_size_continuous(name = "Count") +
  labs(
    title = paste("KEGG enrichment for", TF, "'s regulon in", cell_type),
    x = "Gene Ratio", y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(family = "Times New Roman", face = "bold", size = 18),
    axis.text.y = element_text(family = "Times New Roman", size = 16, color = "black"),
    axis.text.x = element_text(family = "Times New Roman", size = 16),
    axis.title = element_text(family = "Times New Roman", size = 16),
    plot.title = element_text(hjust = 1, family = "Times New Roman", size = 16, face = "bold"),
    legend.title = element_text(family = "Times New Roman", size = 16),
    legend.text = element_text(family = "Times New Roman", size = 14)
  )
ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","KEGG_enrich_",TF,".png"), plot = p_k, width = 8, height = 7, dpi = 600)




enrich_pathway_result <- enrichPathway(gene = genes_entrez$ENTREZID, 
                                       organism = "mouse", 
                                       pAdjustMethod = "BH", 
                                       pvalueCutoff = 0.05)




top_terms <- enrich_pathway_result %>%
  as.data.frame() %>%
  arrange(desc(GeneRatio), p.adjust) %>%
  slice_head(n = 15)

top_terms$GeneRatio <- sapply(strsplit(as.character(top_terms$GeneRatio), "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]))

p_r<- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                            color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
  scale_size_continuous(name = "Count") +
  labs(
    title = paste("Reactome enrichment for", TF, "'s regulon in", cell_type),
    x = "Gene Ratio", y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(family = "Times New Roman", face = "bold", size = 18),
    axis.text.y = element_text(family = "Times New Roman", size = 16, color="black"),
    axis.text.x = element_text(family = "Times New Roman", size = 16),
    axis.title = element_text(family = "Times New Roman", size = 16),
    plot.title = element_text(hjust = 1, family = "Times New Roman", size = 16,face = "bold"),
    legend.title = element_text(family = "Times New Roman", size = 16),
    legend.text = element_text(family = "Times New Roman", size = 14))

ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","React_enrich_",TF,".png"), plot = p_r, width = 8, height = 7, dpi = 600)


combined_plot <- p_go + p_k + p_r + plot_layout(ncol = 3)

ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","Combined_enrich_",TF,".png"), plot = combined_plot, width = 30, height = 7, dpi = 600)



# MACROPHAGES -------------------------------------------------------------


## Mef2c -------------------------------------------------------------------


line <- "Macrophages"
TF <- "Mef2c"
cell_type <- "Macrophages"

regulon <- c(
  "Arhgap4", "Mgat4a", "Cx3cr1", "Rhoh", "Lair1", "Sesn1", "Ly86", "Plec", "Itga6", "Sgk1",
  "A630001G21Rik", "Il10rb", "Serinc3", "Ms4a6c", "Gm5086", "Col14a1", "Slc16a6", "Ntpcr", "Asap1", "Mbnl1",
  "Fgd2", "Sla", "Fap", "Tmem273", "Wasf2", "Extl3", "Sgms1", "Disc1", "Exoc6", "Ptgs1",
  "Lifr", "Relch", "Clmp", "Lpcat2", "Mvb12b", "Tnfaip8", "Cd38", "Ifnar2", "F630028O10Rik", "Agps",
  "Ctsc", "Gatm", "Rap1gap2", "Snx2", "Gm35154", "Apbb1ip", "Slc7a11", "Dab2", "Chst11", "Mgat5",
  "Gnaq", "Cttnbp2nl", "Phc2", "Rnase4", "Actn1", "Ppp1r21", "Tes", "C1qb", "Phf14", "P2ry12",
  "Xdh", "Pla2g4a", "Slc12a6", "Cd9", "Hck", "Blnk", "Akap10", "Peli2", "Cables1", "Ssh2",
  "Slco2b1", "Sptbn1", "Stk17b", "Phlpp1", "Otulinl", "Tgfbr2", "Abi3", "Rapgef5", "Rab3il1", "Lyn",
  "Cd33", "Fcrls", "Gdi2", "Dock4", "Bmp2k", "2610203C22Rik", "Map3k5", "Shtn1", "Adam10", "Pag1",
  "Ttc28", "Rap1gds1", "Epb41l1", "Gm4258", "Ahnak", "Pde3b", "Vim", "Elmo1", "Itpr2", "Bcl2"
)

genes_entrez <- bitr(regulon, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
enrich_go_result <- enrichGO(gene = genes_entrez$ENTREZID, 
                             OrgDb = org.Mm.eg.db, 
                             ont = "BP",  
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.05)

top_terms <- enrich_go_result %>%
  as.data.frame() %>%
  arrange(desc(GeneRatio), p.adjust) %>%
  slice_head(n = 15)


top_terms$GeneRatio <- sapply(strsplit(as.character(top_terms$GeneRatio), "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]))

p_go <- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                              color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
  scale_size_continuous(name = "Count") +
  labs(
    title = paste("GO BP enrichment for", TF, "'s regulon in", cell_type),
    x = "Gene Ratio", y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(family = "Times New Roman", face = "bold", size = 18),
    axis.text.y = element_text(family = "Times New Roman", size = 16, color="black"),
    axis.text.x = element_text(family = "Times New Roman", size = 16),
    axis.title = element_text(family = "Times New Roman", size = 16),
    plot.title = element_text(hjust = 1, family = "Times New Roman", size = 16,face = "bold"),
    legend.title = element_text(family = "Times New Roman", size = 16),
    legend.text = element_text(family = "Times New Roman", size = 14)
  )

ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","GO_enrich_",TF,".png"), plot = p_go, width = 10, height = 7, dpi = 600)



enrich_KEGG <- enrichKEGG(gene = genes_entrez$ENTREZID, 
                          organism = "mouse", 
                          pAdjustMethod = "BH", 
                          pvalueCutoff = 0.05)

top_terms <- enrich_KEGG %>%
  as.data.frame() %>%
  arrange(desc(GeneRatio), p.adjust) %>%
  slice_head(n = 15)

top_terms$GeneRatio <- sapply(strsplit(as.character(top_terms$GeneRatio), "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]))

p_k<- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                            color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
  scale_size_continuous(name = "Count") +
  labs(
    title = paste("KEGG enrichment for", TF, "'s regulon in", cell_type),
    x = "Gene Ratio", y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(family = "Times New Roman", face = "bold", size = 18),
    axis.text.y = element_text(family = "Times New Roman", size = 16, color = "black"),
    axis.text.x = element_text(family = "Times New Roman", size = 16),
    axis.title = element_text(family = "Times New Roman", size = 16),
    plot.title = element_text(hjust = 1, family = "Times New Roman", size = 16, face = "bold"),
    legend.title = element_text(family = "Times New Roman", size = 16),
    legend.text = element_text(family = "Times New Roman", size = 14)
  )
ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","KEGG_enrich_",TF,".png"), plot = p_k, width = 8, height = 7, dpi = 600)




enrich_pathway_result <- enrichPathway(gene = genes_entrez$ENTREZID, 
                                       organism = "mouse", 
                                       pAdjustMethod = "BH", 
                                       pvalueCutoff = 0.1)




top_terms <- enrich_pathway_result %>%
  as.data.frame() %>%
  arrange(desc(GeneRatio), p.adjust) %>%
  slice_head(n = 15)

top_terms$GeneRatio <- sapply(strsplit(as.character(top_terms$GeneRatio), "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]))

p_r<- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                            color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
  scale_size_continuous(name = "Count") +
  labs(
    title = paste("Reactome enrichment for", TF, "'s regulon in", cell_type),
    x = "Gene Ratio", y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(family = "Times New Roman", face = "bold", size = 18),
    axis.text.y = element_text(family = "Times New Roman", size = 16, color="black"),
    axis.text.x = element_text(family = "Times New Roman", size = 16),
    axis.title = element_text(family = "Times New Roman", size = 16),
    plot.title = element_text(hjust = 1, family = "Times New Roman", size = 16,face = "bold"),
    legend.title = element_text(family = "Times New Roman", size = 16),
    legend.text = element_text(family = "Times New Roman", size = 14))

ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","React_enrich_",TF,".png"), plot = p_r, width = 8, height = 7, dpi = 600)


combined_plot <- p_go + p_k + plot_layout(ncol = 2)

ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","Combined_enrich_",TF,".png"), plot = combined_plot, width = 20, height = 7, dpi = 600)


## Pias1 -------------------------------------------------------------------



line <- "Macrophages"
TF <- "Pias1"
cell_type <- "Macrophages"

regulon <-  c(
  "Sesn1", "Tab2", "Serinc3", "Dyrk1a", "Slc12a6", "Aftph", "Dleu2", "Ap3b1", "Pan3", "Rnf169",
  "Relch", "Il6ra", "Mtmr3", "Peli1", "Scarb2", "Nedd4l", "Slc16a6", "Rassf3", "Ankfy1", "Dclre1c",
  "Vps54", "Ncoa3", "Cask", "Ptpro", "Ptbp3", "Agps", "Mon2", "Stk17b", "Hnrnpa2b1", "Galnt1",
  "Gdi2", "Pard3b", "Nek7", "Fnip1", "Itm2b", "Tle4", "Lrp1", "Fgr", "Slc30a7", "Ifngr1",
  "Wdr26", "Mki67", "Smurf2", "Tmem164", "Macf1", "Cltc", "Actr3", "Cenpe", "Capza2", "Lrp6",
  "Anln", "Clasp2", "Kdm6a", "Abi1", "Gapvd1", "Dip2b", "Actr2", "Cfh", "Rnf149", "Uhrf2",
  "Acer3", "Gtdc1", "Atf7ip", "Cbl", "Slc38a6", "Top2a", "Epsti1", "Plcl2", "Ehd4", "Knl1",
  "Kif15", "Fn1", "Malat1", "F13a1", "Mrc1", "Plxdc2", "Frmd4b", "Ophn1", "Psd3", "Tgfbr1",
  "Slc8a1", "Neat1", "Mmp12", "Alcam", "Slc9a9", "Abca1", "Pde7b", "Ms4a7", "Arhgap24", "Nrp1",
  "Tmcc3", "Pdgfc", "Stab1", "Aoah", "Zdhhc14", "Junos", "Pde4d", "Atp11b", "Fmnl2", "Itga9"
)

genes_entrez <- bitr(regulon, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
enrich_go_result <- enrichGO(gene = genes_entrez$ENTREZID, 
                             OrgDb = org.Mm.eg.db, 
                             ont = "BP",  
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.05)

top_terms <- enrich_go_result %>%
  as.data.frame() %>%
  arrange(desc(GeneRatio), p.adjust) %>%
  slice_head(n = 15)


top_terms$GeneRatio <- sapply(strsplit(as.character(top_terms$GeneRatio), "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]))

p_go <- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                              color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
  scale_size_continuous(name = "Count") +
  labs(
    title = paste("GO BP enrichment for", TF, "'s regulon in", cell_type),
    x = "Gene Ratio", y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(family = "Times New Roman", face = "bold", size = 18),
    axis.text.y = element_text(family = "Times New Roman", size = 16, color="black"),
    axis.text.x = element_text(family = "Times New Roman", size = 16),
    axis.title = element_text(family = "Times New Roman", size = 16),
    plot.title = element_text(hjust = 1, family = "Times New Roman", size = 16,face = "bold"),
    legend.title = element_text(family = "Times New Roman", size = 16),
    legend.text = element_text(family = "Times New Roman", size = 14)
  )

ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","GO_enrich_",TF,".png"), plot = p_go, width = 10, height = 7, dpi = 600)



enrich_KEGG <- enrichKEGG(gene = genes_entrez$ENTREZID, 
                          organism = "mouse", 
                          pAdjustMethod = "BH", 
                          pvalueCutoff = 0.05)

top_terms <- enrich_KEGG %>%
  as.data.frame() %>%
  arrange(desc(GeneRatio), p.adjust) %>%
  slice_head(n = 15)

top_terms$GeneRatio <- sapply(strsplit(as.character(top_terms$GeneRatio), "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]))

p_k<- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                            color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
  scale_size_continuous(name = "Count") +
  labs(
    title = paste("KEGG enrichment for", TF, "'s regulon in", cell_type),
    x = "Gene Ratio", y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(family = "Times New Roman", face = "bold", size = 18),
    axis.text.y = element_text(family = "Times New Roman", size = 16, color = "black"),
    axis.text.x = element_text(family = "Times New Roman", size = 16),
    axis.title = element_text(family = "Times New Roman", size = 16),
    plot.title = element_text(hjust = 1, family = "Times New Roman", size = 16, face = "bold"),
    legend.title = element_text(family = "Times New Roman", size = 16),
    legend.text = element_text(family = "Times New Roman", size = 14)
  )
ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","KEGG_enrich_",TF,".png"), plot = p_k, width = 8, height = 7, dpi = 600)




enrich_pathway_result <- enrichPathway(gene = genes_entrez$ENTREZID, 
                                       organism = "mouse", 
                                       pAdjustMethod = "BH", 
                                       pvalueCutoff = 0.05)




top_terms <- enrich_pathway_result %>%
  as.data.frame() %>%
  arrange(desc(GeneRatio), p.adjust) %>%
  slice_head(n = 15)

top_terms$GeneRatio <- sapply(strsplit(as.character(top_terms$GeneRatio), "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]))

p_r<- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                            color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
  scale_size_continuous(name = "Count") +
  labs(
    title = paste("Reactome enrichment for", TF, "'s regulon in", cell_type),
    x = "Gene Ratio", y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(family = "Times New Roman", face = "bold", size = 18),
    axis.text.y = element_text(family = "Times New Roman", size = 16, color="black"),
    axis.text.x = element_text(family = "Times New Roman", size = 16),
    axis.title = element_text(family = "Times New Roman", size = 16),
    plot.title = element_text(hjust = 1, family = "Times New Roman", size = 16,face = "bold"),
    legend.title = element_text(family = "Times New Roman", size = 16),
    legend.text = element_text(family = "Times New Roman", size = 14))

ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","React_enrich_",TF,".png"), plot = p_r, width = 8, height = 7, dpi = 600)


combined_plot <- p_go + p_k + plot_layout(ncol = 2)

ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","Combined_enrich_",TF,".png"), plot = combined_plot, width = 20, height = 7, dpi = 600)



## ZFP292 ------------------------------------------------------------------



line <- "Macrophages"
TF <- "Zfp292"
cell_type <- "Macrophages"

regulon <- c(
  "Setd5", "Hk2", "Npepps", "Pkm", "P4ha1", "Irak2", "Vegfa", "R3hdm1", "Luc7l2", "Kdm3a",
  "Dock7", "Gpi1", "Cep170", "Rnf217", "Rrbp1", "Fto", "Bnip3", "Asph", "Eno1", "Ogt",
  "Map2k1", "Chd7", "Vps13d", "Arrb1", "Epb41", "Pcnx", "P4hb", "Map4k4", "Nisch", "Rbms1",
  "Nampt", "Ecpas", "Fndc3b", "Picalm", "Ranbp2", "Ndrg1", "Maml3", "Birc6", "Malat1", "Ankrd11",
  "Ctnnb1", "Qki", "Gys1", "Plekha2", "Fth1", "Sbf2", "Hipk2", "Spidr", "Cyrib", "Eno2",
  "Cmip", "Trip12", "Iqgap1", "Atp2b1", "Slc9a9", "Arhgap24", "Tpd52", "Neat1", "Epb41l2", "Trrap",
  "Ero1a", "Gsn", "Mthfd1l", "Afg3l2", "Aoah", "Slc2a1", "ENSMUSG00000120385", "Spag9", "Plcl2", "Soat1",
  "Tmem87b", "Ehd4", "Iqgap2", "Pgm1", "Appl2", "Jmjd1c", "Gbe1", "Tbc1d4", "Tes", "Fnbp1",
  "Btg1", "Rap1a", "Pitpnc1", "Lpp", "Nsd3", "Fcgr2b", "Lnpep", "Lrp1", "Pde7b", "Junos",
  "Mtmr3", "Utrn", "Igf1r", "Camk1", "Ric1", "Atp6v0a1", "Tnrc6b", "Itpk1", "Dennd1a", "Reps"
)

genes_entrez <- bitr(regulon, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
enrich_go_result <- enrichGO(gene = genes_entrez$ENTREZID, 
                             OrgDb = org.Mm.eg.db, 
                             ont = "BP",  
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.05)

top_terms <- enrich_go_result %>%
  as.data.frame() %>%
  arrange(desc(GeneRatio), p.adjust) %>%
  slice_head(n = 15)


top_terms$GeneRatio <- sapply(strsplit(as.character(top_terms$GeneRatio), "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]))

p_go <- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                              color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
  scale_size_continuous(name = "Count") +
  labs(
    title = paste("GO BP enrichment for", TF, "'s regulon in", cell_type),
    x = "Gene Ratio", y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(family = "Times New Roman", face = "bold", size = 18),
    axis.text.y = element_text(family = "Times New Roman", size = 16, color="black"),
    axis.text.x = element_text(family = "Times New Roman", size = 16),
    axis.title = element_text(family = "Times New Roman", size = 16),
    plot.title = element_text(hjust = 1, family = "Times New Roman", size = 16,face = "bold"),
    legend.title = element_text(family = "Times New Roman", size = 16),
    legend.text = element_text(family = "Times New Roman", size = 14)
  )

ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","GO_enrich_",TF,".png"), plot = p_go, width = 10, height = 7, dpi = 600)



enrich_KEGG <- enrichKEGG(gene = genes_entrez$ENTREZID, 
                          organism = "mouse", 
                          pAdjustMethod = "BH", 
                          pvalueCutoff = 0.05)

top_terms <- enrich_KEGG %>%
  as.data.frame() %>%
  arrange(desc(GeneRatio), p.adjust) %>%
  slice_head(n = 15)

top_terms$GeneRatio <- sapply(strsplit(as.character(top_terms$GeneRatio), "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]))

p_k<- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                            color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
  scale_size_continuous(name = "Count") +
  labs(
    title = paste("KEGG enrichment for", TF, "'s regulon in", cell_type),
    x = "Gene Ratio", y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(family = "Times New Roman", face = "bold", size = 18),
    axis.text.y = element_text(family = "Times New Roman", size = 16, color = "black"),
    axis.text.x = element_text(family = "Times New Roman", size = 16),
    axis.title = element_text(family = "Times New Roman", size = 16),
    plot.title = element_text(hjust = 1, family = "Times New Roman", size = 16, face = "bold"),
    legend.title = element_text(family = "Times New Roman", size = 16),
    legend.text = element_text(family = "Times New Roman", size = 14)
  )
ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","KEGG_enrich_",TF,".png"), plot = p_k, width = 8, height = 7, dpi = 600)




enrich_pathway_result <- enrichPathway(gene = genes_entrez$ENTREZID, 
                                       organism = "mouse", 
                                       pAdjustMethod = "BH", 
                                       pvalueCutoff = 0.05)




top_terms <- enrich_pathway_result %>%
  as.data.frame() %>%
  arrange(desc(GeneRatio), p.adjust) %>%
  slice_head(n = 15)

top_terms$GeneRatio <- sapply(strsplit(as.character(top_terms$GeneRatio), "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]))

p_r<- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                            color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "lightcoral", high = "darkblue", name = "p.adjust") +
  scale_size_continuous(name = "Count") +
  labs(
    title = paste("Reactome enrichment for", TF, "'s regulon in", cell_type),
    x = "Gene Ratio", y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(family = "Times New Roman", face = "bold", size = 18),
    axis.text.y = element_text(family = "Times New Roman", size = 16, color="black"),
    axis.text.x = element_text(family = "Times New Roman", size = 16),
    axis.title = element_text(family = "Times New Roman", size = 16),
    plot.title = element_text(hjust = 1, family = "Times New Roman", size = 16,face = "bold"),
    legend.title = element_text(family = "Times New Roman", size = 16),
    legend.text = element_text(family = "Times New Roman", size = 14))

ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","React_enrich_",TF,".png"), plot = p_r, width = 8, height = 7, dpi = 600)


combined_plot <- p_go + p_k + p_r + plot_layout(ncol = 3)

ggsave(filename= paste0("/home/workdir/HDAC_tfm/DAC_anti-PDL1_tfm/fig/PEA_SimiC/",line,"/","Combined_enrich_",TF,".png"), plot = combined_plot, width = 30, height = 7, dpi = 600)


