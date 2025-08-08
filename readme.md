# Machine learning models to elucidate the transcriptional dynamics in TNBC at single-cell resolution #


This repository contains the **pipeline** for the analysis of **Triple Negative Breast Cancer (TNBC)** at the single-cell level, combining **machine learning** techniques with transcriptomic analysis to unravel cellular and molecular dynamics.  
The workflow integrates from initial exploratory visualization to predictive modeling based on *single-cell RNA-seq*.

This framework is based on:

![](/fig/md_fig.png)


## Step-by-Step Description

1. **0_UMAP_starting_point** – Initial UMAP of clustered and annotated cells with scType.
2. **1_cell_type_proportions** – Estimation and comparison of cell type proportions across conditions and cell lines.
3. **2_DGEA** – Identification of differentially expressed genes (DEGs).
4. **3_PEA** – Functional pathway enrichment analysis GO BP, KEGG and Reactome.
5. **4_reannotation** – Reannotation of clusters based on marker genes shared between tumor and unannotated cells.
6. **5_subclustering** – Discovery of subpopulations within tumor cells.
7. **6_gene_signatures** – Evaluation of gene signatures associated with TNBC subtypes and immune related ones in TME.
8. **7_Liana** – Modeling of cell-cell communication using Liana.
9. **8_CellChat** – Cell-cell communication network analysis using CellChat.
10. **9_cell_cycle** – Determination of cell cycle phases and their link to tumor phenotype.
**SIMIC** - computational framework designed to infer regulatory dynamics of genes from scRNA-Seq data
---


##  Requirements

Docker installed with client version 23.0.1 or higher, and preferably server version 28.1.1 or higher.

Docker daemon running and user permissions to execute Docker commands.

Ports 8080 and 8787 free or modifiable as needed to publish services.

Git installed to clone the repository.

In configs folder, a dockerfile for liana can be found. See 2. Seting up environments

## Istalation

### 1. Clone repository
```
git clone https://github.com/paudbr/DAC_anti-PDL1_tfm.git
cd DAC_anti-PDL1_tfm
```

### 2. Setting up environments 

#### R DOCKER USED for Seurat and the rest of R packages

```
docker pull pdeblas/r_seurat:5.0.1

docker run -it -p 8889:8787 --name x pdeblas/r_seurat:5.0.1 /bin/bash

```

#### SIMIC DOCKER

```
docker pull guisesanz/simic:1.1.0

docker run -it -p 8080:8080 --name guisesanz/simic:2.1.0 /bin/bash

```

####LIANA DOCKER
```
cd configs/liana-py

docker build -t [your tag] f liana-py .

```

## Contact

For questions or support, please contact Paula de Blas Rioja (pau.dbr2002@gmail.com) or open an issue :)!