# ------------------------------------------------------------------------------
# Title: MACIG.R
# Author: Irene
# Date: 22/01/2025
# Version: v1
# Description: Script to perform RMagic imputation on the raw data from HDAC_PDL1 project
#
# ------------------------------------------------------------------------------
# Script to run MAGIC imputation on the raw data
library(reticulate)
mat_path <- "./workdir/DGE_filtered"

reticulate::use_python("/usr/bin/python3")
py_run_string("import magic")
# reticulate::py_install("anndata")
py_discover_config("magic")
library(Rmagic)
pymagic_is_available()
# Load the raw data
loaded_sparse_matrix <- Matrix::readMM(paste0(mat_path,"/HDAC_PDL1_seurat_FILTERED_metadata_matrix.txt"))
# Attach the row and column names to the loaded matrix
rownames(loaded_sparse_matrix) <- read.table(paste0(mat_path,"/FILTERED_genes.csv"), header = FALSE, sep = ",")$V1
colnames(loaded_sparse_matrix) <- read.table(paste0(mat_path,"/FILTERED_cells.csv"), header = FALSE, sep = ",")$V1
str(loaded_sparse_matrix)


### RUN MAGIC #### #MAGIC works on a cells x genes matrix, seurat gives a genes x cells matrix
norm_counts <-Rmagic::library.size.normalize(t(loaded_sparse_matrix))
norm_counts <- sqrt(norm_counts)
data_MAGIC <- Rmagic::magic(norm_counts,genes='all_genes') 
data_MAGIC2 <- Rmagic::magic(norm_counts,genes='all_genes',solver='approximate') 

data_MAGIC_df <- as.data.frame(data_MAGIC)
data_MAGIC_df2 <- as.data.frame(data_MAGIC22)
saveRDS(data_MAGIC_df, paste0(mat_path,"/Magic_imputed_data.rds"))
saveRDS(data_MAGIC_df2, paste0(mat_path,"/Magic_aprox_imputed_data.rds"))