################################################################################
# Title: Running CellChat2 on our data
# Author: Paula
# Description: Script to run CellChat analysis on our data comparing between treats and cell lines
################################################################################

# Install CellChat2 from GitHub (commented out because you don't usually reinstall every run)
# remotes::install_github("SiYangming/CellChat2")

# Load required libraries
library(CellChat)     # Main package for cell-cell communication inference
library(patchwork)    # For combining ggplot-based visualizations
options(stringsAsFactors = FALSE)  # Avoid automatic factor conversion

# Load the mouse CellChat database (contains ligand-receptor interaction information)
data("CellChatDB.mouse")

# Track start time for performance monitoring
ptm = Sys.time()

# Load preprocessed Seurat object with annotated immune cells
cells <- readRDS("/home/workdir/analysis_augie/20250120_labeled_immune.combined_edit.RDS")

# Merge layers (if data was stored in multiple layers such as RNA/ADT)
cells_joined <- JoinLayers(cells)

# Replace "Unannotated" cells with "Cancer cells" label
cells_joined@meta.data$immuno[cells_joined@meta.data$immuno == "Unannotated"] <- "Cancer cells"

# Split the 'sample' column into separate columns: cell line and treatment
split_info <- do.call(rbind, strsplit(as.character(cells_joined@meta.data$sample), "_"))
cells_joined@meta.data$cell_line <- split_info[, 1]
cells_joined@meta.data$treat     <- split_info[, 2]

# Get unique cell lines and treatments present in the dataset
cell_lines <- unique(cells_joined@meta.data$cell_line)
treats     <- unique(cells_joined@meta.data$treat)

# For this run, focus only on one specific cell line and treatment
cell_lines <- c("KAP25L")
treats     <- c("Combination")

# Loop through chosen combinations of cell line and treatment
for (cell_line in cell_lines) {
  for (treat in treats) {
    
    cat("Processing:", cell_line, "+", treat, "\n")
    
    # Select cells that match current cell line and treatment
    selected_cells <- rownames(cells_joined@meta.data)[
      cells_joined@meta.data$cell_line == cell_line &
      cells_joined@meta.data$treat     == treat
    ]
    
    # Extract normalized RNA expression matrix for selected cells
    data.input <- cells_joined[["RNA"]]$data[, selected_cells]
    
    # Prepare metadata for CellChat (labels = cell type, sample ID)
    meta_subset <- cells_joined@meta.data[selected_cells, c("immuno", "sample")]
    colnames(meta_subset) <- c("labels", "sample")
    
    # Create CellChat object using selected data
    cellchat <- createCellChat(object = data.input, meta = meta_subset, group.by = "labels")
    
    # Use only secreted signaling interactions from the mouse database
    CellChatDB.use <- subsetDB(CellChatDB.mouse, search = "Secreted Signaling", key = "annotation")
    cellchat@DB <- CellChatDB.use
    
    # Pre-filter the expression data for faster processing
    cellchat <- subsetData(cellchat)
    
    # Enable parallel processing with 12 workers
    future::plan("multisession", workers = 12)
    
    # Identify overexpressed genes and ligand-receptor interactions
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    
    # Record time just before projection
    ptm = Sys.time()
    
    # Project mouse interactions into the human protein-protein interaction network
    cellchat <- projectData(cellchat, PPI.human)
    
    # Compute communication probabilities between cell types
    cellchat <- computeCommunProb(cellchat, type = "triMean", population.size = TRUE, raw.use = FALSE)
    
    # Filter communications: keep interactions where each cell type has ≥ 1 cell
    # (Other treatments may have been filtered with min.cells = 10)
    cellchat <- filterCommunication(cellchat, min.cells = 1)
    
    # Compute communication probabilities at the pathway level
    cellchat <- computeCommunProbPathway(cellchat)
    
    # Aggregate the network into a more compact representation
    cellchat <- aggregateNet(cellchat)
    
    # Calculate total execution time for this combination
    execution.time = Sys.time() - ptm
    print(paste("Time for", cell_line, "+", treat, ":",
                round(as.numeric(execution.time, units = "secs"), 2), "seconds"))
    
    # Save results to file
    saveRDS(cellchat, file = paste0("/home/workdir/HDAC_tfm/data/cellchat/cellchat_SS_30_05_",
                                    cell_line, "_", treat, "_min_cells_5.rds"))
  }
}

# -------------------------------------------------------------------------
# SECOND RUN: Full secreted signaling analysis without subset to cell line/treat
# -------------------------------------------------------------------------

# Select secreted signaling interactions from the entire database
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")
cellchat@DB <- CellChatDB.use

# Subset data according to selected database
cellchat <- subsetData(cellchat)

# Parallelize with 4 workers (fewer cores for this run)
future::plan("multisession", workers = 4)

# Identify overexpressed genes and interactions
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# This prints: "The number of highly variable ligand-receptor pairs used for signaling inference is 692"

# Measure execution time for this stage
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

# Project to human PPI network
cellchat <- projectData(cellchat, PPI.human)

# Compute communication probabilities
ptm = Sys.time()
cellchat <- computeCommunProb(cellchat, type = "triMean", population.size = TRUE, raw.use = FALSE)

# Save intermediate results
saveRDS(cellchat, file = "/home/workdir/HDAC_tfm/data/cellchat/cellchat_SS_30_05.rds")

# Filter out low-cell-count interactions (≥ 10 cells per group)
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Compute pathway-level probabilities
cellchat <- computeCommunProbPathway(cellchat)

# Aggregate network into summary form
cellchat <- aggregateNet(cellchat)

# Measure final execution time
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))