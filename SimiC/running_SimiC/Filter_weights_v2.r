# --------------------------------------------------------------------
# Title: Filter_weights.R
# Author: Irene Marín
# Date: 05/02/2025
# Version: v2 (con path como único argumento)
# Description: Script to filter weights from SimicLASSO_op
# Usage: Called from Python with a single path argument
# -----------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Instructions:
# 1. Load libraries
# 2. Set Run paths and parameters
# 3. Run script
# 4. Save output weight matrix
# ------------------------------------------------------------------------------

library(reticulate)

# Definir el intérprete de Python
reticulate::use_python("/usr/bin/python3")

# Leer el único argumento (ruta de la matriz de pesos)
args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 1) {
    stop("Error: Se requiere un argumento (ruta del archivo de pesos)")
}

# Ruta del archivo de pesos (recibido desde Python)
weights_file <- args[1]

# Extraer información desde el path
path_parts <- unlist(strsplit(weights_file, "/"))
file_name <- path_parts[length(path_parts)]
path_dir <- paste(path_parts[1:(length(path_parts)-1)], collapse="/")

# Generar nombres de salida dinámicamente
new_fl_name <- sub(".pickle", "_filtered_BIC.pickle", weights_file)
fig_name <- sub(".pickle", "_Filtered_Weights_BIC.pdf", weights_file)

# Cargar los pesos desde Python
weights <- py_load_object(weights_file)
print(paste("Processing:", weights_file))

# Filtrado de pesos
pdf(fig_name)

for (phenotype in names(weights$weight_dic)) {
    print(phenotype)
    max_l <- c()
    all_data <- as.data.frame(weights$weight_dic[[phenotype]])
    colnames(all_data) <- weights[['query_targets']]
    all_data <- all_data[, colSums(abs(all_data)) > 0]

    bottom <- all_data[101, ]
    all_data <- as.data.frame(all_data[-nrow(all_data), ])
    rownames(all_data) <- weights[['TF_ids']]
    tmp <- scale(all_data, center = FALSE)

    for (target in colnames(all_data)) {
        all_tf <- rownames(tmp)[order(abs(tmp[, target]), decreasing = TRUE)]
        l <- 1
        tf_2_keep <- all_tf[1:l]
        while (sum(tmp[tf_2_keep, target]^2) / sum(tmp[, target]^2) < 0.9) {
            l <- l + 1
            tf_2_keep <- all_tf[1:l]
        }
        max_l <- c(max_l, l)
        tmp[!rownames(tmp) %in% tf_2_keep, target] <- 0
    }

    tmp <- rbind(tmp, bottom)
    tmp <- as.matrix(tmp)
    weights$weight_dic[[phenotype]] <- tmp

    hist(max_l, breaks = 20, main = paste0("Histogram of TFs with filtered weight for ", file_name), col = 'darkgrey')
}

# Guardar la nueva matriz de pesos filtrados
reticulate::py_save_object(weights, filename = new_fl_name)
dev.off()
print(paste("Filtered file saved:", new_fl_name))
