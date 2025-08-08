generar_dotplot <- function(df, linea_celular, cell_type, modo = c("source", "target"), output_dir = "./") {
  modo <- match.arg(modo)
  
  # Filtrado
  datos_filtrados <- df %>%
    filter(linea_celular == linea_celular) %>%
    filter((!!sym(modo)) == cell_type)
  
  # Extraer los 15 primeros por tratamiento
  top15_por_tratamiento <- datos_filtrados %>%
    group_by(tratamiento) %>%
    arrange(aggregate_rank) %>%
    slice_head(n = 15) %>%
    ungroup()
  
  # Crear columna para interacción
  top15_por_tratamiento <- top15_por_tratamiento %>%
    mutate(interaccion_y = paste(ligand.complex, receptor.complex, sep = " - "))
  
  # Reordenar factores para eje Y
  top15_por_tratamiento <- top15_por_tratamiento %>%
    group_by(tratamiento) %>%
    mutate(interaccion_y = fct_reorder(interaccion_y, aggregate_rank)) %>%
    ungroup()
  
  # Reordenar eje X
  eje_x <- if (modo == "source") "target" else "source"
  top15_por_tratamiento <- top15_por_tratamiento %>%
    mutate(!!sym(eje_x) := factor(!!sym(eje_x)))
  
  # Crear el plot
  p <- ggplot(top15_por_tratamiento, aes(x = .data[[eje_x]], y = interaccion_y,
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
    labs(
      x = ifelse(modo == "source", "Target", "Source"),
      y = "Interacción (Ligand complex - Receptor complex)",
      title = paste("Top 15 Interacciones en", linea_celular, "por tratamiento"),
      subtitle = paste0(ucfirst(modo), ": ", cell_type),
      color = "Expresión",
      size = "Especificidad"
    )
  
  # Guardar el plot
  file_name <- paste0(output_dir, "Dotplot_", linea_celular, "_", modo, "_", gsub(" ", "_", cell_type), ".png")
  ggsave(file_name, plot = p, height = 15, width = 12, dpi = 600)
  
  # Mostrar el plot (opcional)
  print(p)
}