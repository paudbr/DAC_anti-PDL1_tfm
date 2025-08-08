generar_dotplot <- function(df, linea_celular, cell_type, modo = c("source", "target"), output_dir = "./") {
  modo <- match.arg(modo)
  
  # Filtrado
  datos_filtrados <- df %>%
    filter(linea_celular == linea_celular) %>%
    filter((!!sym(modo)) == cell_type)
  
  # Crear columna para interacción
  top15_por_tratamiento <- datos_filtrados %>%
    mutate(interaccion_y = paste(ligand_complex, receptor_complex, sep = " - "))
  
  # Reordenar factores para eje Y
  top15_por_tratamiento <- top15_por_tratamiento %>%
    group_by(tratamiento) %>%
    arrange(magnitude_rank, .by_group = TRUE)%>% filter(specificity_rank <= 0.01) %>% mutate(rank_in_group = row_number()) %>%
    filter(rank_in_group <= 15) %>%
    ungroup()
  
  
  # Reordenar eje X
  eje_x <- if (modo == "source") "target" else "source"
  top15_por_tratamiento <- top15_por_tratamiento %>%
    mutate(!!sym(eje_x) := factor(!!sym(eje_x)))
  
  # Crear el plot
  p <- ggplot(top15_por_tratamiento, aes(x = .data[[eje_x]], y = interaccion_y,
                                         color = magnitude_rank, size = specificity_rank)) +
  
    geom_point(alpha = 0.8) +
    scale_color_viridis_c(option = "plasma", name = "Magnitude_rank", trans = "reverse") +
    scale_size(range = c(1, 6), name = "Specificity_rank", trans = "reverse") +
      facet_wrap(~ tratamiento, scales = "free") +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 75, hjust = 1, size = 7),
        axis.text.y = element_text(size = 9),
        strip.text = element_text(face = "bold", size = 10)
      ) +
      labs(
        x = ifelse(modo == "source", "Target", "Source"),
        y = "Interaction (Ligand complex - Receptor complex)",
        title = paste("Top 15 Interacciones en ", linea_celular, " por tratamiento"),
        subtitle = paste0(ucfirst(modo), ": ", cell_type),
        color = "Magnitude_rank",
        size = "Specificity_rank"
      )
  
  # Guardar el plot
  file_name <- paste0(output_dir, "Dotplot_", linea_celular, "_", modo, "_", gsub(" ", "_", cell_type), ".png")
  ggsave(file_name, plot = p, height = 15, width = 12, dpi = 600)
  
  # Mostrar el plot (opcional)
  print(p)
}


