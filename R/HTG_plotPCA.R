#' HTG_plotPCA
#'
#' @description This function generates a PCA plot specifically for genes, along with plots showing explained variability and accumulated explained variability. It highlights the samples that are the most distant from the center of the PCA plot.
#'
#' @param counts_data A data frame containing counts data with genes as rows and samples as columns.
#' @param n_samples The number of samples to label as outliers based on their distance from the center of the PCA plot.
#' @param pattern An optional pattern to filter out rows based on row names. in HTG is normally: "^NC-|^POS-|^GDNA-|^ERCC-"
#'
#' @return A combined plot consisting of the PCA plot, explained variability plot, and accumulated explained variability plot.
#' @export
#'
#' @examples
#' HTG_plotPCA(counts_data_tutorial, n_samples = 3)
#' @name HTG_plotPCA
#'
HTG_plotPCA <- function(counts_data, n_samples = 3, pattern = NULL) {
  library(ggplot2)
  library(ggrepel)
  library(gridExtra)  # Añadimos gridExtra para la combinación de gráficos

  # Filtrado opcional de datos
  if (!is.null(pattern)) {
    counts_data <- subset(counts_data, !grepl(pattern, rownames(counts_data)))
  }

  # Realizar PCA
  pca_result <- prcomp(t(counts_data))
  pca_data <- as.data.frame(pca_result$x[,1:2])
  pca_data$label <- rownames(pca_data)

  # Calcular distancia al centro promedio
  centro_promedio <- colMeans(pca_data[, c("PC1", "PC2")])
  pca_data$distancia_al_centro <- sqrt((pca_data$PC1 - centro_promedio[1])^2 + (pca_data$PC2 - centro_promedio[2])^2)
  pca_data <- pca_data[order(-pca_data$distancia_al_centro), ]

  # Seleccionar las muestras a etiquetar
  muestras_a_etiquetar <- head(pca_data$label, n_samples)
  pca_data$Tag <- ifelse(pca_data$label %in% muestras_a_etiquetar, "Tag", "No Tag")

  # Calcular el porcentaje de variabilidad explicada
  porcentaje_explicado <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 2)
  titulo_x <- paste0("PC1 (", porcentaje_explicado[1], "%)")
  titulo_y <- paste0("PC2 (", porcentaje_explicado[2], "%)")

  # Gráfico 1: PCA plot
  p1 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Tag)) +
    geom_point() +
    ggrepel::geom_text_repel(data = subset(pca_data, Tag == "Tag"), aes(label = label), color = "red") +
    labs(title = "PCA",
         x = titulo_x,
         y = titulo_y) +
    theme(legend.position = "none",
          plot.margin = margin(0.5, 1, 0.5, 0.5, "cm")) +
    scale_color_manual(values = c("Tag" = "red", "No Tag" = "black")) +
    scale_y_continuous(labels = scales::scientific_format())

  # Gráfico 2: Barras de los primeros 10 componentes principales
  porcentaje_explicado_10 <- head(porcentaje_explicado, 10)
  varianza_acumulada <- cumsum(porcentaje_explicado)
  varianza_acumulada_10 <- head(varianza_acumulada, 10)

  bardata <- data.frame(Componente = factor(paste("PC", 1:10), levels = paste("PC", 1:10)),
                        Porcentaje = porcentaje_explicado_10)

  p2 <- ggplot(bardata, aes(x = Componente, y = Porcentaje)) +
    geom_bar(stat = "identity", fill = "#4793AF") +
    geom_text(aes(label = sprintf("%.2f", varianza_acumulada_10)),
              vjust = -0.5, size = 3, color = "black") +
    labs(title = "The First 10 Principal Components",
         x = "Principal Components",
         y = "Explained Variability (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Gráfico 3: Línea de varianza acumulada
  df <- data.frame(varianza_acumulada_10, pc = factor(paste("PC", 1:10, sep=""), levels = paste("PC", 1:10, sep="")))

  p3 <- ggplot(data = df, aes(x = pc, y = varianza_acumulada_10, group = 1)) +
    geom_point() +
    geom_line() +
    theme_bw() +
    labs(x = "Principal Component",
         y = "Explained Variability Accumulated (%)")

  # Combinar gráficos en una sola visualización
  combined_plot <- grid.arrange(p1, p2, p3, nrow = 1)

  # Guardar el gráfico combinado en un archivo PDF
  ggsave("plot_PCA.pdf", combined_plot, width = 14, height = 5)

}


