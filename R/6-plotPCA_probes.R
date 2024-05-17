
#' PCA Plot and Explained Variability
#'
#' This function generates a PCA plot of samples and additional plots showing the explained variability of the principal components.
#'
#' @param counts An object with sample names as column names and genes as row names, containing counts data.
#' @param n_samples The number of samples to label in red, which are farthest from the center.
#'
#' @return A combined plot of PCA samples, the first 10 principal components, and accumulated explained variability.
#' @export
#'
#' @examples
#' plotPCA_probes(counts_data, 5)
#' plotPCA_probes(counts_data)
#' @name HTG_plotPCA_probes

HTG_plotPCA_probes <- function(counts, n_samples = 3) {
  library(ggplot2)
  library(ggrepel)
  pca_result <- prcomp(t(counts))
  pca_data <- as.data.frame(pca_result$x[,1:2])
  pca_data$label <- rownames(pca_data)
  centro_promedio <- colMeans(pca_data[, c("PC1", "PC2")])
  pca_data$distancia_al_centro <- sqrt((pca_data$PC1 - centro_promedio[1])^2 + (pca_data$PC2 - centro_promedio[2])^2)
  pca_data <- pca_data[order(-pca_data$distancia_al_centro), ]

  muestras_a_etiquetar <- head(pca_data$label, n_samples)
  pca_data$Tag <- ifelse(pca_data$label %in% muestras_a_etiquetar, "Tag", "No Tag")

  porcentaje_explicado <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 2)
  titulo_x <- paste0("PC1 (", porcentaje_explicado[1], "%)")
  titulo_y <- paste0("PC2 (", porcentaje_explicado[2], "%)")

  p1 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Tag)) +
    geom_point() +
    ggrepel::geom_text_repel(data = subset(pca_data, Tag == "Tag"), aes(label = label), color = "red") +
    labs(title = "PCA samples (genes only)",
         x = titulo_x,
         y = titulo_y) +
    theme(legend.position = "none",
          plot.margin = margin(0.5, 1, 0.5, 0.5, "cm")) +
    scale_color_manual(values = c("Tag" = "red", "No Tag" = "black")) +
    scale_y_continuous(labels = scales::scientific_format())

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

  varianza_acumulada_10 <- head(varianza_acumulada, 10)
  df <- data.frame(varianza_acumulada_10, pc = factor(paste("PC", 1:10, sep=""), levels = paste("PC", 1:10, sep="")))

  p3 <- ggplot(data = df, aes(x = pc, y = varianza_acumulada_10, group = 1)) +
    geom_point() +
    geom_line() +
    theme_bw() +
    labs(x = "Principal Component",
         y = "Explained Variability Accumulated (%)")

  # Return combined plots
  cowplot::plot_grid(p1, p2, p3, nrow = 1, align = "hv")
}
