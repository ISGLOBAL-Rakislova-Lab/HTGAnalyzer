#' Plot PCA for Genes Only
#'
#' This function generates a PCA plot for genes only, along with the explained variability and accumulated explained variability plots.
#'
#' @param counts_filtered A data frame containing filtered counts data with genes as rows and samples as columns.
#' @param n_samples Number of samples to label as outliers.
#'
#' @return A combined plot consisting of PCA, explained variability, and accumulated explained variability plots.
#' @export
#'
#' @examples
#' plotPCA_genes(counts_filtered, n_samples = 3)
#' @name HTG_plotPCA_genes
HTG_plotPCA_genes <- function(counts_filtered, n_samples = 3) {
  library(ggplot2)
  library(ggrepel)
  pca_result <- prcomp(t(counts_filtered))
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
