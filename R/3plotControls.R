#' Plot Controls
#'
#' This function generates multiple plots to visualize various quality control (QC) metrics.
#'
#' @param counts_filtered data frame without probes.
#' @param ratios A data frame containing QC ratios and other relevant metrics.
#' @param threshold_superior_pos Threshold for positive control upper limit.
#' @param threshold_inferior_pos Threshold for positive control lower limit.
#' @param threshold_line_pos Threshold line for positive control.
#' @param threshold_superior_gdna Threshold for genomic DNA upper limit.
#' @param threshold_inferior_gdna Threshold for genomic DNA lower limit.
#' @param threshold_line_gdna Threshold line for genomic DNA.
#' @param threshold_superior_nc Threshold for negative control upper limit.
#' @param threshold_inferior_nc Threshold for negative control lower limit.
#' @param threshold_line_nc Threshold line for negative control.
#' @param threshold_superior_median Threshold for median upper limit.
#' @param threshold_inferior_median Threshold for median lower limit.
#' @param threshold_line_median Threshold line for median.
#' @param threshold_superior_ercc Threshold for ERCC control upper limit.
#' @param threshold_inferior_ercc Threshold for ERCC control lower limit.
#' @param threshold_line_ercc Threshold line for ERCC control.
#' @param n_samples Number of samples to label as outliers.
#'
#' @return Multiple plots displaying different QC metrics.
#' @export
#'
#' @examples
#' plotControls(counts_filtered,ratios)
#' @name HTG_calculate_ratios
#'
HTG_plotControls <- function(counts_filtered,ratios,
                         threshold_superior_pos = 5,
                         threshold_inferior_pos = 3,
                         threshold_line_pos = 4,
                         threshold_superior_lib = 8e+06,
                         threshold_inferior_lib = 5e+06,
                         threshold_lib = 7e+06,
                         threshold_superior_nc = 0.05,
                         threshold_inferior_nc = 0.035,
                         threshold_line_nc = 0.045,
                         threshold_superior_gdna = 0.025,
                         threshold_inferior_gdna = 0.015,
                         threshold_line_gdna = 0.02,
                         threshold_superior_ercc = 0.03,
                         threshold_inferior_ercc = 0.015,
                         threshold_line_ercc = 0.025,
                         threshold_superior_median = 7,
                         threshold_inferior_median = 3,
                         threshold_line_median = 5,
                         n_samples = 3) {

  # Calculate library size
  library_size <- colSums(counts_filtered)

  # Create dataframe of library size
  lib_s <- data.frame(Size = library_size)

  # Density plot
  plot(density(lib_s$Size),
       col = "#4793AF",
       main = "Library Size per Sample",
       xlab = "Counts")


  # Controles positivos
  colores_pos <- ifelse(ratios$`pos/gens` < threshold_inferior_pos, "#4793AF",
                        ifelse(ratios$`pos/gens` > threshold_superior_pos, "red", "#FFC470"))
  plot(ratios$`pos/gens`, xlab = "", ylab = "pos/gens", col = colores_pos,
       xaxt = "n", pch=19, main = "Positive control 4% (QC0)")
  axis(1, at = 1:nrow(ratios), labels = rownames(ratios),
       las = 2, cex.axis = 0.8)
  abline(h = threshold_line_pos, col = "red")

  # library size
  lib_s2 <- data.frame(Sample = colnames(counts_filtered), Size = library_size)
  colores <- ifelse(lib_s2$Size < threshold_inferior_lib, "#4793AF",
                    ifelse(lib_s2$Size > threshold_superior_lib, "red", "#FFC470") )
  plot(lib_s2$Size, xlab = "", ylab = "Library Size", col = colores,
       xaxt = "n", pch=19, main = "Library Size per Sample (QC1)", cex.axis = 0.8)
  abline(h = threshold_lib, col = "red")
  axis(1, at = 1:length(lib_s2$Sample), labels = lib_s2$Sample, las = 2, cex.axis = 0.8)

  # Controles Negativos
  colores_nc <- ifelse(ratios$`nc/gens` < threshold_inferior_nc, "#4793AF",
                       ifelse(ratios$`nc/gens` > threshold_superior_nc, "red", "#FFC470"))
  plot(ratios$`nc/gens`, xlab = "", ylab = "nc/gens", col = colores_nc,
       xaxt = "n", pch=19,main = "Negative Control (QC2)")
  axis(1, at = 1:nrow(ratios), labels = rownames(ratios),
       las = 2, cex.axis = 0.8)
  abline(h = threshold_line_nc, col = "red")

  # Controles gDNA
  colores_gdna <- ifelse(ratios$`gdna/gens` < threshold_inferior_gdna, "#4793AF",
                         ifelse(ratios$`gdna/gens` > threshold_superior_gdna, "red","#FFC470"))
  plot(ratios$`gdna/gens`, xlab = "", ylab = "gdna/gens", col = colores_gdna,
       xaxt = "n", pch=19,main = "Genomic DNA (QC3)")
  axis(1, at = 1:nrow(ratios), labels = rownames(ratios),
       las = 2, cex.axis = 0.8)
  abline(h = threshold_line_gdna, col = "red")

  # ERCC
  colores_ercc <- ifelse(ratios$"ERCC/gens" < threshold_inferior_ercc, "#4793AF",
                         ifelse(ratios$"ERCC/gens" > threshold_superior_ercc, "red", "#FFC470"))
  plot(ratios$"ERCC/gens", xlab = "", ylab = "ERCC", col = colores_ercc,
       xaxt = "n", pch=19, main = "ERCC Control (QC4)")
  axis(1, at = 1:nrow(ratios), labels = rownames(ratios),
       las = 2, cex.axis = 0.8)
  abline(h = threshold_line_ercc, col = "red")

  # Median
  colores_median <- ifelse(ratios$median < threshold_inferior_median, "#4793AF",
                           ifelse(ratios$median > threshold_superior_median, "red", "#FFC470"))
  plot(ratios$median, xlab = "", ylab = "Median", col = colores_median,
       xaxt = "n", pch=19, main= "Median (QC5)")
  axis(1, at = 1:nrow(ratios), labels = rownames(ratios),
       las = 2, cex.axis = 0.8)
  abline(h = threshold_line_median, col = "red")

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

