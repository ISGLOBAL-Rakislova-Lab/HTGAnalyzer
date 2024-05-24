#' HTG_HeatmapQC
#'
#' @description Generate a heatmap based on QC ratios
#'
#' @param ratios A data frame containing QC ratios and other relevant metrics.
#' @param counts_filtered A matrix of counts data with samples as columns and features as rows.
#' @param n_samples Number of samples to label as outliers in PCA.
#' @param threshold_pos Threshold for positive QC value.
#' @param threshold_gdna Threshold for gdna QC value.
#' @param threshold_neg Threshold for negative QC value.
#' @param threshold_size Threshold for library size.
#' @param threshold_median Threshold for median.
#' @param threshold_pca Threshold for PCA genes.
#' @param threshold_ERCC Threshold for ERCC QC value.
#'
#' @return A heatmap plot displaying QC ratios.
#' @export
#'
#' @examples
#' HTG_HeatmapQC(ratios, counts_filtered, n_samples = 3)
#' @name HTG_HeatmapQC

HTG_HeatmapQC <- function(ratios, counts_filtered, n_samples = 3,
                      threshold_pos = 4, threshold_gdna = 0.02,
                      threshold_neg = 0.05, threshold_size = 7e+06,
                      threshold_median = 5, threshold_pca = 1,
                      threshold_ERCC = 0.05) {
  library(ggplot2)
  pca_result <- prcomp(t(counts_filtered))
  pca_data <- as.data.frame(pca_result$x[,1:2])
  pca_data$label <- rownames(pca_data)
  centro_promedio <- colMeans(pca_data[, c("PC1", "PC2")])
  pca_data$distancia_al_centro <- sqrt((pca_data$PC1 - centro_promedio[1])^2 + (pca_data$PC2 - centro_promedio[2])^2)
  pca_data <- pca_data[order(-pca_data$distancia_al_centro), ]

  muestras_a_etiquetar <- head(pca_data$label, n_samples)
  library_size <- colSums(counts_filtered)
  lib_s2 <- data.frame(Sample = colnames(counts_filtered), Size = library_size)


  ratios$PCA_genes <- ifelse(rownames(ratios) %in% muestras_a_etiquetar, "2", "0")
  ratios_heat <- ratios
  ratios_heat <- as.matrix(ratios_heat)

  # Add a fourth column to ratios_heat with library sizes from lib_s2
  ratios_heat <- cbind(ratios_heat, Size = "")
  ratios_heat[,"Size"] <- lib_s2[match(rownames(ratios_heat), rownames(lib_s2)),"Size"]
  ratios_heat <- as.data.frame(ratios_heat)

  # Convert values of ratios_heat to numeric
  ratios_heat[] <- as.data.frame(lapply(ratios_heat, as.numeric))

  # Function to assign 0 or 1 according to QC value
  assign_01_QC <- function(valor, threshold){
    ifelse(valor < threshold, 0, 1)
  }

  # Function to assign 0 or 1 according to library size value
  assign_01_size <- function(valor, threshold){
    ifelse(valor > threshold, 0, 1)
  }

  # Create binary matrix for the heatmap
  bin_matrix <- matrix(0, nrow = nrow(ratios_heat), ncol = 7)
  for (i in 1:nrow(ratios_heat)) {
    bin_matrix[i, 1] <- assign_01_QC(ratios_heat[i, "pos/gens"], threshold_pos)
    bin_matrix[i, 2] <- assign_01_size(ratios_heat[i,"Size"], threshold_size)
    bin_matrix[i, 3] <- assign_01_QC(ratios_heat[i, "nc/gens"], threshold_neg)
    bin_matrix[i, 4] <- assign_01_QC(ratios_heat[i, "gdna/gens"], threshold_gdna)
    bin_matrix[i, 5] <- assign_01_QC(ratios_heat[i,"ERCC/gens"], threshold_ERCC)
    bin_matrix[i, 6] <- assign_01_size(ratios_heat[i,"median"], threshold_median)
    bin_matrix[i, 7] <- assign_01_QC(ratios_heat[i,"PCA_genes"], threshold_pca)
  }

  # Row and column names
  rownames(bin_matrix) <- rownames(ratios_heat)
  colnames(bin_matrix) <- c("QC0","QC1","QC2","QC3","QC4","Median", "PCA_Genes")

  # Plot heatmap
  pheatmap::pheatmap(bin_matrix,
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     fontsize_row = 7,
                     fontsize = 8,
                     legend_breaks = c(0, 1),
                     legend_labels = c("OK", "Possible Outlier"),
                     cellwidth = 40,
                     name = "QC",
                     color = c("#FFF9D0","red")
  )

  rows_with_1 <- suppressWarnings(rownames(bin_matrix)[apply(bin_matrix, 1, any)])
  return(rows_with_1)
  print()
  print(paste("These are the samples plotted at least once in the heatmap:"))
  print(rows_with_1)
  print()
  print(paste("The number of samples that are outliers are:", length(sorted_rows)))
}
