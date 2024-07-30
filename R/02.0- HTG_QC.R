#' HTG_QC
#'
#' @description
#' This function performs various quality control (QC). The QC checks include:
#' QC0: Percentage of positive values greater than 4%.
#' QC1: Library size greater than 7e+06.
#' QC2: Negative control threshold greater than 0.045.
#' QC3: Genomic DNA threshold greater than 0.02.
#' QC4: ERCC threshold greater than 0.025.
#' Calculation of median values.
#' Principal Component Analysis (PCA).
#' These thresholds are tailored for the HTG EdgeSeq transcriptomic panel but can be adjusted as needed. The function generates various plots and calculates all ratios, saving them in .pdf and .xlsx files. Additionally, it includes an optional heatmap to highlight potential outlier samples.
#'
#' @param countsdata A data frame containing the HTG count data. The data must include probes that start with "^NC-|^POS-|^GDNA-|^ERCC-" for the function to work correctly.
#' @param threshold_superior_pos Threshold for upper limit of positive control ratio.
#' @param threshold_inferior_pos Threshold for lower limit of positive control ratio.
#' @param threshold_line_pos Threshold line for positive control ratio.
#' @param threshold_superior_lib Threshold for upper limit of library size.
#' @param threshold_inferior_lib Threshold for lower limit of library size.
#' @param threshold_lib Threshold line for library size.
#' @param threshold_superior_nc Threshold for upper limit of negative control ratio.
#' @param threshold_inferior_nc Threshold for lower limit of negative control ratio.
#' @param threshold_line_nc Threshold line for negative control ratio.
#' @param threshold_superior_gdna Threshold for upper limit of genomic DNA ratio.
#' @param threshold_inferior_gdna Threshold for lower limit of genomic DNA ratio.
#' @param threshold_line_gdna Threshold line for genomic DNA ratio.
#' @param threshold_superior_ercc Threshold for upper limit of ERCC control ratio.
#' @param threshold_inferior_ercc Threshold for lower limit of ERCC control ratio.
#' @param threshold_line_ercc Threshold line for ERCC control ratio.
#' @param threshold_superior_median Threshold for upper limit of median ratio.
#' @param threshold_inferior_median Threshold for lower limit of median ratio.
#' @param threshold_line_median Threshold line for median ratio.
#' @param n_samples Number of samples to label as outliers in plots.
#' @param save_csv Logical, whether to save the ratios as a CSV file. Default is FALSE.
#' @param csv_file The name of the CSV file to save the ratios if save_csv is TRUE. Default is "QC_results.csv".
#' @param show_heatmap Logical, whether to show the heatmap of potential outliers. Default is FALSE.
#'
#' @return This function generates multiple plots displaying various QC metrics, saves an Excel file with all the ratios, and optionally creates a heatmap highlighting potential outlier samples. Additionally, it identifies and returns the most probable outliers based on the QC analysis.
#'
#' @export
#'
#' @examples
#' # Run the function with example data
#' HTG_QC(counts_data, save_csv = TRUE, show_heatmap = TRUE)
#'

HTG_QC <- function(counts_data,
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
                             n_samples = 3,
                             save_csv = FALSE,
                             csv_file = "QC_results.csv",
                             show_heatmap = FALSE) {
  # Load required libraries
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
  library(pheatmap)

  # Filter counts_data data
  counts_filtered <- subset(counts_data, !grepl("^NC-|^POS-|^GDNA-|^ERCC-", rownames(counts_data)))

  # Subsets
  cts_ERCC <- as.data.frame(subset(counts_data, grepl("^ERCC-", rownames(counts_data))))
  cts_NC <- as.data.frame(subset(counts_data, grepl("^NC-", rownames(counts_data))))
  cts_POS <- as.data.frame(subset(counts_data, grepl("^POS-", rownames(counts_data))))
  cts_GDNA <- as.data.frame(subset(counts_data, grepl("^GDNA-", rownames(counts_data))))

  # Ratios
  total_gens <- colSums(counts_filtered)
  total_POS <- colSums(cts_POS)
  total_NC <- colSums(cts_NC)
  total_GDNA <- colSums(cts_GDNA)
  total_ERCC <- colSums(cts_ERCC)

  # Calculate ratios
  ratios <- data.frame(
    total_POS = total_POS,
    total_GDNA = total_GDNA,
    total_gens = total_gens,
    total_NC = total_NC,
    total_ERCC = total_ERCC
  )
  ratios$`pos/gens` <- (ratios$total_POS / ratios$total_gens) * 100
  ratios$`gdna/gens` <- (ratios$total_GDNA / ratios$total_gens) * 100
  ratios$`nc/gens` <- (ratios$total_NC / ratios$total_gens) * 100
  ratios$`ERCC/gens` <- (ratios$total_ERCC / ratios$total_gens) * 100

  # Add median column
  summary_stats <- HTG_calculate_summary_stats(counts_filtered)
  ratios$median <- summary_stats$Median
  ratiosb<-ratios
  ratiosb$min<- summary_stats$Min
  ratiosb$max<- summary_stats$Max
  ratiosb$mean<- summary_stats$Mean

  # Add sample names as a factor column
  ratios$samples <- factor(rownames(ratios))

  # Optionally save as CSV
  if (save_csv) {
    write.csv(ratiosb, csv_file, row.names = TRUE)
    cat("Ratios saved as", csv_file, "\n")
  }

  # Calculate library size
  library_size <- colSums(counts_filtered)

  # Create dataframe of library size
  lib_s <- data.frame(Size = library_size)

  pdf("QC_plots.pdf")
  # Density plot
  plot(density(lib_s$Size),
       col = "#4793AF",
       main = "Library Size per Sample",
       xlab = "Counts")

  # Positive controls
  colores_pos <- ifelse(ratios$`pos/gens` < threshold_inferior_pos, "#4793AF",
                        ifelse(ratios$`pos/gens` > threshold_superior_pos, "red", "#FFC470"))
  plot(ratios$`pos/gens`, xlab = "", ylab = "pos/gens", col = colores_pos,
       xaxt = "n", pch = 19, main = "Positive control 4% (QC0)")
  axis(1, at = 1:nrow(ratios), labels = rownames(ratios),
       las = 2, cex.axis = 0.8)
  abline(h = threshold_line_pos, col = "red")

  # Library size
  lib_s2 <- data.frame(Sample = colnames(counts_filtered), Size = library_size)
  colores <- ifelse(lib_s2$Size < threshold_inferior_lib, "#4793AF",
                    ifelse(lib_s2$Size > threshold_superior_lib, "red", "#FFC470"))
  plot(lib_s2$Size, xlab = "", ylab = "Library Size", col = colores,
       xaxt = "n", pch = 19, main = "Library Size per Sample (QC1)", cex.axis = 0.8)
  abline(h = threshold_lib, col = "red")
  axis(1, at = 1:length(lib_s2$Sample), labels = lib_s2$Sample, las = 2, cex.axis = 0.8)

  # Negative controls
  colores_nc <- ifelse(ratios$`nc/gens` < threshold_inferior_nc, "#4793AF",
                       ifelse(ratios$`nc/gens` > threshold_superior_nc, "red", "#FFC470"))
  plot(ratios$`nc/gens`, xlab = "", ylab = "nc/gens", col = colores_nc,
       xaxt = "n", pch = 19, main = "Negative Control (QC2)")
  axis(1, at = 1:nrow(ratios), labels = rownames(ratios),
       las = 2, cex.axis = 0.8)
  abline(h = threshold_line_nc, col = "red")

  # Genomic DNA
  colores_gdna <- ifelse(ratios$`gdna/gens` < threshold_inferior_gdna, "#4793AF",
                         ifelse(ratios$`gdna/gens` > threshold_superior_gdna, "red", "#FFC470"))
  plot(ratios$`gdna/gens`, xlab = "", ylab = "gdna/gens", col = colores_gdna,
       xaxt = "n", pch = 19, main = "Genomic DNA (QC3)")
  axis(1, at = 1:nrow(ratios), labels = rownames(ratios),
       las = 2, cex.axis = 0.8)
  abline(h = threshold_line_gdna, col = "red")

  # ERCC
  colores_ercc <- ifelse(ratios$`ERCC/gens` < threshold_inferior_ercc, "#4793AF",
                         ifelse(ratios$`ERCC/gens` > threshold_superior_ercc, "red", "#FFC470"))
  plot(ratios$`ERCC/gens`, xlab = "", ylab = "ERCC", col = colores_ercc,
       xaxt = "n", pch = 19, main = "ERCC (QC4)")
  axis(1, at = 1:nrow(ratios), labels = rownames(ratios),
       las = 2, cex.axis = 0.8)
  abline(h = threshold_line_ercc, col = "red")

  # Median
  colores_med <- ifelse(ratios$median < threshold_inferior_median, "#4793AF",
                        ifelse(ratios$median > threshold_superior_median, "red", "#FFC470"))
  plot(ratios$median, xlab = "", ylab = "Median", col = colores_med,
       xaxt = "n", pch = 19, main = "Median")
  axis(1, at = 1:nrow(ratios), labels = rownames(ratios),
       las = 2, cex.axis = 0.8)
  abline(h = threshold_line_median, col = "red")
  dev.off()

  pca_result <- prcomp(t(counts_filtered))
  pca_data <- as.data.frame(pca_result$x[,1:2])
  pca_data$label <- rownames(pca_data)
  centro_promedio <- colMeans(pca_data[, c("PC1", "PC2")])
  pca_data$distancia_al_centro <- sqrt((pca_data$PC1 - centro_promedio[1])^2 + (pca_data$PC2 - centro_promedio[2])^2)
  pca_data <- pca_data[order(-pca_data$distancia_al_centro), ]

  samples_a_etiquetar <- head(pca_data$label, n_samples)
  library_size <- colSums(counts_filtered)
  lib_s2 <- data.frame(Sample = colnames(counts_filtered), Size = library_size)


  ratios$PCA_genes <- ifelse(rownames(ratios) %in% samples_a_etiquetar, "2", "0")
  ratios_heat <- ratios
  ratios_heat <- as.matrix(ratios_heat)

  # Add a fourth column to ratios_heat with library sizes from lib_s2
  ratios_heat <- cbind(ratios_heat, Size = "")
  ratios_heat[,"Size"] <- lib_s2[match(rownames(ratios_heat), rownames(lib_s2)),"Size"]
  ratios_heat <- as.data.frame(ratios_heat)

  # Convert values of ratios_heat to numeric
  cols_to_convert <- c("total_POS", "total_GDNA", "total_gens", "total_NC", "total_ERCC",
                       "pos/gens", "gdna/gens", "nc/gens", "ERCC/gens", "median", "PCA_genes", "Size")
  ratios_heat[cols_to_convert] <- lapply(ratios_heat[cols_to_convert], as.numeric)
  str(ratios_heat)

  if (show_heatmap) {
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
      bin_matrix[i, 1] <- assign_01_QC(ratios_heat[i, "pos/gens"], threshold_line_pos)
      bin_matrix[i, 2] <- assign_01_size(ratios_heat[i,"Size"], threshold_lib)
      bin_matrix[i, 3] <- assign_01_QC(ratios_heat[i, "nc/gens"], threshold_line_nc)
      bin_matrix[i, 4] <- assign_01_QC(ratios_heat[i, "gdna/gens"], threshold_line_gdna)
      bin_matrix[i, 5] <- assign_01_QC(ratios_heat[i,"ERCC/gens"], threshold_line_ercc)
      bin_matrix[i, 6] <- assign_01_size(ratios_heat[i,"median"], threshold_line_median)
      bin_matrix[i, 7] <- assign_01_QC(ratios_heat[i,"PCA_genes"], 1)
    }

    # Row and column names
    rownames(bin_matrix) <- rownames(ratios_heat)
    colnames(bin_matrix) <- c("QC0","QC1","QC2","QC3","QC4","Median", "PCA_Genes")

    # Plot heatmap
    pdf("heatmap_outliers.pdf")
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
    dev.off()
    rows_with_1 <- suppressWarnings(rownames(bin_matrix)[apply(bin_matrix, 1, any)])
    cat("\033[32m                              ***\033[0m\n")
    cat(paste("\033[32mThese are the samples plotted at least once in the heatmap:\033[0m\n"))
    print(rows_with_1)
    cat(paste("The number of samples that are outliers are:", length(rows_with_1)))


    return(rows_with_1)
  }
}


