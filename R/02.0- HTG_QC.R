#' HTG_QC
#'
#' @description
#' This function performs various quality control (QC) checks. The QC checks include:
#' - QC0: Percentage of positive values greater than 4%.
#' - QC1: Library size greater than 7e+06.
#' - QC2: Negative control threshold greater than 0.045.
#' - QC3: Genomic DNA threshold greater than 0.02.
#' - QC4: ERCC threshold greater than 0.025.
#'
#' These thresholds are tailored for the HTG EdgeSeq transcriptomic panel but can be adjusted as needed.
#'
#' In addition to these plots, this function generates a data frame that can be saved as a .csv file (a preview of it will be shown). This data frame includes:
#' - The sum of each probe for each sample (total genes, positive, negative, gdna, and ercc)
#' - The ratio for each sample
#' - The size of each sample
#' - A column called PCA which indicates the samples furthest from the center (you can specify how many samples to highlight as the furthest from the center)
#'
#' This function also includes an optional heatmap to highlight potential outlier samples, which will be saved in a vector.
#'
#' The plots will be shown in the console and saved in the current working directory.
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
#'
#' @return This function generates multiple plots displaying various QC metrics, saves an Excel file with all the ratios, and optionally creates a heatmap highlighting potential outlier samples. Additionally, it identifies and returns the most probable outliers based on the QC analysis.
#'
#' @export
#'
#' @examples
#' # Run the function with example data
#' HTG_QC(counts_data, save_csv = TRUE)
#'

HTG_QC <- function(counts_data, pattern = "^NC-|^POS-|^GDNA-|^ERCC-",
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
                             save_csv = TRUE,
                             csv_file = "QC_results.csv") {
  # Load required libraries
  suppressMessages({
    library(ggplot2)
    library(ggrepel)
    library(cowplot)
    library(reshape2)
  })

  # Filter counts_data data
  counts_filtered <- subset(counts_data, !grepl(pattern, rownames(counts_data)))

  min_values <- apply(counts_filtered, 2, min)
  max_values <- apply(counts_filtered, 2, max)
  mean_values <- apply(counts_filtered, 2, mean)
  median_values <- apply(counts_filtered, 2, median)
  mode_values <- apply(counts_filtered, 2, function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  })
  sd_values <- apply(counts_filtered, 2, sd)
  var_values <- apply(counts_filtered, 2, var)
  range_values <- apply(counts_filtered, 2, function(x) max(x) - min(x))
  quartile_1 <- apply(counts_filtered, 2, function(x) quantile(x, 0.25))
  quartile_3 <- apply(counts_filtered, 2, function(x) quantile(x, 0.75))
  iqr_values <- quartile_3 - quartile_1
  skewness_values <- apply(counts_filtered, 2, function(x) {
    n <- length(x)
    mean_x <- mean(x)
    sd_x <- sd(x)
    sum((x - mean_x)^3) / ((n - 1) * (sd_x^3))
  })
  kurtosis_values <- apply(counts_filtered, 2, function(x) {
    n <- length(x)
    mean_x <- mean(x)
    sd_x <- sd(x)
    sum((x - mean_x)^4) / ((n - 1) * (sd_x^4)) - 3
  })
  missing_values <- apply(counts_filtered, 2, function(x) sum(is.na(x)))
  cv_values <- sd_values / mean_values

  summary_stats <- data.frame(
    Min = min_values,
    Max = max_values,
    Mean = mean_values,
    Median = median_values,
    Mode = mode_values,
    SD = sd_values,
    Variance = var_values,
    Range = range_values,
    Q1 = quartile_1,
    Q3 = quartile_3,
    IQR = iqr_values,
    Skewness = skewness_values,
    Kurtosis = kurtosis_values,
    Missing = missing_values,
    CV = cv_values
  )
  write.csv(summary_stats, file = "summary_stats.csv")
  ####################### PCA
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
    labs(title = "PCA",
         x = titulo_x,
         y = titulo_y) +
    theme(legend.position = "none",
          plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
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
         y = "Explained Variability Accumulated (%)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Return combined plots
  pdf("plot_PCA.pdf", width = 14, height = 5) # Adjust these values as needed
  print(cowplot::plot_grid(p1, p2, p3, nrow = 1, align = "hv"))
  dev.off()


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

  #####################

  library_size <- colSums(counts_filtered)

  # Create dataframe of library size
  lib_s2 <- data.frame(Sample = colnames(counts_filtered), Size = library_size)
  pca_result <- prcomp(t(counts_filtered))
  pca_data <- as.data.frame(pca_result$x[, 1:2])
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
  ratios_heat[, "Size"] <- lib_s2[match(rownames(ratios_heat), rownames(lib_s2)), "Size"]
  ratios_heat <- as.data.frame(ratios_heat)

  # Convert values of ratios_heat to numeric
  cols_to_convert <- c("total_POS", "total_GDNA", "total_gens", "total_NC", "total_ERCC",
                       "pos/gens", "gdna/gens", "nc/gens", "ERCC/gens", "median", "PCA_genes", "Size")
  ratios_heat[cols_to_convert] <- lapply(ratios_heat[cols_to_convert], as.numeric)
  str(ratios_heat)

  assign_01_QC <- function(valor, threshold) {
    ifelse(valor < threshold, 0, 1)
  }

  # Function to assign 0 or 1 according to library size value
  assign_01_size <- function(valor, threshold) {
    ifelse(valor > threshold, 0, 1)
  }

  # Create binary matrix for the heatmap
  bin_matrix <- matrix(0, nrow = nrow(ratios_heat), ncol = 7)
  for (i in 1:nrow(ratios_heat)) {
    bin_matrix[i, 1] <- assign_01_QC(ratios_heat[i, "pos/gens"], threshold_line_pos)
    bin_matrix[i, 2] <- assign_01_size(ratios_heat[i, "Size"], threshold_lib)
    bin_matrix[i, 3] <- assign_01_QC(ratios_heat[i, "nc/gens"], threshold_line_nc)
    bin_matrix[i, 4] <- assign_01_QC(ratios_heat[i, "gdna/gens"], threshold_line_gdna)
    bin_matrix[i, 5] <- assign_01_QC(ratios_heat[i, "ERCC/gens"], threshold_line_ercc)
    bin_matrix[i, 6] <- assign_01_size(ratios_heat[i, "median"], threshold_line_median)
    bin_matrix[i, 7] <- assign_01_QC(ratios_heat[i, "PCA_genes"], 1)
  }

  # Row and column names
  rownames(bin_matrix) <- rownames(ratios_heat)
  colnames(bin_matrix) <- c("QC0", "QC1", "QC2", "QC3", "QC4", "Median", "PCA_Genes")

  # Convert the matrix to a data frame for ggplot2
  bin_df <- as.data.frame(bin_matrix)
  bin_df$Sample <- rownames(bin_df)

  # Melt the data frame
  bin_df_melted <- melt(bin_df, id.vars = "Sample")

  pdf("QC_plots.pdf")

  # Violin plot for each sample
  b<- ggplot(lib_s2, aes(x = Sample, y = Size)) +
    geom_violin(fill = "#4793AF", adjust = 1.5) +  # Adjust the `adjust` parameter for smoothing
    geom_jitter(height = 0, width = 0.2, color = "red") +  # Add jittered points to show individual sizes
    labs(title = "Violin Plot of Library Sizes", x = "Sample", y = "Library Size") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels if needed
    theme_minimal()  # Optional: Use a minimal theme for cleaner look
  print(b)


  # Positive controls
  max_value <- max(ratios$`pos/gens`, threshold_line_pos)
  colores_pos <- ifelse(ratios$`pos/gens` < threshold_inferior_pos, "#4793AF",
                        ifelse(ratios$`pos/gens` > threshold_superior_pos, "red", "#FFC470"))
  plot(ratios$`pos/gens`, xlab = "", ylab = "pos/gens", col = colores_pos,
       xaxt = "n", pch = 19, main = "Positive control 4% (QC0)", ylim = c(0, max_value))
  axis(1, at = 1:nrow(ratios), labels = rownames(ratios), las = 2, cex.axis = 0.8)
  abline(h = threshold_line_pos, col = "red")


  # Library size
  max_size <- max(lib_s2$Size, threshold_lib)
  min_size <- min(lib_s2$Size, threshold_lib)
  colores <- ifelse(lib_s2$Size > threshold_inferior_lib, "#4793AF",
                    ifelse(lib_s2$Size < threshold_superior_lib, "red", "#FFC470"))
  plot(lib_s2$Size, xlab = "", ylab = "Library Size", col = colores,
       xaxt = "n", pch = 19, main = "Library Size per Sample (QC1)", cex.axis = 0.8,
       ylim = c(min_size, max_size))
  axis(1, at = 1:length(lib_s2$Sample), labels = lib_s2$Sample, las = 2, cex.axis = 0.8)
  abline(h = threshold_lib, col = "red")

  # Negative controls
  max_value <- max(ratios$`nc/gens`, threshold_line_nc)
  colores_nc <- ifelse(ratios$`nc/gens` < threshold_inferior_nc, "#4793AF",
                       ifelse(ratios$`nc/gens` > threshold_superior_nc, "red", "#FFC470"))
  plot(ratios$`nc/gens`, xlab = "", ylab = "nc/gens", col = colores_nc,
       xaxt = "n", pch = 19, main = "Negative Control (QC2)", ylim = c(0, max_value))
  axis(1, at = 1:nrow(ratios), labels = rownames(ratios), las = 2, cex.axis = 0.8)
  abline(h = threshold_line_nc, col = "red")


  # Genomic DNA
  max_value <- max(ratios$`gdna/gens`, threshold_line_gdna)
  colores_gdna <- ifelse(ratios$`gdna/gens` < threshold_inferior_gdna, "#4793AF",
                         ifelse(ratios$`gdna/gens` > threshold_superior_gdna, "red", "#FFC470"))
  plot(ratios$`gdna/gens`, xlab = "", ylab = "gdna/gens", col = colores_gdna,
       xaxt = "n", pch = 19, main = "Genomic DNA (QC3)", ylim = c(0, max_value))
  axis(1, at = 1:nrow(ratios), labels = rownames(ratios), las = 2, cex.axis = 0.8)
  abline(h = threshold_line_gdna, col = "red")


  # ERCC
  max_value <- max(ratios$`ERCC/gens`, threshold_line_ercc)
  colores_ercc <- ifelse(ratios$`ERCC/gens` < threshold_inferior_ercc, "#4793AF",
                         ifelse(ratios$`ERCC/gens` > threshold_superior_ercc, "red", "#FFC470"))
  plot(ratios$`ERCC/gens`, xlab = "", ylab = "ERCC", col = colores_ercc,
       xaxt = "n", pch = 19, main = "ERCC (QC4)", ylim = c(0, max_value))
  axis(1, at = 1:nrow(ratios), labels = rownames(ratios), las = 2, cex.axis = 0.8)
  abline(h = threshold_line_ercc, col = "red")


  ## Median
  max_value <- max(ratios$median, threshold_line_median)
  colores_med <- ifelse(ratios$median < threshold_inferior_median, "#4793AF",
                        ifelse(ratios$median > threshold_superior_median, "red", "#FFC470"))
  plot(ratios$median, xlab = "", ylab = "Median", col = colores_med,
       xaxt = "n", pch = 19, main = "Median", ylim = c(0, max_value))
  axis(1, at = 1:nrow(ratios), labels = rownames(ratios), las = 2, cex.axis = 0.8)
  abline(h = threshold_line_median, col = "red")


  # Create the heatmap with ggplot2
  a<- ggplot(bin_df_melted, aes(x = variable, y = Sample, fill = factor(value))) +
    geom_tile(color = "white") +
    scale_fill_manual(values = c("0" = "#FFF9D0", "1" = "red"), labels = c("OK", "Possible Outlier")) +
    labs(x = "QC Metrics", y = "Samples", fill = "QC Status") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 7),
          legend.position = "bottom")
  print(a)
    dev.off()

    # Violin plot
    b<- ggplot(lib_s2, aes(x = Sample, y = Size)) +
      geom_violin(fill = "#4793AF") +
      geom_jitter(height = 0, width = 0.2, color = "red") +
      labs(title = "Violin Plot of Library Sizes", x = "Sample", y = "Library Size") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    print(b)

    # Positive controls
    max_value <- max(ratios$`pos/gens`, threshold_line_pos)
    colores_pos <- ifelse(ratios$`pos/gens` < threshold_inferior_pos, "#4793AF",
                          ifelse(ratios$`pos/gens` > threshold_superior_pos, "red", "#FFC470"))
    plot(ratios$`pos/gens`, xlab = "", ylab = "pos/gens", col = colores_pos,
         xaxt = "n", pch = 19, main = "Positive control 4% (QC0)", ylim = c(0, max_value))
    axis(1, at = 1:nrow(ratios), labels = rownames(ratios), las = 2, cex.axis = 0.8)
    abline(h = threshold_line_pos, col = "red")


    # Library size
    max_size <- max(lib_s2$Size, threshold_lib)
    min_size <- min(lib_s2$Size, threshold_lib)
    colores <- ifelse(lib_s2$Size > threshold_inferior_lib, "#4793AF",
                      ifelse(lib_s2$Size < threshold_superior_lib, "red", "#FFC470"))
    plot(lib_s2$Size, xlab = "", ylab = "Library Size", col = colores,
         xaxt = "n", pch = 19, main = "Library Size per Sample (QC1)", cex.axis = 0.8,
         ylim = c(min_size, max_size))
    axis(1, at = 1:length(lib_s2$Sample), labels = lib_s2$Sample, las = 2, cex.axis = 0.8)
    abline(h = threshold_lib, col = "red")

    # Negative controls
    max_value <- max(ratios$`nc/gens`, threshold_line_nc)
    colores_nc <- ifelse(ratios$`nc/gens` < threshold_inferior_nc, "#4793AF",
                         ifelse(ratios$`nc/gens` > threshold_superior_nc, "red", "#FFC470"))
    plot(ratios$`nc/gens`, xlab = "", ylab = "nc/gens", col = colores_nc,
         xaxt = "n", pch = 19, main = "Negative Control (QC2)", ylim = c(0, max_value))
    axis(1, at = 1:nrow(ratios), labels = rownames(ratios), las = 2, cex.axis = 0.8)
    abline(h = threshold_line_nc, col = "red")


    # Genomic DNA
    max_value <- max(ratios$`gdna/gens`, threshold_line_gdna)
    colores_gdna <- ifelse(ratios$`gdna/gens` < threshold_inferior_gdna, "#4793AF",
                           ifelse(ratios$`gdna/gens` > threshold_superior_gdna, "red", "#FFC470"))
    plot(ratios$`gdna/gens`, xlab = "", ylab = "gdna/gens", col = colores_gdna,
         xaxt = "n", pch = 19, main = "Genomic DNA (QC3)", ylim = c(0, max_value))
    axis(1, at = 1:nrow(ratios), labels = rownames(ratios), las = 2, cex.axis = 0.8)
    abline(h = threshold_line_gdna, col = "red")


    # ERCC
    max_value <- max(ratios$`ERCC/gens`, threshold_line_ercc)
    colores_ercc <- ifelse(ratios$`ERCC/gens` < threshold_inferior_ercc, "#4793AF",
                           ifelse(ratios$`ERCC/gens` > threshold_superior_ercc, "red", "#FFC470"))
    plot(ratios$`ERCC/gens`, xlab = "", ylab = "ERCC", col = colores_ercc,
         xaxt = "n", pch = 19, main = "ERCC (QC4)", ylim = c(0, max_value))
    axis(1, at = 1:nrow(ratios), labels = rownames(ratios), las = 2, cex.axis = 0.8)
    abline(h = threshold_line_ercc, col = "red")


    ## Median
    max_value <- max(ratios$median, threshold_line_median)
    colores_med <- ifelse(ratios$median < threshold_inferior_median, "#4793AF",
                          ifelse(ratios$median > threshold_superior_median, "red", "#FFC470"))
    plot(ratios$median, xlab = "", ylab = "Median", col = colores_med,
         xaxt = "n", pch = 19, main = "Median", ylim = c(0, max_value))
    axis(1, at = 1:nrow(ratios), labels = rownames(ratios), las = 2, cex.axis = 0.8)
    abline(h = threshold_line_median, col = "red")

    # Create the heatmap with ggplot2
    a<- ggplot(bin_df_melted, aes(x = variable, y = Sample, fill = factor(value))) +
      geom_tile(color = "white") +
      scale_fill_manual(values = c("0" = "#FFF9D0", "1" = "red"), labels = c("OK", "Possible Outlier")) +
      labs(x = "QC Metrics", y = "Samples", fill = "QC Status") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = 7),
            legend.position = "bottom")
    print(a)
    print(cowplot::plot_grid(p1, p2, p3, nrow = 1, align = "hv"))
    rows_with_1 <- suppressWarnings(rownames(bin_matrix)[apply(bin_matrix, 1, any)])
    cat("\033[32m                              ***\033[0m\n")
    cat(paste("\033[32mThese are the samples plotted at least once in the heatmap:\033[0m\n"))
    print(rows_with_1)
    cat(paste("The number of samples that are outliers are:", length(rows_with_1)))
    return(rows_with_1)
dev.off()
dev.off()

  }
