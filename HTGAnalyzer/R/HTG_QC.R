#' HTG_QC: Quality Control for HTG EdgeSeq Data
#'
#' @description
#' This function performs various quality control (QC) checks for HTG EdgeSeq transcriptomic panel data. The QC checks include:
#'
#' QC0: Percentage of positive values < than 4%; QC1: Library size > than 7e+06; QC2: Negative control threshold < than 0.045; QC3: Genomic DNA threshold < than 0.02; QC4: ERCC threshold < than 0.025; Median:threshold > 5.
#'
#' In addition to these plots, highlighting potential outlier samples. The function also creates a data frame including the sum of each probe for each sample (total genes, positive, negative, gdna, and ercc)
#' the ratio for each sample and the size of each sample. Additionally, a statistical .csv is generated with columns for Min, Max, Mean, Median, Mode, SD, Variance, Range, Q1, Q3, IQR, Skewness, Kurtosis, Missing, and CV.
#'
#' This function also includes an optional heatmap to highlight potential outlier samples, which will be saved in a vector.
#'
#' The plots will saved in the current working directory.
#'
#' @param counts_data A data frame containing the HTG count data. The data must include probes that start with "^NC-|^POS-|^GDNA-|^ERCC-" for the function to work correctly.
#' @param pattern A regular expression pattern to identify control probes in the count data. For HTG data, this could be "^NC-|^POS-|^GDNA-|^ERCC-". If NULL, the pattern will not be applied.
#' @param threshold_superior_pos Threshold for upper limit of positive control ratio.
#' @param threshold_line_pos Threshold line for positive control ratio.
#' @param threshold_inferior_lib Threshold for lower limit of library size.
#' @param threshold_lib Threshold line for library size.
#' @param threshold_superior_nc Threshold for upper limit of negative control ratio.
#' @param threshold_line_nc Threshold line for negative control ratio.
#' @param threshold_superior_gdna Threshold for upper limit of genomic DNA ratio.
#' @param threshold_line_gdna Threshold line for genomic DNA ratio.
#' @param threshold_superior_ercc Threshold for upper limit of ERCC control ratio.
#' @param threshold_line_ercc Threshold line for ERCC control ratio.
#' @param threshold_inferior_median Threshold for lower limit of median ratio.
#' @param threshold_line_median Threshold line for median ratio.
#' @param save_csv Logical, whether to save the ratios as a CSV file. Default is FALSE.
#' @param csv_file The name of the CSV file to save the ratios if save_csv is TRUE. Default is "QC_results.csv".
#'
#' @return This function generates multiple plots displaying various QC metrics, including a violin plots, and saves an Excel file with all the ratios. Additionally, it identifies and returns the most probable outliers based on the QC analysis.
#'
#' @export
#'
#'
#' @examples
#' # Run the function with example data
#' HTG_QC(counts_data = counts_data_tutorial, pattern = "^NC-|^POS-|^GDNA-|^ERCC-", save_csv = TRUE)
#'
#' @name HTG_QC
#'
utils::globalVariables(c("PC1", "PC2", "Tag", "label", "Componente", "Porcentaje", "pc", "Sample", "LogTPM", "variable", "value"))
HTG_QC <- function(counts_data, pattern = "^NC-|^POS-|^GDNA-|^ERCC-",
                             threshold_superior_pos = 5,
                             threshold_line_pos = 4,
                             threshold_inferior_lib = 5e+06,
                             threshold_lib = 7e+06,
                             threshold_superior_nc = 0.05,
                             threshold_line_nc = 0.045,
                             threshold_superior_gdna = 0.025,
                             threshold_line_gdna = 0.02,
                             threshold_superior_ercc = 0.03,
                             threshold_line_ercc = 0.025,
                             threshold_inferior_median = 3,
                             threshold_line_median = 5,
                             save_csv = TRUE,
                             csv_file = "QC_results.csv") {

  # Filter counts_data data
  cat("\033[33mINITIATING DATA FILTERING...\033[0m\n")
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
  cat("\033[32mSummary statistics saved as 'summary_stats.csv'\033[0m\n")

  cat("\033[33mINITIATING QC PLOTS...\033[0m\n")
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
  summary_stats <- HTG_calculate_summary_stats(counts_filtered, pattern= pattern)
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
    cat("\033[32mQC DATA SAVED AS '", csv_file, "'\033[0m\n")
  }
###
  library_size <- colSums(counts_filtered)

  # Create dataframe of library size
  lib_s2 <- data.frame(Sample = colnames(counts_filtered), Size = library_size)
  ratios_heat <- as.matrix(ratios)

  # Add a fourth column to ratios_heat with library sizes from lib_s2
  ratios_heat <- cbind(ratios_heat, Size = "")
  ratios_heat[, "Size"] <- lib_s2[match(rownames(ratios_heat), rownames(lib_s2)), "Size"]
  ratios_heat <- as.data.frame(ratios_heat)

    # Convert values of ratios_heat to numeric
  cols_to_convert <- c("total_POS", "total_GDNA", "total_gens", "total_NC", "total_ERCC",
                       "pos/gens", "gdna/gens", "nc/gens", "ERCC/gens", "median", "Size")
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
  bin_matrix <- matrix(0, nrow = nrow(ratios_heat), ncol = 6)
  for (i in 1:nrow(ratios_heat)) {
    bin_matrix[i, 1] <- assign_01_QC(ratios_heat[i, "pos/gens"], threshold_line_pos)
    bin_matrix[i, 2] <- assign_01_size(ratios_heat[i, "Size"], threshold_lib)
    bin_matrix[i, 3] <- assign_01_QC(ratios_heat[i, "nc/gens"], threshold_line_nc)
    bin_matrix[i, 4] <- assign_01_QC(ratios_heat[i, "gdna/gens"], threshold_line_gdna)
    bin_matrix[i, 5] <- assign_01_QC(ratios_heat[i, "ERCC/gens"], threshold_line_ercc)
    bin_matrix[i, 6] <- assign_01_size(ratios_heat[i, "median"], threshold_line_median)
  }

  # Row and column names
  rownames(bin_matrix) <- rownames(ratios_heat)
  colnames(bin_matrix) <- c("QC0", "QC1", "QC2", "QC3", "QC4", "QC5")

  # Convert the matrix to a data frame for ggplot2
  bin_df <- as.data.frame(bin_matrix)
  bin_df$Sample <- rownames(bin_df)
  # Melt the data frame
  bin_df_melted <- reshape2::melt(bin_df, id.vars = "Sample")

  # VIOLIN PLOT
  # Convert raw counts to TPM
  tpm_counts <- IOBR::count2tpm(counts_data,
                          idType = "Symbol",
                          org = "hsa",
                          source = "biomart")

  # TPM data formatting
  tpm_counts$Gene <- rownames(tpm_counts)
  tpm_long <- reshape2::melt(tpm_counts, id.vars = "Gene", variable.name = "Sample", value.name = "TPM")
  tpm_long$LogTPM <- log1p(tpm_long$TPM)

  # Calculate the 95th percentile threshold for each sample
  #percentile_95 <- aggregate(LogTPM ~ Sample, data = tpm_long, FUN = function(x) quantile(x, 0.95))
  percentile_95 <- dplyr::summarise(dplyr::group_by(tpm_long, Sample),
    percentile_95 = quantile(LogTPM, 0.95))


  colnames(percentile_95)[2] <- "Threshold"
  tpm_with_threshold <- merge(tpm_long, percentile_95, by = "Sample")

  # Filter data based on the 95th percentile threshold
  tpm_filtered <- tpm_with_threshold[tpm_with_threshold$LogTPM < tpm_with_threshold$Threshold, ]

  # Function to create violin plot
  create_violin_plot <- function(data, title) {
    ggplot2::ggplot(data, ggplot2::aes(x = Sample, y = LogTPM)) +
      ggplot2::geom_violin(trim = FALSE, fill = "#4793AF", color = "black") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                     panel.background = ggplot2::element_rect(fill = "white"),
                     plot.background = ggplot2::element_rect(fill = "white"),
                     panel.grid.major = ggplot2::element_line(color = "gray80"),
                     panel.grid.minor = ggplot2::element_line(color = "gray90")) +
      ggplot2::labs(title = title,
                    x = "Sample",
                    y = "Log-Transformed TPM")
  }

  # Check if there are more than 50 unique samples and split if needed
  if (length(unique(tpm_filtered$Sample)) > 50) {
    samples <- unique(tpm_filtered$Sample)
    half <- ceiling(length(samples) / 2)
    subset1 <- samples[1:half]
    subset2 <- samples[(half + 1):length(samples)]
    data1 <- tpm_filtered[tpm_filtered$Sample %in% subset1, ]
    data2 <- tpm_filtered[tpm_filtered$Sample %in% subset2, ]
    p4 <- create_violin_plot(data1, "Distribution of Log-Transformed TPM (Up to 95th Percentile) - Part 1")
    p5 <- create_violin_plot(data2, "Distribution of Log-Transformed TPM (Up to 95th Percentile) - Part 2")
    combined_plot <- ggpubr::ggarrange(p4, p5, ncol = 1, nrow = 2)
  } else {
    p6 <- create_violin_plot(tpm_filtered, "Distribution of Log-Transformed TPM (Up to 95th Percentile)")
    print(summary(tpm_filtered))
    combined_plot <- p6
  }

  # Save the violin plot to a PDF
  pdf("QC_plots_violin_plot.pdf", width = 14, height = 10)
  print(combined_plot)
  dev.off()


  cat("\033[32mViolin plot saved as 'QC_plots_violin_plot.pdf'\033[0m\n")


  pdf("QC_plots.pdf")

  # Positive controls
  max_value <- max(ratios$`pos/gens`, threshold_line_pos)
  min_size <- min(ratios$`pos/gens`, threshold_line_pos)

  colores_pos <- ifelse(ratios$`pos/gens` <= threshold_line_pos, "#4793AF",
                        ifelse(ratios$`pos/gens` <= threshold_superior_pos, "#FFC470", "red"))
  plot(ratios$`pos/gens`, xlab = "", ylab = "pos/gens", col = colores_pos,
       xaxt = "n", pch = 19, main = "Positive control 4% (QC0)", ylim = c(0, max_value))
  axis(1, at = 1:nrow(ratios), labels = rownames(ratios), las = 2, cex.axis = 0.8)
  abline(h = threshold_line_pos, col = "red")


  # Library size
  max_size <- max(lib_s2$Size, threshold_lib)
  min_size <- min(lib_s2$Size, threshold_lib)
    colores <- ifelse(lib_s2$Size < threshold_inferior_lib, "red",
                    ifelse(lib_s2$Size <= threshold_lib, "#FFC470", "#4793AF"))
  plot(lib_s2$Size, xlab = "", ylab = "Library Size", col = colores,
       xaxt = "n", pch = 19, main = "Library Size per Sample (QC1)", cex.axis = 0.8,
       ylim = c(min_size, max_size))
  axis(1, at = 1:length(lib_s2$Sample), labels = lib_s2$Sample, las = 2, cex.axis = 0.8)
  abline(h = threshold_lib, col = "red")


  # Negative controls
  max_value <- max(ratios$`nc/gens`, threshold_line_nc)
  min_value <- min(ratios$`nc/gens`, threshold_line_nc)
  colores_nc <- ifelse(ratios$`nc/gens` <= threshold_line_nc, "#4793AF",
                       ifelse(ratios$`nc/gens` <= threshold_superior_nc, "#FFC470", "red"))
  plot(ratios$`nc/gens`, xlab = "", ylab = "nc/gens", col = colores_nc,
       xaxt = "n", pch = 19, main = "Negative Control (QC2)", ylim = c(min_value, max_value))
  axis(1, at = 1:nrow(ratios), labels = rownames(ratios), las = 2, cex.axis = 0.8)
  abline(h = threshold_line_nc, col = "red")


  # Genomic DNA
  max_value <- max(ratios$`gdna/gens`, threshold_line_gdna)
  min_size <- min(ratios$`gdna/gens`, threshold_line_gdna)
  colores_gdna <- ifelse(ratios$`gdna/gens` <= threshold_line_gdna, "#4793AF",
                         ifelse(ratios$`gdna/gens` <= threshold_superior_gdna, "#FFC470", "red"))
  plot(ratios$`gdna/gens`, xlab = "", ylab = "gdna/gens", col = colores_gdna,
       xaxt = "n", pch = 19, main = "Genomic DNA (QC3)", ylim = c(0, max_value))
  axis(1, at = 1:nrow(ratios), labels = rownames(ratios), las = 2, cex.axis = 0.8)
  abline(h = threshold_line_gdna, col = "red")

  # ERCC
  max_value <- max(ratios$`ERCC/gens`, threshold_line_ercc)
  min_size <- min(ratios$`ERCC/gens`, threshold_line_ercc)

  colores_ercc <- ifelse(ratios$`ERCC/gens` <= threshold_line_ercc, "#4793AF",
                         ifelse(ratios$`ERCC/gens` <= threshold_superior_ercc, "#FFC470", "red"))
  plot(ratios$`ERCC/gens`, xlab = "", ylab = "ERCC", col = colores_ercc,
       xaxt = "n", pch = 19, main = "ERCC (QC4)", ylim = c(0, max_value))
  axis(1, at = 1:nrow(ratios), labels = rownames(ratios), las = 2, cex.axis = 0.8)
  abline(h = threshold_line_ercc, col = "red")


  ## Median
  max_value <- max(ratios$median, threshold_line_median)
  min_size <- min(ratios$median, threshold_line_median)
  colores_med <- ifelse(ratios$median < threshold_inferior_median, "red",
                        ifelse(ratios$median <= threshold_line_median, "#FFC470", "#4793AF"))
  plot(ratios$median, xlab = "", ylab = "Median", col = colores_med,
       xaxt = "n", pch = 19, main = "Median (QC5)", ylim = c(0, max_value))
  axis(1, at = 1:nrow(ratios), labels = rownames(ratios), las = 2, cex.axis = 0.8)
  abline(h = threshold_line_median, col = "red")
    dev.off()

    cat("\033[32mQC plots saved as 'plot_QC.pdf'\033[0m\n")

  # Create the heatmap with ggplot2
  a<- ggplot2::ggplot(bin_df_melted, ggplot2::aes(x = variable, y = Sample, fill = factor(value))) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_manual(values = c("0" = "#FFF9D0", "1" = "red"), labels = c("OK", "Possible Outlier")) +
    ggplot2::labs(x = "QC Metrics", y = "Samples", fill = "QC Status") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          axis.text.y = ggplot2::element_text(size = 7),
          legend.position = "bottom")

  pdf("HTG_QC_heatmap.pdf", width = 10, height = 14)
  print(a)
      dev.off()

  cat("\033[32mHeatmap PLOTS SAVED AS 'HTG_QC_heatmap.pdf'\033[0m\n")

  rows_with_1 <- suppressWarnings(rownames(bin_matrix)[apply(bin_matrix, 1, any)])
  cat("\033[32m                              ***\033[0m\n")
  cat(paste("\033[32mThese are the samples plotted at least once in the heatmap:\033[0m\n"))
  cat(paste("The number of samples that are outliers are:", length(rows_with_1)))
  cat(rows_with_1)
  return(rows_with_1)
}
