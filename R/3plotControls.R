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
                         threshold_line_median = 5) {

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
}

