#' HTG_calculate_ratios
#'
#' @description  This function calculates ratios based on counts data for different categories such as positive controls, genomic DNA, etc.
#'
#' @param counts_filtered A data frame containing filtered counts data (Counts data without probes.)
#' @param cts_POS A data frame containing counts data for positive controls.
#' @param cts_NC A data frame containing counts data for negative controls.
#' @param cts_GDNA A data frame containing counts data for genomic DNA.
#' @param cts_ERCC A data frame containing counts data for ERCC controls.
#' @param save_csv Logical, indicating whether to save the calculated ratios as a CSV file. Default is FALSE.
#' @param csv_file Character string specifying the file name if `save_csv` is TRUE. Default is "ratios.csv".
#'
#' @return A data frame containing calculated ratios for different categories.
#' @export
#'
#' @examples
#' HTG_calculate_ratios(counts_filtered, cts_POS, cts_NC, cts_GDNA, cts_ERCC,
#'                      save_csv = TRUE, csv_file = "ratios.csv")
#' @name HTG_calculate_ratios

HTG_calculate_ratios <- function(counts_filtered, cts_POS, cts_NC, cts_GDNA, cts_ERCC,
                                 save_csv = FALSE, csv_file = "ratios.csv") {
  # Calculate total counts for each category
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
  ratios$"pos/gens" <- (ratios$total_POS / ratios$total_gens) * 100
  ratios$"gdna/gens" <- (ratios$total_GDNA / ratios$total_gens) * 100
  ratios$"nc/gens" <- (ratios$total_NC / ratios$total_gens) * 100
  ratios$"ERCC/gens" <- (ratios$total_ERCC / ratios$total_gens) * 100

  # Add median column
  summary_stats <- HTG_calculate_summary_stats(counts_filtered)
  ratios$median <- summary_stats$Median
  ratios$min<- summary_stats$Min
  ratios$max<- summary_stats$Max
  ratios$mean<- summary_stats$Mean

  # Add sample names as a factor column
  ratios$samples <- factor(rownames(ratios))

  # Optionally save as CSV
  if (save_csv) {
    write.csv(ratios, csv_file, row.names = TRUE)
    cat("Ratios saved as", csv_file, "\n")
  }

  return(ratios)
}
