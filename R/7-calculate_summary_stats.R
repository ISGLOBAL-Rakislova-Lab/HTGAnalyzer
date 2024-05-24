#' Calculate Summary Statistics
#'
#' @description This function calculates summary statistics including minimum, maximum, mean, and median for each column of the input counts_filtered.
#'
#' @param counts_filtered A counts_filtered frame containing the input counts_filtered.
#'
#' @return A counts_filtered frame containing summary statistics for each column of the input counts_filtered.
#' @export
#'
#' @examples
#' HTG_calculate_summary_stats(counts_filtered)
#' @name HTG_calculate_summary_stats
HTG_calculate_summary_stats <- function(counts_filtered) {
  min_values <- apply(counts_filtered, 2, min)
  max_values <- apply(counts_filtered, 2, max)
  mean_values <- apply(counts_filtered, 2, mean)
  median_values <- apply(counts_filtered, 2, median)

  summary_stats <- data.frame(
    Min = min_values,
    Max = max_values,
    Mean = mean_values,
    Median = median_values
  )

  return(summary_stats)
}
