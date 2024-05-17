#' Calculate Summary Statistics
#'
#' This function calculates summary statistics including minimum, maximum, mean, and median for each column of the input data.
#'
#' @param data A data frame containing the input data.
#'
#' @return A data frame containing summary statistics for each column of the input data.
#' @export
#'
#' @examples
#' calculate_summary_stats(data)
#' @name HTG_calculate_summary_stats
HTG_calculate_summary_stats <- function(data) {
  min_values <- apply(data, 2, min)
  max_values <- apply(data, 2, max)
  mean_values <- apply(data, 2, mean)
  median_values <- apply(data, 2, median)

  summary_stats <- data.frame(
    Min = min_values,
    Max = max_values,
    Mean = mean_values,
    Median = median_values
  )

  return(summary_stats)
}
