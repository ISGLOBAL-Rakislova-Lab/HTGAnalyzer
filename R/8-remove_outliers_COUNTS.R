#' Remove outliers from counts data based on specified sample IDs
#'
#' This function removes specified sample IDs from a data frame containing counts data.
#'
#' @param counts_data A data frame containing counts data.
#' @param outlier_rows A vector of sample IDs (column names) to be removed from the counts data frame.
#'
#' @return Returns a data frame with the outlier columns removed.
#'
#' @export
#'
#' @examples
#' outlier_rows <- c("B04-23456-B", "B16-39493-1")
#' cleaned_counts <- remove_outliers_COUNTS(counts_data, outlier_rows)
#' @name HTG_remove_outliers_COUNTS

HTG_remove_outliers_COUNTS <- function(counts_data, outlier_rows) {
  counts_filtered <- counts_data[, !colnames(counts_data) %in% outlier_rows]
}
