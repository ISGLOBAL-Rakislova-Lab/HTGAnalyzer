#' Subset Counts Data Frame
#'
#' This function subsets a counts data frame based on a specified prefix in row names.
#'
#' @param counts A data frame containing counts data.
#' @param prefix A prefix string used to subset rows from the counts data frame.
#'
#' @return A subset of the counts data frame containing rows with row names that match the specified prefix.
#' @export
#'
#' @examples
#' cts_ERCC <- subset_counts(counts, "ERCC")
HTG_subset_counts <- function(counts, prefix) {
  subset_df <- as.data.frame(subset(counts, grepl(paste0("^", prefix, "-"), rownames(counts))))
  cat("Dimensions:", dim(subset_df), "\n")
  return(subset_df)
}
