#' HTG_subset
#'
#' @description This function subsets a data frame based on a specified prefix in the row names. It allows you to extract specific rows that match a given pattern.
#'
#' @param data_frame A data frame from which rows need to be subsetted.
#' @param prefix A string specifying the prefix used to match and subset rows from the data frame.
#'
#' @return A subset of the data frame containing rows with row names that match the specified prefix.
#' @export
#'
#' @examples
#' # Subset rows with the prefix "ERCC"
#' cts_ERCC <- HTG_subset(counts, "ERCC")
#'
#' # Subset rows with the prefix "POS"
#' cts_POS <- HTG_subset(counts, "POS")
#'
#' # Subset rows with the prefix "NC"
#' cts_NC <- HTG_subset(counts, "NC")
#'
#' # Subset rows with the prefix "GDNA"
#' cts_GDNA <- HTG_subset(counts, "GDNA")
#' @name HTG_subset
#'
HTG_subset <- function(data_frame, prefix) {
  subset_df <- as.data.frame(subset(data_frame, grepl(paste0("^", prefix), rownames(data_frame))))
  cat("Dimensions:", dim(subset_df), "\n")
  return(subset_df)
}
