#' Filter Counts Data
#'
#' @description This function filters counts data to remove rows with specific prefixes such as "NC-", "POS-", "GDNA-", and "ERCC-".
#'
#' @param counts A data frame containing counts data with genes as rows and samples as columns.
#'
#' @return A filtered counts data frame with specified prefixes removed.
#' @export
#'
#' @examples
#' counts_filtered <- HTG_filterCounts(counts)
#'@name HTG_filterCounts
#'
HTG_filterCounts <- function(counts) {
  counts_filtered <- subset(counts, !grepl("^NC-|^POS-|^GDNA-|^ERCC-", rownames(counts)))
  return(counts_filtered)
}
