#' HTG_filter
#'
#' @description This function filters out rows in a data frame that have specified prefixes in their row names.
#'
#' @param data_frame A data frame containing the data to be filtered.
#' @param pattern A string specifying the pattern of prefixes to filter out. The default pattern is "^NC-|^POS-|^GDNA-|^ERCC-", which will remove rows starting with these prefixes.
#'
#' @return A filtered data frame with rows that match the specified prefixes removed.
#' @export
#'
#' @examples
#' # Filter out rows with default prefixes
#' counts_filtered <- HTG_filter(counts_data)
#'
#' # Filter out rows with specific prefixes
#' counts_filtered <- HTG_filter(counts_data, "^NC-|^POS-")
#'
#' # Filter out rows with a specific pattern
#' counts_filtered <- HTG_filter(counts_data, "A1BG")
#' @name HTG_filter
#'
HTG_filter<- function(data_frame, pattern = "^NC-|^POS-|^GDNA-|^ERCC-") {
  data_frame_filtered <- subset(data_frame, !grepl(pattern, rownames(data_frame)))
  return(data_frame_filtered)
}
