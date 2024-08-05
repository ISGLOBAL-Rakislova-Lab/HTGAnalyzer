#' HTG_subset
#'
#' @description This function subsets a data frame based on a specified prefix in the row names. It allows you to extract specific rows that match a given pattern and optionally normalizes the data using TPM (Transcripts Per Million) before subsetting. The function also provides the dimensions of the resulting data frame.
#'
#' @param data_frame A numeric data frame or matrix containing count data. The row names should include identifiers used for subsetting.
#' @param prefix A string specifying the prefix used to match and subset rows from the data frame.
#' @param normalize A logical value indicating whether to normalize the data using TPM before subsetting. Default is FALSE.
#'
#' @return A data frame containing rows with row names that match the specified prefix. If `normalize` is TRUE, the data will be normalized before subsetting.
#' @export
#'
#' @examples
#' # Subset rows with the prefix "ERCC"
#' cts_ERCC <- HTG_subset(counts_data, "ERCC")
#'
#' # Subset rows with the prefix "POS"
#' cts_POS <- HTG_subset(counts_data, "POS")
#'
#' # Subset rows with the prefix "NC"
#' cts_NC <- HTG_subset(counts_data, "NC")
#'
#'# Subset rows with the prefix "GDNA"
#' cts_GDNA <- HTG_subset(counts_data, "GDNA")
#'
#' # Subset rows with the prefix "ZZZ3"
#' cts_ZZZ3 <- HTG_subset(counts_data, "ZZZ3")
#'
#' #' # Subset rows with the prefix "ZZZ3" and normalize the data
#' cts_ZZZ3 <- HTG_subset(counts_data, "ZZZ3", normalize = TRUE)
#'
#' resAAAS<- HTG_subset(res, "AAAS")
#'
#' @name HTG_subset
#'
HTG_subset <- function(data_frame, prefix, normalize = FALSE) {
  data_frame <- as.data.frame(data_frame)

  if (!all(sapply(data_frame, is.numeric))) {
    stop("Input must be a numeric data frame or matrix.")
  }

  if (normalize) {
    cat("\033[32mNormalizing data using TPM (Transcripts Per Million)...\033[0m\n")

    # Load the IOBR library, handle messages and warnings
    suppressMessages(library(IOBR))

    # Perform normalization and handle possible warnings
    tpm_counts <- suppressWarnings(
      suppressMessages(count2tpm(data_frame, idType = "Symbol", org = "hsa", source = "biomart"))
    )

    # Subset the normalized data based on the prefix
    subset_df <- as.data.frame(base::subset(tpm_counts, grepl(paste0("^", prefix), rownames(tpm_counts))))
  } else {
    cat("\033[32mSubsetting without normalization...\033[0m\n")

    # Perform subsetting without normalization
    subset_df <- as.data.frame(base::subset(data_frame, grepl(paste0("^", prefix), rownames(data_frame))))
  }

  # Print the dimensions of the subsetted data
  cat("Dimensions:", dim(subset_df), "\n")

  # Return the subsetted data
  return(subset_df)
}

