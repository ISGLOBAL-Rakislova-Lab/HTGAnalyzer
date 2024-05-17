#' Remove outliers from annotation data based on specified sample IDs
#'
#' This function removes specified sample IDs from a data frame containing annotation data.
#'
#' @param annot_data A data frame containing annotation data.
#' @param sample_id_col The name of the column containing sample IDs.
#' @param outlier_rows A vector of sample IDs (values in the specified column) to be removed from the annotation data frame.
#'
#' @return Returns a data frame with the rows containing outlier sample IDs removed.
#'
#' @export
#'
#' @examples
#' # Remove outliers from annot_data based on outlier_rows
#' outlier_rows <- c("sample2", "sample3")
#' cleaned_annot <- remove_outliers_ANNOT(annot_data, "SampleID", outlier_rows)
HTG_remove_outliers_ANNOT <- function(annot_data, sample_id_col, outlier_rows) {
  annot_data <- annot_data[!annot_data[[sample_id_col]] %in% outlier_rows, ]
  return(annot_data)
}
