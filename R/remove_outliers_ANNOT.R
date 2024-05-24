#' HTG_remove_outliers_ANNOT
#'
#' @description This function removes specified sample IDs from a data frame containing annotation data.Remove outliers from annotation data based on specified sample IDs
#'
#' @param AnnotData A data frame containing annotation data.
#' @param sample_id_col The name of the column containing sample IDs.
#' @param outlier_rows A vector of sample IDs (values in the specified column) to be removed from the annotation data frame.
#'
#' @return Returns a data frame with the rows containing outlier sample IDs removed.
#'
#' @export
#'
#' @examples
#' # Remove outliers from AnnotData based on outlier_rows
#' outlier_rows <- c("sample2", "sample3")
#' remove_outliers_ANNOT <- HTG_remove_outliers_ANNOT(AnnotData, "id", outlier_rows)
#' @name HTG_remove_outliers_ANNOT
HTG_remove_outliers_ANNOT <- function(AnnotData, sample_id_col, outlier_rows) {
  AnnotData <- AnnotData[!AnnotData[[sample_id_col]] %in% outlier_rows, ]
  return(AnnotData)
}
