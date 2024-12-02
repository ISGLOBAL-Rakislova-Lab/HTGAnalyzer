#' HTG_quant_to_qual
#'
#' @description
#' This function converts quantitative columns in a data frame into qualitative ones based on a specified threshold.
#'
#' @param data A data frame containing the quantitative column to be transformed.
#' @param column The name of the quantitative column to be transformed.
#' @param threshold The threshold value to classify the quantitative values into qualitative categories.
#' @param above_label The label for values above the threshold.
#' @param below_label The label for values below or equal to the threshold.
#'
#' @return Returns the modified data frame with the specified column transformed into qualitative values. The name of the new column will be the same as the column that will be transformed but with a 2 at the end.
#'
#' @export
#'
#' @examples
#' # Transform the "Ciclina_D1" column in the "data" data frame with threshold 60
#' AnnotData_tutorial <- HTG_quant_to_qual(AnnotData_tutorial, "Cyclin_D1", 60, "high", "low")
#'
#'@name HTG_quant_to_qual


HTG_quant_to_qual <- function(data, column, threshold, above_label = "yes", below_label = "no") {
  if (!is.data.frame(data)) {
    stop("This is not a data.table. Transform into data.table")
  }
    if (!column %in% colnames(data)) {
    stop(paste("Column", column, "doen't exist"))
  }

  if (!is.numeric(data[[column]])) {
    stop(paste("Column", column, "is not numeric."))
  }
    new_column <- paste0(column, "2")
    data[[new_column]] <- ifelse(data[[column]] > threshold, above_label, below_label)
  return(data)
}
