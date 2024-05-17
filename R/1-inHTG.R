#' Import HTG counts data from an Excel file
#'
#' @param file_path Path to the Excel file containing HTG counts
#' @return A data frame with HTG counts data
#' @export
#' @examples
#' htg_data <- inHTG("path/to/your/excel/file.xlsx")
#'
#' @name HTG_import
HTG_import <- function(file_path) {
  library(readxl)
  htg_db <- readxl::read_excel(file_path)
  htg_db <- as.data.frame(htg_db)
  rownames_db <- htg_db$`Sample Name`
  htg_db <- htg_db[, -1]
  rownames(htg_db) <- rownames_db

  cts <- htg_db
  cts <- cts[-c(1:4), ]
  cts <- apply(cts, 2, as.numeric)
  rownames(cts) <- rownames(htg_db)[-c(1:4)]
  cts <- as.data.frame(cts)
  return(cts)
}

