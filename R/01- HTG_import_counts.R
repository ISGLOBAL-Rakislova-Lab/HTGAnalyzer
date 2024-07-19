#' HTG_import_counts
#'
#' @description Import counts data from an Excel file. It can be either HTG excel or RNAseq.
#' @param file_path Path to the Excel file containing HTG counts. All the columns have to had a name in excel file.
#' @param file_type Type of file being imported, either "HTG" or "RNAseq".
#' @return A data frame with HTG counts data.
#' @export
#' @examples
#' \dontrun{
#' # Replace "path/to/your/excel/file.xlsx" with the actual path to your Excel file
#' path <- "path/to/your/excel/file.xlsx"
#' htg_data <- HTG_import_counts(path, "HTG")
#' rna_data <- HTG_import_counts(path, "RNAseq")
#' }
#' @name HTG_import_counts
#'
#'
HTG_import_counts <- function(file_path, file_type) {
  library(readxl)

  if (file_type == "HTG") {
    htg_db <- readxl::read_excel(file_path)
    htg_db <- as.data.frame(htg_db)
    rownames_db <- htg_db[[1]]
    htg_db <- htg_db[, -1]
    rownames(htg_db) <- rownames_db

    cts <- htg_db
    cts <- cts[-c(1:4), ]
    cts <- apply(cts, 2, as.numeric)
    rownames(cts) <- rownames(htg_db)[-c(1:4)]
    cts <- as.data.frame(cts)
    return(cts)
  } else if (file_type == "RNAseq") {
    rna_db <- readxl::read_excel(file_path)
    rna_db <- as.data.frame(rna_db)
    return(rna_db)
  } else {
    stop("file_type has to be 'HTG' or 'RNAseq'")
  }
}

