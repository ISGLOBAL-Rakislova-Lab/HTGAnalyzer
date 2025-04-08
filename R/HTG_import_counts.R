#' HTG_import_counts
#'
#' @description  Import counts data from an Excel file. It can be either HTG excel or RNAseq.
#' The first row in the Excel file must contain the column headers, with "id" as the first column header followed by the names of each sample.
#' Please note that you might need to modify the Excel file to ensure it is in the correct format for importing data into R. Ensure column names are free of special characters and spaces to avoid import issues.
#'
#' @param file_path Path to the Excel file containing HTG counts. All the columns have to had a name in excel file.
#' @param file_type Type of file being imported, either "HTG" or "RNAseq".
#' @return A data frame with counts data.
#' @export
#' @importFrom readxl read_excel
#'
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

if (file_type == "HTG") {
  # Leer y procesar archivo HTG
  htg_db <- readxl::read_excel(file_path)
  htg_db <- as.data.frame(htg_db)
  htg_db <- gsub("[[:space:]-]", "_", rownames(htg_db))
  rownames_db <- htg_db[[1]]
  htg_db <- htg_db[, -1]
  rownames(htg_db) <- rownames_db

  cts <- htg_db
  if (grepl("^Sample ID", rownames(cts)[1])) {
    cts <- cts[-c(1:4), ]
    cts <- apply(cts, 2, as.numeric)
    rownames(cts) <- rownames(htg_db)[-c(1:4)]
    cts <- as.data.frame(cts)
  }
  counts_data <- cts
} else if (file_type == "RNAseq") {
  # Leer y procesar archivo RNAseq
  rna_db <- readxl::read_excel(file_path)
  rna_db <- as.data.frame(rna_db)
  rna_db <- gsub("[[:space:]-]", "_", rownames(rna_db))
  cat("First column will become rownames")
  rownames(rna_db) <- rna_db[[1]]
  rna_db <- rna_db[,-1]
  counts_data <- rna_db
  
} else {
  stop("File_type has to be 'HTG' or 'RNAseq'")
}
}
