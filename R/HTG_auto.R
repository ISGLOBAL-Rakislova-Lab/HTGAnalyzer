#' HTG_auto: Automated Pipeline for HTG Data Analysis
#'
#' @description
#' This function automates the process of analyzing HTG data by performing quality control, differential expression analysis, and generating plots. It integrates multiple steps including importing counts data, performing quality control (QC), and conducting differential expression analysis using DESeq2. The function also includes options for generating volcano plots, heatmaps, and performing GSEA, TME, and survival analysis.
#'
#' The pipeline serves to conduct Differential Expression Analysis (DEA), Gene Set Enrichment Analysis (GSEA), Tumor Microenvironment (TME) analysis, and survival analysis of the top 10 genes specified. If the data type is RNAseq, the function will skip the QC and outlier removal steps. The `genes_to_use` parameter is used in survival analysis. If DEA is performed (`DEA = TRUE`), the survival analysis will focus on the top 10 differentially expressed genes. If DEA is not performed, the analysis will use the genes specified in `genes_to_use`.
#'
#' @param counts_file_path Character. Path to the file containing the HTG counts data in Excel format.
#' @param file_type Character. Type of file being imported, either "HTG" or "RNAseq".
#' @param AnnotData_file_path Character. Path to the file containing the annotation data in Excel format. Sample id column have to be named "id"
#' @param pattern Character. A regular expression pattern to identify control probes in the count data. Default is "^NC-|^POS-|^GDNA-|^ERCC-" but for RNAseq it should be NULL.
#' @param QC Logical. Indicates whether to perform quality control on the data. Default is TRUE. If set to FALSE, the QC step is skipped (for RNAseq data).
#' @param threshold_superior_pos Numeric. Upper limit threshold for positive control ratio. Default is 5.
#' @param threshold_line_pos Numeric. Threshold line for positive control ratio. Default is 4.
#' @param threshold_inferior_lib Numeric. Lower limit threshold for library size. Default is 5e+06.
#' @param threshold_lib Numeric. Threshold line for library size. Default is 7e+06.
#' @param threshold_superior_nc Numeric. Upper limit threshold for negative control ratio. Default is 0.05.
#' @param threshold_line_nc Numeric. Threshold line for negative control ratio. Default is 0.045.
#' @param threshold_superior_gdna Numeric. Upper limit threshold for genomic DNA ratio. Default is 0.025.
#' @param threshold_line_gdna Numeric. Threshold line for genomic DNA ratio. Default is 0.02.
#' @param threshold_superior_ercc Numeric. Upper limit threshold for ERCC control ratio. Default is 0.03.
#' @param threshold_line_ercc Numeric. Threshold line for ERCC control ratio. Default is 0.025.
#' @param threshold_inferior_median Numeric. Lower limit threshold for median ratio. Default is 3.
#' @param threshold_line_median Numeric. Threshold line for median ratio. Default is 5.
#' @param save_csv Logical. Whether to save the QC results as a CSV file. Default is TRUE.
#' @param csv_file Character. The name of the CSV file to save the QC results if `save_csv` is TRUE. Default is "QC_results.csv".
#' @param design_formula Character. The design formula for DESeq2 analysis, specified as a string without the tilde (~).
#' @param percentage_gene Numeric. Minimum fraction of samples in which a gene must be expressed to be retained. Default is 0.2.
#' @param threshold_gene Numeric. Minimum count threshold per gene. Default is 200.
#' @param threshold_subject Numeric. Minimum count threshold per subject. Default is 10.
#' @param genes_to_use Character vector. Specifies top genes for analysis and survival analysis. Default is c("CCND1", "MMP10", "CTTN").
#' @param heatmap_columns Character vector. Specifies columns to be used for annotations in the heatmap.
#' @param contrast Character vector. Specifies the contrast for differential expression analysis. Default is NULL.
#' @param pCutoff Numeric. The p-value cutoff for generating the volcano plot. Default is 0.05.
#' @param variable_01 Character. Name of the survival event variable (e.g., "Recurrence_01"). Required for TME. Set NULL if you don't have it and turn FALSE TME
#' @param time Character. Name of the time variable (e.g., "Time_to_death_surv"). Required for survival analysis. Set NULL if you don't have it and turn FALSE survival_analysis
#' @param DEA Logical. Indicates whether to perform DESeq2 analysis with filtering and without Lfc shrinkage. Default is TRUE.
#' @param remove_outliers Logical. Indicates whether to remove outliers. Default is TRUE. Turn FALSE for RNAseq data.
#' @param GSEA Logical. Indicates whether to perform Gene Set Enrichment Analysis (GSEA). Default is TRUE.
#' @param generate_heatmap Logical. Indicates whether to generate a heatmap. Default is TRUE.
#' @param TME Logical. Indicates whether to perform Tumor Microenvironment (TME) analysis. Default is TRUE. Parameter variable_01 will be needed.
#' @param survival_analysis Logical. Indicates whether to perform survival analysis. Default is TRUE. Parameter time will be needed.
#'
#' @return This function performs the complete analysis pipeline, including:
#' \item{Quality Control (QC):}{Conducts QC on the imported data and optionally saves results as a CSV file.}
#' \item{Differential Expression Analysis:}{Performs DESeq2 analysis based on the provided design formula.}
#' \item{Gene Set Enrichment Analysis (GSEA):}{Performs GSEA and generates relevant plots.}
#' \item{Tumor Microenvironment (TME):}{Conducts TME analysis if specified.}
#' \item{Survival Analysis:}{Generates Kaplan-Meier survival plots if requested. The analysis will focus on the top 10 differentially expressed genes if DEA is performed. If DEA is not performed, the analysis will use the genes specified in `genes_to_use`.}
#' \item{Plots and Files:}{Saves plots (volcano, heatmap) and analysis results (Excel files) in the current working directory.}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- HTG_auto(
#'   counts_file_path = "path/to/your/file.xlsx",
#'   file_type = "HTG",
#'   AnnotData_file_path = "path/to/your/annot_file.xlsx",
#'   pattern = "^NC-|^POS-|^GDNA-|^ERCC-",
#'   design_formula = "Group",
#'   percentage_gene = 0.2,
#'   threshold_gene = 200,
#'   threshold_subject = 10,
#'   genes_to_use = c("CCND1", "MMP10", "CTTN"),
#'   heatmap_columns = c("Group", "Smoker"),
#'   contrast = c("Group", "high", "low"),
#'   pCutoff = 0.05,
#'   variable_01 = "Recurrence_01",
#'   time = "Time_to_event",
#'   DEA = TRUE,
#'   remove_outliers = TRUE,
#'   GSEA = FALSE,
#'   generate_heatmap = TRUE,
#'   TME = TRUE,
#'   survival_analysis = TRUE
#' )
#' }
#' @name HTG_auto


HTG_auto <- function(counts_file_path, file_type = "HTG",
                     AnnotData_file_path,
                     pattern = "^NC-|^POS-|^GDNA-|^ERCC-",
                     QC = TRUE,
                     threshold_superior_pos = 5,
                     threshold_line_pos = 4,
                     threshold_inferior_lib = 5e+06,
                     threshold_lib = 7e+06,
                     threshold_superior_nc = 0.05,
                     threshold_line_nc = 0.045,
                     threshold_superior_gdna = 0.025,
                     threshold_line_gdna = 0.02,
                     threshold_superior_ercc = 0.03,
                     threshold_line_ercc = 0.025,
                     threshold_inferior_median = 3,
                     threshold_line_median = 5,
                     save_csv = TRUE,
                     csv_file = "QC_results.csv",
                     design_formula = NULL,
                     percentage_gene = 0.2,
                     threshold_gene = 200,
                     threshold_subject = 10,
                     genes_to_use = c("CCND1", "MMP10", "CTTN"),
                     heatmap_columns = NULL,
                     contrast = NULL,
                     pCutoff = 5e-2,
                     variable_01 = NULL,
                     time = NULL,
                     DEA = TRUE,
                     remove_outliers = TRUE,
                     GSEA = TRUE,
                     generate_heatmap = TRUE,
                     TME = TRUE,
                     survival_analysis = TRUE){

  # Importar los datos
  cat("\033[33mINICIATING IMPORTATION OF DATA...\033[0m\n")
  counts_data <- HTG_import_counts(counts_file_path, file_type)
  counts_data <- as.data.frame(counts_data)
  col_data <- readxl::read_excel(AnnotData_file_path)
  col_data <- as.data.frame(col_data)

  # Normalizar los nombres de las columnas en col_data para que sean válidos en R
  colnames(col_data) <- gsub(" ", "_", colnames(col_data))

  # Obtener los nombres de las columnas en counts_data
  counts_colnames <- colnames(counts_data)
  rownames(col_data) <- col_data$id
  col_data_rownames <- rownames(col_data)

  common_ids <- base::intersect(counts_colnames, col_data_rownames)

  counts_data <- counts_data[, counts_colnames %in% common_ids]
  col_data <- col_data[common_ids, , drop = FALSE]

  cat("\033[32mDimensions of filtered counts_data:\033[0m\n")
  print(dim(counts_data))

  cat("\033[32mDimensions of filtered col_data:\033[0m\n")
  print(dim(col_data))

  if (ncol(counts_data) == nrow(col_data)) {
    cat("\033[32mThe filtered data sets are now aligned.\033[0m\n")
  } else {
    cat("\033[31mThere is a mismatch in the number of columns in counts_data_filtered and rows in col_data_filtered.\033[0m\n")
  }

  # Perform quality control if QC is TRUE
  if (QC) {
    cat("\033[33mPERFORMING QUALITY CONTROL...\033[0m\n")

    outliers <- HTG_QC(
      counts_data = counts_data,
      pattern = pattern,
      threshold_superior_pos = threshold_superior_pos,
      threshold_line_pos = threshold_line_pos,
      threshold_inferior_lib = threshold_inferior_lib,
      threshold_lib = threshold_lib,
      threshold_superior_nc = threshold_superior_nc,
      threshold_line_nc = threshold_line_nc,
      threshold_superior_gdna = threshold_superior_gdna,
      threshold_line_gdna = threshold_line_gdna,
      threshold_superior_ercc = threshold_superior_ercc,
      threshold_line_ercc = threshold_line_ercc,
      threshold_inferior_median = threshold_inferior_median,
      threshold_line_median = threshold_line_median,
      save_csv = save_csv,
      csv_file = csv_file
    )
  } else {
    outliers <- NULL
  }



  if (!is.null(design_formula)) {
    if (!all(all.vars(design_formula) %in% colnames(col_data))) {
      stop("Variables in design_formula must be columns in col_data.")
    }
  }

  cat("\n")
  cat("\033[33mPERFORMING ANALYSIS...\033[0m\n")

  result <- HTG_analysis(
    outliers = outliers,
    pattern = pattern,
    counts_data = counts_data,
    col_data = col_data,
    design_formula = design_formula,
    percentage_gene = percentage_gene,
    threshold_gene = threshold_gene,
    threshold_subject = threshold_subject,
    heatmap_columns = heatmap_columns,
    contrast = contrast,
    pCutoff = pCutoff,
    variable_01 = variable_01,
    time = time,
    DEA = DEA,
    remove_outliers = remove_outliers,
    GSEA = GSEA,
    generate_heatmap = generate_heatmap,
    TME = TME,
    survival_analysis = survival_analysis
  )

  return(result)
}
