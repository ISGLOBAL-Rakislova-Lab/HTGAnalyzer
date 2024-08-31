#' @title Data Documentation for Tutorial Package
#'
#' @description This file contains documentation for the data sets included in the tutorial package.
#' The datasets provided are used for demonstrating various data analysis techniques.
#'
#' @details
#' The following datasets are included:
#' \itemize{
#'   \item \code{counts_data_tutorial}: A matrix of gene expression counts.
#'   \item \code{AnnotData_tutorial}: Annotation data for the samples.
#'   \item \code{outliers_tutorial}: A character vector of outlier sample IDs.
#'   \item \code{res_tutorial}: Results from differential expression analysis.
#'   \item \code{TME_data_tutorial}: Tumor microenvironment data including different cell types and scores.
#' }
#'
#' @name data_documentation
#' @docType data
#' @keywords datasets
#'
#' @examples
#' # Load data
#' data(counts_data_tutorial)
#' data(AnnotData_tutorial)
#' data(outliers_tutorial)
#' data(res_tutorial)
#' data(TME_data_tutorial)
#'
#' # Inspect data
#' head(counts_data_tutorial)
#' head(AnnotData_tutorial)
#' head(outliers_tutorial)
#' head(res_tutorial)
#' head(TME_data_tutorial$EPIC)
#'
NULL

#' @name counts_data_tutorial
#' @title Gene Expression Counts Data
#' @description A matrix of gene expression counts with genes as rows and samples as columns.
#' @format A matrix with dimensions 19612 x 11.
#' @source Generated for tutorial purposes.
NULL

#' @name AnnotData_tutorial
#' @title Sample Annotation Data
#' @description A data frame containing annotation data for samples, including HPV status, FIGO stage, and recurrence information.
#' @format A data frame with 11 rows and 7 columns:
#' \itemize{
#'   \item \code{id}: Sample ID
#'   \item \code{HPV_status}: HPV status of the sample (Positive/Negative)
#'   \item \code{Ciclina_D1}: Cyclin D1 expression
#'   \item \code{FIGO_2021_STAGE}: FIGO stage of the sample
#'   \item \code{Recurrence}: Recurrence status (yes/no)
#'   \item \code{Recurrence_01}: Recurrence status as a binary variable (1/0)
#'   \item \code{Time_to_death_surv}: Time to death or survival
#' }
#' @source Generated for tutorial purposes.
NULL

#' @name outliers_tutorial
#' @title Outlier Sample IDs
#' @description A character vector containing the IDs of samples identified as outliers.
#' @format A character vector with sample IDs.
#' @source Generated for tutorial purposes.
NULL

#' @name res_tutorial
#' @title Differential Expression Results
#' @description A data frame containing results from differential expression analysis.
#' @format A data frame with 6 rows and 5 columns:
#' \itemize{
#'   \item \code{baseMean}: Average expression level
#'   \item \code{log2FoldChange}: Log2 fold change
#'   \item \code{lfcSE}: Standard error of log2 fold change
#'   \item \code{pvalue}: Wald test p-value
#'   \item \code{padj}: Adjusted p-value
#' }
#' @source Generated for tutorial purposes.
NULL

#' @name TME_data_tutorial
#' @title Tumor Microenvironment Data
#' @description A list containing different tumor microenvironment scores and cell type abundances.
#' @format A list with the following components:
#' \itemize{
#'   \item \code{EPIC}: Data frame with cell type abundances estimated by EPIC.
#'   \item \code{QTI}: Data frame with cell type abundances estimated by QTI.
#'   \item \code{XCELL}: Data frame with cell type abundances estimated by XCELL.
#' }
#' @source Generated for tutorial purposes.
NULL
