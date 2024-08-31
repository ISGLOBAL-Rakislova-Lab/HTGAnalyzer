#' @title HTGAnalyzer Package
#' @description This package performs HTG analysis and various associated tasks.
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix counts DESeq vst results resultsNames lfcShrink plotMA
#' @importFrom SummarizedExperiment assay
#' @importFrom pheatmap pheatmap
#' @importFrom EnhancedVolcano EnhancedVolcano
#' @importFrom ggplot2 ggplot aes geom_point labs theme element_text element_blank scale_color_manual scale_x_continuous scale_y_continuous
#' @importFrom grDevices pdf dev.off colorRampPalette
#' @importFrom graphics boxplot abline axis legend par text plot
#' @importFrom RColorBrewer brewer.pal
#' @importFrom utils write.csv capture.output head globalVariables data
#' @importFrom maxstat maxstat.test
#' @importFrom survival Surv survfit survdiff
#' @importFrom stats as.formula t.test lm prcomp na.omit quantile reorder sd median pchisq aggregate aov cor dist var
#' @importFrom purrr map
#' @importFrom readxl read_excel
#' @importFrom reshape2 melt dcast
#' @importFrom rlang .data
#' @importFrom tibble tibble
#' @importFrom ComplexHeatmap HeatmapAnnotation Heatmap
#' @import immunedeconv
#' @importFrom IOBR count2tpm
#' @import apeglm
#' @importFrom scales percent_format
#' @import ggridges
#' @import ggupset
#' @import EPIC
#' @import xCell
#' @importFrom data.table .SD :=
#'
#'
#' @details
#' This file centralizes the import statements for the HTGAnalyzer package.
#' It ensures that functions from the imported packages are available throughout the package.
#'
#' @name HTGAnalyzer-package
NULL
