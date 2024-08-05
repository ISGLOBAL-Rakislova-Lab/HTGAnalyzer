#' HTG_DEA
#'
#' Perform differential expression analysis using DESeq2, with options for filtering and log fold change (LFC) shrinkage.
#'
#' @param outliers A character vector of sample IDs to be removed as outliers.
#' @param counts_data A matrix or data frame of raw count data, with genes as rows and samples as columns.
#' @param col_data A data frame containing the sample metadata. Must include an 'id' column that matches the column names of `counts_data`.
#' @param design_formula A character string representing the design formula for DESeq2. Must specify two groups for comparison.
#' @param heatmap_columns A character vector of columns in `col_data` to be used for heatmap annotation.
#' @param contrast A character vector specifying the contrast to be used for differential expression analysis. Should include two levels for comparison.
#' @param pattern A regular expression pattern to filter out unwanted genes from the count data.
#' @param remove_outliers A logical value indicating whether to remove outliers specified in `outliers`.
#' @param percentage_gene A numeric value specifying the minimum percentage of samples a gene must be expressed in to be kept.
#' @param percentage_zero A numeric value specifying the maximum percentage of samples a gene can be zero in to be kept.
#' @param threshold_gene A numeric value specifying the minimum count for a gene to be considered expressed.
#' @param threshold_subject A numeric value specifying the minimum number of subjects a gene must be expressed in to be kept.
#' @param pCutoff A numeric value specifying the p-value cutoff for significance in the volcano plot.
#' @param apply_filtering A logical value indicating whether to apply filtering of genes based on expression thresholds.
#' @param apply_lfc_shrinkage A logical value indicating whether to apply log fold change shrinkage.
#' @param extract_shrinkage A logical value indicating whether to extract shrinkage results.
#'
#' @return A DESeqResults object containing the differential expression results. If `apply_lfc_shrinkage` is TRUE, the function returns the results with log fold change shrinkage applied.
#' @export
#'
#' @examples
#' res1 <- HTG_DEA(outliers, counts_data, AnnotData, design_formula = "Ciclina2", heatmap_columns = c("Ciclina2", "Smoker"),
#'                 contrast = c("Ciclina2", "high", "low"), pattern = "^NC-|^POS-|^GDNA-|^ERCC-", remove_outliers = TRUE,
#'                 percentage_gene = 0.2, percentage_zero = 0.2, threshold_gene = 200, threshold_subject = 10,
#'                 pCutoff = 5e-2, apply_filtering = TRUE, apply_lfc_shrinkage = TRUE, extract_shrinkage = FALSE)
#'
#' @name HTG_DEA

HTG_DEA <- function(outliers, counts_data, col_data, design_formula, heatmap_columns, contrast,
                    pattern, remove_outliers = TRUE, percentage_gene = 0.2, percentage_zero = 0.2, threshold_gene = 200,
                    threshold_subject = 10, pCutoff = 5e-2, apply_filtering = TRUE, apply_lfc_shrinkage = FALSE, extract_shrinkage = FALSE) {
  library(DESeq2)
  library(EnhancedVolcano)
  library(ggplot2)
  library(apeglm)

  if (remove_outliers) {
    filtered <- base::subset(counts_data, !grepl(pattern, rownames(counts_data)))
    counts_filtered <- filtered[, !colnames(filtered) %in% outliers]
    AnnotData <- col_data[!col_data[["id"]] %in% outliers, ]
  } else {
    counts_filtered <- counts_data
    AnnotData <- col_data
  }

  cat("\033[32mDifferential Expression Analysis design formula\033[0m\n")
  design_formul <- as.formula(paste("~ ", design_formula))
  print(design_formul)

  ### Variables should not have spaces
  cat("\033[32mSpaces will be changed for _\033[0m\n")
  colnames(AnnotData) <- gsub(" ", "_", colnames(AnnotData))

  col_data <- AnnotData[order(AnnotData$id), ]
  counts_filtered <- counts_filtered[, order(colnames(counts_filtered))]

  # Create DESeqDataSet object
  rownames(col_data) <- col_data$id

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts_filtered, colData = col_data, design = design_formul)
  cat("\033[32mDESeqDataSet before filtering\033[0m\n")
  print(dds)

  if (apply_filtering) {
    cat("\033[32mApplying filtering\033[0m\n")
    n_genes <- nrow(DESeq2::counts(dds))
    n_subj <- ncol(DESeq2::counts(dds))
    zero_threshold <- ceiling(n_subj * percentage_zero)  # 20% threshold for zeros
    keep_genes <- rowSums(DESeq2::counts(dds) == 0) <= zero_threshold
    # Filter genes that are present in at least 20% of the cases
    smallest_group_size <- ceiling(n_subj * percentage_gene)
    keep_genes <- keep_genes & (rowSums(DESeq2::counts(dds) >= threshold_gene) >= smallest_group_size)
    # Apply gene filters to the DESeq2 object
    dds <- dds[keep_genes, ]

    # Mostrar el resultado
    cat("\033[32mDESeqDataSet after filtering\033[0m\n")
    print(dds)
  }

  # Perform DESeq2 analysis
  cat("\033[32m Variance Stabilized Data of DESeqDataSet\033[0m\n")
  vsd <- DESeq2::vst(dds, blind = FALSE)

  dds <- DESeq2::DESeq(dds)

  # Results contrast
  res <- results(dds, contrast = contrast, cooksCutoff = TRUE)
  cat("\033[32mResults of the contrast will be stored in results_HTG_DEA.csv\033[0m\n")
  write.csv(res, "results_HTG_DEA.csv", row.names = TRUE)
  cat("\033[32mTOP10 results will be shown\033[0m\n")
  top_genes <- head(res[order(res$padj), ], 10)
  print(top_genes)

  # Check available coefficients
  print(resultsNames(dds))

  if (apply_lfc_shrinkage) {
    coef_name <- resultsNames(dds)[2]
    resLFC <- lfcShrink(dds, coef = coef_name, type = "apeglm")
    cat("\033[32mTOP10 resLFC will be shown\033[0m\n")
    top_genes <- head(resLFC[order(resLFC$padj), ], 10)
    print(top_genes)
  }

  # Plots and return result
  if (extract_shrinkage) {
    if (apply_lfc_shrinkage) {
      cat("\033[32mExtracting shrinkage results\033[0m\n")
      print(head(resLFC[order(resLFC$padj), ], 10))
      plot <- EnhancedVolcano(resLFC,
                              lab = rownames(resLFC),
                              x = 'log2FoldChange',
                              y = 'padj',
                              pCutoff = pCutoff,
                              title = "Volcano Plot (LFC Shrinkage)")
      print(plot)
      pdf("deseq2_analysis_plot_resLFC.pdf")
      print(plot)
      dev.off()
      return(resLFC)
    } else {
      cat("\033[32mReturning results without shrinkage\033[0m\n")
      print(head(res[order(res$padj), ], 10))
      plot <- EnhancedVolcano(res,
                              lab = rownames(res),
                              x = 'log2FoldChange',
                              y = 'padj',
                              pCutoff = pCutoff,
                              title = "Volcano Plot")
      print(plot)
      pdf("deseq2_analysis_plot_res.pdf")
      print(plot)
      dev.off()
      return(res)
    }
  } else {
    return(res)
  }
}
