#' HTG_DEA
#'
#' Perform differential expression analysis using DESeq2, with options for filtering and log fold change (LFC) shrinkage.
#'
#' @param outliers A character vector of sample IDs to be removed as outliers.
#' @param count_data A matrix or data frame of raw count data, with genes as rows and samples as columns.
#' @param col_data A data frame containing the sample metadata. Must include an 'id' column that matches the column names of `count_data`.
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
#'
#' @return A DESeqResults object containing the differential expression results.
#' @export
#'
#' @examples
#' res<- HTG_DEA(outliers =outliers,count_data =  counts,col_data= AnnotData, design_formula = "Ciclina2",
#' heatmap_columns = c("Ciclina2", "Smoker"),
#' contrast = c("Ciclina2", "high", "low"),
#' pattern = "^NC-|^POS-|^GDNA-|^ERCC-", remove_outliers = TRUE, percentage_gene = 0.2,
#' percentage_zero = 0.2, threshold_gene = 200, threshold_subject = 10, pCutoff = 5e-2,
#' apply_filtering = TRUE, apply_lfc_shrinkage = TRUE)
#'
#' @name HTG_DEA

HTG_DEA <- function(outliers, count_data, col_data, design_formula, heatmap_columns, contrast,
                            pattern, remove_outliers = TRUE, percentage_gene, percentage_zero, threshold_gene,
                            threshold_subject, pCutoff, apply_filtering = TRUE, apply_lfc_shrinkage = TRUE) {
  library(DESeq2)
  library(EnhancedVolcano)
  library(ggplot2)

  if (remove_outliers) {
    filtered <- subset(count_data, !grepl(pattern, rownames(count_data)))
    counts_filtered <- filtered[, !colnames(filtered) %in% outliers]
    AnnotData <- col_data[!col_data[["id"]] %in% outliers, ]
  } else {
    counts_filtered <- filtered
    AnnotData <- col_data
  }
  cat("\033[32mDEA desing formula\033[0m\n")
  design_formul <- as.formula(paste("~ ", design_formula))
  print(design_formul)

  ### Variables should not have spaces
  cat("\033[32mSpaces will be changed for _\033[0m\n")
  colnames(AnnotData) <- gsub(" ", "_", colnames(AnnotData))

  ### CHECK IF COLNAMES OF SAMPLE IDs IN col_data ARE THE NAMES OF COL IN COUNT DATA
  col_data <- AnnotData[order(AnnotData$id), ]
  counts_filtered <- counts_filtered[, order(colnames(counts_filtered))]

  # Create DESeqDataSet object
  rownames(col_data) <- col_data$id

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts_filtered, colData = col_data, design = design_formul)
  cat("\033[32mDDS before filtering\033[0m\n")
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
    dds_filtered <- dds[keep_genes, ]

    # Mostrar el resultado
    cat("\033[32mDDS after filtering\033[0m\n")
    print(dds_filtered)
  }

  # Perform DESeq2 analysis
  cat("\033[32mVSD of dds\033[0m\n")
  vsd <- DESeq2::vst(dds, blind = FALSE)
  if (apply_filtering) {
    cat("\033[32mVSD of dds filtered\033[0m\n")
    vsd_filtered <- DESeq2::vst(dds_filtered, blind = FALSE)

  }

  dds <- DESeq2::DESeq(dds)
  if (apply_filtering) {
    cat("\033[32mnormalization of dds filtered\033[0m\n")
    dds_filtered <- DESeq2::DESeq(dds_filtered)
  }

  # Results contrast
  res <- results(dds, contrast = contrast, cooksCutoff = TRUE)
  cat("\033[32mresults of the contrast will be stored in results_HTG_DEA.csv\033[0m\n")
  write.csv(res, "results_HTG_DEA.csv", row.names = TRUE)
  cat("\033[32mTOP10 results will be shown\033[0m\n")
  top_genes <- head(res[order(res$padj), ], 10)
  print(top_genes)


  if (apply_filtering) {
    cat("\033[32mresults of the contrast with dds filtered will be stored in results_HTG_DEA.csv\033[0m\n")
    res_filtered <- results(dds_filtered, contrast = contrast, cooksCutoff = TRUE)
    write.csv(res_filtered, "results_HTG_DEA.csv", row.names = TRUE)
    cat("\033[32mTOP10 results will be shown\033[0m\n")
    top_genes <- head(res_filtered[order(res_filtered$padj), ], 10)
    print(top_genes)

  }

  # Check available coefficients
  print(resultsNames(dds))
  if (apply_filtering) {
    print(resultsNames(dds_filtered))
  }

  # Use the correct coefficient name
  coef_name <- resultsNames(dds)[2]  # Adjust this index based on the printed resultsNames

  if (apply_lfc_shrinkage) {
    resLFC <- lfcShrink(dds, coef = coef_name, type = "apeglm")
    cat("\033[32mTOP10 resLFC will be shown\033[0m\n")
    top_genes <- head(resLFC[order(resLFC$padj), ], 10)
    print(top_genes)
    if (apply_filtering) {
      resLFC_filtered <- lfcShrink(dds_filtered, coef = coef_name, type = "apeglm")
      cat("\033[32mTOP10 resLFC off filtered will be shown\033[0m\n")
      top_genes <- head(resLFC_filtered[order(resLFC_filtered$padj), ], 10)
      print(top_genes)
    }
  }

  # Plots and return results
  if (apply_lfc_shrinkage) {
    if (apply_filtering) {
      print(head(resLFC_filtered[order(resLFC_filtered$padj), ], 10))
      plot<- EnhancedVolcano(resLFC_filtered,
                      lab = rownames(resLFC_filtered),
                      x = 'log2FoldChange',
                      y = 'padj',
                      pCutoff = pCutoff,
                      title = "Volcano Plot with Filtering (LFC Shrinkage)")
      print(plot)
      return(resLFC_filtered)
    } else {
      print(head(resLFC[order(resLFC$padj), ], 10))
      plot<- EnhancedVolcano(resLFC,
                      lab = rownames(resLFC),
                      x = 'log2FoldChange',
                      y = 'padj',
                      pCutoff = pCutoff,
                      title = "Volcano Plot without Filtering (LFC Shrinkage)")
      print(plot)
      return(resLFC)
    }
  } else {
    if (apply_filtering) {
      print(head(res_filtered[order(res_filtered$padj), ], 10))
      plot<- EnhancedVolcano(res_filtered,
                      lab = rownames(res_filtered),
                      x = 'log2FoldChange',
                      y = 'padj',
                      pCutoff = pCutoff,
                      title = "Volcano Plot with Filtering")
      print(plot)
      return(res_filtered)
    } else {
      print(head(res[order(res$padj), ], 10))
      plot<- EnhancedVolcano(res,
                      lab = rownames(res),
                      x = 'log2FoldChange',
                      y = 'padj',
                      pCutoff = pCutoff,
                      title = "Volcano Plot without Filtering")
      print(plot)
      pdf("deseq2_analysis_plot.pdf")
      plot
      dev.off()
      return(res)
    }
  }
}
