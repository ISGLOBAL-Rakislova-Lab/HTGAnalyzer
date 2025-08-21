#' HTG_DEA
#'
#' Perform differential expression analysis using DESeq2, with options for filtering and log fold change (LFC) shrinkage.
#'
#' @param outliers A character vector of sample IDs to be removed as outliers.
#' @param counts_data A matrix or data frame of raw count data, with genes as rows and samples as columns. Probes will be deleted.
#' @param col_data A data frame containing the sample metadata. Must include an 'id' column that matches the column names of `counts_data`.
#' @param design_formula A character string representing the design formula for DESeq2. Must specify two groups for comparison.
#' @param heatmap_columns A character vector of columns in `col_data` to be used for heatmap annotation. It is recommended to include the column used in `design_formula`.
#' @param contrast A character vector specifying the contrast to be used for differential expression analysis. Should include two levels for comparison.
#' @param pattern A regular expression pattern to filter out unwanted genes from the count data. For HTG, this could be "^NC-|^POS-|^GDNA-|^ERCC-". If NULL, the pattern will not be applied.
#' @param remove_outliers A logical value indicating whether to remove outliers specified in `outliers`.
#' @param percentage_gene A numeric value specifying the minimum percentage of samples a gene must be expressed in to be kept.
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
#' results <- HTG_DEA(
#'   outliers = outliers_tutorial,
#'   counts_data = counts_data_tutorial,
#'   col_data = AnnotData_tutorial,
#'   design_formula = "HPV_status",
#'   heatmap_columns = c("HPV_status", "Recurrence"),
#'    contrast = c("HPV_status", "Associated", "Independent"),
#'   pattern = "^NC-|^POS-|^GDNA-|^ERCC-",
#'   remove_outliers = TRUE,
#'   percentage_gene = 0.2,
#'   threshold_gene = 200,
#'   threshold_subject = 10,
#'   pCutoff = 5e-2,
#'   apply_filtering = TRUE,
#'   apply_lfc_shrinkage = TRUE,
#'   extract_shrinkage = FALSE
#' )
#' @name HTG_DEA
#'
utils::globalVariables(c("SampleID", "status"))
HTG_DEA <- function(outliers = NULL,
                    counts_data,
                    col_data,
                    design_formula,
                    heatmap_columns,
                    contrast,
                    pattern = NULL,
                    remove_outliers = TRUE,
                    percentage_gene = 0.2,
                    threshold_gene = 200,
                    threshold_subject = 10,
                    pCutoff = 5e-2,
                    apply_filtering = TRUE,
                    apply_lfc_shrinkage = FALSE,
                    extract_shrinkage = FALSE) {

  # Filtering counts data based on provided pattern
  if (!is.null(pattern)) {
    cat("\033[33mINITIATING DATA FILTERING...\033[0m\n")
    counts_data <- subset(counts_data, !grepl(pattern, rownames(counts_data)))
  }

  # Removing outliers if specified
  if (remove_outliers && !is.null(outliers) && length(outliers) > 0) {
    cat("\033[33mREMOVING SPECIFIED OUTLIERS...\033[0m\n")
    counts_filtered <- counts_data[, !colnames(counts_data) %in% outliers]
    AnnotData <- col_data[!col_data[["id"]] %in% outliers, ]
  } else {
    counts_filtered <- counts_data
    AnnotData <- col_data
  }

  # Prepare data for DESeq2 analysis
  cat("\033[33mPREPARING DATA FOR DESeq2 ANALYSIS...\033[0m\n")
  design_formul <- as.formula(paste("~ ", design_formula))
  colnames(AnnotData) <- gsub(" ", "_", colnames(AnnotData))
  col_data <- AnnotData[order(AnnotData$id), ]
  counts_filtered <- counts_filtered[, order(colnames(counts_filtered))]
  rownames(col_data) <- col_data$id

  # Convert design variables to factors if needed
  for (var in all.vars(design_formul)) {
    if (var %in% colnames(col_data)) {
      col_data[[var]] <- as.factor(col_data[[var]])
    }
  }

  # Create DESeqDataSet object
  cat("\033[33mCREATING DESeqDataSet OBJECT...\033[0m\n")
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts_filtered, colData = col_data, design = design_formul)

  # Filtering genes if specified
  if (apply_filtering) {
    cat("\033[33mAPPLYING GENE FILTERING BASED ON THRESHOLDS...\033[0m\n")
    n_genes <- nrow(DESeq2::counts(dds))
    n_subj <- ncol(DESeq2::counts(dds))
    zero_threshold <- ceiling(n_subj * 0.2)
    keep_genes <- rowSums(DESeq2::counts(dds) == 0) <= zero_threshold
    smallest_group_size <- ceiling(n_subj * percentage_gene)
    keep_genes <- keep_genes & (rowSums(DESeq2::counts(dds) >= threshold_gene) >= smallest_group_size)
    dds <- dds[keep_genes, ]
  }

  # Perform DESeq2 analysis
  cat("\033[33mPERFORMING DESeq2 ANALYSIS...\033[0m\n")
  vsd <- DESeq2::vst(dds, blind = FALSE)
  dds <- DESeq2::DESeq(dds)
  dds
  res <- DESeq2::results(dds, contrast = contrast, cooksCutoff = TRUE)

  # Save results
  write.csv(res, "results_HTG_DEA.csv", row.names = TRUE)
  cat("\033[32mResults saved as 'results_HTG_DEA.csv'\033[0m\n")

  # Print top 10 results
  cat("\033[32mDisplaying top 10 differential expression results...\033[0m\n")
  top_genes <- head(res[order(res$padj), ], 10)
  print(top_genes)

  # Apply LFC shrinkage if specified
  if (apply_lfc_shrinkage) {
    cat("\033[33mAPPLYING LOG-FOLD CHANGE SHRINKAGE...\033[0m\n")
    coef_name <- DESeq2::resultsNames(dds)[2]
    resLFC <- DESeq2::lfcShrink(dds, coef = coef_name, type = "apeglm")
    top_genes <- head(resLFC[order(resLFC$padj), ], 10)
    print(top_genes)
  }

  # Save plots to PDF
  pdf("DEA_plots.pdf", width = 10, height = 8)
  cat("\033[32mSaving plots as 'DEA_plots.pdf'\033[0m\n")

  # Volcano plot
  suppressWarnings(invisible(print(EnhancedVolcano::EnhancedVolcano(res,
                                                                    lab = rownames(res),
                                                                    x = 'log2FoldChange',
                                                                    y = 'padj',
                                                                    pCutoff = pCutoff))))

  # Heatmap of sample-to-sample distances
  vsd_cor <- stats::cor(SummarizedExperiment::assay(vsd))
  rownames(vsd_cor) <- paste(vsd$SampleID)
  colnames(vsd_cor) <- paste(vsd$HTG_RUN)
  suppressWarnings(invisible(print(pheatmap::pheatmap(vsd_cor, main = "Sample-to-Sample Correlation"))))
  suppressWarnings(invisible(print(pheatmap::pheatmap(vsd_cor, main = "Sample-to-Sample Correlation"))))


  invisible(capture.output(
    boxplot(SummarizedExperiment::assay(vsd), las = 2, main = "VST-transformed Data", outline = FALSE, col = "#4793AF")
  ))

  # MA-plot
  suppressWarnings(invisible(print(DESeq2::plotMA(res, main = "MA Plot of Results"))))

  # Individual gene plots for top 10 genes
  top_genes_indices <- rownames(top_genes)
  for (gene_index in top_genes_indices) {
    cat("\033[33mPlotting gene expression for gene:\033[0m", gene_index, "\n")
    gen_a2m <- as.data.frame(SummarizedExperiment::assay(vsd)[gene_index, ])
    colnames(gen_a2m) <- "expression"
    gen_a2m$SampleID <- rownames(gen_a2m)
    gen_a2m$status <- col_data[[design_formula]]

    # Create ggplot for the gene
    levels_design_formula <- unique(col_data[[design_formula]])

    # Define colors based on number of levels
    color_palette <- if (length(levels_design_formula) == 2) {
      c("red", "#4793AF")
    } else {
      colorRampPalette(brewer.pal(min(length(levels_design_formula), 9), "Set1"))(length(levels_design_formula))
    }

    plot_gene <- ggplot2::ggplot(gen_a2m, ggplot2::aes(x = factor(SampleID, levels = SampleID[order(status)]), y = expression, color = status)) +
      ggplot2::geom_point(size = 2) +
      ggplot2::labs(x = "", y = gene_index, title = paste("Gene:", gene_index)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1), legend.title = ggplot2::element_blank()) +
      ggplot2::scale_color_manual(values = color_palette)

    suppressWarnings(invisible(print(plot_gene)))
  }

  dev.off()

  if (extract_shrinkage) {
    if (apply_lfc_shrinkage) {
      return(resLFC)
    } else {
      return(res)
    }
  } else {
    return(res)
  }
}


