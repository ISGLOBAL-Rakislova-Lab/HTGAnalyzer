#' Process DESeq2 analysis
#'
#' This function processes DESeq2 analysis with customizable parameters.
#'
#' @param count_data A matrix or data frame containing count data. It is recommended that the count data frame does not include probes nor outliers.
#' @param col_data A data frame containing sample annotations.
#' @param design_formula The design formula for DESeq2 analysis without the tilde (~).
#' @param threshold_gene Minimum count threshold per gene.
#' @param threshold_subject Minimum count threshold per subject.
#' @param heatmap_columns A character vector specifying the columns to be used for annotations in the heatmap.
#' @param contrast A character vector specifying the contrast for differential expression analysis.
#' @param pCutoff The p-value cutoff for generating the volcano plot.
#'
#' @return Returns an object with the results of contrast demanded and an excel with the results.
#'
#' @export
#'
#' @examples
#' dds_processed <- process_DESeq(counts_filtered, AnnotData, "column_name", threshold_gene = 200, threshold_subject = 10, heatmap_columns = c("GROUP", "RUN"), contrast = c('column_name', 'variable1','variable2'), pCutoff = 5e-2)
#' @name HTG_DESeq
#'
HTG_DESeq <- function(count_data, col_data, design_formula, threshold_gene = 200, threshold_subject = 10, heatmap_columns = c("GROUP", "RUN"), contrast = c('column_name', 'variable1','variable2'), pCutoff = 5e-2) {
  design_formul <- as.formula(paste("~", " ",design_formula))
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(PoiClaClu)
  library(RColorBrewer)
  library(EnhancedVolcano)
  # Create DESeqDataSet object
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design =  design_formul)
  print(" ")
  print(        "BEFORE FILTERING")
  print(dds)
  print(" ")
  # Determine number of genes and subjects
  n_genes <- nrow(count_data)
  n_subj <- nrow(col_data)
  smallest_group_size <- ceiling(n_subj * 0.2)
  keep_rows <- rowSums(DESeq2::counts(dds) >= threshold_gene) >= smallest_group_size
  keep_cols <- colSums(DESeq2::counts(dds) >= threshold_subject) >= ceiling(n_genes * 0.3)
  dds <- dds[keep_rows, keep_cols]
  print(" ")
  print(        "AFTER FILTERING")
  print(dds)
  print(" ")
  # Perform DESeq2 analysis
  vsd <- DESeq2::vst(dds, blind = FALSE)
  print(        "RESULTS OF VST NORMALIZATION")
  print(vsd)
  print(" ")
  dds <- DESeq2::DESeq(dds)
  print("HEATMAP")

  # Generate heatmap
  select <- order(rowMeans(DESeq2::counts(dds, normalized = TRUE)), decreasing = TRUE)[1:500]
  df <- as.data.frame(colData(dds)[, heatmap_columns])
  pheatmap(assay(vsd)[select,], cluster_rows = FALSE, show_rownames = FALSE,
           cluster_cols = TRUE, annotation_col = df, replace = FALSE)

  # Results contrast
  res <- results(dds, contrast = contrast, cooksCutoff = TRUE)
  print("SUMMARY OF RESULT OF THE CONTRAST")
  print(summary(res))
  print(" ")
  print("RESULTS OF THE CONTRAST (TOP 10 p-adj)")
  print(head(res[order(res$padj), ], 10))
  print(" ")
  print("RESULTS WILL BE SAVE IN .CSV IN YOUR CURRENT DIRECTORY")
  write.csv(res, "results.csv", row.names = TRUE)
  GWASTools::qqPlot(res$pvalue)

  # Plots
  dds <- estimateSizeFactors(dds)
  plot(sizeFactors(dds), colSums(counts(dds)))
  text(sizeFactors(dds), colSums(counts(dds)), labels = dds$SampleID, pos = 3, cex = 0.7)
  abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))

  vsd_cor <- cor(assay(vsd))
  rownames(vsd_cor) <- paste(vsd$SampleID)
  colnames(vsd_cor) <- paste(vsd$HTG_RUN)
  pheatmap(vsd_cor)

  poisd <- PoissonDistance(t(counts(dds)))
  samplePoisDistMatrix <- as.matrix(poisd$dd)
  colnames(samplePoisDistMatrix) <- dds$SampleID
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

  pheatmap(samplePoisDistMatrix,
           clustering_distance_rows = poisd$dd,
           clustering_distance_cols = poisd$dd,
           col = colors,
           main = "Poisson distances",
           width = 800,
           height = 600)

  sampleDists <- dist(t(assay(vsd)))

  options(repr.plot.height=7, repr.plot.width=10)
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsd$SampleID, sep = " - ")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors,
           fontsize = 8)

  boxplot(assay(vsd), las = 2, main = "vsd", cex.axis = 0.6)

  # COOK's DISTANCE
  boxplot(log10(assays(dds)[["cooks"]]), range = 0, las = 2, cex.axis = 0.9,main = "COOK'S DISTANCE")

  plotMA(res, main="PLOT MA OF RESULTS")
  res_sorted <- res[order(res$padj), ]
  top_genes_indices <- head(row.names(res_sorted), 10)

  res_sorted <- res[order(res$padj), ]
  top_genes_indices <- head(row.names(res_sorted), 10)

  for (gene_index in top_genes_indices) {
    gen_a2m <- as.data.frame(assay(vsd)[gene_index, ])
    rownames(gen_a2m)
    gen_a2m$status <- col_data[[design_formula]]
    gen_a2m_ordered <- gen_a2m[order(gen_a2m$status), ]

    # Define colors based on status
    levels_design_formula <- unique(col_data[[design_formula]])
    palette <- colorRampPalette(brewer.pal(length(levels_design_formula), "Set1"))
    colors <- palette(length(levels_design_formula))
    group_colors <- colors[as.numeric(factor(col_data[[design_formula]], levels = levels_design_formula))]
    plot(gen_a2m_ordered$`assay(vsd)[gene_index, ]`, xlab = "", ylab = gene_index, col = group_colors,
         xaxt = "n", pch=19)
    axis(1, at = 1:nrow(gen_a2m_ordered), labels = rownames(gen_a2m_ordered), las = 2, cex.axis = 0.6)
    legend("topright", legend = levels_design_formula, fill = colors)

  }

  EnhancedVolcano(res,
                  lab = rownames(res),
                  x = 'log2FoldChange',
                  y = 'padj',
                  pCutoff = pCutoff)
  return(res)
  plotDispEsts(dds, main = "DIPERSION PLOT")
}
