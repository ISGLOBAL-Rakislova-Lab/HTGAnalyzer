#' HTG_analysis: Perform DESeq2 Analysis, GSEA, TME, and Survival Analysis
#'
#' @description This function conducts a comprehensive analysis pipeline including DESeq2 differential expression analysis (DEA), Gene Set Enrichment Analysis (GSEA), TME analysis, and survival analysis. The pipeline supports optional steps for generating volcano plots and heatmaps. The function is suitable for both HTG and RNA-seq data.
#'
#' @param outliers A character vector specifying the IDs of outlier samples to be removed. If you have RNA-seq data or don't have outliers, this parameter can be set to NULL. Note: For RNA-seq data, `remove_outliers` should be set to FALSE.
#' @param pattern A regular expression pattern to identify control probes in the count data. For HTG, this could be "^NC-|^POS-|^GDNA-|^ERCC-". If NULL, the pattern will not be applied.
#' @param counts_data A matrix or data frame containing count data. It is recommended that the count data does not include control probes or outliers if not needed.
#' @param col_data A data frame containing sample annotations. For survival analysis, it must include variables `time` and `variable_01`.
#' @param design_formula The design formula for DESeq2 analysis, specified as a string without the tilde (~).
#' @param threshold_gene Minimum count threshold per gene. Default is 200.
#' @param threshold_subject Minimum count threshold per subject. Default is 10.
#' @param heatmap_columns A character vector specifying the columns to be used for annotations in the heatmap. It is recommended to include the column used in `design_formula`. Default is NULL.
#' @param contrast A character vector specifying the contrast for differential expression analysis. Should include the column name and two levels for comparison.
#' @param pCutoff The p-value cutoff for generating the volcano plot. Default is 0.05.
#' @param remove_outliers A logical value indicating whether to remove outliers. Default is TRUE. Note: For RNA-seq data, this should be set to FALSE.
#' @param GSEA A logical value indicating whether to perform Gene Set Enrichment Analysis. Default is FALSE.
#' @param generate_heatmap A logical value indicating whether to generate a heatmap. Default is TRUE.
#' @param TME A logical value indicating whether to perform Tumor Microenvironment (TME) analysis. Default is TRUE. Parameter variable_01 will be needed.
#' @param survival_analysis A logical value indicating whether to perform survival analysis. Default is FALSE. Parameter time will be needed.
#' @param percentage_gene A numeric value between 0 and 1 indicating the minimum fraction of samples in which a gene must be expressed to be retained. Default is 0.2.
#' @param genes_to_use A character vector specifying top genes for analysis. Default is c("CCND1", "MMP10", "CTTN").
#' @param DEA A logical value indicating whether to perform DESeq2 analysis with filtering and without LFC shrinkage. Default is TRUE.
#' @param variable_01 A character string specifying the event/censoring variable used in the survival analysis. Set NULL if you don't have it and turn FALSE TME
#' @param time A character string specifying the time variable used in the survival analysis. Set NULL if you don't have it and turn FALSE survival_analysis
#'
#' @return Returns an object with the results of the specified contrast and saves an Excel file with the results, along with PDF files of the generated plots.
#' @export
#'
#' @import org.Hs.eg.db
#' @import immunedeconv
#' @importFrom maxstat maxstat.test
#' @importFrom survival Surv survfit survdiff
#' @importFrom stats as.formula t.test lm prcomp na.omit quantile reorder sd median pchisq
#' @importFrom utils write.csv head capture.output str
#' @importFrom grDevices colorRampPalette pdf dev.off
#' @importFrom graphics boxplot abline axis legend par text
#' @importFrom dplyr filter
#' @importFrom DESeq2 DESeqDataSetFromMatrix counts DESeq vst results resultsNames lfcShrink plotMA
#' @importFrom SummarizedExperiment assay
#' @importFrom pheatmap pheatmap
#' @importFrom EnhancedVolcano EnhancedVolcano
#' @importFrom ggplot2 ggplot aes geom_point labs theme element_text element_blank scale_color_manual
#' @importFrom grDevices pdf dev.off colorRampPalette
#' @importFrom graphics boxplot
#' @importFrom RColorBrewer brewer.pal
#' @importFrom utils write.csv capture.output head
#' @importFrom stats aggregate aov cor dist var
#'
#' @examples
#' ALL_done <- HTG_analysis(
#'   outliers = outliers_tutorial,
#'   pattern = "^NC-|^POS-|^GDNA-|^ERCC-",
#'   counts_data = counts_data_tutorial,
#'   col_data = AnnotData_tutorial,
#'   design_formula = "HPV_status",
#'   percentage_gene = 0.2,
#'   threshold_gene = 200,
#'   threshold_subject = 10,
#'   genes_to_use = c("CCND1", "MMP10", "CTTN"),
#'   heatmap_columns = c("HPV_status", "Cyclin_D1"),
#'   contrast = c("HPV_status", "Associated", "Independent"),
#'   pCutoff = 5e-2,
#'   variable_01 = "Recurrence_01",
#'   time = "time_to_recurrence",
#'   DEA = TRUE,
#'   remove_outliers = TRUE,
#'   GSEA = TRUE,
#'   generate_heatmap = TRUE,
#'   TME = TRUE,
#'   survival_analysis = TRUE
#' )
#' @name HTG_analysis
utils::globalVariables(c("PC1", "PC2", "Tag", "Sample", "padj", "Description", "Count", ".data", "Condition_Group", "mean_value", "Cell_Type", "Value", "Average", "Fraction", "shapiro_test", "id"))
HTG_analysis <- function(outliers = NULL,
                         pattern = NULL,
                         counts_data,
                         col_data,
                         design_formula = NULL ,
                         percentage_gene = 0.2,
                         threshold_gene = 200,
                         threshold_subject = 10,
                         genes_to_use = c("CCND1", "MMP10", "CTTN"),
                         heatmap_columns = NULL,
                         pvalueCutoff = 0.05,
                         contrast = NULL,
                         pCutoff = 5e-2,
                         variable_01 = NULL,
                         time = NULL,
                         DEA = TRUE,
                         remove_outliers = TRUE,
                         GSEA = FALSE,
                         generate_heatmap = TRUE,
                         TME = TRUE,
                         survival_analysis = FALSE) {
  if (!requireNamespace("IOBR", quietly = TRUE)) {
    stop("Package 'IOBR' is required but not installed.")
  }
  if (!requireNamespace("immunedeconv", quietly = TRUE)) {
    stop("Package 'immunedeconv' is required but not installed.")
  }

  if (!requireNamespace("EPIC", quietly = TRUE)) {
    stop("Package 'EPIC' is required but not installed.")
  }

  if (!requireNamespace("xCell", quietly = TRUE)) {
    stop("Package 'xCell' is required but not installed.")
  }
  utils::data("xCell.data", package = "xCell", envir = environment())

  if (!is.null(pattern)) {
    cat("\033[33mFILTERING THE COUNT DATA. DELETING THE PROVES...\033[0m\n")
    counts_data <- subset(counts_data, !grepl(pattern, rownames(counts_data)))
  }

  # Removing outliers if specified
  if (remove_outliers && !is.null(outliers) && length(outliers) > 0) {
    cat("\033[33mREMOVING OUTLIERS...\033[0m\n")
    counts_filtered <- counts_data[, !colnames(counts_data) %in% outliers]
    AnnotData <- col_data[!col_data[["id"]] %in% outliers, ]
  } else {
    counts_filtered <- counts_data
    AnnotData <- col_data
  }


  if (!is.null(DEA) && DEA) {
    cat("\033[33mSTARTING THE DIFERENTIAN EXPRESSION ANALYSIS.\033[0m\n")
    if (is.null(contrast)) {
      stop("Contrast is required for DESeq2 analysis. Remember structure: contrast = c('column_name', 'variable1','variable2')")
    }
    if (is.null(design_formula)) {
      stop("Design formula is required for DESeq2 analysis.")
    }
    print(design_formula)
    design_formul <- as.formula(paste("~", " ", design_formula))

  ### Variables should not have spaces
  colnames(AnnotData) <- gsub(" ", "_", colnames(AnnotData))

  ### CHECK IF COLNAMES OF SAMPLE IDs IN col_data ARE THE NAMES OF COL IN COUNT DATA
  col_data <- AnnotData[order(AnnotData$id), ]
  counts_filtered <- counts_filtered[, order(colnames(counts_filtered))]
  if (!identical(colnames(counts_filtered), col_data$id)) {
    stop("Column names of counts_filtered and IDs in col_data do not match.")
  }

  # Create DESeqDataSet object
  rownames(col_data)<- col_data$id  #important that the columns is called id
  col_data[[contrast[1]]] <- as.factor(col_data[[contrast[1]]])
  counts_filtered<- as.data.frame(counts_filtered)

cat("\033[33m[WARNING] DESeq2 expects integer counts. Are you sure it's raw counts?\033[0m\n")
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts_filtered, colData = col_data, design = design_formul)
cat("\033[33m[WARNING] Yes! It was raw counts\033[0m\n")
    
  cat("\033[32m\033[0m\n")
  cat("\033[32mBEFORE FILTERING\033[0m\n")
  print(dds)
  cat("\033[32m\033[0m\n")
  n_genes <- nrow(DESeq2::counts(dds))
  n_subj <- ncol(DESeq2::counts(dds))
  zero_threshold <- ceiling(n_subj * 0.2)
  keep_genes <- rowSums(DESeq2::counts(dds) == 0) <= zero_threshold
  smallest_group_size <- ceiling(n_subj * percentage_gene)
  keep_genes <- keep_genes & (rowSums(DESeq2::counts(dds) >= threshold_gene) >= smallest_group_size)
  # Apply gene filters to the DESeq2 object
  dds <- dds[keep_genes, ]

  # Mostrar el resultado
  cat("\033[32m \033[0m\n")
  cat("\033[32m\033[32mAFTER FILTERING\033[0m\n")
  print(dds)
  cat("\033[32m\033[32m \033[0m\n")

  # Perform DESeq2 analysis
  vsd <- DESeq2::vst(dds, blind = FALSE)
  cat("\033[32m\033[32mRESULTS OF VST NORMALIZATION\033[0m\n")
  print(vsd)
  cat("\033[32m\033[32m \033[0m\n")
  dds <- DESeq2::DESeq(dds)
  cat("\033[32m\033[32mHEATMAP\033[0m\n")

  # Generate heatmap if generate_heatmap is TRUE
  if (generate_heatmap) {
    cat("\033[33mGENERATING HEATMAP\033[0m\n")
    if (is.null(heatmap_columns)) {
      stop("heatmap_columns are required for generating the heatmap.")
    }
    pdf("HEATMAP_analysis.pdf", width = 12, height = 10)

    # Seleccionar las filas con los 500 genes más expresados
    selecto <- order(rowMeans(DESeq2::counts(dds, normalized = TRUE)), decreasing = TRUE)[1:500]

    # Crear un data frame para las anotaciones de columnas
    df <- as.data.frame(SummarizedExperiment::colData(dds)[, heatmap_columns])

    # Generar el heatmap con nombres de columnas (muestras) y nombres de filas si se desea
    a<-pheatmap::pheatmap(SummarizedExperiment::assay(vsd)[selecto,],
                       cluster_rows = FALSE,
                       show_rownames = FALSE,  # set TRUE to show the gene names
                       show_colnames = TRUE,    # Sample names.
                       cluster_cols = TRUE,
                       annotation_col = df)
    print(a)
    dev.off()
  }

  # Results contrast
  cat("\033[33mGENERATING CONTRAST RESULTS\033[0m\n")
  res <- DESeq2::results(dds, contrast = contrast, cooksCutoff = TRUE)
  cat("\033[32mSUMMARY OF RESULT OF THE CONTRAST\033[0m\n")
  print(summary(res))
  cat("\033[32m \033[0m\n")
  cat("\033[32mRESULTS OF THE CONTRAST (TOP 10 p-adj)\033[0m\n")
  print(head(res[order(res$padj), ], 10))
  cat("\033[32m \033[0m\n")
  cat("\033[32mRESULTS of DEA WILL BE SAVED IN .CSV IN YOUR CURRENT DIRECTORY\033[0m\n")
  write.csv(res, "results_DEA.csv", row.names = TRUE)
  pdf("DEA_plots_HTG_analysis.pdf", width = 10, height = 8)

  #  GWASTools::qqPlot(res$pvalue)
  z_scores <- stats::qnorm(res$pvalue)
  stats::qqnorm(z_scores, main = "QQ Plot de p-values")
  stats::qqline(z_scores, col = "red", lwd = 2)

  # Plots
  # dds <- DESeq2::estimateSizeFactors(dds)
  # plot(DESeq2::sizeFactors(dds), colSums(DESeq2::counts(dds)),
  #      xlab = "Size Factors", ylab = "Column Sums of Counts",
  #      main = "Size Factors vs. Column Sums")
  # text(DESeq2::sizeFactors(dds), colSums(DESeq2::counts(dds)), labels = colnames(DESeq2::counts(dds)), pos = 3, cex = 0.7)
  # abline(lm(colSums(DESeq2::counts(dds)) ~ DESeq2::sizeFactors(dds) + 0))
  #
  # vsd_cor <- stats::cor(SummarizedExperiment::assay(vsd))
  # rownames(vsd_cor) <- paste(vsd$SampleID)
  # colnames(vsd_cor) <- paste(vsd$HTG_RUN)
  # # pheatmap::pheatmap(vsd_cor)
  # pheatmap::pheatmap(vsd_cor,
  #                    main = "Sample-to-Sample Correlation Heatmap",
  #                    display_numbers = TRUE)
  dds <- DESeq2::estimateSizeFactors(dds)
  plot(DESeq2::sizeFactors(dds), colSums(DESeq2::counts(dds)),
       xlab = "Size Factors", ylab = "Column Sums of Counts",
       main = "Size Factors vs. Column Sums")
  text(DESeq2::sizeFactors(dds), colSums(DESeq2::counts(dds)), labels = colnames(DESeq2::counts(dds)), pos = 3, cex = 0.7)
  abline(lm(colSums(DESeq2::counts(dds)) ~ DESeq2::sizeFactors(dds) + 0))
  vsd_cor <- stats::cor(SummarizedExperiment::assay(vsd))
  sample_ids <- vsd$id
  annotation_col <- data.frame(SampleID = sample_ids)
  rownames(annotation_col) <- colnames(vsd_cor)  # Asegúrate de que los nombres coincidan
  pheatmap::pheatmap(vsd_cor,
                     main = "Sample-to-Sample Correlation Heatmap",
                     display_numbers = FALSE,
                     annotation_col = annotation_col,
                     annotation_legend = FALSE)

  # pois_distance <- PoiClaClu::PoissonDistance(t(DESeq2::counts(dds, normalized = TRUE)))
  # samplePoisDistMatrix <- as.matrix(pois_distance$dd)
  # sample_ids <- SummarizedExperiment::colData(dds)$SampleID
  # colnames(samplePoisDistMatrix) <- sample_ids
  # rownames(samplePoisDistMatrix) <- sample_ids
  #  colors <- grDevices::colorRampPalette(c("white", "#4793AF", "#013649"))(255)
  #
  # cat("\033[33mGENERATING POISSON DISTANCES PLOT\033[0m\n")
  #
  # pheatmap::pheatmap(samplePoisDistMatrix,
  #          clustering_distance_rows = pois_distance$dd,
  #          clustering_distance_cols = pois_distance$dd,
  #          col = colors,
  #          main = "Poisson Distances Heatmap",
  #          width = 800,
  #          height = 600,
  #          display_numbers = TRUE)
  #
  #


  # Calcular la distancia de Poisson
  pois_distance <- PoiClaClu::PoissonDistance(t(DESeq2::counts(dds, normalized = TRUE)))
  samplePoisDistMatrix <- as.matrix(pois_distance$dd)
  sample_ids <- dds$id
  colnames(samplePoisDistMatrix) <- sample_ids
  rownames(samplePoisDistMatrix) <- sample_ids
  colors <- grDevices::colorRampPalette(c("white", "#4793AF", "#013649"))(255)

  cat("\033[33mGENERATING POISSON DISTANCES PLOT\033[0m\n")

  pheatmap::pheatmap(samplePoisDistMatrix,
                     clustering_distance_rows = pois_distance$dd,
                     clustering_distance_cols = pois_distance$dd,
                     color = colors,
                     main = "Poisson Distances Heatmap",
                     width = 800,
                     height = 600,
                     display_numbers = FALSE,
                     number_format = "%.2f",
                     angle_col = 45)  # Ajuste opcional para rotar etiquetas de columnas



  sampleDists <- stats::dist(t(SummarizedExperiment::assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsd$id, sep = " - ")
  colnames(sampleDistMatrix) <- rownames(sampleDistMatrix)

  colors <- grDevices::colorRampPalette(c("white", "#4793AF"))(255)
  cat("\033[33mGENERATING VSD DISTANCE HEATMAP\033[0m\n")

  pheatmap::pheatmap(
    sampleDistMatrix,
    clustering_distance_rows = sampleDists,
    clustering_distance_cols = sampleDists,
    col = colors,
    fontsize = 8,
    main = "Sample Distance Heatmap",
    display_numbers = FALSE)

  cat("\033[33mGENERATING VSD BOXPLOT\033[0m\n")
  graphics::boxplot(SummarizedExperiment::assay(vsd), las = 2, main = "vsd", cex.axis = 0.6)


  res_sorted <- res[order(res$padj), ]
  top_genes_indices <- head(row.names(res_sorted), 10)

  for (gene_index in top_genes_indices) {
    gen_a2m <- as.data.frame(SummarizedExperiment::assay(vsd)[gene_index, ])
    rownames(gen_a2m)
    gen_a2m$status <- col_data[[design_formula]]
    gen_a2m_ordered <- gen_a2m[order(gen_a2m$status), ]


    # Define colors based on status
     levels_design_formula <- unique(col_data[[design_formula]])

    num_colors <- length(levels_design_formula)
    palette <- colorRampPalette(c("red", "orange", "yellow", "green", "purple", "#4793AF"))(num_colors)
    colors <- palette

    group_colors <- colors[as.numeric(factor(col_data[[design_formula]], levels = levels_design_formula))]

    # Graficar
    plot(gen_a2m_ordered$`SummarizedExperiment::assay(vsd)[gene_index, ]`,
      xlab = "",
      ylab = gene_index,
      col = group_colors,
      xaxt = "n",
      pch = 19)
    axis(1, at = 1:nrow(gen_a2m_ordered), labels = rownames(gen_a2m_ordered), las = 2, cex.axis = 0.6)
        legend("topright", legend = levels_design_formula, fill = colors)
  }

  cat("\033[33mGENERATING PCA PLOT\033[0m\n")

  normalized_counts <- DESeq2::counts(dds, normalized = TRUE)
  pca_result <- prcomp(t(normalized_counts), scale. = TRUE)

  pca_data <- as.data.frame(pca_result$x)
  pca_data$Sample <- SummarizedExperiment::colData(dds)$id
  pca_data$Condition_Group <- SummarizedExperiment::colData(dds)[[design_formula]]

  color_palette <- c("#4793AF", "#E57373")
  names(color_palette) <- unique(pca_data$Condition_Group)

  pca_plot <- ggplot2::ggplot(pca_data, ggplot2::aes(x = PC1, y = PC2, color = Condition_Group)) +
    ggplot2::geom_point(size = 3) +
    ggrepel::geom_text_repel(ggplot2::aes(label = Sample), size = 8, max.overlaps = 15) +  # Aumenta el tamaño de las etiquetas de muestra
    ggplot2::labs(title = "PCA of Normalized Counts",
                  x = "Principal Component 1",
                  y = "Principal Component 2") +
    ggplot2::theme_minimal() +
    ggplot2::scale_color_manual(values = color_palette) +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(size = 25),          # Título del gráfico
      axis.title.x = ggplot2::element_text(size = 25),        # Título del eje X
      axis.title.y = ggplot2::element_text(size = 25),        # Título del eje Y
      axis.text = ggplot2::element_text(size = 20),           # Texto de los ejes
      legend.title = ggplot2::element_text(size = 25),        # Título de la leyenda
      legend.text = ggplot2::element_text(size = 12)          # Texto de los elementos de la leyenda
    )

  print(pca_plot)






    cat("\033[33mGENERATING VOLCANO PLOT\033[0m\n")
    suppressWarnings(invisible(print(EnhancedVolcano::EnhancedVolcano(res,
                                                                      lab = rownames(res),
                                                                      x = 'log2FoldChange',
                                                                      y = 'padj',
                                                                      pCutoff = pCutoff))))

    cat("\033[33mGENERATING MA PLOT\033[0m\n")

    DESeq2::plotMA(res, main = "PLOT MA OF RESULTS")

    cat("\033[33mGENERATING DISPERSION PLOT\033[0m\n")

  DESeq2::plotDispEsts(dds, main = "DIPERSION PLOT")
  dev.off()
  } else {cat("\033[32mSkipping Diferential expresion analysis\033[0m\n")}

  if (GSEA) {
    if (!exists("res")) {
      cat("\033[31mWarning: The variable 'res' does not exist. GSEA analysis cannot proceed without it.\033[0m\n")
    }
    else if (!inherits(res, "DESeqResults")) {
      cat("\033[31mWarning: The variable 'res' exists but is not of class 'DESeqResults'. GSEA analysis cannot proceed.\033[0m\n")
    }
    else {
    cat("\033[33mSTARTING GSEA\033[0m\n")
  # Prepare gene list for gseGO
  cat("\033[32mPreparing gene list for GSEA\033[0m\n")
  original_gene_list <- res$log2FoldChange
  names(original_gene_list) <- rownames(res)
  gene_list <- na.omit(original_gene_list)
  gene_list <- sort(gene_list, decreasing = TRUE)

  # gseGO  Analysis
  cat("\033[32mBEFORE STARTING...\033[0m\n")
  cat("\033[34mEnsure your gene_list is named, numeric, and sorted decreasingly (important for GSEA).\033[0m\n")
  cat("\033[34mIf your gene_list contains mostly zeros or non-positive values, gseGO might not work properly.\033[0m\n")
  cat("\033[34mChecking gene_list structure...\033[0m\n")


  # Checks útils
  print(head(gene_list, 5))
  cat("Gene list length: ", length(gene_list), "\n")
  cat("Any negative values? ", any(gene_list < 0), "\n")
  cat("Any NA values? ", any(is.na(gene_list)), "\n")
  cat("Length gene_list:", length(gene_list),"\n")

  # Comprovació de noms
  if (is.null(names(gene_list))) {
    stop("Gene_list must be a named numeric vector. Please provide gene symbols as names.")
  }

  cat("\033[32mPerforming gseGO Analysis\033[0m\n")
      gse2 <- clusterProfiler::gseGO(geneList = gene_list, ont = "BP", keyType = "SYMBOL", nPermSimple = 500000,
                                     minGSSize = 3, maxGSSize = 800,  pvalueCutoff = pvalueCutoff, verbose = TRUE, eps = 0,
                                     OrgDb = org.Hs.eg.db::org.Hs.eg.db, pAdjustMethod = "bonferroni")

      # Dotplot for gseGO
      cat("\033[32mCreating Dotplot for gseGO \033[0m\n")
      dotplot1 <- clusterProfiler::dotplot(gse2, showCategory = 10, split = ".sign", font.size = 9, label_format = 40,
                                           title = "gseGO Enrichment Results: Pathways", color = "p.adjust", size = "Count")

      dotplot2 <- clusterProfiler::dotplot(gse2, showCategory = 10, split = ".sign", font.size = 9, label_format = 40,
                                           title = "gseGO Enrichment Results: Pathways", color = "p.adjust", size = "Count") + ggplot2::facet_grid(.~.sign)

      # Emaplot for gseGO
      cat("\033[32mCreating Emaplot for gseGO \033[0m\n")
      x2 <- enrichplot::pairwise_termsim(gse2)
      #emapplot1 <- enrichplot::emapplot(x2, max.overlaps = 70, min.segment.length = 0.3, point_size = 0.3, font.size = 5) +   ggplot2::ggtitle("Enrichment Map gseGO ")
      emapplot1 <- enrichplot::emapplot(x2) + ggplot2::ggtitle("Enrichment Map gseGO ")


      # Ridgeplot for gseGO
      cat("\033[32mCreating Ridgeplot for gseGO \033[0m\n")
      ridgeplot1 <- enrichplot::ridgeplot(gse2)  +  ggplot2::labs(x = "gseGO enrichment distribution", font.size = 7) +  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9))

      # Heatplot for gseGO
      cat("\033[32mCreating Heatplot for gseGO \033[0m\n")
      heatplot1 <- enrichplot::heatplot(gse2, showCategory = 10) + ggplot2::ggtitle("gseGO Heatplot")

      # Treeplot for gseGO
      cat("\033[32mCreating Treeplot for gseGO \033[0m\n")
      treeplot1 <- suppressWarnings(enrichplot::treeplot(x2) + ggplot2::ggtitle("gseGO Treeplot"))

      # Create gseaplot2 plots with titles
      a <- enrichplot::gseaplot2(gse2, geneSetID = 1, title = paste("GSEA Plot:", gse2$Description[1]))
      b <- enrichplot::gseaplot2(gse2, geneSetID = 1:5, pvalue_table = TRUE, title = "GSEA: Top 5 Gene Sets")

      # KEGG Analysis
      cat("\033[32mPerforming KEGG Analysis\033[0m\n")
      ids <- clusterProfiler::bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb =  org.Hs.eg.db::org.Hs.eg.db)
      dedup_ids <- ids[!duplicated(ids[c("SYMBOL")]), ]
      df2 <- res[rownames(res) %in% dedup_ids$SYMBOL, ]
      df2$Y <- dedup_ids$ENTREZID
      kegg_gene_list <- df2$log2FoldChange
      names(kegg_gene_list) <- df2$Y
      kegg_gene_list <- na.omit(kegg_gene_list)
      kegg_gene_list <- sort(kegg_gene_list, decreasing = TRUE)

      kk2 <- clusterProfiler::gseKEGG(geneList = kegg_gene_list, organism = "hsa", minGSSize = 3, maxGSSize = 800,
                                       pvalueCutoff = pvalueCutoff, pAdjustMethod = "none", keyType = "ncbi-geneid", nPermSimple = 100000)

      # Dotplot for KEGG
      cat("\033[32mCreating Dotplot for KEGG\033[0m\n")
      dotplot3 <- clusterProfiler::dotplot(kk2, showCategory = 10, title = "Enriched Pathways for KEGG", split = ".sign", font.size = 9) + ggplot2::facet_grid(.~.sign)

      # Emaplot for KEGG
      cat("\033[32mCreating Emaplot for KEGG\033[0m\n")
      x3 <- enrichplot::pairwise_termsim(kk2)
      #emapplot2 <- enrichplot::emapplot(x3, font.size = 8 +  ggplot2::ggtitle("KEGG Enrichment Map"))
      emapplot2 <- enrichplot::emapplot(x3) + ggplot2::ggtitle("KEGG Enrichment Map")


      # Ridgeplot for KEGG
      cat("\033[32mCreating Ridgeplot for KEGG\033[0m\n")
      ridgeplot2 <- enrichplot::ridgeplot(kk2) +  ggplot2::labs(x = "KEGG enrichment distribution", font.size = 6) +  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9))

      # Heatplot for KEGG
      cat("\033[32mCreating Heatplot for KEGG \033[0m\n")
      heatplot2 <- enrichplot::heatplot(kk2, showCategory = 10) + ggplot2::ggtitle("KEGG Heatplot")

      # Treeplot for KEGG
      cat("\033[32mCreating Treeplot for KEGG \033[0m\n")
      treeplot2 <- suppressWarnings(enrichplot::treeplot(x3) + ggplot2::ggtitle("KEGG Treeplot"))

      upset_plot <- enrichplot::upsetplot(kk2) + ggplot2::labs(title = "Up set plot for KEGG")

      # enrichGO Analysis
      cat("\033[32mPerforming GO Enrichment Analysis\033[0m\n")
      sig_genes_df <- subset(res, padj < 0.05)

      if (nrow(sig_genes_df) > 0) {
        # Prepare the gene list for enrichment analysis
        genes <- sig_genes_df$log2FoldChange
        names(genes) <- rownames(sig_genes_df)

        # Perform GO enrichment analysis
        cat("\033[32mPerforming GO Enrichment Analysis\033[0m\n")
        go_enrich <- clusterProfiler::enrichGO(
          gene = names(genes),
          universe = names(gene_list),
          OrgDb =  org.Hs.eg.db::org.Hs.eg.db,
          keyType = 'SYMBOL',
          readable = TRUE,
          ont = "BP",
           pvalueCutoff = pvalueCutoff,
          qvalueCutoff = 0.10
        )

        # Create a bar plot for significant GO terms
        cat("\033[32mCreating Bar Plot for GO Enrichment\033[0m\n")
        go_results <- go_enrich@result
        significant_terms <- go_results[go_results$qvalue < 0.05, ]
        significant_terms <- significant_terms[order(significant_terms$qvalue), ]

        # Check if there are significant terms to plot
        if (nrow(significant_terms) > 0) {
          bar_plot <- ggplot2::ggplot(significant_terms, ggplot2::aes(x = reorder(Description, -Count), y = Count)) +
            ggplot2::geom_bar(stat = "identity", fill = "skyblue") +
            ggplot2::labs(title = "Significant GO Terms", x = "GO Term", y = "Count") +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))

        } else {
          cat("\033[31mNo significant GO terms found to plot.\033[0m\n")
        }

        # Save the GO enrichment results to a CSV file
        write.csv(go_results, "go_results.csv")
        cat("\033[32mGO enrichment results saved to 'go_results.csv'\033[0m\n")
      } else {
        cat("\033[31mNo significant genes found for GO enrichment analysis.\033[0m\n")
      }

      pdf("GSEA_analysis_plots1_of_2.pdf", width = 11, height = 14)
      print(dotplot1)
      print(dotplot2)
      print(emapplot1)
      print(ridgeplot1)
      print(treeplot1)
      print(a)
      print(b)
      print(dotplot3)
      print(emapplot2)
      print(ridgeplot2)
      print(treeplot2)
      print(upset_plot)
      if (nrow(sig_genes_df) > 0) {
        genes <- sig_genes_df$log2FoldChange
        names(genes) <- rownames(sig_genes_df)
        go_enrich <- clusterProfiler::enrichGO(gene = names(genes), universe = names(gene_list), OrgDb =  org.Hs.eg.db::org.Hs.eg.db,
                                               keyType = 'SYMBOL', readable = TRUE, ont = "BP",
                                                pvalueCutoff = pvalueCutoff, qvalueCutoff = 0.10)
        go_results <- go_enrich@result
        significant_terms <- go_results[go_results$qvalue < 0.05, ]
        significant_terms <- significant_terms[order(significant_terms$qvalue), ]
        bar_plot <- ggplot2::ggplot(significant_terms, ggplot2::aes(x = reorder(Description, -Count), y = Count)) +
          ggplot2::geom_bar(stat = "identity", fill = "skyblue") +
          ggplot2::labs(title = "Significant GO Terms", x = "GO Term", y = "Count") +
          ggplot2::theme_minimal() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
        dev.off()
        pdf("GSEA_analysis_plots2_of_2.pdf", width = 35, height = 15)
        print(heatplot1)
        print(heatplot2)
        dev.off()


        # Save tables
        write.csv(gene_list, "gene_list.csv")
        write.csv(kegg_gene_list, "kegg_gene_list.csv")
  }}}
  #####
  ####
  ####
  if (TME) {
    cat("\033[33mSTARTING TME\033[0m\n")

    if (!is.null(pattern)) {
      # Remove outliers based on pattern
      filtered <- subset(counts_data, !grepl(pattern, rownames(counts_data)))
    } else {
      # If no pattern is provided, do not filter based on pattern
      filtered <- counts_data
    }

    if (remove_outliers) {
      # Remove columns corresponding to outliers
      counts_data <- filtered[, !colnames(filtered) %in% outliers]
      AnnotData <- AnnotData[!AnnotData[["id"]] %in% outliers, ]
    } else {
      counts_data <- counts_data
      AnnotData <- AnnotData
    }

    print(paste0("Inicial gene number: ", dim(counts_data)[1]))
    ## Normalización TPM
    cat("\033[32mTPM normalization performed and stored on tpm_counts.csv\033[0m\n")
    count2tpm <- function(countMat, idType = "Ensembl", org = "hsa",  source = "local", effLength = NULL, id = "id", gene_symbol = "symbol", length = "eff_length", check_data = FALSE){


  if(!org%in%c("hsa", "mmus")) stop(">>>== `org` must be hsa or mmus...")
  # requireNamespace("biomaRt")
  if(!is.matrix(countMat)){
    countMat<-as.matrix(countMat)
    countMat<-matrix(as.numeric(countMat), dim(countMat), dimnames = dimnames(countMat))
  }

  if(sum(is.na(countMat))>0|check_data){
    message(paste0("There are ", sum(is.na(countMat)) ," missing value in count matrix, these genes will be removed."))
    feas<-feature_manipulation(data = countMat, feature = rownames(countMat), is_matrix = T)
    countMat<-countMat[rownames(countMat)%in%feas,]

  }

  if(is.null(effLength) & source == "biomart"){
    datasets = paste0(c("hsapiens", "mmusculus", "btaurus", "cfamiliaris",
                        "ptroglodytes", "rnorvegicus", "sscrofa"), "_gene_ensembl")
    type = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", "start_position", "end_position")
    if(org =="mmu") type[3] = "mgi_symbol"
    # listEnsemblArchives()
    # listMarts()
    # listAttributes()
    ds <- datasets[grepl(org, datasets)]
    mart <- biomaRt::useMart(host = "https://www.ensembl.org", biomart = 'ENSEMBL_MART_ENSEMBL', dataset = ds)
    ensembl <- biomaRt::getBM(attributes=type, mart = mart)
    #######################################

    ensembl$Length <- abs(ensembl$end_position - ensembl$start_position)

    message(">>>--- This function is being optimised and we strongly recommend that you should set `source` as `local`....")
    #######################################
    if(toupper(idType) == "ENSEMBL"){

      len <- ensembl[match(rownames(countMat),ensembl$ensembl_gene_id), "Length"]
      rownames(countMat) = ensembl[match(rownames(countMat),ensembl$ensembl_gene_id), 3]
    }else if(toupper(idType) == "SYMBOL"){
      len <- ensembl[match(rownames(countMat), ensembl[,3]), "Length"]
    }else if(toupper(idType) == "ENTREZ"){
      len <- ensembl[match(rownames(countMat), ensembl[,2]), "Length"]
    }else{
      stop("Please input right type of gene name, such as `ensembl`, `entrez`, or `symbol` ...")
    }
  }


  if(source == "local" & tolower(idType) == "ensembl" & org == "hsa") {

    rownames(countMat) <- substring(rownames(countMat), 1, 15)
    data("anno_grch38", package = "IOBR")
    message(">>>--- Using variables (anno_grch38) and gene lengths (eff_length)  built into the IOBR package to perform TPM transformation")
    message(">>>--- The gene lengths (eff_length) was estimated by function `getGeneLengthAndGCContent` from EDASeq package with default parameters at 2023-02-10")

    length_ensembl<-anno_grch38[,c("id", "eff_length", "symbol")]
    length_ensembl<-length_ensembl[order(length_ensembl$eff_length, decreasing = T), ]

    countMat<- countMat[rownames(countMat)%in%length_ensembl$id,]

    if(dim(countMat)[1]==0) stop("Identifier of matrix is not match to references.")
    length_ensembl<-length_ensembl[length_ensembl$id%in%rownames(countMat),]
    len<- length_ensembl[match(rownames(countMat), length_ensembl$id), "eff_length"]
    rownames(countMat)<- length_ensembl[match(rownames(countMat),length_ensembl$id), 3]

    countMat <- matrix(as.numeric(countMat), dim(countMat), dimnames = dimnames(countMat))
  }else if(source == "local" & tolower(idType) == "entrez"  & org == "hsa"){


    message(">>>--- This is a fuzzy calculation. We recommend that users provide expression matrices with ENSEMBL as row names")
    message(">>>--- Using variables (anno_grch38) and gene lengths (eff_length)  built into the IOBR package to perform TPM transformation")
    message(">>>--- The gene lengths (eff_length) was estimated by function `getGeneLengthAndGCContent` from EDASeq package with default parameters at 2023-02-10")

    length_ensembl<-anno_grch38[,c("entrez", "eff_length", "symbol")]
    length_ensembl<-length_ensembl[order(length_ensembl$eff_length, decreasing = T), ]
    colnames(length_ensembl)[1]<-"id"
    length_ensembl<-length_ensembl[!duplicated(length_ensembl$id), ]

    countMat<-countMat[rownames(countMat)%in%length_ensembl$id, ]
    if(dim(countMat)[1]==0) stop("Identifier of matrix is not match to references.")
    length_ensembl<-length_ensembl[length_ensembl$id%in%rownames(countMat), ]
    len<- length_ensembl[match(rownames(countMat), length_ensembl$id), "eff_length"]
    rownames(countMat) <- length_ensembl[match(rownames(countMat),length_ensembl$id), 3]
    countMat<-matrix(as.numeric(countMat), dim(countMat), dimnames = dimnames(countMat))

  }else if(source == "local" & tolower(idType) == "symbol"  & org == "hsa"){

    message(">>>--- This is a fuzzy calculation. We recommend that users provide expression matrices with ENSEMBL as row names")
    message(">>>--- Using variables (anno_grch38) and gene lengths (eff_length)  built into the IOBR package to perform TPM transformation")
    message(">>>--- The gene lengths (eff_length) was estimated by function `getGeneLengthAndGCContent` from EDASeq package with default parameters at 2023-02-10")
    length_ensembl<-anno_grch38[, c("symbol", "eff_length", "gc")]
    length_ensembl<-length_ensembl[order(length_ensembl$eff_length, decreasing = T), ]

    colnames(length_ensembl)[1] <-"id"

    # print(head(length_ensembl))
    # print(head(countMat))

    length_ensembl <- length_ensembl[!duplicated(length_ensembl$id), ]
    countMat <- countMat[rownames(countMat)%in%length_ensembl$id, ]

    if(dim(countMat)[1]==0) stop("Identifier of matrix is not match to references.")
    length_ensembl<-length_ensembl[length_ensembl$id%in%rownames(countMat),]
    len<- length_ensembl[match(rownames(countMat), length_ensembl$id), "eff_length"]

    rownames(countMat)<- length_ensembl[match(rownames(countMat),length_ensembl$id), 1]
    countMat<-matrix(as.numeric(countMat), dim(countMat), dimnames = dimnames(countMat))
  }

  #######################################################################
  if(source == "local" & tolower(idType) == "ensembl" & org == "mmus") {

    message(">>>--- Using variables (anno_gc_vm32) and gene lengths (eff_length)  built into the IOBR package to perform TPM transformation")
    message(">>>--- The gene lengths (eff_length) was estimated by function `getGeneLengthAndGCContent` from EDASeq package with default parameters at 2023-02-10")

    length_ensembl<-anno_gc_vm32[,c("id", "eff_length", "symbol")]
    length_ensembl<-length_ensembl[order(length_ensembl$eff_length, decreasing = T), ]

    countMat<-countMat[rownames(countMat)%in%length_ensembl$id,]
    if(dim(countMat)[1]==0) stop("Identifier of matrix is not match to references.")
    length_ensembl<-length_ensembl[length_ensembl$id%in%rownames(countMat),]
    len<- length_ensembl[match(rownames(countMat), length_ensembl$id), "eff_length"]
    rownames(countMat)<- length_ensembl[match(rownames(countMat),length_ensembl$id), 3]
    countMat<-matrix(as.numeric(countMat), dim(countMat), dimnames = dimnames(countMat))
  }else if(source == "local" & tolower(idType) == "mgi"  & org == "mmus"){

    message(">>>--- This is a fuzzy calculation. We recommend that users provide expression matrices with ENSEMBL ID as row names")
    message(">>>--- Using variables (anno_gc_vm32) and gene lengths (eff_length)  built into the IOBR package to perform TPM transformation")
    message(">>>--- The gene lengths (eff_length) was estimated by function `getGeneLengthAndGCContent` from EDASeq package with default parameters at 2023-02-10")

    length_ensembl<-anno_gc_vm32[,c("mgi_id", "eff_length", "symbol")]
    length_ensembl<-length_ensembl[order(length_ensembl$eff_length, decreasing = T), ]
    colnames(length_ensembl)[1]<-"id"
    length_ensembl<-length_ensembl[!duplicated(length_ensembl$id), ]

    countMat<-countMat[rownames(countMat)%in%length_ensembl$id,]
    if(dim(countMat)[1]==0) stop("Identifier of matrix is not match to references.")
    length_ensembl<-length_ensembl[length_ensembl$id%in%rownames(countMat),]
    len<- length_ensembl[match(rownames(countMat), length_ensembl$id), "eff_length"]
    rownames(countMat) <- length_ensembl[match(rownames(countMat),length_ensembl$id), 3]
    countMat<-matrix(as.numeric(countMat), dim(countMat), dimnames = dimnames(countMat))

  }else if(source == "local" & tolower(idType) == "symbol"  & org == "mmus"){

    message(">>>--- This is a fuzzy calculation. We recommend that users provide expression matrices with ENSEMBL ID as row names")
    message(">>>--- Using variables (anno_gc_vm32) and gene lengths (eff_length)  built into the IOBR package to perform TPM transformation")
    message(">>>--- The gene lengths (eff_length) was estimated by function `getGeneLengthAndGCContent` from EDASeq package with default parameters at 2023-02-10")
    length_ensembl<-anno_gc_vm32[,c("symbol", "eff_length", "gc")]
    length_ensembl<-length_ensembl[order(length_ensembl$eff_length, decreasing = T), ]

    colnames(length_ensembl)[1]<-"id"
    length_ensembl<-length_ensembl[!duplicated(length_ensembl$id), ]

    countMat<-countMat[rownames(countMat)%in%length_ensembl$id,]
    if(dim(countMat)[1]==0) stop("Identifier of matrix is not match to references.")
    length_ensembl<-length_ensembl[length_ensembl$id%in%rownames(countMat),]
    len<- length_ensembl[match(rownames(countMat), length_ensembl$id), "eff_length"]
    rownames(countMat)<- length_ensembl[match(rownames(countMat),length_ensembl$id), 1]
    countMat<-matrix(as.numeric(countMat), dim(countMat), dimnames = dimnames(countMat))
  }

  #########################################################################
  if(!is.null(effLength)){
    effLength<-as.data.frame(effLength)
    colnames(effLength)[which(colnames(effLength)==id)]<-"id"
    colnames(effLength)[which(colnames(effLength)==length)]<-"eff_length"
    effLength<-effLength[!duplicated(effLength$id),]

    countMat<-as.matrix(countMat)
    countMat<-countMat[rownames(countMat)%in%effLength$id, ]
    effLength<-effLength[effLength$id%in%rownames(countMat), ]

    if(id!= gene_symbol){
      # countMat<-as.matrix(countMat)
      colnames(effLength)[which(colnames(effLength)==gene_symbol)]<-"gene_symbol"
      rownames(countMat)<- effLength[match(rownames(countMat),effLength$id), "gene_symbol"]

    }else{
      # countMat<-as.matrix(countMat)
      effLength$gene_symbol<-effLength$id
      # colnames(effLength)[which(colnames(effLength)==gene_symbol)]<-"gene_symbol"
      rownames(countMat)<- effLength[match(rownames(countMat),effLength$id), "gene_symbol"]
    }

    len<- effLength[match(rownames(countMat), effLength[,"gene_symbol"]), "eff_length"]

  }

  na_idx <- which(is.na(len))
  if(length(na_idx)>0){
    warning(paste0(">>>--- Omit ", length(na_idx), " genes of which length is not available !"))
    countMat <- countMat[!is.na(len), ]
    len = len[!is.na(len)]
  }
  #####################################
  tmp <- countMat / c(len/1000) # (`per million` scaling factor)
  TPM <- 1e6 * t(t(tmp) / colSums(tmp))
  TPM <- TPM[!is.na(rownames(TPM)),]
  TPM <- TPM[!rownames(TPM)==" ",]

  # TPM <- rownames_to_column(as.data.frame(TPM), var = "symbol")
  symbol.id = rownames(TPM)
  TPM = as.data.frame(TPM)
  TPM$symbol = symbol.id

  TPM <- remove_duplicate_genes(eset = TPM, column_of_symbol = "symbol")
  # TPM <- TPM[,!is.na(colnames(TPM))]
  # TPM <- TPM[,!colnames(TPM)==" "]
  return(TPM)
}
    suppressWarnings({
      tpm_counts <- count2tpm(counts_data,
                              idType = "Symbol",
                              org = "hsa",
                              source = "biomart")
      write.csv(tpm_counts, "tpm_counts.csv")
    })

    # Se almacenan los genes omitidos en la normalización TPM
    genes_omitidos <- base::setdiff(rownames(counts_data), rownames(tpm_counts))
    print(paste0("Number of genes omitted during TPM normalization due to their length not being available in Biomart:   ", dim(counts_data)[1]-dim(tpm_counts)[1]))
    tpm_counts <- as.data.frame(tpm_counts)

    if (!is.null(DEA) && DEA) {
      cat("\033[32mWe are going to use information from DEA\033[0m\n")
      dds <- DESeq2::estimateSizeFactors(dds)
      normalized_counts <- DESeq2::counts(dds, normalized = TRUE)
    } else {
      cat("Performing normalization.\n")
      design_formul <- as.formula("~ 1")
      colnames(AnnotData) <- gsub(" ", "_", colnames(AnnotData))
      AnnotData <- AnnotData[order(AnnotData$id), ]
      counts_data <- counts_data[, order(colnames(counts_data))]
      if (!identical(colnames(counts_data), AnnotData$id)) {
        stop("Column names of counts_data and IDs in AnnotData do not match.")
      }
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts_data, colData = AnnotData, design = design_formul)
      dds <- DESeq2::estimateSizeFactors(dds)
      normalized_counts <- DESeq2::counts(dds, normalized = TRUE)
    }

    ## Deconvolución
    imm_epic <- immunedeconv::deconvolute(tpm_counts, method = "epic")
    imm_qti <- immunedeconv::deconvolute(tpm_counts, method = "quantiseq")
    imm_xcell <- immunedeconv::deconvolute(tpm_counts, method = "xcell")
    cat("\033[32mresults of the TME will be stored in imm_epic.csv, imm_qti.csv and imm_xcell.csv \033[0m\n")
    write.csv(imm_epic, file = "imm_epic.csv")
    write.csv(imm_qti, file = "imm_qti.csv")
    write.csv(imm_xcell, file = "imm_xcell.csv")

    # Transponer y cambiar colnames
    std.im.df <- function(imm_df){
      imm_df <- as.data.frame(t(imm_df))
      celltype_imm <- imm_df[1,]
      imm_df <- imm_df[-1,]
      colnames(imm_df) <- celltype_imm
      imm_df[-1,]
      rn_imm <- rownames(imm_df)
      imm_df <- as.data.frame(sapply(imm_df, as.numeric))
      rownames(imm_df) <- rn_imm
      return(imm_df)
    }

    imm_epic <- std.im.df(imm_epic)
    imm_qti <- std.im.df(imm_qti)
    imm_xcell <- std.im.df(imm_xcell)

    # Se incluye la variable AnnotData$design_formula en los dataframes
    cat("\033[32mAre they in the same order?\033[0m\n")

    rownames(AnnotData)<- AnnotData$id
    AnnotData <- AnnotData[order(rownames(AnnotData)), ]
    imm_epic <- imm_epic[order(rownames(imm_epic)), ]
    imm_qti <- imm_qti[order(rownames(imm_qti)), ]
    imm_xcell <- imm_xcell[order(rownames(imm_xcell)), ]
    AnnotData$id
    rownames(imm_epic)
    rownames(imm_qti)
    rownames(imm_xcell)

    cat("\033[32mHave to be true.\033[0m\n")
    cat("\033[32mEPIC\033[0m\n")
    print(all(rownames(AnnotData)==rownames(imm_epic)))
    cat("\033[32mqti\033[0m\n")
    print(all(rownames(AnnotData)==rownames(imm_qti)))
    cat("\033[32mxcell\033[0m\n")
    print(all(rownames(AnnotData)==rownames(imm_xcell)))

    imm_epic[[design_formula]] <- factor(AnnotData[[design_formula]])
    imm_qti[[design_formula]] <- factor(AnnotData[[design_formula]])
    imm_xcell[[design_formula]] <- factor(AnnotData[[design_formula]])

    replace_space_with_underscore <- function(df) {
      colnames(df) <- gsub(" ", "_", colnames(df))
      return(df)
    }
    imm_epic <- replace_space_with_underscore(imm_epic)
    imm_qti <- replace_space_with_underscore(imm_qti)
    imm_xcell <- replace_space_with_underscore(imm_xcell)

    ############# imm_epic
    design_formula_sym <- rlang::sym(design_formula)
    # Vector with column names except the last one
    column_names <- names(imm_epic)[-ncol(imm_epic)]
    # Loop to generate plots for each column

    pdf("plots_imm_EPIC.pdf")
    for (col_name in column_names) {
      # Calculate means by design formula
      grouped_data <- dplyr::group_by(imm_epic, !!design_formula_sym)
      means_df <- dplyr::summarize(grouped_data, mean_value = mean(.data[[col_name]] * 100, na.rm = TRUE), .groups = "drop")

      # Perform t-test or ANOVA
      formula <- as.formula(paste0("`", col_name, "` ~ ", rlang::as_label(design_formula_sym)))
      p_value <- tryCatch({
        if (dplyr::n_distinct(imm_epic[[rlang::as_label(design_formula_sym)]]) == 2) {
          t.test(formula, data = imm_epic)$p.value
        } else {
          summary(stats::aov(formula, data = imm_epic))[[1]]$`Pr(>F)`[1]
        }
      }, error = function(e) {
        NA  # In case of error, return NA for p-value
      })
      # Generate plot
      plot <- ggplot2::ggplot(imm_epic, ggplot2::aes(x = !!design_formula_sym , y = .data[[col_name]] * 100,
                                   fill = !!design_formula_sym, color = !!design_formula_sym)) +
        ggplot2::geom_jitter(alpha = 1, width = 0.3, height = 0) +
        ggplot2::geom_boxplot(fill = "white", alpha = 0.5, outlier.alpha = 1) +
        ggplot2::geom_point(data = means_df, ggplot2::aes(x = !!design_formula_sym, y = mean_value, fill = !!design_formula_sym),
                   shape = 22, color = "black", size = 3, stroke = 1.5,
                   show.legend = F) +
        ggplot2::labs(x = NULL, y = "Abundance (%)", title = col_name, subtitle = "imm_epic") +
        ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1),
                           limits = c(0, NA)) +
        ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
        ggplot2::annotate("text", x = Inf, y = Inf, label = paste("p-value:", format(p_value, digits = 3)),
                 hjust = 1.1, vjust = 1.1, size = 5, color = "red")
      print(plot)
    }
    dev.off()

    ############# imm_qti
    # Vector with column names except the last one
    column_names <- names(imm_qti)[-ncol(imm_qti)]

    pdf("plots_imm_qti.pdf")
    for (col_name in column_names) {
      # Calculate means by design formula
      grouped_data <- dplyr::group_by(imm_qti, !!design_formula_sym)
      means_df <- dplyr::summarize(grouped_data, mean_value = mean(.data[[col_name]] * 100, na.rm = TRUE), .groups = "drop")

      # Perform t-test or ANOVA
      formula <- as.formula(paste0("`", col_name, "` ~ ", rlang::as_label(design_formula_sym)))
      p_value <- tryCatch({
        if (dplyr::n_distinct(imm_qti[[rlang::as_label(design_formula_sym)]]) == 2) {
          t.test(formula, data = imm_qti)$p.value
        } else {
          summary(stats::aov(formula, data = imm_qti))[[1]]$`Pr(>F)`[1]
        }
      }, error = function(e) {
        NA  # In case of error, return NA for p-value
      })

      # Generate plot
      plot <- ggplot2::ggplot(imm_qti, ggplot2::aes(x = !!design_formula_sym , y = .data[[col_name]] * 100,
                                                    fill = !!design_formula_sym, color = !!design_formula_sym)) +
        ggplot2::geom_jitter(alpha = 1, width = 0.3, height = 0) +
        ggplot2::geom_boxplot(fill = "white", alpha = 0.5, outlier.alpha = 1) +
        ggplot2::geom_point(data = means_df, ggplot2::aes(x = !!design_formula_sym, y = mean_value, fill = !!design_formula_sym),
                            shape = 22, color = "black", size = 3, stroke = 1.5,
                            show.legend = F) +
        ggplot2::labs(x = NULL, y = "Abundance (%)", title = col_name, subtitle = "imm_qti") +
        ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1),
                                    limits = c(0, NA)) +
        ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
        ggplot2::annotate("text", x = Inf, y = Inf, label = paste("p-value:", format(p_value, digits = 3)),
                          hjust = 1.1, vjust = 1.1, size = 5, color = "red")
      print(plot)
    }
    dev.off()

    ############# imm_xcell
    # Vector with column names except the last one
    column_names <- names(imm_xcell)[-ncol(imm_xcell)]
    pdf("plots_imm_xcell.pdf")
    for (col_name in column_names) {
      # Calculate means by design formula
      grouped_data <- dplyr::group_by(imm_xcell, !!design_formula_sym)
      means_df <- dplyr::summarize(grouped_data, mean_value = mean(.data[[col_name]] * 100, na.rm = TRUE), .groups = "drop")

      # Perform t-test or ANOVA
      formula <- as.formula(paste0("`", col_name, "` ~ ", rlang::as_label(design_formula_sym)))
      p_value <- tryCatch({
        if (dplyr::n_distinct(imm_xcell[[rlang::as_label(design_formula_sym)]]) == 2) {
          t.test(formula, data = imm_xcell)$p.value
        } else {
          summary(stats::aov(formula, data = imm_xcell))[[1]]$`Pr(>F)`[1]
        }
      }, error = function(e) {
        NA  # In case of error, return NA for p-value
      })
      # Generate plot
      plot <- ggplot2::ggplot(imm_xcell, ggplot2::aes(x = !!design_formula_sym , y = .data[[col_name]] * 100,
                                                      fill = !!design_formula_sym, color = !!design_formula_sym)) +
        ggplot2::geom_jitter(alpha = 1, width = 0.3, height = 0) +
        ggplot2::geom_boxplot(fill = "white", alpha = 0.5, outlier.alpha = 1) +
        ggplot2::geom_point(data = means_df, ggplot2::aes(x = !!design_formula_sym, y = mean_value, fill = !!design_formula_sym),
                            shape = 22, color = "black", size = 3, stroke = 1.5,
                            show.legend = F) +
        ggplot2::labs(x = NULL, y = "Abundance (%)", title = col_name, subtitle = "imm_xcell") +
        ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1),
                                    limits = c(0, NA)) +
        ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
        ggplot2::annotate("text", x = Inf, y = Inf, label = paste("p-value:", format(p_value, digits = 3)),
                          hjust = 1.1, vjust = 1.1, size = 5, color = "red")
      print(plot)
      }
    dev.off()
    # Composición celular del TME i por grupo

    # Function to generate the bar plot for each sample
    plot_bar <- function(df, paleta, titulo, legend.position) {
      # Convertir las filas en una columna llamada "Sample"
      df <- tibble::rownames_to_column(df, var = "Sample")
      df <- tidyr::pivot_longer(
        df, cols = colnames(df)[2:(ncol(df) - 1)],
        names_to = "Cell_Type", values_to = "Value")
      df <- dplyr::mutate(df,
        Sample = factor(Sample, levels = rev(unique(Sample))),
        Cell_Type = factor(Cell_Type, levels = rev(unique(Cell_Type))))

      p <- ggplot2::ggplot(df, ggplot2::aes(x = Sample, y = Value, fill = Cell_Type)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::labs(title = titulo,
                      x = "Samples",
                      y = "Cell Fraction (%)") +
        ggplot2::coord_flip() +
        ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE)) +
        ggplot2::scale_fill_manual(values = paleta) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = legend.position,
                       axis.text.y = ggplot2::element_text(size = 5)) +
        ggplot2::scale_y_continuous(labels = scales::percent)

      return(p)
    }

    plot_bar_group <- function(df, paleta, titulo, design_formula, legend_position = "right") {
      suppressWarnings({
        design_formula_sym <- rlang::sym(design_formula)
        niveles_tipo_cel <- colnames(df)[1:(ncol(df) - 1)]
        df_rownames <- tibble::rownames_to_column(df, var = "Sample")
        df_long <- tidyr::pivot_longer(
          df_rownames,
          cols = niveles_tipo_cel,
          names_to = "Cell_Type",
          values_to = "Value")
        df_grouped <- dplyr::group_by(df_long, !!design_formula_sym, Cell_Type)
        df_summarised <- dplyr::summarise(df_grouped,
          Average = mean(Value, na.rm = TRUE),
          .groups = "drop")
        df_ungrouped <- dplyr::ungroup(df_summarised)
        promedios <- dplyr::mutate(df_ungrouped,
          !!design_formula_sym := factor(!!design_formula_sym, levels = rev(unique(!!design_formula_sym))),
          Cell_Type = factor(Cell_Type, levels = rev(niveles_tipo_cel)))
          p <- ggplot2::ggplot(promedios, ggplot2::aes(x = !!design_formula_sym, y = Average, fill = Cell_Type)) +
          ggplot2::geom_bar(stat = "identity") +
          ggplot2::labs(title = titulo,
                        x = "Samples",
                        y = "Cell Fraction (%)") +
          ggplot2::coord_flip() +
          ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE)) +
          ggplot2::scale_fill_manual(values = paleta) +
          ggplot2::theme_minimal() +
          ggplot2::theme(legend.position = legend_position,
                         axis.text.y = ggplot2::element_text(size = 5)) +
          ggplot2::scale_y_continuous(labels = scales::percent)
        return(p)
      })
    }

  # Function to combine both plots into one
  plot_combined <- function(df, paleta, titulo_individual, titulo_grupo, design_formula, legend_position = "right") {
    p1 <- plot_bar(df, paleta, titulo_individual, legend_position)
    p2 <- plot_bar_group(df, paleta, titulo_grupo, design_formula, legend_position)

    combined_plot <- ggpubr::ggarrange(p1, p2, ncol = 1, nrow = 2, heights = c(1, 1))
    return(combined_plot)
  }

    # Example usage with extended palette for larger datasets
    paleta_imm <- c("grey95","#FB8072","#FFED6F","#6F6C87","#94DFD1","#FDB462", "#B3DE69", "#FFB1D9")
    paleta_qti <- c("grey95","#B3DE69","#6F6C87","#94DFD1","#FDB462", "#FB8072","#FFFFB3","#8BB07A","#FFED6F","#80B1D3","#FFB1D9")
    paleta_extendida <- c("grey95","#2CA02C","#E6F8E0","#7F7F7F","#FF8000","#DF7401","#6F6C80","#F5A9E1",
                          "#F7D358","#6F6C99","#FB8072","#FFFFB3","#F5A9F2",
                          "#8BB07A","#E377C2","#FE2E2E","#FFED6F","#1F77B4","#80B1D3","#0489B1","#FF0000","#8BB07A",
                          "#0000FF","#D6616B","#58FA82","#98DF8A",
                          "#8C564B", "#C49C94","#CEE3F6","#E0F8EC","#58FAF4", "#ADD8E6","#9EDAE5","#94DFD1","#FFA500",
                          "#FF7F0E","#FDB462","#FFB1D9","#B3DE69")
    # Assuming 'imm_epic', 'imm_qti', and 'imm_xcell' dataframes are already loaded
    combined_plot_EPIC <- plot_combined(imm_epic, paleta_imm, "EPIC Individual", "EPIC Average", design_formula, "right")
    combined_plot_quanTIseq <- plot_combined(imm_qti, paleta_qti, "quanTIseq Individual", "quanTIseq Average", design_formula, "right")
    plot_bar <- function(df, paleta, titulo, legend.position) {
      # Convertir las filas en una columna llamada "Sample"
      df <- tibble::rownames_to_column(df, var = "Sample")
      df <- tidyr::pivot_longer(
        df, cols = colnames(df)[2:(ncol(df) - 1)],
        names_to = "Cell_Type", values_to = "Value")
      df <- dplyr::mutate(df,
                          Sample = factor(Sample, levels = rev(unique(Sample))),
                          Cell_Type = factor(Cell_Type, levels = rev(unique(Cell_Type))))

      p <- ggplot2::ggplot(df, ggplot2::aes(x = Sample, y = Value, fill = Cell_Type)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::labs(title = titulo,
                      x = "Samples",
                      y = "Enrichment Scores") +
        ggplot2::coord_flip() +
        ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE)) +
        ggplot2::scale_fill_manual(values = paleta) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = legend.position,
                       axis.text.y = ggplot2::element_text(size = 5)) +
        ggplot2::scale_y_continuous(labels = scales::percent)

      return(p)
    }

    plot_bar_group <- function(df, paleta, titulo, design_formula, legend_position = "right") {
      suppressWarnings({
        design_formula_sym <- rlang::sym(design_formula)
        niveles_tipo_cel <- colnames(df)[1:(ncol(df) - 1)]
        df_rownames <- tibble::rownames_to_column(df, var = "Sample")
        df_long <- tidyr::pivot_longer(
          df_rownames,
          cols = niveles_tipo_cel,
          names_to = "Cell_Type",
          values_to = "Value")
        df_grouped <- dplyr::group_by(df_long, !!design_formula_sym, Cell_Type)
        df_summarised <- dplyr::summarise(df_grouped,
                                          Average = mean(Value, na.rm = TRUE),
                                          .groups = "drop")
        df_ungrouped <- dplyr::ungroup(df_summarised)
        promedios <- dplyr::mutate(df_ungrouped,
                                   !!design_formula_sym := factor(!!design_formula_sym, levels = rev(unique(!!design_formula_sym))),
                                   Cell_Type = factor(Cell_Type, levels = rev(niveles_tipo_cel)))
        p <- ggplot2::ggplot(promedios, ggplot2::aes(x = !!design_formula_sym, y = Average, fill = Cell_Type)) +
          ggplot2::geom_bar(stat = "identity") +
          ggplot2::labs(title = titulo,
                        x = "Samples",
                        y = "Enrichment Scores") +
          ggplot2::coord_flip() +
          ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE)) +
          ggplot2::scale_fill_manual(values = paleta) +
          ggplot2::theme_minimal() +
          ggplot2::theme(legend.position = legend_position,
                         axis.text.y = ggplot2::element_text(size = 5)) +
          ggplot2::scale_y_continuous(labels = scales::percent)
        return(p)
      })
    }

    # Function to combine both plots into one
    plot_combined <- function(df, paleta, titulo_individual, titulo_grupo, design_formula, legend_position = "right") {
      p1 <- plot_bar(df, paleta, titulo_individual, legend_position)
      p2 <- plot_bar_group(df, paleta, titulo_grupo, design_formula, legend_position)

      combined_plot <- ggpubr::ggarrange(p1, p2, ncol = 1, nrow = 2, heights = c(1, 1))
      return(combined_plot)
    }
    combined_plot_xCell <- plot_combined(imm_xcell, paleta_extendida, "xCell Individual", "xCell Average", design_formula, "right")
    combined_plot_xCell <- plot_combined(imm_xcell, paleta_extendida, "xCell Individual", "xCell Average", design_formula, "right")

    pdf("plot_cell_fraction_Average_cell_fraction_EPIC.pdf", width = 11, height = 14)
    print(combined_plot_EPIC)
    dev.off()

    pdf("plot_cell_fraction_Average_cell_fraction_quanTIseq.pdf", width = 11, height = 14)
    print(combined_plot_quanTIseq)
    dev.off()

    pdf("plot_cell_fraction_Average_cell_fraction_xCell.pdf", width = 11, height = 14)
    print(combined_plot_xCell)
    dev.off()

    ########################
    trans_formato_largo <- function(df, design_formula_sym) {
      df_largo <- tidyr::pivot_longer(df, cols = -dplyr::all_of(design_formula_sym), names_to = "Cell_Type",
                                      values_to = "Fraction")
      return(df_largo)
    }

    # Function to perform Shapiro-Wilk normality test
    prueba_norm <- function(df, design_formula_sym) {
      df_largo <- trans_formato_largo(df, design_formula_sym)

      resultados_normalidad <- dplyr::group_by(df_largo, Cell_Type)
      resultados_normalidad <- dplyr::summarise(resultados_normalidad,
                                                shapiro_test = list(stats::shapiro.test(Fraction)))
      resultados_normalidad <- dplyr::mutate(resultados_normalidad,
                                             p.value = purrr::map_dbl(shapiro_test, "p.value"))


      return(resultados_normalidad)
    }

    # Function to check if all values in a column are the same
    check_column_equal <- function(column) {
      all_equal <- length(unique(column)) == 1
      return(all_equal)
    }

    # Function to filter valid columns
    filter_valid_columns <- function(df, design_formula_sym, dataset_name) {
      cat("\033[32mChecking if All Values in Each Column Are Equal in dataset:", dataset_name, "\033[0m\n")

      equality_results <- base::sapply(df, check_column_equal)
      unequal_columns_count <- base::sum(!equality_results)

      if (unequal_columns_count > 0) {
        cat("\033[33mIn dataset", dataset_name, "the following columns have the same value in all rows and will be excluded from the analysis:\033[0m\n")
        print(names(df)[equality_results])
      } else {
        cat("\033[32mAll columns in dataset", dataset_name, "have variability. No columns to exclude.\033[0m\n")
      }

      valid_columns <- !equality_results
      df_filtered <- df[, valid_columns]

      return(df_filtered)
    }

    # Function to perform parametric tests (t-test and ANOVA)
    perform_parametric_tests <- function(df, design_formula_sym) {
      df_largo <- trans_formato_largo(df, design_formula_sym)
      cell_types <- unique(df_largo$Cell_Type)

      results <- list()

      for (cell_type in cell_types) {

        df_cell_type <- dplyr::filter(df_largo, Cell_Type == cell_type)
        groups <- unique(df_cell_type[[design_formula_sym]])


        if (length(groups) == 2) {
          # Perform t-test
          t_test_result <- t.test(Fraction ~ df_cell_type[[design_formula_sym]], data = df_cell_type)
          results[[cell_type]] <- list(Test = "t-test", p.value = t_test_result$p.value)
        } else {
          # Perform ANOVA
          anova_result <- stats::aov(Fraction ~ df_cell_type[[design_formula_sym]], data = df_cell_type)
          p_value <- summary(anova_result)[[1]]$`Pr(>F)`[1]
          results[[cell_type]] <- list(Test = "ANOVA", p.value = p_value)
        }
      }

      results_df <- dplyr::bind_rows(lapply(names(results), function(cell_type) {
        result <- results[[cell_type]]
        tibble::tibble(Cell_Type = cell_type, Test = result$Test, p.value = result$p.value)
      }))

      return(results_df)
    }

    # Apply the filtering and normality test for each dataset
    cat("\033[33mApply the filtering and normality test for each dataset\033[0m\n")
    imm_qti_filtered <- filter_valid_columns(imm_qti, design_formula_sym, "quanTIseq")
    imm_epic_filtered <- filter_valid_columns(imm_epic, design_formula_sym, "EPIC")
    imm_xcell_filtered <- filter_valid_columns(imm_xcell, design_formula_sym, "xcell")

    # Perform the normality test
    if (ncol(imm_qti_filtered) > 1) {
      cat("\033[32mPerforming Shapiro-Wilk test for quanTIseq\033[0m\n")
      resultados_norm_imm_qti <- prueba_norm(imm_qti_filtered, design_formula_sym)
      print(resultados_norm_imm_qti)
      # Perform parametric tests if normality is satisfied
      cat("\033[32mPerforming parametric tests for quanTIseq\033[0m\n")
      parametric_results_imm_qti <- perform_parametric_tests(imm_qti_filtered, design_formula_sym)
      print(parametric_results_imm_qti)
    } else {
      cat("\033[31mCan't perform Shapiro-Wilk test for quanTIseq. No valid columns available.\033[0m\n")
    }

    if (ncol(imm_epic_filtered) > 1) {
      cat("\033[32mPerforming Shapiro-Wilk test for EPIC\033[0m\n")
      resultados_norm_epic <- prueba_norm(imm_epic_filtered, design_formula_sym)
      print(resultados_norm_epic)
      # Perform parametric tests if normality is satisfied
      cat("\033[32mPerforming parametric tests for EPIC\033[0m\n")
      parametric_results_epic <- perform_parametric_tests(imm_epic_filtered, design_formula_sym)
      print(parametric_results_epic)
    } else {
      cat("\033[31mCan't perform Shapiro-Wilk test for EPIC. No valid columns available.\033[0m\n")
    }

    if (ncol(imm_xcell_filtered) > 1) {
      cat("\033[32mPerforming Shapiro-Wilk test for xcell\033[0m\n")
      resultados_norm_xcell <- prueba_norm(imm_xcell_filtered, design_formula_sym)
      print(resultados_norm_xcell)
      # Perform parametric tests if normality is satisfied
      cat("\033[32mPerforming parametric tests for xcell\033[0m\n")
      parametric_results_xcell <- perform_parametric_tests(imm_xcell_filtered, design_formula_sym)
      print(parametric_results_xcell)
    } else {
      cat("\033[31mCan't perform Shapiro-Wilk test for xcell. No valid columns available.\033[0m\n")
    }

    ############# Heatmaps
    # Transponer y estandarizar por filas
    h_imm_epic <- as.data.frame(t(imm_epic))
    ## Delete the row uncharacterized cell
    h_imm_epic <- head(h_imm_epic, -1)
    rn_himmepic <- rownames(h_imm_epic)
    # "Cancer associated fibroblast" will become CAFs
    rn_himmepic[2] <- "CAFs"
    cl_himmepic <- colnames(h_imm_epic)
    h_imm_epic <- apply(h_imm_epic, 2, as.numeric)
    h_imm_epic <- t(apply(h_imm_epic,1,scale))
    rownames(h_imm_epic) <- rn_himmepic
    colnames(h_imm_epic) <- cl_himmepic
    h_imm_epic <- as.data.frame(h_imm_epic)


    # Transponer y estandarizar por filas
    h_imm_qti <- as.data.frame(t(imm_qti))
    h_imm_qti <- head(h_imm_qti, -1)
    rn_himmqti <- rownames(h_imm_qti)
    rn_himmqti[7] <- "T cell CD4+"
    rn_himmqti[9] <- "T cell regulatory"
    cl_himmqti <- colnames(h_imm_qti)
    h_imm_qti <- apply(h_imm_qti, 2, as.numeric)
    h_imm_qti <- t(apply(h_imm_qti,1,scale))
    rownames(h_imm_qti) <- rn_himmqti
    colnames(h_imm_qti) <- cl_himmqti
    h_imm_qti <- as.data.frame(h_imm_qti)

    # Crear dataframe para heatmap
    h_imm_xcell <- as.data.frame(imm_xcell)
    # Transponer y estandarizar por filas
    h_imm_xcell <- as.data.frame(t(h_imm_xcell))
    h_imm_xcell <- head(h_imm_xcell, -4)

    # Poblaciones celulares no interesantes
    xcell_row_delete <- c('Common lymphoid progenitor', 'Common myeloid progenitor',
                          'Granulocyte-monocyte progenitor', 'Hematopoietic stem cell')

    h_imm_xcell <- dplyr::filter(h_imm_xcell, !rownames(h_imm_xcell) %in% xcell_row_delete)

    rn_himmxcell <- rownames(h_imm_xcell)
    cl_himmxcell <- colnames(h_imm_xcell)
    h_imm_xcell <- apply(h_imm_xcell, 2, as.numeric)
    h_imm_xcell <- t(apply(h_imm_xcell,1,scale))
    rownames(h_imm_xcell) <- rn_himmxcell
    colnames(h_imm_xcell) <- cl_himmxcell
    h_imm_xcell <- as.data.frame(h_imm_xcell)

    # Factorización de variables y renombrado en AnnotData
    AnnotData[[design_formula]] <- as.factor(AnnotData[[design_formula]])

    ################################### HEATMAP
    cat("\033[33mGENERATING HEATMAPS\033[0m\n")
    remove_nan_rows <- function(df) {
      # Remove rows with any NaN values
      df_cleaned <- df[!apply(df, 1, function(row) any(is.nan(row))), ]
      return(df_cleaned)
    }

    # Generar el heatmap para h_imm_epic
    combined_data <- remove_nan_rows(h_imm_epic)

    if (ncol(combined_data) < 3) {
      cat("Cannot generate the heatmap because the dataset has less than 3 columns after removing rows with NaN values.\n")
    } else {
      if (any(is.na(combined_data))) {
        cat("Cannot generate the heatmap because there are NaN values in the data.\n")
      } else {
        col_ann_data <- AnnotData[colnames(combined_data), , drop = FALSE]

        if (ncol(combined_data) == nrow(col_ann_data)) {
          col_ann_data <- as.data.frame(col_ann_data)
        } else {
          stop("Dimensions of col_ann_data and combined_data do not match")
        }

        annotation_col <- as.data.frame(col_ann_data[, design_formula, drop = FALSE])
        rownames(annotation_col) <- colnames(combined_data)

        HeatmapEPIC <- pheatmap::pheatmap(
          as.matrix(combined_data),
          cluster_rows = TRUE,
          cluster_cols = TRUE,
          annotation_col = annotation_col,
          fontsize = 9,
          color = colorRampPalette(c("#4793AF", "white", "#013649"))(50),
          border_color = "grey60",
          main = "Heatmap EPIC",
          legend = TRUE,
          angle_col = 45,
          silent = TRUE
        )
      }
    }

    # Generar el heatmap para h_imm_qti
    combined_data <- remove_nan_rows(h_imm_qti)

    if (ncol(combined_data) < 3) {
      cat("Cannot generate the heatmap because the dataset has less than 3 columns after removing rows with NaN values.\n")
    } else {
      if (any(is.na(combined_data))) {
        cat("Cannot generate the heatmap because there are NaN values in the data.\n")
      } else {
        col_ann_data <- AnnotData[colnames(combined_data), , drop = FALSE]

        if (ncol(combined_data) == nrow(col_ann_data)) {
          col_ann_data <- as.data.frame(col_ann_data)
        } else {
          stop("Dimensions of col_ann_data and combined_data do not match")
        }

        annotation_col <- as.data.frame(col_ann_data[, design_formula, drop = FALSE])
        rownames(annotation_col) <- colnames(combined_data)

        Heatmap_qti <- pheatmap::pheatmap(
          as.matrix(combined_data),
          cluster_rows = TRUE,
          cluster_cols = TRUE,
          annotation_col = annotation_col,
          fontsize = 9,
          color = colorRampPalette(c("#4793AF", "white", "#013649"))(50),
          border_color = "grey60",
          main = "Heatmap QTI",
          legend = TRUE,
          angle_col = 45,
          silent = FALSE
        )
      }
    }

    # Generar el heatmap para h_imm_xcell
    combined_data <- remove_nan_rows(h_imm_xcell)

    if (ncol(combined_data) < 3) {
      cat("Cannot generate the heatmap because the dataset has less than 3 columns after removing rows with NaN values.\n")
    } else {
      if (any(is.na(combined_data))) {
        cat("Cannot generate the heatmap because there are NaN values in the data.\n")
      } else {
        col_ann_data <- AnnotData[colnames(combined_data), , drop = FALSE]

        if (ncol(combined_data) == nrow(col_ann_data)) {
          col_ann_data <- as.data.frame(col_ann_data)
        } else {
          stop("Dimensions of col_ann_data and combined_data do not match")
        }

        annotation_col <- as.data.frame(col_ann_data[, design_formula, drop = FALSE])
        rownames(annotation_col) <- colnames(combined_data)

        Heatmap_xcell <- pheatmap::pheatmap(
          as.matrix(combined_data),
          cluster_rows = TRUE,
          cluster_cols = TRUE,
          annotation_col = annotation_col,
          fontsize = 9,
          color = colorRampPalette(c("#4793AF", "white", "#013649"))(50),
          border_color = "grey60",
          main = "Heatmap XCELL",
          legend = TRUE,
          angle_col = 45,
          silent = TRUE
        )
      }
    }

    cat("\033[32mHeatmap of qti, EPIC and xcell will be stored on plots_TME_heatmap.pdf\033[0m\n")

    if (exists("Heatmap_qti")) {
      pdf("plots_TME_Heatmap_qti.pdf", width = 15, height = 14)
      print(Heatmap_qti)
      dev.off()
    } else {
      cat("Heatmap_qti object does not exist.\n")
    }

    if (exists("HeatmapEPIC")) {
      pdf("plots_TME_HeatmapEPIC.pdf", width = 15, height = 14)
      print(HeatmapEPIC)
      dev.off()
    } else {
      cat("HeatmapEPIC object does not exist.\n")
    }
    if (exists("Heatmap_xcell")) {
      pdf("plots_TME_Heatmap_xcell.pdf", width = 15, height = 14)
      print(Heatmap_xcell)
      dev.off()
    } else {
      cat("Heatmap_xcell object does not exist.\n")
    }


    }else {cat("\033[32mTME analysis skipped.\033[0m\n")}

  #####
  #####
  #####
    if (survival_analysis) {
    cat("\033[33mSTARTING SURVIVAL ANALYSIS\033[0m\n")
    if (is.null(col_data[[time]]) || is.null(col_data[[variable_01]])) {
      stop("Variables for survival analysis are required.")
    }

    if (!is.null(pattern)) {
      # Remove outliers based on pattern
      filtered <- subset(counts_data, !grepl(pattern, rownames(counts_data)))
    } else {
      filtered <- counts_data
    }

    if (remove_outliers) {
      counts_filtered <- filtered[, !colnames(filtered) %in% outliers]
      AnnotData <- col_data[!col_data[["id"]] %in% outliers, ]
    } else {
      counts_filtered <- counts_data
      AnnotData <- col_data
    }

    cat("\033[32mNormalizing data...\033[0m\n")
    colnames(AnnotData) <- gsub(" ", "_", colnames(AnnotData))
    col_data <- AnnotData[order(AnnotData$id), ]
    col_data<- as.data.frame(col_data)
    counts_filtered <- counts_filtered[, order(colnames(counts_filtered))]
    if (!identical(colnames(counts_filtered), col_data$id)) {
      stop("Column names of counts_filtered and IDs in col_data do not match.")
    }

    # Clean column names to avoid issues with special characters
    clean_column_names <- function(names) {
      gsub("[^[:alnum:]_]", "_", names)
    }

    # Create DESeqDataSet object
    if (!exists("dds")) {
      rownames(col_data) <- col_data$id
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts_filtered, colData = col_data, design = ~ 1)
      dds <- DESeq2::estimateSizeFactors(dds)
    }
    normalized_counts <- DESeq2::counts(dds, normalized = TRUE)
    df_t <- t(normalized_counts)

    # Checking for duplicate columns
    cat("\033[32mChecking for duplicate columns\033[0m\n")
    duplicated_columns <- colnames(df_t)[duplicated(colnames(df_t))]
    if (length(duplicated_columns) > 0) {
      print(paste("Duplicate columns found:", paste(duplicated_columns, collapse = ", ")))
    } else {
      print("No duplicate columns.")
    }

    rownames(col_data) <- col_data$id
    ids_data <- rownames(df_t)

    # Filter data
    cat("\033[32mFiltering data.\033[0m\n")
    subset_data <- dplyr::filter(col_data, id %in% ids_data)
    df_ta <- as.data.frame(df_t)
    df_ta$id <- rownames(df_ta)

    # Selecting genes
    if (!is.null(genes_to_use)) {
      cat("\033[32mUsing provided genes\033[0m\n")
      top_genes <- genes_to_use
    } else if (exists("res")) {
      cat("\033[32mSelecting TOP 10 genes with the lowest padj\033[0m\n")
      top_genes <- rownames(head(res[order(res$padj), ], 10))
    } else {
      stop("No genes provided or found in res object.")
    }
    
    selected_df_t <- df_ta[, top_genes, drop = FALSE]
    selected_df_t <- as.data.frame(selected_df_t)
    selected_df_t$id <- rownames(selected_df_t)
    col_data$id <- rownames(col_data)
    merged_data <- merge(col_data, selected_df_t, by = "id")
    merged_data[[time]] <- as.numeric(merged_data[[time]])

    cat("\033[32mStarting survival analysis.\033[0m\n")

    pdf("survival_analysis_plots_HTG_analyzer.pdf")
    top_genes_clean <- clean_column_names(top_genes)

    # Perform survival analysis for each gene
    for (i in top_genes_clean) {
      if (!is.numeric(merged_data[[i]])) {
        next
      }
      cat("\n")
      cat("\033[32mPerforming analysis for column:\033[0m ", i, "\n")
      # Perform MAXSTAT test
      merged_data$time <- merged_data[[time]]
      merged_data$variable_01 <- merged_data[[variable_01]]
      merged_data$variable_01<- as.numeric(merged_data$variable_01)
      gene_column <- merged_data[[i]]    #get(i, merged_data)
      MAXSTAT <- maxstat::maxstat.test(survival::Surv(time, variable_01) ~ gene_column, data = merged_data,
                                       smethod = "LogRank", pmethod = "Lau92", iscores = TRUE, minprop = 0.1, maxprop = 0.9)
      cut.off <- MAXSTAT$estimate
      cat("\033[32mCUT OFF\033[0m\n")
      print(cut.off)

      # Create a new variable based on the cutoff
      new_column_name <- paste0(i, "_mRNA_expression")
      merged_data[[new_column_name]] <- ifelse(gene_column > cut.off, "High", "Low")
      merged_data[[new_column_name]] <- factor(merged_data[[new_column_name]])

      # Fit survival model
      cat("\033[32mFitting survival model\033[0m\n")
      surv_object <- survival::Surv(merged_data$time, merged_data$variable_01)
      surv_formula <- as.formula(paste("surv_object ~", new_column_name))

      fit1 <- survival::survfit(surv_formula, data = merged_data)

      # Summary of the fit
      cat("\033[32mSummary of the fit\033[0m\n")
      print(summary(fit1))
      fit1_df <- data.frame(
        time = fit1$time,
        n_risk = fit1$n.risk,
        n_event = fit1$n.event,
        n_censor = fit1$n.censor,
        survival = fit1$surv,
        std_err = fit1$std.err,
        cumhaz = fit1$cumhaz,
        std_chaz = fit1$std.chaz,
        lower = fit1$lower,
        upper = fit1$upper
      )
      fit1_df$strata <- rep(names(fit1$strata), times = fit1$strata)
      merged_data_subset <- data.frame(
        time = merged_data$time,
        n_event = merged_data$variable_01,
        id = merged_data$id
      )
      fit1_df <- merge(fit1_df, merged_data_subset, by = c("time", "n_event"), all.x = TRUE)
      csv_filename <- paste0("Summary_of_the_fit_", new_column_name, ".csv")
      write.csv(fit1_df, file = csv_filename, row.names = FALSE)


      # Log-rank test and p-value
      cat("\033[32mPerforming log-rank test and obtaining p-value\033[0m\n")
      surv_diff <- survival::survdiff(surv_formula, data = merged_data)
      p_value <- 1 - stats::pchisq(surv_diff$chisq, length(surv_diff$n) - 1)

      print(surv_diff)
      surv_diff_df <- data.frame(
        group = attr(surv_diff$n, "dimnames")[[1]],
        N = as.vector(surv_diff$n),
        Observed = surv_diff$obs,
        Expected = surv_diff$exp,
        `O-E^2/E` = (surv_diff$obs - surv_diff$exp)^2 / surv_diff$exp,
        `Chisq` = rep(surv_diff$chisq, 2),  # chi-cuadrado repetido dos veces
        p_value = rep(surv_diff$pvalue, 2)    # p-valor repetido dos veces
      )
      csv_filename <- paste0("surv_diff_summary_", new_column_name, ".csv")
      write.csv(surv_diff_df, file = csv_filename, row.names = FALSE)
      cat("\033[32mP-value\033[0m\n")
      print(p_value)

      # Generate Kaplan-Meier plot
      cat("\033[32mGenerating Kaplan-Meier plot\033[0m\n")
      palette <- c("#9A3449", "#D4A8B1")
      plot(fit1, lty = 1, col = palette, lwd = 4, main = paste("Survival analysis for", i, "\n", "p-value =", format(p_value, digits = 3)))

      # Add a legend
      legend("topright",
             legend = c("High", "Low"),
             lty = 1,
             col = palette,
             lwd = 4)
      cat("\033[32mPlots saved in survival_analysis_plots_top10.pdf\033[0m\n")
    }
    dev.off()


    } else {
      cat("\033[32mSkipping survival analysis.\033[0m\n")
    }
  }
