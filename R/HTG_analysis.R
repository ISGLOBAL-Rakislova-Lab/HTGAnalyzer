#' HTG_analysis: Perform DESeq2 Analysis, GSEA, TME, and Survival Analysis
#'
#' @description This function conducts a comprehensive analysis pipeline including DESeq2 differential expression analysis (DEA), Gene Set Enrichment Analysis (GSEA), TME analysis, and survival analysis. The pipeline supports optional steps for generating volcano plots and heatmaps. The function is suitable for both HTG and RNA-seq data.
#'
#' @param outliers A character vector specifying the IDs of outlier samples to be removed. Outliers can also be identified using the HTG_QC function.
#' @param pattern (Optional) A regular expression pattern to identify control probes in the count data. For HTG, this could be "^NC-|^POS-|^GDNA-|^ERCC-". If NULL, the pattern will not be applied.
#' @param counts_data A matrix or data frame containing count data. It is recommended that the count data does not include control probes or outliers if not needed.
#' @param col_data A data frame containing sample annotations. For survival analysis, it must include variables `time` and `variable_01`.
#' @param design_formula The design formula for DESeq2 analysis, specified as a string without the tilde (~).
#' @param threshold_gene Minimum count threshold per gene. Default is 200.
#' @param threshold_subject Minimum count threshold per subject. Default is 10.
#' @param heatmap_columns A character vector specifying the columns to be used for annotations in the heatmap. Default is c("design_formula", "Smoker").
#' @param contrast A character vector specifying the contrast for differential expression analysis. Default is c('column_name', 'variable1', 'variable2').
#' @param pCutoff The p-value cutoff for generating the volcano plot. Default is 0.05.
#' @param generate_volcano A logical value indicating whether to generate a volcano plot. Default is TRUE.
#' @param remove_outliers A logical value indicating whether to remove outliers. Default is TRUE.
#' @param GSEA A logical value indicating whether to perform GSEA analysis. Default is FALSE.
#' @param generate_heatmap A logical value indicating whether to generate a heatmap. Default is TRUE.
#' @param TME A logical value indicating whether to perform TME analysis. Default is TRUE.
#' @param survival_analysis A logical value indicating whether to perform survival analysis. Default is FALSE.
#' @param percentage_gene A numeric value between 0 and 1 indicating the minimum fraction of samples in which a gene must be expressed to be retained. Default is 0.2.
#' @param percentage_zero A numeric value between 0 and 1 indicating the maximum fraction of samples in which a gene can be zero to be retained. Default is 0.2.
#' @param top_genes A character vector specifying top genes for analysis. Default is c("CCND1", "MMP10", "CTTN").
#' @param DEA A logical value indicating whether to perform DESeq2 analysis with filtering and without Lfc shrinkage. Default is TRUE.
#'
#' @return Returns an object with the results of the specified contrast and saves an Excel file with the results, along with PDF files of the generated plots.
#'
#' @export
#'
#' @examples
#'ALL_done<- HTG_analysis(outliers = outliers, pattern = "^NC-|^POS-|^GDNA-|^ERCC-", counts_data, col_data = AnnotData,
#'design_formula = "Ciclina2" , percentage_gene = 0.2, percentage_zero = 0.2,
#'threshold_gene = 200, threshold_subject = 10, top_genes = c("CCND1", "MMP10", "CTTN"), heatmap_columns = c("Ciclina2", "Smoker"),
#'contrast = c("Ciclina2", "high", "low"), pCutoff = 5e-2,variable_01 = "smoke_01", time = "time",
#'DEA = TRUE, generate_volcano = TRUE, remove_outliers = TRUE, GSEA = TRUE, generate_heatmap = TRUE, TME = TRUE,
#'survival_analysis = TRUE)
#'
#'
#' @name HTG_analysis


HTG_analysis <- function(outliers, pattern = NULL, counts_data, col_data, design_formula = NULL , percentage_gene = 0.2, percentage_zero = 0.2,
                        threshold_gene = 200, threshold_subject = 10, top_genes = c("CCND1", "MMP10", "CTTN"), heatmap_columns = NULL,
                        contrast = NULL, pCutoff = 5e-2,variable_01 = NULL, time = NULL,
                        DEA = TRUE, generate_volcano = TRUE, remove_outliers = TRUE, GSEA = FALSE, generate_heatmap = TRUE, TME = TRUE,
                        survival_analysis = FALSE) {


  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(PoiClaClu)
  library(RColorBrewer)
  library(EnhancedVolcano)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(msigdbr)
  library(fgsea)
  library(DOSE)
  library(enrichplot)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(ggupset)
  library(grid)
  library(survival)
  library(survminer)
  library(dplyr)
  library(IOBR)
  library(immunedeconv)
  library(tidyverse)
  library(ggpubr)
  library(tidyr)
  library(maxstat)

  if (!is.null(pattern)) {
    cat("\033[33mFILTERING THE COUNT DATA. DELETING THE PROVES.\033[0m\n")
    counts_data <- subset(counts_data, !grepl(pattern, rownames(counts_data)))
  }

  if (remove_outliers) {
    cat("\033[33mREMOVING OUTLIERS\033[0m\n")
    if (!is.null(pattern)) {
      # Remove outliers based on pattern
      filtered <- base::subset(counts_data, !grepl(pattern, rownames(counts_data)))
    } else {
      # If no pattern is provided, do not filter based on pattern
      filtered <- counts_data
    }
    # Remove columns corresponding to outliers
    counts_filtered <- filtered[, !colnames(filtered) %in% outliers]
    AnnotData <- col_data[!col_data[["id"]] %in% outliers, ]
  } else {
    counts_filtered <- counts_data
    AnnotData <- col_data
  }
  if (DEA) {
    cat("\033[33mSTARTING THE DIFERENTIAN EXPRESSION ANALYSIS.\033[0m\n")
    if (is.null(contrast)) {
      stop("Contrast is required for DESeq2 analysis. Remember structure: contrast = c('column_name', 'variable1','variable2')")
    }
    if (is.null(design_formula)) {
      stop("Design formula is required for DESeq2 analysis.")
    }
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

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts_filtered, colData = col_data, design = design_formul)
  cat("\033[32m\033[0m\n")
  cat("\033[32mBEFORE FILTERING\033[0m\n")
  print(dds)
  cat("\033[32m\033[0m\n")
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
      stop("heatmap_columns are required for genereting the heatmap.")
    }
    select <- order(rowMeans(DESeq2::counts(dds, normalized = TRUE)), decreasing = TRUE)[1:500]
    df <- as.data.frame(colData(dds)[, heatmap_columns])
    pheatmap(assay(vsd)[select,], cluster_rows = FALSE, show_rownames = FALSE,
             cluster_cols = TRUE, annotation_col = df)
  }

  # Results contrast
  cat("\033[33mGENERATING CONTRAST RESULTS\033[0m\n")
  res <- results(dds, contrast = contrast, cooksCutoff = TRUE)
  cat("\033[32mSUMMARY OF RESULT OF THE CONTRAST\033[0m\n")
  print(summary(res))
  cat("\033[32m \033[0m\n")
  cat("\033[32mRESULTS OF THE CONTRAST (TOP 10 p-adj)\033[0m\n")
  print(head(res[order(res$padj), ], 10))
  cat("\033[32m \033[0m\n")
  cat("\033[32mRESULTS of DEA WILL BE SAVED IN .CSV IN YOUR CURRENT DIRECTORY\033[0m\n")
  write.csv(res, "results_DEA.csv", row.names = TRUE)
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
  cat("\033[33mGENERATING POISSON DISTANCES PLOT\033[0m\n")

  pheatmap(samplePoisDistMatrix,
           clustering_distance_rows = poisd$dd,
           clustering_distance_cols = poisd$dd,
           col = colors,
           main = "Poisson distances",
           width = 800,
           height = 600)

  sampleDists <- dist(t(assay(vsd)))

  options(repr.plot.height = 7, repr.plot.width = 10)
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsd$SampleID, sep = " - ")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors,
           fontsize = 8)
  cat("\033[33mGENERATING VSD BOXPLOT\033[0m\n")
  boxplot(assay(vsd), las = 2, main = "vsd", cex.axis = 0.6)

  # COOK's DISTANCE
  cat("\033[33mGENERATING COOKS DISTANCE PLOT\033[0m\n")
  boxplot(log10(assays(dds)[["cooks"]]), range = 0, las = 2, cex.axis = 0.9, main = "COOK'S DISTANCE")

  plotMA(res, main = "PLOT MA OF RESULTS")
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
         xaxt = "n", pch = 19)
    axis(1, at = 1:nrow(gen_a2m_ordered), labels = rownames(gen_a2m_ordered), las = 2, cex.axis = 0.6)
    legend("topright", legend = levels_design_formula, fill = colors)
  }

  if (generate_volcano) {
    cat("\033[33mGENERATING VOLCANO PLOT\033[0m\n")
    EnhancedVolcano(res,
                    lab = rownames(res),
                    x = 'log2FoldChange',
                    y = 'padj',
                    pCutoff = pCutoff)
  }

  plotDispEsts(dds, main = "DIPERSION PLOT")
  } else {cat("\033[32mSkipping Diferential expresion analysis\033[0m\n")}

  if (GSEA) {
    cat("\033[33mSTARTING GSEA\033[0m\n")
    suppressMessages(library(clusterProfiler))
    suppressMessages(library(dplyr))
    suppressMessages(library(msigdbr))
    suppressMessages(library(enrichplot))
    suppressMessages(library(org.Hs.eg.db))
    suppressMessages(library(fgsea))
    suppressMessages(library(DOSE))
    suppressMessages(library(ggplot2))
    suppressMessages(library(ggupset))
    suppressMessages(library(grid))

    cat("\033[32mPerforming gseGO analysis\033[0m\n")

    # Prepare gene list for gseGO
    cat("\033[32mPreparing gene list for GSEA\033[0m\n")
    original_gene_list <- res$log2FoldChange
    names(original_gene_list) <- rownames(res)
    gene_list <- na.omit(original_gene_list)
    gene_list <- sort(gene_list, decreasing = TRUE)

    # gseGO Analysis
    cat("\033[32mPerforming gseGO Analysis\033[0m\n")
    gse2 <- gseGO(geneList = gene_list, ont = "BP", keyType = "SYMBOL", nPermSimple = 500000,
                  minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, eps = 0,
                  OrgDb = "org.Hs.eg.db", pAdjustMethod = "bonferroni")

    # Create and print plots for gseGO
    cat("\033[32mCreating Plots for gseGO\033[0m\n")
    dotplot1 <- dotplot(gse2, showCategory = 10, split = ".sign", font.size = 9, label_format = 40,
                        title = "gseGO Enrichment Results: Pathways", color = "p.adjust", size = "Count")
    dotplot2 <- dotplot(gse2, showCategory = 10, split = ".sign", font.size = 9, label_format = 40,
                        title = "gseGO Enrichment Results: Pathways", color = "p.adjust", size = "Count") + facet_grid(.~.sign)
    x2 <- pairwise_termsim(gse2)
    emapplot1 <- emapplot(x2, max.overlaps = 70, min.segment.length = 0.3, point_size = 0.3, font.size = 5) + ggtitle("Enrichment Map gseGO")
    ridgeplot1 <- ridgeplot(gse2) + labs(x = "gseGO enrichment distribution", font.size = 7) + theme(axis.text.y = element_text(size = 9))
    heatplot1 <- heatplot(gse2, showCategory = 10) + ggtitle("gseGO Heatplot")
    treeplot1 <- suppressWarnings(treeplot(x2)) + ggtitle("gseGO Treeplot")
    a <- gseaplot2(gse2, geneSetID = 1, title = paste("GSEA Plot:", gse2$Description[1]))
    b <- gseaplot2(gse2, geneSetID = 1:5, pvalue_table = TRUE, title = "GSEA: Top 5 Gene Sets")

    # KEGG Analysis
    cat("\033[32mPerforming KEGG Analysis\033[0m\n")
    ids <- bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    dedup_ids <- ids[!duplicated(ids[c("SYMBOL")]), ]
    df2 <- res[rownames(res) %in% dedup_ids$SYMBOL, ]
    df2$Y <- dedup_ids$ENTREZID
    kegg_gene_list <- df2$log2FoldChange
    names(kegg_gene_list) <- df2$Y
    kegg_gene_list <- na.omit(kegg_gene_list)
    kegg_gene_list <- sort(kegg_gene_list, decreasing = TRUE)

    kk2 <- gseKEGG(geneList = kegg_gene_list, organism = "hsa", minGSSize = 3, maxGSSize = 800,
                   pvalueCutoff = 0.05, pAdjustMethod = "none", keyType = "ncbi-geneid", nPermSimple = 100000)

    # Create and print plots for KEGG
    cat("\033[32mCreating Plots for KEGG\033[0m\n")
    dotplot3 <- dotplot(kk2, showCategory = 10, title = "Enriched Pathways for KEGG", split = ".sign", font.size = 9) + facet_grid(.~.sign)
    x3 <- pairwise_termsim(kk2)
    emapplot2 <- emapplot(x3, font.size = 8) + ggtitle("KEGG Enrichment Map")
    ridgeplot2 <- ridgeplot(kk2) + labs(x = "KEGG enrichment distribution", font.size = 6) + theme(axis.text.y = element_text(size = 9))
    heatplot2 <- heatplot(kk2, showCategory = 10) + ggtitle("KEGG Heatplot")
    treeplot1 <- suppressWarnings(treeplot(x3)) + ggtitle("KEGG Treeplot")
    upset_plot <- upsetplot(kk2) + labs(title = "UpSet Plot for KEGG")

    # enrichGO Analysis
    cat("\033[32mPerforming GO Enrichment Analysis\033[0m\n")
    sig_genes_df <- subset(res, padj < 0.05)

    if (nrow(sig_genes_df) > 0) {
      genes <- sig_genes_df$log2FoldChange
      names(genes) <- rownames(sig_genes_df)

      go_enrich <- enrichGO(
        gene = names(genes),
        universe = names(gene_list),
        OrgDb = org.Hs.eg.db,
        keyType = 'SYMBOL',
        readable = TRUE,
        ont = "BP",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.10
      )

      go_results <- go_enrich@result
      significant_terms <- go_results[go_results$qvalue < 0.05, ]
      significant_terms <- significant_terms[order(significant_terms$qvalue), ]

      if (nrow(significant_terms) > 0) {
        bar_plot <- ggplot(significant_terms, aes(x = reorder(Description, -Count), y = Count)) +
          geom_bar(stat = "identity", fill = "skyblue") +
          labs(title = "Significant GO Terms", x = "GO Term", y = "Count") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
        print(bar_plot)
      } else {
        cat("\033[31mNo significant GO terms found to plot.\033[0m\n")
      }

      # Save GO enrichment results to CSV
      write.csv(go_results, "go_results.csv")
      cat("\033[32mGO enrichment results saved to 'go_results.csv'\033[0m\n")
    } else {
      cat("\033[31mNo significant genes found for GO enrichment analysis.\033[0m\n")
    }

    # Save plots to PDF
    cat("\033[33mGENERATING GSEA PLOTS ANALYSIS\033[0m\n")
    pdf("GSEA_analysis_plots.pdf", width = 11, height = 14)
    print(dotplot1)
    print(dotplot2)
    print(emapplot1)
    print(ridgeplot1)
    print(heatplot1)
    print(treeplot1)
    print(a)
    print(b)
    print(dotplot3)
    print(emapplot2)
    print(ridgeplot2)
    print(heatplot2)
    print(treeplot2)
    print(upset_plot)
    if (nrow(sig_genes_df) > 0 && nrow(significant_terms) > 0) {
      print(bar_plot)
    }
    dev.off()

    # Save tables
    cat("\033[33mGENERATING .CSV OF GSEA RESULTS033[0m\n")
    write.csv(gene_list, "gene_list.csv")
    write.csv(kegg_gene_list, "kegg_gene_list.csv")
  }


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
    suppressWarnings({
      tpm_counts <- count2tpm(counts_data,
                              idType = "Symbol",
                              org = "hsa",
                              source = "biomart")
      write.csv(tpm_counts, "tpm_counts.csv")
    })

    # Se almacenan los genes omitidos en la normalización TPM
    genes_omitidos <- setdiff(rownames(counts_data), rownames(tpm_counts))
    print(paste0("Number of genes omitted during TPM normalization due to their length not being available in Biomart:   ", dim(counts_data)[1]-dim(tpm_counts)[1]))
    tpm_counts <- as.data.frame(tpm_counts)

    if (!is.null(dds)) {
      cat("\033[32mWe are going to use information from DEA\033[0m\n")
      dds <- estimateSizeFactors(dds)
      normalized_counts <- counts(dds, normalized = TRUE)
    } else {
      cat("Performing normalization.\n")
      design_formul <- as.formula("~ 1")
      colnames(AnnotData) <- gsub(" ", "_", colnames(AnnotData))
      AnnotData <- AnnotData[order(AnnotData$id), ]
      counts_data <- counts_data[, order(colnames(counts_data))]
      if (!identical(colnames(counts_data), AnnotData$id)) {
        stop("Column names of counts_data and IDs in AnnotData do not match.")
      }
    }

    ## Deconvolución
    imm_epic <- deconvolute(tpm_counts, method = "epic")
    imm_qti <- deconvolute(tpm_counts, method = "quantiseq")
    imm_xcell <- deconvolute(tpm_counts, method = "xcell")
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
    design_formula_sym <- sym(design_formula)
    # Vector with column names except the last one
    column_names <- names(imm_epic)[-ncol(imm_epic)]
    # Loop to generate plots for each column
    for (col_name in column_names) {
      # Calculate means by design formula
      means_df <- imm_epic %>%
        group_by(!!design_formula_sym) %>%
        summarize(mean_value = mean(.data[[col_name]] * 100, na.rm = TRUE), .groups = "drop")

      # Perform t-test or ANOVA
      formula <- as.formula(paste0("`", col_name, "` ~ ", as_label(design_formula_sym)))
      p_value <- tryCatch({
        if (n_distinct(imm_epic[[as_label(design_formula_sym)]]) == 2) {
          t.test(formula, data = imm_epic)$p.value
        } else {
          summary(aov(formula, data = imm_epic))[[1]]$`Pr(>F)`[1]
        }
      }, error = function(e) {
        NA  # In case of error, return NA for p-value
      })

      # Generate plot
      plot <- ggplot(imm_epic, aes(x = !!design_formula_sym , y = .data[[col_name]] * 100,
                                   fill = !!design_formula_sym, color = !!design_formula_sym)) +
        geom_jitter(alpha = 1, width = 0.3, height = 0) +
        geom_boxplot(fill = "white", alpha = 0.5, outlier.alpha = 1) +
        geom_point(data = means_df, aes(x = !!design_formula_sym, y = mean_value, fill = !!design_formula_sym),
                   shape = 22, color = "black", size = 3, stroke = 1.5,
                   show.legend = F) +
        labs(x = NULL, y = "Abundance (%)", title = col_name, subtitle = "imm_epic") +
        scale_y_continuous(labels = scales::percent_format(scale = 1),
                           limits = c(0, NA)) +
        theme(axis.text.x = element_blank()) +
        annotate("text", x = Inf, y = Inf, label = paste("p-value:", format(p_value, digits = 3)),
                 hjust = 1.1, vjust = 1.1, size = 5, color = "red")

      print(plot)
    }
    pdf("plots_imm_EPIC.pdf")
    for (col_name in column_names) {
      # Calculate means by design formula
      means_df <- imm_epic %>%
        group_by(!!design_formula_sym) %>%
        summarize(mean_value = mean(.data[[col_name]] * 100, na.rm = TRUE), .groups = "drop")

      # Perform t-test or ANOVA
      formula <- as.formula(paste0("`", col_name, "` ~ ", as_label(design_formula_sym)))
      p_value <- tryCatch({
        if (n_distinct(imm_epic[[as_label(design_formula_sym)]]) == 2) {
          t.test(formula, data = imm_epic)$p.value
        } else {
          summary(aov(formula, data = imm_epic))[[1]]$`Pr(>F)`[1]
        }
      }, error = function(e) {
        NA  # In case of error, return NA for p-value
      })

      # Generate plot
      plot <- ggplot(imm_epic, aes(x = !!design_formula_sym , y = .data[[col_name]] * 100,
                                   fill = !!design_formula_sym, color = !!design_formula_sym)) +
        geom_jitter(alpha = 1, width = 0.3, height = 0) +
        geom_boxplot(fill = "white", alpha = 0.5, outlier.alpha = 1) +
        geom_point(data = means_df, aes(x = !!design_formula_sym, y = mean_value, fill = !!design_formula_sym),
                   shape = 22, color = "black", size = 3, stroke = 1.5,
                   show.legend = F) +
        labs(x = NULL, y = "Abundance (%)", title = col_name, subtitle = "imm_epic") +
        scale_y_continuous(labels = scales::percent_format(scale = 1),
                           limits = c(0, NA)) +
        theme(axis.text.x = element_blank()) +
        annotate("text", x = Inf, y = Inf, label = paste("p-value:", format(p_value, digits = 3)),
                 hjust = 1.1, vjust = 1.1, size = 5, color = "red")

      print(plot)
    }
    dev.off()
    while (!is.null(dev.list())) dev.off()


    ############# imm_qti
    # Vector with column names except the last one
    column_names <- names(imm_qti)[-ncol(imm_qti)]
    # Loop to generate plots for each column
    for (col_name in column_names) {
      # Calculate means by design formula
      means_df <- imm_qti %>%
        group_by(!!design_formula_sym) %>%
        summarize(mean_value = mean(.data[[col_name]] * 100, na.rm = TRUE), .groups = "drop")

      # Perform t-test or ANOVA
      formula <- as.formula(paste0("`", col_name, "` ~ ", as_label(design_formula_sym)))
      p_value <- tryCatch({
        if (n_distinct(imm_qti[[as_label(design_formula_sym)]]) == 2) {
          t.test(formula, data = imm_qti)$p.value
        } else {
          summary(aov(formula, data = imm_qti))[[1]]$`Pr(>F)`[1]
        }
      }, error = function(e) {
        NA  # In case of error, return NA for p-value
      })

      # Generate plot
      plot <- ggplot(imm_qti, aes(x = !!design_formula_sym , y = .data[[col_name]] * 100,
                                  fill = !!design_formula_sym, color = !!design_formula_sym)) +
        geom_jitter(alpha = 1, width = 0.3, height = 0) +
        geom_boxplot(fill = "white", alpha = 0.5, outlier.alpha = 1) +
        geom_point(data = means_df, aes(x = !!design_formula_sym, y = mean_value, fill = !!design_formula_sym),
                   shape = 22, color = "black", size = 3, stroke = 1.5,
                   show.legend = F) +
        labs(x = NULL, y = "Abundance (%)", title = col_name, subtitle = "imm_qti") +
        scale_y_continuous(labels = scales::percent_format(scale = 1),
                           limits = c(0, NA)) +
        theme(axis.text.x = element_blank()) +
        annotate("text", x = Inf, y = Inf, label = paste("p-value:", format(p_value, digits = 3)),
                 hjust = 1.1, vjust = 1.1, size = 5, color = "red")

      print(plot)
    }
    pdf("plots_imm_qti.pdf")
    for (col_name in column_names) {
      # Calculate means by design formula
      means_df <- imm_qti %>%
        group_by(!!design_formula_sym) %>%
        summarize(mean_value = mean(.data[[col_name]] * 100, na.rm = TRUE), .groups = "drop")

      # Perform t-test or ANOVA
      formula <- as.formula(paste0("`", col_name, "` ~ ", as_label(design_formula_sym)))
      p_value <- tryCatch({
        if (n_distinct(imm_qti[[as_label(design_formula_sym)]]) == 2) {
          t.test(formula, data = imm_qti)$p.value
        } else {
          summary(aov(formula, data = imm_qti))[[1]]$`Pr(>F)`[1]
        }
      }, error = function(e) {
        NA  # In case of error, return NA for p-value
      })

      # Generate plot
      plot <- ggplot(imm_qti, aes(x = !!design_formula_sym , y = .data[[col_name]] * 100,
                                  fill = !!design_formula_sym, color = !!design_formula_sym)) +
        geom_jitter(alpha = 1, width = 0.3, height = 0) +
        geom_boxplot(fill = "white", alpha = 0.5, outlier.alpha = 1) +
        geom_point(data = means_df, aes(x = !!design_formula_sym, y = mean_value, fill = !!design_formula_sym),
                   shape = 22, color = "black", size = 3, stroke = 1.5,
                   show.legend = F) +
        labs(x = NULL, y = "Abundance (%)", title = col_name, subtitle = "imm_qti") +
        scale_y_continuous(labels = scales::percent_format(scale = 1),
                           limits = c(0, NA)) +
        theme(axis.text.x = element_blank()) +
        annotate("text", x = Inf, y = Inf, label = paste("p-value:", format(p_value, digits = 3)),
                 hjust = 1.1, vjust = 1.1, size = 5, color = "red")

      print(plot)
    }
    dev.off()
    while (!is.null(dev.list())) dev.off()


    ############# imm_xcell
    # Vector with column names except the last one
    column_names <- names(imm_xcell)[-ncol(imm_xcell)]
    # Loop to generate plots for each column
    for (col_name in column_names) {
      # Calculate means by design formula
      means_df <- imm_xcell %>%
        group_by(!!design_formula_sym) %>%
        summarize(mean_value = mean(.data[[col_name]] * 100, na.rm = TRUE), .groups = "drop")

      # Perform t-test or ANOVA
      formula <- as.formula(paste0("`", col_name, "` ~ ", as_label(design_formula_sym)))
      p_value <- tryCatch({
        if (n_distinct(imm_xcell[[as_label(design_formula_sym)]]) == 2) {
          t.test(formula, data = imm_xcell)$p.value
        } else {
          summary(aov(formula, data = imm_xcell))[[1]]$`Pr(>F)`[1]
        }
      }, error = function(e) {
        NA  # In case of error, return NA for p-value
      })

      # Generate plot
      plot <- ggplot(imm_xcell, aes(x = !!design_formula_sym , y = .data[[col_name]] * 100,
                                    fill = !!design_formula_sym, color = !!design_formula_sym)) +
        geom_jitter(alpha = 1, width = 0.3, height = 0) +
        geom_boxplot(fill = "white", alpha = 0.5, outlier.alpha = 1) +
        geom_point(data = means_df, aes(x = !!design_formula_sym, y = mean_value, fill = !!design_formula_sym),
                   shape = 22, color = "black", size = 3, stroke = 1.5,
                   show.legend = F) +
        labs(x = NULL, y = "Abundance (%)", title = col_name, subtitle = "imm_xcell") +
        scale_y_continuous(labels = scales::percent_format(scale = 1),
                           limits = c(0, NA)) +
        theme(axis.text.x = element_blank()) +
        annotate("text", x = Inf, y = Inf, label = paste("p-value:", format(p_value, digits = 3)),
                 hjust = 1.1, vjust = 1.1, size = 5, color = "red")

      print(plot)
    }
    pdf("plots_imm_xcell.pdf")
    for (col_name in column_names) {
      # Calculate means by design formula
      means_df <- imm_xcell %>%
        group_by(!!design_formula_sym) %>%
        summarize(mean_value = mean(.data[[col_name]] * 100, na.rm = TRUE), .groups = "drop")

      # Perform t-test or ANOVA
      formula <- as.formula(paste0("`", col_name, "` ~ ", as_label(design_formula_sym)))
      p_value <- tryCatch({
        if (n_distinct(imm_xcell[[as_label(design_formula_sym)]]) == 2) {
          t.test(formula, data = imm_xcell)$p.value
        } else {
          summary(aov(formula, data = imm_xcell))[[1]]$`Pr(>F)`[1]
        }
      }, error = function(e) {
        NA  # In case of error, return NA for p-value
      })

      # Generate plot
      plot <- ggplot(imm_xcell, aes(x = !!design_formula_sym , y = .data[[col_name]] * 100,
                                    fill = !!design_formula_sym, color = !!design_formula_sym)) +
        geom_jitter(alpha = 1, width = 0.3, height = 0) +
        geom_boxplot(fill = "white", alpha = 0.5, outlier.alpha = 1) +
        geom_point(data = means_df, aes(x = !!design_formula_sym, y = mean_value, fill = !!design_formula_sym),
                   shape = 22, color = "black", size = 3, stroke = 1.5,
                   show.legend = F) +
        labs(x = NULL, y = "Abundance (%)", title = col_name, subtitle = "imm_xcell") +
        scale_y_continuous(labels = scales::percent_format(scale = 1),
                           limits = c(0, NA)) +
        theme(axis.text.x = element_blank()) +
        annotate("text", x = Inf, y = Inf, label = paste("p-value:", format(p_value, digits = 3)),
                 hjust = 1.1, vjust = 1.1, size = 5, color = "red")

      print(plot)
    }
    dev.off()
    while (!is.null(dev.list())) dev.off()


    # Composición celular del TME i por grupo

    # Function to generate the bar plot for each sample
    plot_bar <- function(df, paleta, titulo, legend.position) {
      df <- df %>%
        rownames_to_column(var = "Sample") %>%
        pivot_longer(cols = colnames(df)[2:(ncol(df) - 1)],
                     names_to = "Cell_Type", values_to = "Value") %>%
        mutate(Sample = factor(Sample, levels = rev(unique(Sample))),
               Cell_Type = factor(Cell_Type, levels = rev(unique(Cell_Type))))

      p <- ggplot(df, aes(x = Sample, y = Value, fill = Cell_Type)) +
        geom_bar(stat = "identity") +
        labs(title = titulo,
             x = "Samples",
             y = "Cell Fraction (%)") +
        coord_flip() +
        guides(fill = guide_legend(reverse = TRUE)) +
        scale_fill_manual(values = paleta) +
        theme_minimal() +
        theme(legend.position = legend.position,
              axis.text.y = element_text(size = 5)) +
        scale_y_continuous(labels = scales::percent)

      return(p)
    }

    # Function to generate the average bar plot by group (e.g., HPV status)
    plot_bar_group <- function(df, paleta, titulo, design_formula, legend_position = "right") {
      suppressWarnings({
        design_formula_sym <- sym(design_formula)
        niveles_tipo_cel <- colnames(df)[1:(ncol(df) - 1)]

        promedios <- df %>%
          rownames_to_column(var = "Sample") %>%
          pivot_longer(cols = niveles_tipo_cel,
                       names_to = "Cell_Type", values_to = "Value") %>%
          group_by(!!design_formula_sym, Cell_Type) %>%
          summarise(Average = mean(Value, na.rm = TRUE), .groups = "drop") %>%
          ungroup() %>%
          mutate(!!design_formula_sym := factor(!!design_formula_sym, levels = rev(unique(!!design_formula_sym))),
                 Cell_Type = factor(Cell_Type, levels = rev(niveles_tipo_cel)))

        p <- ggplot(promedios, aes(x = !!design_formula_sym, y = Average, fill = Cell_Type)) +
          geom_bar(stat = "identity") +
          labs(title = titulo,
               x = design_formula,
               y = "Average Cell Fraction (%)") +
          coord_flip() +
          guides(fill = guide_legend(reverse = TRUE)) +
          scale_fill_manual(values = paleta) +
          theme_minimal() +
          theme(legend.position = legend_position,
                axis.text.y = element_text(size = 7)) +
          scale_y_continuous(labels = scales::percent)

        return(p)
      })
    }

    # Function to combine both plots into one
    plot_combined <- function(df, paleta, titulo_individual, titulo_grupo, design_formula, legend_position = "right") {
      p1 <- plot_bar(df, paleta, titulo_individual, legend_position)
      p2 <- plot_bar_group(df, paleta, titulo_grupo, design_formula, legend_position)


      combined_plot <- ggarrange(p1, p2, ncol = 1, nrow = 2, heights = c(1, 1))
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
    print(combined_plot_EPIC)

    combined_plot_quanTIseq <- plot_combined(imm_qti, paleta_qti, "quanTIseq Individual", "quanTIseq Average", design_formula, "right")
    print(combined_plot_quanTIseq)

    combined_plot_xCell <- plot_combined(imm_xcell, paleta_extendida, "xCell Individual", "xCell Average", design_formula, "right")
    print(combined_plot_xCell)

    pdf("plot_cell_fraction_Average_cell_fraction_EPIC.pdf", width = 11, height = 14)
    print(combined_plot_EPIC)
    dev.off()

    pdf("plot_cell_fraction_Average_cell_fraction_quanTIseq.pdf", width = 11, height = 14)
    print(combined_plot_quanTIseq)
    dev.off()

    pdf("plot_cell_fraction_Average_cell_fraction_xCell.pdf", width = 11, height = 14)
    print(combined_plot_xCell)
    dev.off()
    while (!is.null(dev.list())) dev.off()


    ########################
    trans_formato_largo <- function(df, design_formula_sym) {
      df_largo <- df %>%
        pivot_longer(-all_of(design_formula_sym), names_to = "Cell_Type", values_to = "Fraction")
      return(df_largo)
    }

    # Function to perform Shapiro-Wilk normality test
    prueba_norm <- function(df, design_formula_sym) {
      df_largo <- trans_formato_largo(df, design_formula_sym)

      resultados_normalidad <- df_largo %>%
        group_by(Cell_Type) %>%
        summarise(shapiro_test = list(shapiro.test(Fraction))) %>%
        mutate(p.value = map_dbl(shapiro_test, "p.value"))

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

      equality_results <- sapply(df, check_column_equal)
      unequal_columns_count <- sum(!equality_results)

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

        df_cell_type <- df_largo %>% filter(Cell_Type == cell_type)
        groups <- unique(df_cell_type[[design_formula_sym]])

        if (length(groups) == 2) {
          # Perform t-test
          t_test_result <- t.test(Fraction ~ df_cell_type[[design_formula_sym]], data = df_cell_type)
          results[[cell_type]] <- list(Test = "t-test", p.value = t_test_result$p.value)
        } else {
          # Perform ANOVA
          anova_result <- aov(Fraction ~ df_cell_type[[design_formula_sym]], data = df_cell_type)
          p_value <- summary(anova_result)[[1]]$`Pr(>F)`[1]
          results[[cell_type]] <- list(Test = "ANOVA", p.value = p_value)
        }
      }

      results_df <- bind_rows(lapply(names(results), function(cell_type) {
        result <- results[[cell_type]]
        tibble(Cell_Type = cell_type, Test = result$Test, p.value = result$p.value)
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

    h_imm_xcell <- h_imm_xcell %>%
      filter(!row.names(.) %in% xcell_row_delete)  # Eliminar filas

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

    #### EPIC
    combined_data <- remove_nan_rows(h_imm_epic)

    if (ncol(combined_data) < 3) {
      cat("Cannot generate the heatmap because the dataset has less than 3 columns after removing rows with NaN values.\n")
    } else {
      # Check if there are still any NaN values
      if (any(is.na(combined_data))) {
        cat("Cannot generate the heatmap because there are NaN values in the data.\n")
      } else {
        # Ensure col_ann_data has the same dimensions as combined_data
        col_ann_data <- AnnotData[colnames(combined_data), , drop = FALSE]

        if (ncol(combined_data) == nrow(col_ann_data)) {
          col_ann_data <- as.data.frame(col_ann_data)  # Convert to data frame
        } else {
          stop("Dimensions of col_ann_data and combined_data do not match")
        }

        # Annotations
        col_ann <- HeatmapAnnotation(
          df = col_ann_data[, design_formula, drop = FALSE],
          gp = gpar(col = "grey60", lwd = 1),
          annotation_name_gp = gpar(fontsize = 9),
          show_legend = TRUE
        )
        combined_data <- as.matrix(combined_data)

        # Generar Heatmap
        HeatmapEPIC <- Heatmap(combined_data,
                               column_title = "Heatmap EPIC",
                               cluster_rows = TRUE,
                               cluster_columns = TRUE,
                               column_names_gp = gpar(fontsize = 9),
                               row_names_gp = gpar(fontsize = 9),
                               rect_gp = gpar(col = "grey60", lwd = 1),
                               name = 'Z-score',
                               top_annotation = col_ann)

        print(HeatmapEPIC)
      }
    }
    ###### qti
    combined_data <- remove_nan_rows(h_imm_qti)

    # Verificar si hay NaN en combined_data
    if (ncol(combined_data) < 3) {
      cat("Cannot generate the heatmap because the dataset has less than 3 columns after removing rows with NaN values.\n")
    } else {
      # Check if there are still any NaN values
      if (any(is.na(combined_data))) {
        cat("Cannot generate the heatmap because there are NaN values in the data.\n")
      } else {
        # Ensure col_ann_data has the same dimensions as combined_data
        col_ann_data <- AnnotData[colnames(combined_data), , drop = FALSE]

        if (ncol(combined_data) == nrow(col_ann_data)) {
          col_ann_data <- as.data.frame(col_ann_data)  # Convert to data frame
        } else {
          stop("Dimensions of col_ann_data and combined_data do not match")
        }

        # Annotations
        col_ann <- HeatmapAnnotation(
          df = col_ann_data[, design_formula, drop = FALSE],
          gp = gpar(col = "grey60", lwd = 1),
          annotation_name_gp = gpar(fontsize = 9),
          show_legend = TRUE
        )
        combined_data <- as.matrix(combined_data)

        # Generar Heatmap
        Heatmap_qti <- Heatmap(combined_data,
                               column_title = "Heatmap qti",
                               cluster_rows = TRUE,
                               cluster_columns = TRUE,
                               column_names_gp = gpar(fontsize = 9),
                               row_names_gp = gpar(fontsize = 9),
                               rect_gp = gpar(col = "grey60", lwd = 1),
                               name = 'Z-score',
                               top_annotation = col_ann)
        print(Heatmap_qti)
      }
    }
    ####### h_imm_xcell
    combined_data <- remove_nan_rows(h_imm_xcell)

    # Verificar si hay NaN en combined_data
    if (ncol(combined_data) < 3) {
      cat("Cannot generate the heatmap because the dataset has less than 3 columns after removing rows with NaN values.\n")
    } else {
      # Check if there are still any NaN values
      if (any(is.na(combined_data))) {
        cat("Cannot generate the heatmap because there are NaN values in the data.\n")
      } else {
        # Ensure col_ann_data has the same dimensions as combined_data
        col_ann_data <- AnnotData[colnames(combined_data), , drop = FALSE]

        if (ncol(combined_data) == nrow(col_ann_data)) {
          col_ann_data <- as.data.frame(col_ann_data)  # Convert to data frame
        } else {
          stop("Dimensions of col_ann_data and combined_data do not match")
        }

        # Annotations
        col_ann <- HeatmapAnnotation(
          df = col_ann_data[, design_formula, drop = FALSE],
          gp = gpar(col = "grey60", lwd = 1),
          annotation_name_gp = gpar(fontsize = 9),
          show_legend = TRUE
        )
        combined_data <- as.matrix(combined_data)

        # Generar Heatmap
        Heatmap_xcell <- Heatmap(combined_data,
                                 column_title = "Heatmap xcell",
                                 cluster_rows = TRUE,
                                 cluster_columns = TRUE,
                                 column_names_gp = gpar(fontsize = 9),
                                 row_names_gp = gpar(fontsize = 9),
                                 rect_gp = gpar(col = "grey60", lwd = 1),
                                 name = 'Z-score',
                                 top_annotation = col_ann)
        print(Heatmap_xcell)
      }
    }

    cat("\033[32mHeatmap of qti, EPIC and xcell will be stored on plots_TME_heatmap.pdf\033[0m\n")
    pdf("plots_TME_heatmap.pdf", width = 11, height = 14)
    if (exists("Heatmap_qti")) {
      print(Heatmap_qti)
    } else {
      cat("Heatmap_qti object does not exist.\n")
    }
    if (exists("HeatmapEPIC")) {
      print(HeatmapEPIC)
    } else {
      cat("HeatmapEPIC object does not exist.\n")
    }

    if (exists("Heatmap_xcell")) {
      print(Heatmap_xcell)
    } else {
      cat("Heatmapxcell object does not exist.\n")
    }
    dev.off()
    while (!is.null(dev.list())) dev.off()

    data_table_list<- list(EPIC= imm_epic, QTI = imm_qti, XCELL= imm_xcell)

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
    design_formul <- as.formula("~ 1")
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
    if (is.null(dds)) {
      # Create DESeqDataSet object
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
    subset_data <- col_data %>% dplyr::filter(id %in% ids_data)

    df_ta <- as.data.frame(df_t)
    df_ta$id <- rownames(df_ta)

    # Selecting genes
    if (!is.null(res)) {
      cat("\033[32mSelecting TOP 10 genes with the lowest padj\033[0m\n")
      top_genes <- rownames(head(res[order(res$padj), ], 10))
      selected_df_t <- df_t[, top_genes, drop = FALSE]
      selected_df_t <- as.data.frame(selected_df_t)
      selected_df_t$id <- rownames(selected_df_t)
      col_data$id <- rownames(col_data)
      merged_data <- merge(col_data, selected_df_t, by = "id")
      merged_data[[time]] <- as.numeric(merged_data[[time]])

      cat("\033[32mStarting survival analysis.\033[0m\n")

      pdf("survival_analysis_plots.pdf")

      # Replace special characters in top_genes
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
        gene_column <- get(i, merged_data)  # Access the column dynamically
        MAXSTAT <- maxstat.test(Surv(time, variable_01) ~ gene_column, data = merged_data,
                                smethod = "LogRank", pmethod = "Lau92", iscores = TRUE, minprop = 0.45, maxprop = 0.55)
        cut.off <- MAXSTAT$estimate
        cat("\033[32mCUT OFF\033[0m\n")
        print(cut.off)

        # Create a new variable based on the cutoff
        new_column_name <- paste0(i, "_mRNA_expression")
        merged_data[[new_column_name]] <- ifelse(gene_column > cut.off, "High", "Low")
        merged_data[[new_column_name]] <- factor(merged_data[[new_column_name]])

        # Fit survival model
        cat("\033[32mFitting survival model\033[0m\n")
        surv_object <- Surv(merged_data$time, merged_data$variable_01)
        surv_formula <- as.formula(paste("surv_object ~", new_column_name))

        fit1 <- survfit(surv_formula, data = merged_data)

        # Summary of the fit
        cat("\033[32mSummary of the fit\033[0m\n")
        print(summary(fit1))

        # Log-rank test and p-value
        cat("\033[32mPerforming log-rank test and obtaining p-value\033[0m\n")
        surv_diff <- survdiff(surv_formula, data = merged_data)
        p_value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)

        print(surv_diff)
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
        cat("\033[32mPlots saved in survival_analysis_plots.pdf\033[0m\n")
      }
      dev.off()
      while (!is.null(dev.list())) dev.off()

    } else if (!is.null(genes_to_use)) {
      cat("\033[32mUsing provided genes\033[0m\n")
      top_genes <- genes_to_use
      selected_df_t <- df_t[, top_genes, drop = FALSE]
      selected_df_t <- as.data.frame(selected_df_t)
      selected_df_t$id <- rownames(selected_df_t)
      col_data$id <- rownames(col_data)
      merged_data <- merge(col_data, selected_df_t, by = "id")
      merged_data[[time]] <- as.numeric(merged_data[[time]])


      cat("\033[32mStarting survival analysis.\033[0m\n")

      pdf("survival_analysis_plots_SELECTED_GENES.pdf")

      # Replace special characters in top_genes
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
        gene_column <- get(i, merged_data)  # Access the column dynamically
        MAXSTAT <- maxstat.test(Surv(time, variable_01) ~ gene_column, data = merged_data,
                                smethod = "LogRank", pmethod = "Lau92", iscores = TRUE, minprop = 0.45, maxprop = 0.55)
        cut.off <- MAXSTAT$estimate
        cat("\033[32mCUT OFF\033[0m\n")
        print(cut.off)

        # Create a new variable based on the cutoff
        new_column_name <- paste0(i, "_mRNA_expression")
        merged_data[[new_column_name]] <- ifelse(gene_column > cut.off, "High", "Low")
        merged_data[[new_column_name]] <- factor(merged_data[[new_column_name]])

        # Fit survival model
        cat("\033[32mFitting survival model\033[0m\n")
        surv_object <- Surv(merged_data$time, merged_data$variable_01)
        surv_formula <- as.formula(paste("surv_object ~", new_column_name))

        fit1 <- survfit(surv_formula, data = merged_data)

        # Summary of the fit
        cat("\033[32mSummary of the fit\033[0m\n")
        print(summary(fit1))

        # Log-rank test and p-value
        cat("\033[32mPerforming log-rank test and obtaining p-value\033[0m\n")
        surv_diff <- survdiff(surv_formula, data = merged_data)
        p_value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)

        print(surv_diff)
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
        cat("\033[32mPlots saved in survival_analysis_plots.pdf\033[0m\n")
      }
      dev.off()
      while (!is.null(dev.list())) dev.off()

    } else if (!is.null(TME)) {
      cat("\033[32mUsing rownames from TME\033[0m\n")
      TME <- TME[, -ncol(TME)]
      colnames(TME) <- gsub(" ", "_", colnames(TME))
      col_data$id <- rownames(col_data)
      TME$id <- rownames(TME)
      merged_data <- merge(col_data, TME, by = "id")
      merged_data[[time]] <- as.numeric(merged_data[[time]])

      cat("\033[32mStarting survival analysis.\033[0m\n")
      colnames(merged_data) <- gsub(" ", "_", colnames(merged_data))
      colnames(merged_data) <- gsub("\\+", "_", colnames(merged_data))


      pdf("survival_analysis_plots_TME.pdf")

      for (i in colnames(TME)) {
        if (!is.numeric(merged_data[[i]])) {
          next
        }

        cat("\n")
        cat("\033[32mPerforming analysis for column:\033[0m ", i, "\n")

        # Perform MAXSTAT test
        merged_data$time<- merged_data[[time]]
        merged_data$variable_01<- merged_data[[variable_01]]
        MAXSTAT <- maxstat.test(Surv(time, variable_01) ~ merged_data[[i]], data = merged_data,
                                smethod = "LogRank", pmethod = "Lau92", iscores = TRUE, minprop = 0.45, maxprop = 0.55)
        cut.off <- MAXSTAT$estimate
        cat("\033[32mCUT OFF\033[0m\n")
        print(cut.off)

        # Create a new variable based on the cutoff

        merged_data[[paste0(i)]] <- ifelse(merged_data[[i]] > cut.off, "High", "Low")
        merged_data[[paste0(i)]] <- factor(merged_data[[paste0(i)]])


        # Fit survival model
        cat("\033[32mFitting survival model\033[0m\n")
        column_name <- paste0(i)
        surv_object <- Surv( merged_data$time , merged_data$variable_01)
        surv_formula <- as.formula(paste("surv_object ~", column_name))

        fit1 <- survfit(surv_formula, data = merged_data)
        # Summary of the fit

        cat("\033[32mSummary of the fit\033[0m\n")
        print(summary(fit1))

        # Log-rank test and p-value
        cat("\033[32mPerforming log-rank test and obtaining p-value\033[0m\n")
        surv_diff <- survdiff(surv_formula, data = merged_data)
        p_value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)

        print(surv_diff)
        cat("\033[32mP-value\033[0m\n")
        print(p_value)

        # Generate Kaplan-Meier plot
        cat("\033[32mGenerating Kaplan-Meier plot\033[0m\n")
        palette <- c("#9A3449", "#D4A8B1")
        plot(fit1, lty = 1, col = palette, lwd = 4, main = paste("Survival analysis for", i, "\n", "p-value =", format(p_value, digits = 3)))

        # Añadir una leyenda
        legend("topright",
               legend = c("High", "Low"),
               lty = 1,
               col = palette,
               lwd = 4)
        cat("\033[32mPlots saved in survival_analysis_plots.pdf\033[0m\n")
      }
      dev.off()
      while (!is.null(dev.list())) dev.off()


    } else {
      stop("Either 'res', 'genes_to_use', or 'TME' must be provided.")
    }

    } else {
      cat("\033[31mSkipping survival analysis.\033[0m\n")
    }
}
