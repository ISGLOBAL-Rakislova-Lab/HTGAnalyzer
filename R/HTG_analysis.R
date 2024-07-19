#' HTG_analysis: Perform DESeq2 Analysis, GSEA, Convolution, and Survival Analysis
#'
#' @description This function conducts a comprehensive analysis pipeline including DESeq2 differential expression analysis (DEA), Gene Set Enrichment Analysis (GSEA), convolution analysis, and survival analysis. The pipeline supports optional steps for generating volcano plots and heatmaps. The function is suitable for both HTG and RNA-seq data.
#'
#' @param outliers A character vector specifying the IDs of outlier samples to be removed. Outliers can also be identified using the HTG_QC function.
#' @param pattern (Optional) A regular expression pattern to identify control probes in the count data. For HTG, this could be "^NC-|^POS-|^GDNA-|^ERCC-". If NULL, the pattern will not be applied.
#' @param count_data A matrix or data frame containing count data. It is recommended that the count data does not include control probes or outliers if not needed.
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
#' @param Convolution A logical value indicating whether to perform convolution analysis. Default is TRUE.
#' @param grupos A list specifying groups to be compared. Default is c("brain", "lymph node").
#' @param survival_analysis A logical value indicating whether to perform survival analysis. Default is FALSE.
#' @param percentage_gene A numeric value between 0 and 1 indicating the minimum fraction of samples in which a gene must be expressed to be retained. Default is 0.2.
#' @param percentage_zero A numeric value between 0 and 1 indicating the maximum fraction of samples in which a gene can be zero to be retained. Default is 0.2.
#' @param top_genes A character vector specifying top genes for analysis. Default is c("CCND1", "MMP10", "CTTN").
#' @param correlation A logical value indicating whether to perform correlation analysis within the convolution analysis. Default is FALSE.
#' @param dds A logical value indicating whether to perform DESeq2 analysis with filtering and without Lfc shrinkage. Default is TRUE.
#'
#' @return Returns an object with the results of the specified contrast and saves an Excel file with the results, along with PDF files of the generated plots.
#'
#' @export
#'
#' @examples
#' # Perform HTG analysis
#' results <- HTG_analysis(outliers, count_data = counts, col_data = AnnotData, design_formula = "site")
#'
#' # Perform analysis without generating a volcano plot
#' results <- HTG_analysis(outliers, count_data = counts, col_data = AnnotData, design_formula = "site", generate_volcano = FALSE)
#'
#' # Perform analysis with additional options for contrast and heatmap columns
#' results <- HTG_analysis(outliers, count_data = counts, col_data = AnnotData, design_formula = "site",
#'                         contrast = c("Smoker", "yes", "no"), heatmap_columns = c("site", "Smoker"), grupos = c("brain", "lymph node"))
#'
#' @name HTG_analysis


HTG_analysis <- function(outliers, pattern = NULL, count_data, col_data, design_formula = NULL , percentage_gene = 0.2, percentage_zero = 0.2,
                        threshold_gene = 200, threshold_subject = 10, top_genes = c("CCND1", "MMP10", "CTTN"), heatmap_columns = NULL,
                        contrast = NULL, pCutoff = 5e-2,variable_01 = NULL, time = NULL, correlation = FALSE,
                        dds = TRUE, generate_volcano = TRUE, remove_outliers = TRUE, GSEA = FALSE, generate_heatmap = TRUE, Convolution = TRUE,
                        survival_analysis = FALSE, grupos = NULL) {
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

  if (remove_outliers) {
    if (!is.null(pattern)) {
      # Remove outliers based on pattern
      filtered <- subset(count_data, !grepl(pattern, rownames(count_data)))
    } else {
      # If no pattern is provided, do not filter based on pattern
      filtered <- count_data
    }
    # Remove columns corresponding to outliers
    counts_filtered <- filtered[, !colnames(filtered) %in% outliers]
    AnnotData <- col_data[!col_data[["id"]] %in% outliers, ]
  } else {
    counts_filtered <- count_data
    AnnotData <- col_data
  }
  if (dds) {
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
    if (is.null(heatmap_columns)) {
      stop("heatmap_columns are required for genereting the heatmap.")
    }
    select <- order(rowMeans(DESeq2::counts(dds, normalized = TRUE)), decreasing = TRUE)[1:500]
    df <- as.data.frame(colData(dds)[, heatmap_columns])
    pheatmap(assay(vsd)[select,], cluster_rows = FALSE, show_rownames = FALSE,
             cluster_cols = TRUE, annotation_col = df)
  }

  # Results contrast
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

  boxplot(assay(vsd), las = 2, main = "vsd", cex.axis = 0.6)

  # COOK's DISTANCE
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
    EnhancedVolcano(res,
                    lab = rownames(res),
                    x = 'log2FoldChange',
                    y = 'padj',
                    pCutoff = pCutoff)
  }

  plotDispEsts(dds, main = "DIPERSION PLOT")
  } else {cat("\033[32mSkipping Diferential expresion analysis\033[0m\n")}

  if (GSEA) {
    cat("\033[32mPerforming GSEA analysis\033[0m\n")

    # Function to add title and description to each plot
    add_plot <- function(plot, title, description) {
      cat("Adding plot:", title, "\n")
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(2, 1)))
      grid.text(title, vp = viewport(layout.pos.row = 1, layout.pos.col = 1), gp = gpar(fontsize = 14, fontface = "bold"))
      grid.text(description, vp = viewport(layout.pos.row = 2, layout.pos.col = 1), gp = gpar(fontsize = 12))
      print(plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
    }

    # Prepare gene list for GSEA
    cat("Preparing gene list for GSEA\n")
    original_gene_list <- res$log2FoldChange
    names(original_gene_list) <- rownames(res)
    gene_list <- na.omit(original_gene_list)
    gene_list <- sort(gene_list, decreasing = TRUE)

    # GSEA Analysis
    cat("Performing GSEA Analysis\n")
    gse2 <- gseGO(geneList = gene_list, ont = "ALL", keyType = "SYMBOL", nPermSimple = 500000,
                  minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, eps = 0,
                  OrgDb = "org.Hs.eg.db", pAdjustMethod = "bonferroni")

    # Dotplot for GSEA
    cat("Creating Dotplot for GSEA\n")
    dotplot1 <- dotplot(gse2, showCategory = 10, split = ".sign", font.size = 9, label_format = 40,
                        title = "Enrichment Results: Pathways", color = "p.adjust", size = "Count")
    add_plot(dotplot1, "GSEA Dotplot", "This plot shows the results of Gene Set Enrichment Analysis (GSEA).")

    dotplot2 <- dotplot(gse2, showCategory = 10, split = ".sign", font.size = 9, label_format = 40,
                        title = "Enrichment Results: Pathways", color = "p.adjust", size = "Count") +
      facet_grid(.~.sign)
    add_plot(dotplot2, "GSEA Dotplot with Facet", "This plot shows the results of GSEA with facets.")

    # Emaplot for GSEA
    cat("Creating Emaplot for GSEA\n")
    x2 <- pairwise_termsim(gse2)
    emapplot1 <- emapplot(x2, max.overlaps = 50, min.segment.length = 0.3, point_size = 0.5, font.size = 8)
    add_plot(emapplot1, "GSEA Emaplot", "This plot shows the enriched terms and their relationships.")

    # Ridgeplot for GSEA
    cat("Creating Ridgeplot for GSEA\n")
    ridgeplot1 <- ridgeplot(gse2) + labs(x = "enrichment distribution", font.size = 7)
    add_plot(ridgeplot1, "GSEA Ridgeplot", "This plot shows the enrichment distribution of gene sets.")

    # KEGG Analysis
    cat("Performing KEGG Analysis\n")
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

    # Dotplot for KEGG
    cat("\033[32mCreating Dotplot for KEGG\033[0m\n")
    dotplot3 <- dotplot(kk2, showCategory = 10, title = "Enriched Pathways", split = ".sign", font.size = 9) +
      facet_grid(.~.sign)
    add_plot(dotplot3, "KEGG Dotplot with Facet", "This plot shows the results of KEGG pathway enrichment analysis.")

    # Emaplot for KEGG
    cat("\033[32mCreating Emaplot for KEGG\033[0m\n")
    x3 <- pairwise_termsim(kk2)
    emapplot2 <- emapplot(x3, font.size = 8)
    add_plot(emapplot2, "KEGG Emaplot", "This plot shows the enriched KEGG pathways and their relationships.")

    # Ridgeplot for KEGG
    cat("\033[32mCreating Ridgeplot for KEGG\033[0m\n")
    ridgeplot2 <- ridgeplot(kk2) + labs(x = "enrichment distribution", font.size = 8)
    add_plot(ridgeplot2, "KEGG Ridgeplot", "This plot shows the enrichment distribution of KEGG pathways.")

    # GO Enrichment Analysis
    cat("\033[32mPerforming GO Enrichment Analysis\033[0m\n")
    sig_genes_df <- subset(res, padj < 0.05)

    if (nrow(sig_genes_df) > 0) {
      genes <- sig_genes_df$log2FoldChange
      names(genes) <- rownames(sig_genes_df)

      go_enrich <- enrichGO(gene = names(genes), universe = names(gene_list), OrgDb = org.Hs.eg.db,
                            keyType = 'SYMBOL', readable = TRUE, ont = "BP",
                            pvalueCutoff = 0.05, qvalueCutoff = 0.10)

      # Simplified bar plot for GO enrichment
      cat("\033[32mCreating Bar Plot for GO Enrichment\033[0m\n")
      go_results <- go_enrich@result
      significant_terms <- go_results[go_results$qvalue < 0.05, ]
      significant_terms <- significant_terms[order(significant_terms$qvalue), ]

      bar_plot <- ggplot(significant_terms, aes(x = reorder(Description, -Count), y = Count)) +
        geom_bar(stat = "identity", fill = "skyblue") +
        labs(title = "Significant GO Terms", x = "GO Term", y = "Count") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

      add_plot(bar_plot, "GO Enrichment Bar Plot", "This plot shows the significant GO terms enriched in the dataset.")
      write.csv(go_results, "go_results.csv")

      pdf("GSEA_analysis_plots.pdf")
      dotplot1
      add_plot(dotplot1, "GSEA Dotplot", "This plot shows the results of Gene Set Enrichment Analysis (GSEA).")
      dotplot2
      add_plot(dotplot2, "GSEA Dotplot with Facet", "This plot shows the results of GSEA with facets.")
      emapplot1
      add_plot(emapplot1, "GSEA Emaplot", "This plot shows the enriched terms and their relationships.")
      ridgeplot1
      add_plot(ridgeplot1, "GSEA Ridgeplot", "This plot shows the enrichment distribution of gene sets.")
      dotplot3
      add_plot(dotplot3, "KEGG Dotplot with Facet", "This plot shows the results of KEGG pathway enrichment analysis.")
      emapplot2
      add_plot(emapplot2, "KEGG Emaplot", "This plot shows the enriched KEGG pathways and their relationships.")
      if (nrow(sig_genes_df) > 0) {
        genes <- sig_genes_df$log2FoldChange
        names(genes) <- rownames(sig_genes_df)
        go_enrich <- enrichGO(gene = names(genes), universe = names(gene_list), OrgDb = org.Hs.eg.db,
                              keyType = 'SYMBOL', readable = TRUE, ont = "BP",
                              pvalueCutoff = 0.05, qvalueCutoff = 0.10)
        go_results <- go_enrich@result
        significant_terms <- go_results[go_results$qvalue < 0.05, ]
        significant_terms <- significant_terms[order(significant_terms$qvalue), ]
        bar_plot <- ggplot(significant_terms, aes(x = reorder(Description, -Count), y = Count)) +
          geom_bar(stat = "identity", fill = "skyblue") +
          labs(title = "Significant GO Terms", x = "GO Term", y = "Count") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
        add_plot(bar_plot, "GO Enrichment Bar Plot", "This plot shows the significant GO terms enriched in the dataset.")
        dev.off()

        # Save tables
        write.csv(gene_list, "gene_list.csv")
        write.csv(kegg_gene_list, "kegg_gene_list.csv")
      }
    }

    } else {cat("\033[31m Skipping GO Enrichment Analysis.\033[0m\n")}
  #####
  ####
  ####
  if (Convolution) {
    cat("\033[32mCOnvolution analysis skipped.\033[0m\n")
    print(paste0("Número inicial de genes: ", dim(counts_filtered)[1]))
    ## Normalización TPM
    cat("\033[32mTPM normalization performed and stored on tpm_counts.csv\033[0m\n")
    tpm_counts <- count2tpm(counts_filtered,
                            idType = "Symbol",
                            org = "hsa",
                            source = "biomart")
    write.csv(tpm_counts, "tpm_counts.csv")

    # Se almacenan los genes omitidos en la normalización TPM
    genes_omitidos <- setdiff(rownames(counts_filtered), rownames(tpm_counts))
    print(paste0("Number of genes omitted during TPM normalization due to their length not being available in Biomart:   ", dim(counts_filtered)[1]-dim(tpm_counts)[1]))
    tpm_counts <- as.data.frame(tpm_counts)

    if (exists("dds")) {
      cat("We are going to use dds from DEA.\n")
      dds <- estimateSizeFactors(dds)
      normalized_counts <- counts(dds, normalized = TRUE)
    } else {
      cat("Performing dds normalization.\n")
      design_formul <- as.formula("~ 1")
      colnames(AnnotData) <- gsub(" ", "_", colnames(AnnotData))
      col_data <- AnnotData[order(AnnotData$id), ]
      counts_filtered <- counts_filtered[, order(colnames(counts_filtered))]
      if (!identical(colnames(counts_filtered), col_data$id)) {
        stop("Column names of counts_filtered and IDs in col_data do not match.")
      }
    }

    ## Deconvolución
    imm_epic <- deconvolute(tpm_counts, method = "epic")
    imm_qti <- deconvolute(tpm_counts, method = "quantiseq")
    imm_xcell <- deconvolute(tpm_counts, method = "xcell")
    cat("\033[32mresults of the devonvolution will be stored in imm_epic.csv, imm_qti.csv and imm_xcell.csv \033[0m\n")


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

    # Se incluye la variable col_data$design_formula en los dataframes
    cat("\033[32mAre they in the same order?\033[0m\n")
    dim(col_data)
    dim(imm_epic)

    col_data_ids <- col_data$id
    imm_epic_ids <- as.character(rownames(imm_epic))
    common_ids <- intersect(col_data_ids, imm_epic_ids)
    # Filtrar los datos originales según los IDs comunes
    col_data <- col_data[col_data_ids %in% common_ids, ]
    dim(col_data)
    col_data<- as.data.frame(col_data)
    rownames(col_data)<- col_data$id
    col_data <- col_data[order(rownames(col_data)), ]
    imm_epic <- imm_epic[order(rownames(imm_epic)), ]
    imm_qti <- imm_qti[order(rownames(imm_qti)), ]
    imm_xcell <- imm_xcell[order(rownames(imm_xcell)), ]

    cat("\033[32mhave to be true.\033[0m\n")
    cat("\033[32mEPIC\033[0m\n")
    print(all(rownames(col_data)==rownames(imm_epic)))
    cat("\033[32mqti\033[0m\n")
    print(all(rownames(col_data)==rownames(imm_qti)))
    cat("\033[32mxcell\033[0m\n")
    print(all(rownames(col_data)==rownames(imm_xcell)))

    imm_epic[[design_formula]] <- factor(col_data[[design_formula]])
    imm_qti[[design_formula]] <- factor(col_data[[design_formula]])
    imm_xcell[[design_formula]] <- factor(col_data[[design_formula]])

    ############# imm_epic
    design_formula_sym <- sym(design_formula)
    # Vector with column names except the last one
    column_names <- names(imm_epic)[-ncol(imm_epic)]
    # Loop to generate plots for each column
    for (col_name in column_names) {
      # Calculate means by design_formula
      means_df <- imm_epic %>%
        group_by(!!design_formula_sym) %>%
        summarize(mean_value = mean(.data[[col_name]] * 100, na.rm = TRUE))
      # Generate plot
      print(
        ggplot(imm_epic, aes(x = !!design_formula_sym , y = .data[[col_name]] * 100,
                             fill = !!design_formula_sym, color = !!design_formula_sym)) +
          geom_jitter(alpha = 1, width = 0.3, height = 0) +
          geom_boxplot(fill = "white", alpha = 0.5, outlier.alpha = 1) +
          geom_point(data = means_df, aes(x = !!design_formula_sym, y = mean_value, fill = !!design_formula_sym),
                     shape = 22, color = "black", size = 3, stroke = 1.5,
                     show.legend = F) +
          labs(x = NULL, y = "Abundance (%)", title = col_name, subtitle = "imm_epic") +
          scale_y_continuous(labels = scales::percent_format(scale = 1),
                             limits = c(0, NA)) +
          theme(axis.text.x = element_blank())
      )
    }
    ############# imm_qti
    # Vector with column names except the last one
    column_names <- names(imm_qti)[-ncol(imm_qti)]
    # Loop to generate plots for each column
    pdf("plots_imm_qti.pdf")
    for (col_name in column_names) {
      # Calculate means by design_formula
      means_df <- imm_qti %>%
        group_by(!!design_formula_sym) %>%
        summarize(mean_value = mean(.data[[col_name]] * 100, na.rm = TRUE))
      # Generate plot
      print(
        ggplot(imm_qti, aes(x = !!design_formula_sym , y = .data[[col_name]] * 100,
                            fill = !!design_formula_sym, color = !!design_formula_sym)) +
          geom_jitter(alpha = 1, width = 0.3, height = 0) +
          geom_boxplot(fill = "white", alpha = 0.5, outlier.alpha = 1) +
          geom_point(data = means_df, aes(x = !!design_formula_sym, y = mean_value, fill = !!design_formula_sym),
                     shape = 22, color = "black", size = 3, stroke = 1.5,
                     show.legend = F) +
          labs(x = NULL, y = "Abundance (%)",  title = col_name, subtitle = "imm_qti") +
          scale_y_continuous(labels = scales::percent_format(scale = 1),
                             limits = c(0, NA)) +
          theme(axis.text.x = element_blank())
      )
    }
    dev.off()

    ############# imm_xcell
    # Vector con los nombres de las columnas excepto la última
    column_names <- names(imm_xcell)[-ncol(imm_xcell)]
    # Bucle for para generar gráficos para cada columna
    for (col_name in column_names) {
      # Calcular las medias por design_formula
      means_df <- imm_xcell %>%
        group_by(!!design_formula_sym) %>%
        summarize(mean_value = mean(.data[[col_name]] * 100, na.rm = TRUE))

      print(
        ggplot(imm_xcell, aes(x = !!design_formula_sym , y = .data[[col_name]],
                              fill = !!design_formula_sym, color = !!design_formula_sym)) +
          geom_jitter(alpha = 1, width = 0.3, height = 0) +
          geom_boxplot(fill = "white", alpha = 0.5, outlier.alpha = 1) +
          geom_point(data = means_df, aes(x = !!design_formula_sym, y = mean_value, fill = !!design_formula_sym),
                     shape = 22, color = "black", size = 3, stroke = 1.5,
                     show.legend = F) +
          labs(x = NULL, y = "Abundance (arbitrary value)", title = col_name, subtitle = "imm_xcell") +
          scale_y_continuous(expand = expansion(add = c(0, 0.1)))  +
          theme(axis.text.x = element_blank())
      )
    }


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
              axis.text.y = element_text(size = 4)) +
        scale_y_continuous(labels = scales::percent)

      return(p)
    }

    # Function to generate the average bar plot by group (e.g., HPV status)
    plot_bar_group <- function(df, paleta, titulo, design_formula, legend_position = "right") {
      design_formula_sym <- sym(design_formula)
      niveles_tipo_cel <- colnames(df)[1:(ncol(df) - 1)]

      promedios <- df %>%
        rownames_to_column(var = "Sample") %>%
        pivot_longer(cols = niveles_tipo_cel,
                     names_to = "Cell_Type", values_to = "Value") %>%
        group_by(!!design_formula_sym, Cell_Type) %>%
        summarise(Average = mean(Value, na.rm = TRUE)) %>%
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
              axis.text.y = element_text(size = 8)) +
        scale_y_continuous(labels = scales::percent)

      return(p)
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

    pdf("plot_cell_fraction_Average_cell_fraction.pdf")
    print(combined_plot_EPIC)
    print(combined_plot_quanTIseq)
    print(combined_plot_xCell)
    dev.off()

    # # Contrastes de Hipótesis
    prueba_norm <- function(df, design_formula_sym) {
      # Transform the dataframe to long format
      df_largo <- df %>%
        pivot_longer(-design_formula_sym,
                     names_to = "Celular type", values_to = "Fraction")

      # Perform Shapiro-Wilk test for each cell type
      resultados_normalidad <- df_largo %>%
        group_by(`Celular type`) %>%
        summarise(shapiro_test = list(shapiro.test(Fraction))) %>%
        mutate(p.value = map_dbl(shapiro_test, "p.value"))

      return(resultados_normalidad)
    }
    check_column_equal <- function(column) {
      all_equal <- length(unique(column)) == 1
      return(all_equal)
    }

    # Verificación de igualdad para cada columna
    equality_results_qti <- sapply(imm_qti, check_column_equal)
    unequal_columns_qti <-sum(equality_results_qti)

    equality_results_epic <- sapply(imm_epic, check_column_equal)
    unequal_columns_epic<- sum(equality_results_epic)

    equality_results_xcell <- sapply(imm_xcell, check_column_equal)
    print(equality_results_xcell)
    unequal_columns_xcell<-sum(equality_results_xcell)


    if (unequal_columns_qti == 0) {
      resultados_norm_imm_qti <- prueba_norm(imm_qti, design_formula_sym)
      cat("\033[32mTest for normality: Shapiro-Wilk test of quanTIseq\033[0m\n")
      print(resultados_norm_imm_qti)
    }else{
      cat("\033[32mCan't perform Shapiro-Wilk test for . Some columns have same value\033[0m\n")
      print(equality_results_qti)
    }


    if (unequal_columns_epic == 0) {
      resultados_norm_epic <- prueba_norm(imm_epic, design_formula_sym)
      cat("\033[32mTest for normality: Shapiro-Wilk test of EPIC\033[0m\n")
      print(resultados_norm_epic)
    }else{
      cat("\033[32mCan't perform Shapiro-Wilk test fro EPIC. Some columns have same value")
      print(equality_results_epic)
    }
    if (unequal_columns_xcell == 0) {
      resultados_norm_xcell <- prueba_norm(imm_xcell, design_formula_sym)
      cat("\033[32mTest for normality: Shapiro-Wilk test of xcell\033[0m\n")
      print(resultados_norm_xcell)
    }else{
      cat("\033[32mCan't perform Shapiro-Wilk test fro xcell. Some columns have same value\033[0m\n")
      print(equality_results_epic)
    }



    # Heatmaps
    # Function to transform dataframe to long format
    trans_formato_largo <- function(df, design_formula_sym) {
      df_largo <- df %>%
        pivot_longer(-design_formula_sym, names_to = "Celular type", values_to = "Fraction")
      return(df_largo)
    }
    # Example usage to transform dataframes to long format
    epic_largo <- trans_formato_largo(imm_epic, design_formula_sym)
    qti_largo <- trans_formato_largo(imm_qti, design_formula_sym)
    xcell_largo <- trans_formato_largo(imm_xcell, design_formula_sym)

    ## Comparación entre Grupos
    perform_parametric_test <- function(df, variable, groups) {
      resultados <- list()
      if (length(groups) != 2) {
        stop("groups should be a list with exactly two elements.")
      }
      # Iterate over column names except the last one (assuming the last one  or a similar grouping variable)
      for (col_name in colnames(df)[1:(ncol(df) - 1)]) {
        # Perform t-test for group 1
        if (!is.null(groups[[1]]) && length(groups[[1]]) > 1) {
          resultado <- t.test(df[[col_name]] ~ df[[variable]], subset = df[[variable]] %in% groups[[1]])
          resultados[[paste0(col_name, "_", groups[[1]][1], "_vs_", groups[[1]][2])]] <- c(col_name, paste0(groups[[1]][1], " vs ", groups[[1]][2]), resultado$p.value)
        } else {
          warning("Group 1 does not have enough samples.")
        }

        # Perform t-test for group 2
        if (!is.null(groups[[2]]) && length(groups[[2]]) > 1) {
          resultado <- t.test(df[[col_name]] ~ df[[variable]], subset = df[[variable]] %in% groups[[2]])
          resultados[[paste0(col_name, "_", groups[[2]][1], "_vs_", groups[[2]][2])]] <- c(col_name, paste0(groups[[2]][1], " vs ", groups[[2]][2]), resultado$p.value)
        } else {
          warning("Group 2 does not have enough samples.")
        }
      }
      resultados_df <- do.call(rbind, resultados)
      colnames(resultados_df) <- c("Tipo_Celular", "Comparacion", "P-valor")
      rownames(resultados_df) <- NULL
      resultados_df <- as.data.frame(resultados_df)
      print(resultados_df)
    }

    # Heatmaps
    # Crear dataframe para heatmap
    h_imm_epic <- imm_epic
    # Transponer y estandarizar por filas
    h_imm_epic <- as.data.frame(t(h_imm_epic))
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

    # Crear dataframe para heatmap
    h_imm_qti <- imm_qti
    # Transponer y estandarizar por filas
    h_imm_qti <- as.data.frame(t(h_imm_qti))
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
    # Renombrar nombres de poblaciones celulares antes de transponer
    h_imm_xcell <- h_imm_xcell %>%
      rename('Myeloid DC activated' = 'Myeloid dendritic cell activated') %>%
      rename('T cell CD4+' = 'T cell CD4+ (non-regulatory)') %>%
      rename('Myeloid DC' = 'Myeloid dendritic cell') %>%
      rename('CAFs' = 'Cancer associated fibroblast') %>%
      rename('Plasmacytoid DC' = 'Plasmacytoid dendritic cell') %>%
      rename('T cell regulatory' = 'T cell regulatory (Tregs)')

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

    # Factorización de variables y renombrado en col_data
    col_data[[design_formula]] <- as.factor(col_data[[design_formula]])

    ###################################
    #### EPIC
    combined_data <- h_imm_epic
    # Verificar si hay NaN en combined_data
    if (any(is.na(combined_data))) {
      cat("Cannot generate the EPIc heatmap because there are NaN values in the data.\n")
    } else {
      # Continuar con la generación del heatmap
      combined_data <- h_imm_epic
      col_ann_data <- col_data[colnames(combined_data), , drop = FALSE]

      # Asegúrate de que col_ann_data tenga las mismas dimensiones que combined_data
      if (ncol(combined_data) == nrow(col_ann_data)) {
        col_ann_data <- as.matrix(col_ann_data)  # Convertir a matriz
      } else {
        stop("Dimensions of col_ann_data and combined_data do not match")
      }

      # Anotaciones
      col_ann <- HeatmapAnnotation(
        df = col_ann_data[, design_formula, drop = FALSE],  # Selecciona solo la columna design_formula
        gp = gpar(col = "grey60", lwd = 1),
        annotation_name_gp = gpar(fontsize = 10),
        show_legend = TRUE
      )

      # Generar Heatmap
      HeatmapEPIC <- Heatmap(combined_data,
                             column_title = "Heatmap EPIC",
                             cluster_rows = TRUE,
                             cluster_columns = TRUE,
                             column_names_gp = gpar(fontsize = 6),
                             row_names_gp = gpar(fontsize = 8),
                             rect_gp = gpar(col = "grey60", lwd = 1),
                             name = 'Z-score',
                             top_annotation = col_ann)

      print(HeatmapEPIC)
    }
    ###### qti
    combined_data <- h_imm_qti
    # Verificar si hay NaN en combined_data
    if (any(is.na(combined_data))) {
      cat("Cannot generate the qti heatmap because there are NaN values in the data.\n")
    } else {
      # Continuar con la generación del heatmap
      combined_data <- h_imm_qti
      col_ann_data <- col_data[colnames(combined_data), , drop = FALSE]

      # Asegúrate de que col_ann_data tenga las mismas dimensiones que combined_data
      if (ncol(combined_data) == nrow(col_ann_data)) {
        col_ann_data <- as.matrix(col_ann_data)  # Convertir a matriz
      } else {
        stop("Dimensions of col_ann_data and combined_data do not match")
      }

      # Anotaciones
      col_ann <- HeatmapAnnotation(
        df = col_ann_data[, design_formula, drop = FALSE],  # Selecciona solo la columna design_formula
        gp = gpar(col = "grey60", lwd = 1),
        annotation_name_gp = gpar(fontsize = 10),
        show_legend = TRUE
      )

      # Generar Heatmap
      Heatmap_qti <- Heatmap(combined_data,
                             column_title = "Heatmap qti",
                             cluster_rows = TRUE,
                             cluster_columns = TRUE,
                             column_names_gp = gpar(fontsize = 6),
                             row_names_gp = gpar(fontsize = 8),
                             rect_gp = gpar(col = "grey60", lwd = 1),
                             name = 'Z-score',
                             top_annotation = col_ann)
      print(Heatmap_qti)
    }
    ####### h_imm_xcell
    combined_data <- h_imm_xcell
    # Verificar si hay NaN en combined_data
    if (any(is.na(combined_data))) {
      cat("Cannot generate the xcell heatmap because there are NaN values in the data.\n")
    } else {
      # Continuar con la generación del heatmap
      combined_data <- h_imm_xcell
      col_ann_data <- col_data[colnames(combined_data), , drop = FALSE]

      # Asegúrate de que col_ann_data tenga las mismas dimensiones que combined_data
      if (ncol(combined_data) == nrow(col_ann_data)) {
        col_ann_data <- as.matrix(col_ann_data)  # Convertir a matriz
      } else {
        stop("Dimensions of col_ann_data and combined_data do not match")
      }

      # Anotaciones
      col_ann <- HeatmapAnnotation(
        df = col_ann_data[, design_formula, drop = FALSE],
        gp = gpar(col = "grey60", lwd = 1),
        annotation_name_gp = gpar(fontsize = 10),
        show_legend = TRUE
      )

      # Generar Heatmap
      Heatmap_xcell <- Heatmap(combined_data,
                             column_title = "Heatmap xcell",
                             cluster_rows = TRUE,
                             cluster_columns = TRUE,
                             column_names_gp = gpar(fontsize = 6),
                             row_names_gp = gpar(fontsize = 8),
                             rect_gp = gpar(col = "grey60", lwd = 1),
                             name = 'Z-score',
                             top_annotation = col_ann)
      print(Heatmap_qti)
    }

    cat("\033[32mHeatmap of qti, EPIC and xcell will be stored on plots_convolution_heatmap.pdf\033[0m\n")
    pdf("plots_convolution_heatmap.pdf")
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

    if (exists("Heatmapxcell")) {
      print(Heatmapxcell)
    } else {
      cat("Heatmapxcell object does not exist.\n")
    }
    dev.off()

    # Correlación genes \~ tipos celulares
    ## EPIC
    c_imm_epic <- imm_epic
    c_imm_epic <- as.data.frame(t(c_imm_epic))
    c_imm_epic <- head(c_imm_epic, -1)
    rn_cimmepic <- rownames(c_imm_epic)
    rn_cimmepic[2] <- "CAFs"
    cl_cimmepic <- colnames(c_imm_epic)
    c_imm_epic <- apply(c_imm_epic, 2, as.numeric)
    rownames(c_imm_epic) <- rn_cimmepic
    colnames(c_imm_epic) <- cl_cimmepic
    c_imm_epic <- as.data.frame(c_imm_epic)

    ## Quantiseq
    c_imm_qti <- imm_qti
    c_imm_qti <- as.data.frame(t(c_imm_qti))
    c_imm_qti <- head(c_imm_qti, -1)
    rn_cimmqti <- rownames(c_imm_qti)
    rn_cimmqti[7] <- "T cell CD4+"
    rn_cimmqti[9] <- "T cell regulatory"
    cl_cimmqti <- colnames(c_imm_qti)
    c_imm_qti <- apply(c_imm_qti, 2, as.numeric)
    rownames(c_imm_qti) <- rn_cimmqti
    colnames(c_imm_qti) <- cl_cimmqti
    c_imm_qti <- as.data.frame(c_imm_qti)

    ## xCell
    c_imm_xcell <- as.data.frame(imm_xcell)

    # Renombrar nombres de poblaciones celulares antes de transponer
    c_imm_xcell <- c_imm_xcell %>%
      rename('Myeloid DC activated' = 'Myeloid dendritic cell activated') %>%
      rename('T cell CD4+' = 'T cell CD4+ (non-regulatory)') %>%
      rename('Myeloid DC' = 'Myeloid dendritic cell') %>%
      rename('CAFs' = 'Cancer associated fibroblast') %>%
      rename('Plasmacytoid DC' = 'Plasmacytoid dendritic cell') %>%
      rename('T cell regulatory' = 'T cell regulatory (Tregs)')

    # Transponer
    c_imm_xcell <- as.data.frame(t(c_imm_xcell))
    c_imm_xcell <- head(c_imm_xcell, -4)

    # Poblaciones celulares no interesantes
    xcell_row_delete <- c('Common lymphoid progenitor', 'Common myeloid progenitor',
                          'Granulocyte-monocyte progenitor', 'Hematopoietic stem cell')

    c_imm_xcell <- c_imm_xcell %>%
      filter(!row.names(.) %in% xcell_row_delete)  # Eliminar filas

    rn_cimmxcell <- rownames(c_imm_xcell)
    cl_cimmxcell <- colnames(c_imm_xcell)
    c_imm_xcell <- apply(c_imm_xcell, 2, as.numeric)
    rownames(c_imm_xcell) <- rn_cimmxcell
    colnames(c_imm_xcell) <- cl_cimmxcell
    c_imm_xcell <- as.data.frame(c_imm_xcell)

    decounts_correlation <- function(count_data1, deconv_data, gene_list = NULL, cell_type_list = NULL) {
      # Verificar si las columnas son numéricas
      cat("\033[32mAre numeric?\033[0m\n")
      numeric_columns_check <- sapply(deconv_data, is.numeric)
      print(numeric_columns_check)
      print(table(numeric_columns_check))

      # Verificar si hay valores NA
      cat("\033[32mAre any NA?\033[0m\n")
      na_check <- sapply(deconv_data, function(x) any(is.na(x)))
      print(na_check)
      print(table(na_check))

      if (is.null(gene_list)) {
        gene_list <- rownames(count_data1)
      }

      if (is.null(cell_type_list)) {
        cell_type_list <- rownames(deconv_data)
      }

      cat("\033[32mIniciating the correlation...\033[0m\n")
      correlation_df <- data.frame(matrix(NA, nrow = length(gene_list), ncol = length(cell_type_list), dimnames = list(gene_list, cell_type_list)), check.names = FALSE)
      pvalue_df <- data.frame(matrix(NA, nrow = length(gene_list), ncol = length(cell_type_list), dimnames = list(gene_list, cell_type_list)), check.names = FALSE)

      total_genes <- length(gene_list)
      total_cell_types <- length(cell_type_list)
      for (i in seq_along(gene_list)) {
        gene <- gene_list[i]
        start_time_gene <- Sys.time()
        for (j in seq_along(cell_type_list)) {
          cell_type <- cell_type_list[j]
          # Obtener los valores de expresión del gen y la abundancia del tipo celular
          gene_expression <- as.numeric(count_data1[gene, ])
          cell_type_abundance <- as.numeric(deconv_data[cell_type, ])
          # Realizar la prueba de correlación
          correlation_test <- cor.test(gene_expression, cell_type_abundance, method = "spearman")
          # Almacenar los resultados en los dataframes
          correlation_df[gene, cell_type] <- correlation_test$estimate
          pvalue_df[gene, cell_type] <- correlation_test$p.value
        }
        end_time_gene <- Sys.time()
        cat(sprintf("Done with gene %s (%d/%d). Time taken: %f seconds\n", gene, i,
                    total_genes, as.numeric(difftime(end_time_gene, start_time_gene, units = "secs"))))
      }
      # Retornar los resultados como una lista
      return(list(correlation = correlation_df, pvalues = pvalue_df))
    }
    cat("aaaaaaaaaaa")


    # Ejecución de la función
    if (correlation){
    cat("\033[32mResults of correlation between gene expression data and cell type abundance data. This data will be imported in .csv\033[0m\n")
    cor_epic <- decounts_correlation(tpm_counts, c_imm_epic)
    cat("aaaaaaaaaaa")
    write.csv(cor_epic, "correlation_epic.csv")
    cat("BBBBBBBBBBBBBBBB")
    cor_qti <- decounts_correlation(tpm_counts, c_imm_qti)
    cat("cccccccccc")
    write.csv(cor_qti, "correlation_qti.csv")
    cat("dddddddddddddd")
    cor_xcell <-decounts_correlation(tpm_counts, c_imm_xcell)
    cat("fffffffffff")
    write.csv(cor_xcell, "correlation_xcell.csv")}
    cat("gggggggggggggg")
    }else {cat("\033[32mConvolution analysis skipped.\033[0m\n")}

  #####
  #####
  #####
  if (survival_analysis) {
    if (is.null(variable_01) || is.null(time)) {
      stop("Variables for survival analysis are required.")
    }

    if (remove_outliers) {
      if (!is.null(pattern)) {
        # Remove outliers based on pattern
        filtered <- subset(count_data, !grepl(pattern, rownames(count_data)))
      } else {
        filtered <- count_data
      }
      counts_filtered <- filtered[, !colnames(filtered) %in% outliers]
      AnnotData <- col_data[!col_data[["id"]] %in% outliers, ]
    } else {
      counts_filtered <- count_data
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

    # Create DESeqDataSet object
    rownames(col_data) <- col_data$id
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts_filtered, colData = col_data, design = design_formul)
    dds <- DESeq2::estimateSizeFactors(dds)
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
    } else if (!is.null(genes_to_use)) {
      cat("\033[32mUsing provided genes\033[0m\n")
      top_genes <- genes_to_use
    } else {
      stop("Either 'res' or 'genes_to_use' must be provided.")
    }

    selected_df_t <- df_t[, top_genes, drop = FALSE]
    selected_df_t <- as.data.frame(selected_df_t)
    selected_df_t$id <- rownames(selected_df_t)
    col_data$id <- rownames(col_data)
    merged_data <- merge(col_data, selected_df_t, by = "id")
    merged_data[[time]] <- as.numeric(merged_data[[time]])

    cat("\033[32mStarting survival analysis.\033[0m\n")
    print(dim(merged_data))
    print(head(merged_data))
    print(length(merged_data$time))
    print(length(merged_data$variable_01))

    pdf("survival_analysis_plots.pdf")

    for (i in top_genes) {
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

      merged_data[[paste0(i, "_mRNA_expression")]] <- ifelse(merged_data[[i]] > cut.off, "High", "Low")
      merged_data[[paste0(i, "_mRNA_expression")]] <- factor(merged_data[[paste0(i, "_mRNA_expression")]])


      # Fit survival model
      cat("\033[32mFitting survival model\033[0m\n")
      column_name <- paste0(i, "_mRNA_expression")
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

  } else {cat("\033[32msurvival_analysis skipped.\033[0m\n")}

}
