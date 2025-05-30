#' HTG_GSEA
#'
#' Perform Gene Set Enrichment Analysis (GSEA) using DESeq2 results, generating various plots and saving results.
#'
#' @description This function conducts GSEA to identify enriched gene sets based on DESeq2 log2 fold changes.
#'
#' @param res The DESeq2 results object containing log2 fold changes and statistical information.
#'
#' @details
#' This function performs GSEA using DESeq2 results (`res`) to identify enriched gene sets. It generates several plots including Dotplot, Emaplot, and Ridgeplot for both general enrichment and KEGG pathways. The plots are saved in a PDF file named "GSEA_analysis_plots.pdf" in the working directory. If `res` does not contain significant results (padj < 0.05), the function will not generate certain plots and will display an error message.
#'
#' @return This function does not explicitly return a value. It saves plots and tables related to GSEA results in the "GSE_results" folder and generates two Excel files containing the results of `gene_list` and `kegg_gene_list`.
#'
#' @export
#'
#' @examples
#' HTG_GSEA(res_tutorial)
#' @name HTG_GSEA
#'
#'
utils::globalVariables(c("padj", "Description", "Count"))
HTG_GSEA <- function(res, pvalueCutoff = 0.05) {
  cat("\033[32mPerforming gseGO  analysis\033[0m\n")

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
  gse2 <- clusterProfiler::gseGO(geneList = gene_list, ont = "BP", keyType = "SYMBOL",# nPermSimple = 1000,
                                 minGSSize = 3, maxGSSize = 800, pvalueCutoff = pvalueCutoff, verbose = TRUE, eps = 0,
                                 OrgDb = org.Hs.eg.db::org.Hs.eg.db, pAdjustMethod = "bonferroni")
  print(gse2)


  if (nrow(gse2@result) == 0) {
    cat("\033[31mNo enriched terms found.\033[0m\n")
    cat("\033[33mℹ️ Possible reasons:\n")
    cat(" - gene_list has too many zero or low values\n")
    cat(" - Too strict pvalueCutoff or pAdjustMethod\n")
    cat(" - Genes not annotated in the selected ontology\n")
    cat("\033[34m TIP: You can try increasing the pvalueCutoff (e.g., to 0.1), or use the online GSEA tool:\n")
    cat(" https://www.gsea-msigdb.org/gsea/index.jsp\n\033[0m")
    stop("No GSEA results to plot.")
  }

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
  num_terms <- nrow(x2@result)
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
  num_terms <- nrow(x3@result)
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

    pdf("GSEA_analysis_plots1_of_2.pdf", width = 15, height = 14)
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
      pdf("GSEA_analysis_plots2_of_2.pdf", width = 39, height = 15)
      print(heatplot1)
      print(heatplot2)
        dev.off()


      # Save tables
      write.csv(gene_list, "gene_list.csv")
      write.csv(kegg_gene_list, "kegg_gene_list.csv")
    }
}
