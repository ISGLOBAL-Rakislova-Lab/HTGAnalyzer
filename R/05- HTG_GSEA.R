#' HTG_GSEA
#'
#' Perform Gene Set Enrichment Analysis (GSEA) using DESeq2 results, generating various plots and saving results.
#'
#' @description This function conducts GSEA to identify enriched gene sets based on DESeq2 log2 fold changes.
#'
#' @param res The DESeq2 results object containing log2 fold changes and statistical information.
#' @param gene_list A vector of log2 fold changes representing gene expression changes used for GSEA.
#' @param kegg_gene_list A vector of log2 fold changes representing KEGG gene expression changes used for KEGG analysis.
#'
#' @details
#' This function performs GSEA using DESeq2 results (`res`) to identify enriched gene sets. It generates several plots including Dotplot, Emaplot, and Ridgeplot for both general enrichment and KEGG pathways. The plots are saved in a PDF file named "GSEA_analysis_plots.pdf" in the working directory. If `res` does not contain significant results (padj < 0.05), the function will not generate certain plots and will display an error message.
#'
#' @return This function does not explicitly return a value. It saves plots and tables related to GSEA results in the "GSE_results" folder and generates two Excel files containing the results of `gene_list` and `kegg_gene_list`.
#'
#' @export
#'
#' @examples
#' HTG_GSEA(res)
#' @name HTG_GSEA
#'
HTG_GSEA <- function(res) {
  suppressMessages(library(clusterProfiler))
  suppressMessages(library(dplyr))
  suppressMessages(library(msigdbr))
  suppressMessages(library(enrichplot))
  suppressMessages(library(org.Hs.eg.db))
  suppressMessages(library(fgsea))
  suppressMessages(library(DOSE))
  suppressMessages(library(enrichplot))
  suppressMessages(library(ggplot2))
  suppressMessages(library(ggupset))
  suppressMessages(library(grid))

  cat("\033[32mPerforming gseGO  analysis\033[0m\n")

  # Prepare gene list for gseGO
  cat("\033[32mPreparing gene list for GSEA\033[0m\n")
  original_gene_list <- res$log2FoldChange
  names(original_gene_list) <- rownames(res)
  gene_list <- na.omit(original_gene_list)
  gene_list <- sort(gene_list, decreasing = TRUE)

  # gseGO  Analysis
  cat("\033[32mPerforming gseGO Analysis\033[0m\n")
  gse2 <- gseGO(geneList = gene_list, ont = "BP", keyType = "SYMBOL", nPermSimple = 500000,
                minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, eps = 0,
                OrgDb = "org.Hs.eg.db", pAdjustMethod = "bonferroni")

  # Dotplot for gseGO
   cat("\033[32mCreating Dotplot for gseGO \033[0m\n")
  dotplot1 <- dotplot(gse2, showCategory = 10, split = ".sign", font.size = 9, label_format = 40,
                      title = "gseGO Enrichment Results: Pathways", color = "p.adjust", size = "Count")
  print(dotplot1)

  dotplot2 <- dotplot(gse2, showCategory = 10, split = ".sign", font.size = 9, label_format = 40,
                      title = "gseGO Enrichment Results: Pathways", color = "p.adjust", size = "Count") + facet_grid(.~.sign)
  print(dotplot2)

  # Emaplot for gseGO
  cat("\033[32mCreating Emaplot for gseGO \033[0m\n")
  x2 <- pairwise_termsim(gse2)
  emapplot1 <- emapplot(x2, max.overlaps = 70, min.segment.length = 0.3, point_size = 0.3, font.size = 5) +   ggtitle("Enrichment Map gseGO ")
  print(emapplot1)
  # Ridgeplot for gseGO
  cat("\033[32mCreating Ridgeplot for gseGO \033[0m\n")
  ridgeplot1 <- ridgeplot(gse2)  +  labs(x = "gseGO enrichment distribution", font.size = 7) +  theme(axis.text.y = element_text(size = 9))
  print(ridgeplot1)

  # Heatplot for gseGO
  cat("\033[32mCreating Heatplot for gseGO \033[0m\n")
  heatplot1 <- heatplot(gse2, showCategory = 10) + ggtitle("gseGO Heatplot")
  print(heatplot1)

  # Treeplot for gseGO
  cat("\033[32mCreating Treeplot for gseGO \033[0m\n")
  treeplot1 <- treeplot(x2) + ggtitle("gseGO Treeplot")
  print(treeplot1)

  # Create gseaplot2 plots with titles
  a <- gseaplot2(gse2, geneSetID = 1, title = paste("GSEA Plot:", gse2$Description[1]))
  print(a)

  b <- gseaplot2(gse2, geneSetID = 1:5, pvalue_table = TRUE, title = "GSEA: Top 5 Gene Sets")
  print(b)


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

  # Dotplot for KEGG
  cat("\033[32mCreating Dotplot for KEGG\033[0m\n")
  dotplot3 <- dotplot(kk2, showCategory = 10, title = "Enriched Pathways for KEGG", split = ".sign", font.size = 9) + facet_grid(.~.sign)
  print(dotplot3)

  # Emaplot for KEGG
  cat("\033[32mCreating Emaplot for KEGG\033[0m\n")
  x3 <- pairwise_termsim(kk2)
  emapplot2 <- emapplot(x3, font.size = 8 +  ggtitle("KEGG Enrichment Map"))
  print(emapplot2)

  # Ridgeplot for KEGG
  cat("\033[32mCreating Ridgeplot for KEGG\033[0m\n")
  ridgeplot2 <- ridgeplot(kk2) +  labs(x = "KEGG enrichment distribution", font.size = 6) +  theme(axis.text.y = element_text(size = 9))
  print(ridgeplot2)

  # Heatplot for KEGG
  cat("\033[32mCreating Heatplot for KEGG \033[0m\n")
  heatplot2 <- heatplot(kk2, showCategory = 10) + ggtitle("KEGG Heatplot")
  print(heatplot2)

  # Treeplot for KEGG
  cat("\033[32mCreating Treeplot for KEGG \033[0m\n")
  treeplot2 <- treeplot(x3) + ggtitle("KEGG Treeplot")
  print(treeplot2)

  upset_plot <- upsetplot(kk2) + labs(title = "Up set plot for KEGG")
  print(upset_plot)

  # enrichGO Analysis
  cat("\033[32mPerforming GO Enrichment Analysis\033[0m\n")
  sig_genes_df <- subset(res, padj < 0.05)

  if (nrow(sig_genes_df) > 0) {
    # Prepare the gene list for enrichment analysis
    genes <- sig_genes_df$log2FoldChange
    names(genes) <- rownames(sig_genes_df)

    # Perform GO enrichment analysis
    cat("\033[32mPerforming GO Enrichment Analysis\033[0m\n")
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

    # Create a bar plot for significant GO terms
    cat("\033[32mCreating Bar Plot for GO Enrichment\033[0m\n")
    go_results <- go_enrich@result
    significant_terms <- go_results[go_results$qvalue < 0.05, ]
    significant_terms <- significant_terms[order(significant_terms$qvalue), ]

    # Check if there are significant terms to plot
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

    # Save the GO enrichment results to a CSV file
    write.csv(go_results, "go_results.csv")
    cat("\033[32mGO enrichment results saved to 'go_results.csv'\033[0m\n")
  } else {
    cat("\033[31mNo significant genes found for GO enrichment analysis.\033[0m\n")
  }

    pdf("GSEA_analysis_plots.pdf")
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
      dev.off()

      # Save tables
      write.csv(gene_list, "gene_list.csv")
      write.csv(kegg_gene_list, "kegg_gene_list.csv")
    }
}

