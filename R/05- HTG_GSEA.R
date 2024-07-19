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
  library(clusterProfiler)
  library(dplyr)
  library(msigdbr)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(fgsea)
  library(DOSE)
  library(enrichplot)
  library(ggplot2)
  library(ggupset)
  library(grid)

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
   cat("\033[32mPreparing gene list for GSEA\033[0m\n")
  original_gene_list <- res$log2FoldChange
  names(original_gene_list) <- rownames(res)
  gene_list <- na.omit(original_gene_list)
  gene_list <- sort(gene_list, decreasing = TRUE)

  # GSEA Analysis
   cat("\033[32mPerforming GSEA Analysis\033[0m\n")
  gse2 <- gseGO(geneList = gene_list, ont = "ALL", keyType = "SYMBOL", nPermSimple = 500000,
                minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, eps = 0,
                OrgDb = "org.Hs.eg.db", pAdjustMethod = "bonferroni")

  # Dotplot for GSEA
   cat("\033[32mCreating Dotplot for GSEA\033[0m\n")
  dotplot1 <- dotplot(gse2, showCategory = 10, split = ".sign", font.size = 9, label_format = 40,
                      title = "Enrichment Results: Pathways", color = "p.adjust", size = "Count")
  add_plot(dotplot1, "GSEA Dotplot", "This plot shows the results of Gene Set Enrichment Analysis (GSEA).")

  dotplot2 <- dotplot(gse2, showCategory = 10, split = ".sign", font.size = 9, label_format = 40,
                      title = "Enrichment Results: Pathways", color = "p.adjust", size = "Count") +
    facet_grid(.~.sign)
  add_plot(dotplot2, "GSEA Dotplot with Facet", "This plot shows the results of GSEA with facets.")

  # Emaplot for GSEA
   cat("\033[32mCreating Emaplot for GSEA\033[0m\n")
  x2 <- pairwise_termsim(gse2)
  emapplot1 <- emapplot(x2, max.overlaps = 50, min.segment.length = 0.3, point_size = 0.5, font.size = 8)
  add_plot(emapplot1, "GSEA Emaplot", "This plot shows the enriched terms and their relationships.")

  # Ridgeplot for GSEA
   cat("\033[32mCreating Ridgeplot for GSEA\033[0m\n")
  ridgeplot1 <- ridgeplot(gse2) + labs(x = "enrichment distribution", font.size = 7)
  add_plot(ridgeplot1, "GSEA Ridgeplot", "This plot shows the enrichment distribution of gene sets.")

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
}
