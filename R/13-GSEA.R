#' Process GSEA Results
#'
#' @description This function processes GSE (Gene Set Enrichment) results obtained from DESeq2 analysis.
#'
#' @param res The DESeq2 results object containing log2 fold changes.
#' @param gene_list A vector of log2 fold changes representing gene expression changes.
#' @param kegg_gene_list A vector of log2 fold changes representing KEGG gene expression changes.
#'
#' @return This function doesn't return any value explicitly. It saves plots and tables related to GSEA results in the "GSE_results" folder. Additionally, it generates two Excel files containing the results of gene_list and KEGG gene list.
#'
#' @export
#'
#' @examples
#' HTG_GSEAresults(res)
#' @name HTG_GSEAresults
HTG_GSEAresults <- function(res) {
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

  ##### GENE LIST
  original_gene_list <- res$log2FoldChange
  names(original_gene_list) <- rownames(res)
  gene_list <- na.omit(original_gene_list)
  gene_list <- sort(gene_list, decreasing = TRUE)

  gse2 <- gseGO(geneList = gene_list,
                ont = "ALL",
                keyType = "SYMBOL",
                nPermSimple = 500000,
                minGSSize = 3,
                maxGSSize = 800,
                pvalueCutoff = 0.05,
                verbose = TRUE,
                eps = 0,
                OrgDb = "org.Hs.eg.db",
                pAdjustMethod = "bonferroni")

  print(dotplot(gse2,
                showCategory = 10,
                split = ".sign",
                font.size = 8,
                label_format = 40,
                title = "Enrichment Results: Pathways",
                color = "p.adjust",
                size = "Count"))

  print(dotplot(gse2,
                showCategory = 10,
                split = ".sign",
                font.size = 8,
                label_format = 40,
                title = "Enrichment Results: Pathways",
                color = "p.adjust",
                size = "Count") + facet_grid(.~.sign))

  print(dotplot(gse2, showCategory = 8, split = ".sign", font.size = 8) + facet_grid(.~.sign))
  print(dotplot(gse2, showCategory = 5, split = ".sign", font.size = 8) + facet_grid(.~.sign))

  x2 <- pairwise_termsim(gse2)

  print(emapplot(x2, max.overlaps = 50, min.segment.length = 0.3, point_size = 0.5, font.size = 8))
  print(cnetplot(gse2, categorySize = "pvalue", foldChange = gene_list, showCategory = 3, font.size = 8))

  print(ridgeplot(gse2) + labs(x = "enrichment distribution", font.size = 8))

  ids <- bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  dedup_ids <- ids[!duplicated(ids[c("SYMBOL")]), ]

  df2 <- res[rownames(res) %in% dedup_ids$SYMBOL, ]
  df2$Y <- dedup_ids$ENTREZID
  kegg_gene_list <- df2$log2FoldChange
  names(kegg_gene_list) <- df2$Y
  kegg_gene_list <- na.omit(kegg_gene_list)
  kegg_gene_list <- sort(kegg_gene_list, decreasing = TRUE)

  kk2 <- gseKEGG(geneList = kegg_gene_list,
                 organism = "hsa",
                 minGSSize = 3,
                 maxGSSize = 800,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "none",
                 keyType = "ncbi-geneid",
                 nPermSimple = 100000)

  print(dotplot(kk2, showCategory = 10, title = "Enriched Pathways", split = ".sign", font.size = 8) + facet_grid(.~.sign))

  x3 <- pairwise_termsim(kk2)
  print(emapplot(x3, font.size = 8))

  print(ridgeplot(kk2) + labs(x = "enrichment distribution", font.size = 8))

  sig_genes_df <- subset(res, padj < 0.05)
  genes <- sig_genes_df$log2FoldChange
  names(genes) <- rownames(sig_genes_df)

  go_enrich <- enrichGO(gene = names(genes),
                        universe = names(gene_list),
                        OrgDb = org.Hs.eg.db,
                        keyType = 'SYMBOL',
                        readable = TRUE,
                        ont = "BP",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.10)

  #####################
  # Genes y pathways de humanos
  hs_hallmark_sets <- msigdbr(
    species = "Homo sapiens", # Replace with species name relevant to your data
    category = "H"
  )

  # Enrichment
  gsea_results <- GSEA(
    geneList = gene_list, # Ordered ranked gene list
    minGSSize = 5, # Minimum gene set size
    maxGSSize = 500, # Maximum gene set set
    pvalueCutoff = 0.1, # p-value cutoff
    eps = 0, # Boundary for calculating the p value
    seed = TRUE, # Set seed to make results reproducible
    pAdjustMethod = "BH", # Benjamini-Hochberg correction
    nPermSimple = 10000,
    TERM2GENE = dplyr::select(
      hs_hallmark_sets,
      gs_name,
      human_gene_symbol
    )
  )

  # Visualize results
  print(barplot(go_enrich,
                drop = TRUE,
                showCategory = 10,
                title = "GO Biological Pathways",
                font.size = 8))

  print(upsetplot(go_enrich, show.numbers = "yes", font.size = 8))

  # Save tables
  write.csv(gene_list, "gene_list.csv")
  write.csv(kegg_gene_list, "kegg_gene_list.csv")
}
