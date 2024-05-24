\name{res}
\docType{data}
\title{DESeq Results}
\format{
  A `DESeqResults` object, which is a DataFrame with the following columns:
  \describe{
    \item{baseMean}{The average normalized count for each gene.}
    \item{log2FoldChange}{The log2 fold change for each gene.}
    \item{lfcSE}{The standard error of the log2 fold change.}
    \item{stat}{The Wald statistic.}
    \item{pvalue}{The p-value for the statistical test.}
    \item{padj}{The adjusted p-value (FDR) for multiple testing correction.}
  }
}
\source{
  This object is generated using the `DESeq2` package.
}
\examples{
# Example of how to generate `res` using DESeq2
\dontrun{
  library(DESeq2)
  dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = col_data, design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds)
}
}
