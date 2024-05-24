\name{cts_GDNA}
\docType{data}
\title{Genomic DNA Data}
\description{
  This dataset contains genomic DNA controls for the HTGAnalyzer package.
}
\usage{cts_GDNA}
\format{
  A data frame with several variables:
    \describe{
      \item{\code{Sample_ID}}{Character: Sample identifier.}
      \item{\code{Reads}}{Integer: Number of reads for each sample.}
    }
}
\examples{
  data(cts_GDNA)
  head(cts_GDNA)
}
