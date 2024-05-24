\name{cts_POS}
\docType{data}
\title{cts_POSitive Control Data}
\description{
  This dataset contains cts_POSitive control samples for the HTGAnalyzer package.
}
\usage{cts_POS}
\format{
  A data frame with several variables:
    \describe{
      \item{\code{Sample_ID}}{Character: Sample identifier.}
      \item{\code{Reads}}{Integer: Number of reads for each sample.}
    }
}
\examples{
  data(cts_POS)
  head(cts_POS)
}
