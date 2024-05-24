\name{cts_NC}
\docType{data}
\title{Negative Control Data}
\description{
  This dataset contains negative control samples for the HTGAnalyzer package.
}
\usage{cts_NC}
\format{
  A data frame with several variables:
    \describe{
      \item{\code{Sample_ID}}{Character: Sample identifier.}
      \item{\code{Reads}}{Integer: Number of reads for each sample.}
    }
}
\examples{
  data(cts_NC)
  head(cts_NC)
}
