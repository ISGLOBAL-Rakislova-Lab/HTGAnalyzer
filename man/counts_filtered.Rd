\name{counts_filtered}
\docType{data}
\title{Filtered Counts Data}
\description{
  This dataset contains filtered counts data for the HTGAnalyzer package.
}
\usage{counts_filtered}
\format{
  A data frame with several variables:
    \describe{
      \item{\code{Gene}}{Character: Gene names.}
      \item{\code{Sample1}}{Integer: Filtered counts for sample 1.}
      \item{\code{Sample2}}{Integer: Filtered counts for sample 2.}
      \item{\code{...}}{Integer: Filtered counts for additional samples.}
    }
}
\examples{
  data(counts_filtered)
  head(counts_filtered)
}
