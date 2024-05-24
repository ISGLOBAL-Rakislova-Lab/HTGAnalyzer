\name{counts}
\docType{data}
\title{Raw Counts Data}
\description{
  This dataset contains raw counts data for the HTGAnalyzer package.
}
\usage{counts}
\format{
  A data frame with several variables:
    \describe{
      \item{\code{Gene}}{Character: Gene names.}
      \item{\code{Sample1}}{Integer: Counts for sample 1.}
      \item{\code{Sample2}}{Integer: Counts for sample 2.}
      \item{\code{...}}{Integer: Counts for additional samples.}
    }
}
\examples{
  data(counts)
  head(counts)
}
