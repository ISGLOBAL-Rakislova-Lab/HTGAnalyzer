\name{AnnotData}
\docType{data}
\title{Annotation Data}
\description{
  This dataset contains annotation information for the HTGAnalyzer package.
}
\usage{AnnotData}
\format{
  A data frame with several variables:
  \describe{
    \item{\code{Gene}}{Character: Gene names.}
    \item{\code{Description}}{Character: Descriptions of the genes.}
    \item{\code{Chromosome}}{Character: Chromosomal locations.}
    \item{\code{Start}}{Integer: Start positions on the chromosome.}
    \item{\code{End}}{Integer: End positions on the chromosome.}
  }
}
\examples{
data(AnnotData)
head(AnnotData)
}
