\name{CopyNumberSet-methods}
\docType{methods}
\alias{CopyNumerSet-methods}
\alias{cnEmission,CopyNumberSet,ANY-method}
\alias{order,CopyNumberSet-method}
\alias{sd,CopyNumberSet-method}

%\alias{sd,oligoSnpSet-method}
\title{Methods for estimating copy number standard deviations.}

\description{Estimate the standard deviation for \code{CopyNumberSet}
  and \code{oligoSnpSet} objects.}

\usage{

sd(x, na.rm=TRUE)

}

\arguments{

  \item{object}{A \code{CopyNumberSet} or \code{oligoSnpSet}}

  \item{na.rm}{Logical. }

}

\value{

  A matrix.

}

\details{

  The sd method for \code{CopyNumberSet} and \code{oligoSnpSet} objects
  retrieves the copy number confidence scores from the
  \code{cnConfidence} assay data element.  The confidence matrix is a R
  x C matrix for an object with R features and C samples.  Valid
  confidence estimates must be positive and not missing (not
  \code{NA}). If any elements in the confidence matrix are invalid, a
  robust estimate of the standard deviation is computed (described
  below).  If all elements are valid, the standard deviation matrix is
  returned as 1 / confidence.

  If any elements in the confidence matrix are invalid, the standard
  deviation for each marker and sample is calculated as follows.  If
  autosomal markers are present, the standard deviation is estimated as
  the median absolute deviation across autosomal markers for each
  sample.  This gives a vector of length C.  The R x C standard
  deviation matrix is populated by row from the vector of length C (the
  standard deviation for each marker in a sample is given the same
  standard deviation).  If autosomal markers are not present, the median
  absolute deviation across X-chromosome markers and Y-chromosome
  markers are estimated independently, providing to vectors of length
  C. The matrix of standard deviations for the X chromosome is populated
  by the C-length vector for the X-chromosome (by-row) and likewise for
  the Y chromosome.

}

\seealso{
  \code{\link{mad}}
}

\examples{
  data(oligoSetExample)
  sds <- sd(oligoSet)
}

\keyword{methods}