\name{CopyNumberSet-methods}
\docType{methods}
\alias{CopyNumerSet-methods}
\alias{order,CopyNumberSet-method}
\alias{sd,CopyNumberSet-method}
\title{Methods for CopyNumberSet objects}

\description{Methods for objects of class \code{CopyNumberSet}.}

\usage{

%cnEmission(object, stdev, k=5, cnStates=0:4, is.log=FALSE, verbose=TRUE, ...)
order(..., na.last=TRUE, decreasing=FALSE)
sd(x, na.rm=TRUE)



}

\arguments{

%  \item{object}{A \code{CopyNumberSet} or \code{SnpSet} object}

  \item{x}{A \code{CopyNumberSet} or \code{SnpSet} object}

  \item{\dots}{For the method \code{order}, if the first element of the
  \dots list is a \code{SnpSet} or \code{CopyNumberSet} the object is
  reordered by chromosome and physical position.}

  \item{na.rm}{Logical. Ignored}

  \item{na.last}{Ignored}

  \item{decreasing}{Ignored}

}

\value{

  The method \code{order} returns and object of the same class.

%  The method \code{cnEmission} returns an array (features x samples x state).

  The method \code{sd} return a matrix (features x samples)

}

\details{

  The \code{sd} method for \code{CopyNumberSet} and \code{oligoSnpSet}
  objects retrieves the copy number confidence scores from the
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
  \code{\link{mad}}, \code{\link{cnEmission}}
}

\examples{
  data(oligoSetExample)
  sds <- sd(oligoSet)
  oligoSet2 <- order(oligoSet)
}

\keyword{methods}