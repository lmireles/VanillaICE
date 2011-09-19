\name{order}
\alias{order,CopyNumberSet-method}
\alias{order,SnpSet-method}
\alias{order}
\title{Methods for CopyNumberSet objects}

\description{Methods for objects of class \code{CopyNumberSet}.}

\usage{order(..., na.last=TRUE, decreasing=FALSE)}

\arguments{
  \item{\dots}{blah blah}
%  \item{\dots}{For the method \code{order}, if the first element of the
%  \dots list is a \code{SnpSet} or \code{CopyNumberSet} the object is
%  reordered by chromosome and physical position.}

\item{na.last}{Ignored}

\item{decreasing}{Ignored}

}

\value{

An object of the same class as the first element in \dots.

}

\details{

Reorders the object by chromosome and physical position.

}

\seealso{
  \code{\link{mad}}, \code{\link{cnEmission}}, \code{\link{chromosome}}, \code{\link{position}}
}

\examples{
  data(oligoSetExample)
  oligoSet2 <- order(oligoSet)
}

\keyword{methods}