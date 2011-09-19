\name{oligoSnpSet-methods}
\Rdversion{1.1}
\docType{class}
\alias{cnEmission,oligoSnpSet-method}
\alias{order,oligoSnpSet-method}
\alias{sd,oligoSnpSet-method}

\title{Methods for oligoSnpSet class}
\description{Methods for the \code{oligoSnpSet} class}

\section{Methods}{
  Standard generic methods:
  \describe{
  \item{\code{cnEmission(object, cnStates=0:4)}}{Estimate emission
    probabilities for copy number. See \code{\link{cnEmission}}.}

  \item{\code{order(x)}}{Order object by chromosome and physical position.}

  \item{\code{sd(x)}}{Estimate copy number standard deviation. See
    \code{\link{sd}}.}
}

}

\author{
R. Scharpf
}

\seealso{
  \code{\link{sd}}, \code{\link{cnEmission}}, \code{\link{order}}
}
\examples{
showClass("oligoSnpSet")
}
\keyword{classes}
