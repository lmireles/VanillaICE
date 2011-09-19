\name{gtEmission}
\alias{gtEmission}
\title{Calculate emission probabilities for genotypes}
\description{
  Calculate emission probabilities for genotypes diallelic genotypes.
}
\usage{
gtEmission(object, hmmOptions, gt.conf, is.snp, cdfName, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    A object of class \code{SnpSet}, \code{oligoSnpSet}, or \code{matrix}.
  }

  \item{hmmOptions}{
    A \code{HmmOptionList}.
  }

  \item{gt.conf}{Ignored unless object is of class \code{matrix}.}

  \item{is.snp}{Ignored unless object is of class \code{matrix}.}

  \item{cdfName}{Ignored unless the argument \code{ICE} for the
    \code{HmmOptionList} object is TRUE.}

  \item{\dots}{Ignored}
}

\details{

}
\value{
  An array.  The dimensions are features x samples x states.
}
\author{
R.Scharpf
}
\note{
Only chromosomes 1-24 supported (23=X and 24 = Y).
}

\seealso{
  \code{\link{cnEmission}}
}
\examples{
data(oligoSetExample)
oligoSet <- order(oligoSet)
hmmOpts <- HmmOptionList(oligoSet, is.log=TRUE)
gt.emit <- gtEmission(oligoSet, hmmOpts)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
\keyword{array}
