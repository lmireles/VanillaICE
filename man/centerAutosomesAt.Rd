\name{centerAutosomesAt}
\alias{centerAutosomesAt}
\title{Center estimates of copy number for autosomes.}
\description{Center estimates of copy number for autosomes.}
\usage{
centerAutosomesAt(x, at, ...)
}
\arguments{
  \item{x}{A \code{oligoSnpSet} or \code{CopyNumberSet} object.}
  \item{at}{numeric.  Value at which to center the copy number estimates (e.g., 2).}
  \item{\dots}{Ignored}
}
\details{
  The function sweeps out the column median of the autosomal copy number
  estimates, and adds back a constant given by \code{at}.
}

\value{ A object of the same class as \code{x} with centered copy number
  estimates for the autosomal chromosomes.  Chromosomes X and Y are not
  centered.  }

\author{
R. Scharpf
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{link{sweep}}
}

\keyword{manip}

