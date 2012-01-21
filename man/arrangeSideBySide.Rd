\name{arrangeSideBySide}
\alias{arrangeSideBySide}

\title{
  Helper function to arrange two trellis objects side by side on a grid.
}

\description{

  For visualizing copy number alterations, it is often helpful to plot
  estimates of copy number along with the corresponding estimate of the
  B allele frequencies.  Creating a trellis object for the copy number
  estimates and a separate trellis object for the B allele frequencies,
  this function can be used to arrange the two trellis objects side by
  side on a grid.

}

\usage{
arrangeSideBySide(object1, object2)
}

\arguments{
  \item{object1}{
    A trellis object (e.g., a trellis object of the copy number estimates).
}
  \item{object2}{
    A trellis object (e.g., a trellis object of the B allele frequencies).
  }
}


\author{
Rob Scharpf
}

\seealso{
  \code{\link{xypanel}}, \code{\link{xyplot}}
}

\keyword{hplot}

