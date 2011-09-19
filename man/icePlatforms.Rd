\name{icePlatforms}
\alias{icePlatforms}
\title{List platforms for which ICE option is supported.}
\description{
}
\usage{
icePlatforms()
}
%- maybe also 'usage' for other objects documented here.
\details{

  When procecessing genotypes with the \pkg{crlmm}, confidence scores
  for the diallelic genotype calls are available.  One can estimate the
  emission probabilities for the crlmm diallelic genotypes using the
  confidence scores by setting the value of \code{ICE} to TRUE in the
  constructor for the \code{HmmOptionList} class.  Currently, only
  certain platforms are supported for this option.

}
\value{
  A character vector of the annotation packages that are supported for
  the ICE option
}

\references{Scharpf, RB et al., 2008, Annals of Applied Statistics}

\author{
R. Scharpf
}


\seealso{
  \code{\link{HmmOptionList}}, \code{\link{gtEmission}}
}
\examples{
icePlatforms()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}

