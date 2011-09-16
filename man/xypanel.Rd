\name{xypanel}
\alias{xypanel}
\title{A panel function for plotting copy number versus physical position}
\description{A panel function for \code{xyplot} for plotting copy number versus physical position.}
\usage{
xypanel(x, y, gt, is.snp, range, col.hom = "grey60", fill.hom = "lightblue", col.het = "grey60", fill.het = "salmon", col.np = "grey20", fill.np = "green3", ..., subscripts)
}
\arguments{
  \item{x}{Physical position in megabases.}
  \item{y}{Copy number estimates.}
  \item{gt}{Genotype calls.}
  \item{is.snp}{Logical. Whether the marker is polymorphic.}
  \item{range}{A \code{RangedData} or \code{IRanges} object. Note that we
    expect the units returned by \code{start} and \code{end} to be basepairs.}
  \item{col.hom}{A specification for the color of plotting symbols
    for homozygous genotypes.}
  \item{fill.hom}{A specification for the fill color of plotting symbols
    for homozygous genotypes.
  }
  \item{col.het}{A specification for the color of plotting symbols
    for heterozygous genotypes.
}
  \item{fill.het}{A specification for the fill color of plotting symbols
    for heterozygous genotypes.
  }
  \item{col.np}{A specification for the color of plotting symbols
    for nonpolymorphic markers.
}
  \item{fill.np}{A specification for the fill color of plotting symbols
    for nonpolymorphic genotypes.
}
  \item{\dots}{
    Additional arguments passed to lattice functions \code{xyplot},
    \code{lpoints}, and \code{lrect}.
}
  \item{subscripts}{
%%     ~~Describe \code{subscripts} here~~
}
}
\details{

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\note{
  To make the drawing of the \code{range} object border invisible, one can use border="white".
}
\author{
R. Scharpf
}

\seealso{
  \code{\link{xyplot}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{dplot}
\keyword{color}

