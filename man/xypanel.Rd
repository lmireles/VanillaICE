\name{xypanel}
\alias{xypanel}
\title{A panel function for plotting copy number versus physical position}
\description{A panel function for \code{xyplot} for plotting copy number versus physical position.}
\usage{
xypanel(x, y, gt, is.snp, range, col.hom = "grey20", fill.hom =
"lightblue", col.het = "grey20", fill.het = "salmon", col.np = "grey20",
fill.np = "grey60", show.state=TRUE, ..., subscripts)
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
  \item{show.state}{Logical.  Whether to display the predicted state in
    each panel.}
  \item{\dots}{
    Additional arguments passed to lattice functions \code{xyplot},
    \code{lpoints}, and \code{lrect}.
}
  \item{subscripts}{
%%     ~~Describe \code{subscripts} here~~
}
}
\details{

  The order of plotting is (1) nonpolymorphic markers, (2), homozygous
  SNPs, and (3) heterozygous SNPs. Stretches of homozygosity should
  appear as blue using the default color scheme.

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
\dontrun{
if(require("crlmm")){
	data(cnSetExample, package="crlmm")
        cnSetExample <- order(cnSetExample)
	oligoSet <- as(cnSetExample, "oligoSnpSet")
        ## uses genotype emission probabilities
        hmmOpts <- HmmOptionList(oligoSet, is.log=TRUE, ICE=FALSE)
        fit2 <- hmm(oligoSet, hmmOpts, use.baf=FALSE)
        xyplot(cn ~ x | range, data=oligoSet, range=fit2[1:10, ], frame=2e6,
               panel=xypanel, cex=0.3, pch=21, border="blue",
               scales=list(x="free"),
               col.hom="lightblue", col.het="salmon", col.np="grey60",
  fill.np="grey60")

library(lattice)
## now the generic for xyplot is masked.  Must do
xyplot <- VanillaICE:::xyplot

}

}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{dplot}
\keyword{color}

