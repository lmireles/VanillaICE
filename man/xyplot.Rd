\name{xyplot}
\alias{xyplot}
\alias{xyplot,formula,SnpSet-method}
\title{Plot copy number and physical position for a set of genomic intervals.}
\description{Plot copy number and physical position given by a
  \code{CNSet} object for a set of genomic intervals stored in a
  \code{RangedDataCVN} object.
}
\usage{
xyplot(x, data, ...)
}
\arguments{
  \item{x}{A \code{formula}.  Currently, the formula must be one of \code{cn~x},
    \code{cn ~ x | id} or \code{cn ~ x | range}.}
  \item{data}{A \code{CNSet} object.}
  \item{\dots}{A \code{RangedDataCNV} object must be passed by the name
    'range'. Additional arguments are passed to functions \code{panel.xyplot},
    \code{lpoints}, and \code{lrect}. See examples.}
}
\details{


  For a given \code{RangedDataCNV} object, this function will plot the
  copy number estimates versus physical position. The function is
  particularly useful for multi-panel displays in which the copy number
  estimates for a single range of the \code{RangedDataCNV} object
  appears in one panel. The size of the multi-panel display depends on
  the number of ranges (rows) in the \code{RangedDataCNV} object.
  Typically, one would want to pass no more than 10 ranges to the
  \code{xyplot} function.

  For genomic intervals of interest in the \code{RangedDataCNV}, it is
  often helpful to 'frame' the interval by plotting the data surrounding
  the interval. To facilitate this process, one may pass an argument
  called \code{frame} (an integer) that indicates the number of
  basepairs to the left and right of the start / stop points in the
  interval.  By default, the first interval in the \code{RangedDataCNV}
  object will be plotted in the lower left panel and the last interval
  \code{RangedDataCNV} object will be plotted in the upper right panel.
  Overplotting the copy number data in each panel is a rectangle that
  indicates the start and stop coordinates in the \code{RangedDataCNV}
  object.


}

\value{

  An object of class \code{trellis}.

}

\note{

  Note that users must pass a \code{RangedDataCNV} object called
  'range'.  As mentioned previously, it can be helpful to pass an
  integer called 'frame' that indicates how much contextual data we
  should plot surrounding each genomic interval.

}
\author{
R. Scharpf
}

\seealso{

  To modify the plot appearance from the default, additional arguments
  can be passed to \code{\link{panel.xyplot}}, \code{\link{lpoints}}, and
  \code{\link{lrect}}.

}
\examples{
	## simulated data
	data(oligoSetExample, package="oligoClasses")
	hmmOpts <- HmmOptionList(oligoSet,
				 copynumberStates=log2(c(0.1,1,2,2,3)))
	## plot the 2nd range, framing the borders on each side by 2e6
	## bases
	fit <- hmm(oligoSet, hmmOpts)
	xyplot(cn ~ x, data=oligoSet, range=fit[2,], frame=2e6,
	       panel=xypanel, cex=0.5, pch=21, border="orange")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{dplot}
\keyword{methods}

