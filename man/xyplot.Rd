\name{xyplot}
\alias{xyplot}
\alias{xyplot2}
\alias{xyplot,formula,SnpSet-method}
\alias{xyplot,formula,BeadStudioSet-method}
\alias{xyplot2,formula,CNSet,RangedDataCNV-method}
\alias{xyplot2,formula,gSet,RangedDataCNV-method}
\alias{xyplot2,formula,SnpSet,RangedDataCNV-method}
\title{Plot copy number and physical position for a set of genomic intervals.}
\description{Plot copy number and physical position given by a
  \code{CNSet} object for a set of genomic intervals stored in a
  \code{RangedDataCVN} object.
}
\usage{
xyplot(x, data, ...)
xyplot2(x, data, range, frame=50e3L, ...)
}
\arguments{
  \item{x}{A \code{formula}.  Currently, the formula must be one of \code{cn~x},
    \code{cn ~ x | id} or \code{cn ~ x | range} when \code{data} is a
    \code{CNSet}. If \code{data} is a \code{BeadStudioSet}, the formula
    has the form \code{lrr ~ x| range} or \code{baf ~ x | range}.}
  \item{data}{A \code{CNSet}, \code{BeadStudioSet}, or \code{SnpSet} object.}
  \item{\dots}{A \code{RangedDataCNV} object must be passed by the name
    'range'. Arguments for \code{xyplot} are passed to
    \code{xyplot2}. Additional arguments are passed to \code{xypanel} and \code{panel.xyplot}.}
  \item{range}{A \code{RangedDataCNV} object.}
  \item{frame}{The genomic distance (basepairs) to the left and right of the start and stop
    coordinates in the \code{range} object.}
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

  If the \code{lattice} package is loaded after loading
  \code{VanillaICE}, the generic definition for \code{xyplot} in
  \code{VanillaICE} will be masked. To unmask the S4 generic in
  \code{VanillaICE}, do

  xyplot <- VanillaICE:::xyplot

}
\author{
R. Scharpf
}

\seealso{

	\code{\link[lattice]{xyplot}}, \code{\link{xypanel}}

  To modify the plot appearance from the default, additional arguments
  can be passed to \code{\link{panel.xyplot}}, \code{\link{lpoints}}, and
  \code{\link{lrect}}.

}
\examples{
	## simulated data
	data(oligoSetExample, package="oligoClasses")
	data(hmmResults)
	## to visualize each range in it's own panel surrounded by a
	## frame of 2,000,000 bases:
	## (here the frames are overlapping, but the method could be
	## applied more generally to a collection of ranges from
	## different chromsomes and samples)
	xyplot(cn ~ x | range, data=oligoSet, range=hmmResults, frame=2e6,
	       panel=xypanel, cex=0.5, pch=21, border="orange",
	       scales=list(x="free"))

   if(require("lattice")){
        ## now the generic for xyplot in VanillaICE is masked.  To unmask,
        xyplot <- VanillaICE:::xyplot
   }

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{dplot}
\keyword{methods}

