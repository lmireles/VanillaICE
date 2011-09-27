\name{hmmResults}
\alias{hmmResults}
\docType{data}
\title{
  Example output from hmm
}
\description{
  Example output from \code{hmm} method applied to simulated data.
}
\usage{data(hmmResults)}
\format{
  A \code{RangedDataHMM} object.
}
\details{
  The results of a 6-state HMM fit to simulated copy number and genotype
  data.
}

\seealso{\code{\link{xyplot}}}

\examples{
	data(oligoSetExample, package="oligoClasses")
	oligoSet <- oligoSet[chromosome(oligoSet) == 1, ]
	hmmOpts <- hmm.setup(oligoSet, is.log=TRUE)
	\dontrun{
		hmmResults <- hmm(oligoSet, hmmOpts, k=3)
		save(hmmResults, file="~/Software/VanillaICE/data/hmmResults.rda")
	}
	data(hmmResults)
	state(hmmResults)
	hmmOpts$states[state(hmmResults)]
	##
	## Plot region of homozygosity
	##  (Note that the formula cn~x is required at this time)
	xyplot(cn~x, oligoSet, range=hmmResults[2, ], frame=10e6,
	       panel=xypanel, pch=21, cex=0.6)
	##
	## Or, plot each range in its own panel with a frame of 2e6 bases.
	##  (Again, the formula is a standard format with cn, x,
	##  range, and id the only allowed terms)
	## Because these are all the ranges from one individual's
	## chromosome, the ranges are overlapping
	##  The range 'in focus' is demarcated by vertical blue lines
	xyplot(cn~x | range, oligoSet, range=hmmResults, frame=2e6,
	       panel=xypanel, pch=21, cex=0.2,
	       scales=list(x="free"),
	       border="blue",
	       col.np="grey",
	       col.het="grey",
	       col.hom="grey")
}
\keyword{datasets}
