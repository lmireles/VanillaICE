\name{hmm}
\alias{hmm}

\title{Application of the Viterbi algorithm to copy number and/or
genotype data.}
\description{
  A wrapper for fitting the HMM.
}
\usage{
hmm(object, hmm.params, use.baf=FALSE, k=5, ...)
}

\arguments{

\item{object}{one of the following classes derived from eSet:
	\code{SnpSet}, \code{oligoSnpSet}, \code{CopyNumberSet},
	\code{CNSet}}

      \item{hmm.params}{\code{HmmOptionList}. See details.}

      \item{use.baf}{Logical. Whether to use the BAFs instead of the
      genotype calls and confidence scores to estimate the emission
      probabilities. See details.}

    \item{k}{\code{integer}: The window size for a running median
    smoother.  When estimating the probability that the copy number
    estimate is an outlier, a running median smoother is fit to the copy
    number estimates.  A two-component Gaussian mixture model is fitted
    to the residuals of the smoother.  The mixture probability for the
    'high variance' distribution is the probability that the observation
    is an outlier.}

\item{\dots}{The argument \code{k} to the function \code{runmed} can be
  passed for assessing the probability of an outlier. See details.}
}

\details{

Emission probabilities are estimated from a Uniform-Gaussian mixture,
  where the Uniform component handles outliers.  The function used for
  the running median smoother when estimating the mixture probability
  for the outlier distribution is \code{runmed}. One can pass an integer
  value of $k$ to \code{runmed}.  Larger values of \code{k} should
  result in fewer segments with low coverage.

  For \code{oligoSnpSet} objects, the emission probability for the HMM
  is calculated as

  log emission = log GT + log CN

  where log GT is the emission probability for the genotype calls (or
  the B allele frequences if \code{use.baf=TRUE} and log CN is the
  emission probability for the copy number estimates.

}

\value{

	An object of class \code{RangedData}.

}

\seealso{\code{\link{hmm.setup}}, \code{\link{runmed}},
  \code{\link[crlmm]{calculateRBaf}},
  \code{\link{hmm3}}
}

\references{

    RB Scharpf et al. (2008) Hidden Markov Models for the assessment of
  chromosomal alterations using high-throughput SNP arrays, Annals of
  Applied Statistics

  }

\examples{
data(locusLevelData, package="oligoClasses")
oligoSet <- new("oligoSnpSet",
		copyNumber=log2(locusLevelData[["copynumber"]]/100),
		call=locusLevelData[["genotypes"]],
		callProbability=locusLevelData[["crlmmConfidence"]],
		annotation=locusLevelData[["platform"]])
oligoSet <- oligoSet[!is.na(chromosome(oligoSet)), ]
oligoSet <- chromosomePositionOrder(oligoSet)
oligoSet <- oligoSet[chromosome(oligoSet) < 3, ]
hmmOpts <- hmm.setup(oligoSet, is.log=TRUE)
fit <- hmm(oligoSet, hmmOpts, k=3)
\dontrun{
xyplot(cn~x, oligoSet, range=fit[4, ], frame=2e6,pch=21, cex=0.5,panel=xypanel)
}
## Useful accessors for RangedData
ranges(fit)

##Log likelihood ratio comparing likelihood of predicted state to the 'normal' state
## for each segment
fit$LLR
## the number of SNPs / nonpolymorphic loci in each segment
coverage2(fit)
sampleNames(fit)
chromosome(fit)
}

\author{R. Scharpf}

\keyword{models}
\keyword{manip}
\keyword{ts}

