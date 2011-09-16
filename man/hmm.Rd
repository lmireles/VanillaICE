\name{hmm}
\alias{hmm}

\title{Application of the Viterbi algorithm to copy number and/or
genotype data.}
\description{
  A wrapper for fitting the HMM.
}
\usage{
hmm(object, hmm.params, ...)
}

\arguments{

\item{object}{one of the following classes derived from eSet:
	\code{SnpSet}, \code{oligoSnpSet}, \code{CopyNumberSet},
	\code{CNSet}}

\item{hmm.params}{\code{HmmOptionList}. See details.}

\item{\dots}{The argument \code{k} to the function \code{runmed} can be
  passed for assessing the probability of an outlier. See details.}
}

\details{

  The probability that a point estimate is an outlier is assessed by
  subtracting the running median of the copy number estimates from the
  point estimates (referred to as 'delta'). For copy number changes that
  span multiple markers, we expect delta to be near zero and to have a
  small variance.  For single point outliers, delta can be large in
  absolute value and have a large variance. We estimate the mixture
  probabilities (the probability that an observation is an outlier)
  using EM.  We plug in the estimate of the outlier probability into the
  emission probabilities.  Specifally, we model the emission
  probabilities as a Uniform-Gaussian mixture, where the Uniform
  component handles outliers.

  The function used for the running median is \code{runmed}. One can
  pass an integer value of $k$ to \code{runmed}.  Larger values of
  \code{k} should result in fewer segments with low coverage.

}

\value{

	An object of class \code{RangedData}.

}

\seealso{\code{\link{hmm.setup}}, \code{\link{runmed}}}

\references{

    RB Scharpf et al. (2008) Hidden Markov Models for the assessment of
  chromosomal alterations using high-throughput SNP arrays, Annals of
  Applied Statistics

  }

\examples{
data(locusLevelData)
oligoSet <- new("oligoSnpSet",
		copyNumber=log2(locusLevelData[["copynumber"]]/100),
		call=locusLevelData[["genotypes"]],
		callProbability=locusLevelData[["crlmmConfidence"]],
		annotation=locusLevelData[["platform"]])
oligoSet <- oligoSet[!is.na(chromosome(oligoSet)), ]
oligoSet <- oligoSet[order(chromosome(oligoSet), position(oligoSet)),]
oligoSet <- oligoSet[chromosome(oligoSet) < 3, ]
hmmOpts <- hmm.setup(oligoSet,
		     copynumberStates=log2(c(1, 2, 2, 3)),
		     states=c("hem-del", "ROH", "normal", "amp"),
		     normalIndex=3,
		     log.initialP=rep(log(1/4), 4),
		     prGenotypeHomozygous=c(0.99, 0.99, 0.7, 0.7))
fit <- hmm(oligoSet, hmmOpts, k=3)

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

