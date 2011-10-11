\name{hmm}
\alias{hmm}

\title{Application of the Viterbi algorithm to copy number and/or
genotype data.}
\description{
  A wrapper for fitting the HMM.
}
\usage{
hmm(object, hmm.params, use.baf=FALSE, ...)
}

\arguments{

\item{object}{one of the following classes derived from eSet:
	\code{SnpSet}, \code{oligoSnpSet}, \code{CopyNumberSet},
	\code{CNSet}}

      \item{hmm.params}{\code{HmmOptionList}. See details.}

      \item{use.baf}{Logical. Whether to use the BAFs instead of the
      genotype calls and confidence scores to estimate the emission
      probabilities. See details.}

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

  For \code{oligoSnpSet} objects, the emission probability for the HMM
  is calculated as

  log emission = log GT + log CN

  where log GT is the emission probability for the genotype calls and
  log CN is the emission probability for the copy number estimates.
  When \code{use.baf} is TRUE, the emission probabilities are estimated
  from the B allele frequencies instead of the genotype calls /
  confidence scores.  Generally, use of the \code{use.baf} argument
  requires a \code{CNSet} object obtained from preprocessing with the
  \code{crlmm} package.  However, users may supply the B allele
  frequencies from external sources (e.g., BeadStudio software) when
  available.

  See the examples below for using the \code{use.baf} argument and how
  to visualize the predicted states along with the low-level genotype
  and copy number summaries.

}

\value{

	An object of class \code{RangedData}.

}

\seealso{\code{\link{hmm.setup}}, \code{\link{runmed}}, \code{\link[crlmm]{calculateRBaf}}}

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
oligoSet <- oligoSet[order(chromosome(oligoSet), position(oligoSet)),]
oligoSet <- oligoSet[chromosome(oligoSet) < 3, ]
hmmOpts <- hmm.setup(oligoSet, is.log=TRUE)
fit <- hmm(oligoSet, hmmOpts, k=3)
xyplot(cn~x, oligoSet, range=fit[4, ], frame=2e6,pch=21, cex=0.5, panel=xypanel)

## Useful accessors for RangedData
ranges(fit)

##Log likelihood ratio comparing likelihood of predicted state to the 'normal' state
## for each segment
fit$LLR
## the number of SNPs / nonpolymorphic loci in each segment
coverage2(fit)

sampleNames(fit)
chromosome(fit)

##---------------------------------------------------------------------------
##
##  Downstream of CRLMM for genotyping
##
##---------------------------------------------------------------------------
\dontrun{
if(require("crlmm")){
	data(cnSetExample, package="crlmm")
        cnSetExample <- order(cnSetExample)
	oligoSet <- as(cnSetExample, "oligoSnpSet")
        assayDataElement(oligoSet, "baf") <- calculateRBaf(cnSetExample)[["baf"]]
        ## uses genotype emission probabilities
        hmmOpts <- HmmOptionList(oligoSet, is.log=TRUE, ICE=FALSE)
        fit2 <- hmm(oligoSet, hmmOpts, use.baf=FALSE)
        xyplot(cn ~ x | range, data=oligoSet, range=fit2[1:10, ], frame=2e6,
               panel=xypanel, cex=0.3, pch=21, border="blue",
               scales=list(x="free"),
               col.hom="lightblue", col.het="salmon", col.np="grey60", fill.np="grey60")
        xyplot(baf ~ x | range, data=oligoSet, range=fit2[1:10, ], frame=2e6,
               panel=xypanel, cex=0.3, pch=21, border="blue",
               scales=list(x="free"),
               col.hom="lightblue", col.het="salmon", col.np="grey60", fill.np="grey60")
	hmmOpts <- HmmOptionList(oligoSet, is.log=TRUE)

        ## use BAFs for emission probs.
        fit3 <- hmm(oligoSet, hmmOpts, use.baf=TRUE)
        ## plot the copy number for the first 10 ranges
        xyplot(cn ~ x | range, data=oligoSet, range=fit3[1:10, ], frame=2e6,
               panel=xypanel, cex=0.2, pch=21, border="blue",
               scales=list(x="free"),
               col.hom="lightblue", col.het="salmon", col.np="grey60", fill.np="grey60")
        ## plot the BAFs for the first 10 ranges
        xyplot(baf ~ x | range, data=oligoSet, range=fit3[1:10, ], frame=2e6,
               panel=xypanel, cex=0.2, pch=21, border="orange",
               scales=list(x="free"),
               col.hom="lightblue", col.het="salmon", col.np="grey60", fill.np="grey60")
} ## if(require("crlmm"))
} ## \dontrun

}

\author{R. Scharpf}

\keyword{models}
\keyword{manip}
\keyword{ts}

