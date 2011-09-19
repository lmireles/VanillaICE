\name{cnEmission}
\alias{cnEmission}
\title{Calculate emission probabilities for total copy number}
\description{Calculate emission probabilities for total copy number from
a Uniform-Gaussian mixture.}

\usage{
cnEmission(object, stdev, k = 5, cnStates = 0:4, is.log = FALSE, verbose = TRUE, ...)
}

\arguments{
  \item{object}{A \code{CopyNumberSet}, \code{oligoSnpSet}, or \code{matrix}.}
  \item{stdev}{A matrix.  Ignored unless class of object is
    \code{matrix}. See details}
  \item{k}{Integer. Size of window for running median.  A running median
  of the total copy number is used to estimate the probability that a
  copy number estimate is an outlier. }
  \item{cnStates}{Numeric or integer.  The theoretical or expected copy
    number for each hidden state.}
  \item{is.log}{Logical.  TRUE if the copy number estimates in
    \code{object} are on the log scale.}
  \item{verbose}{Logical.}
  \item{\dots}{Ignored}
}
\details{

  We calculate the emission probabilities of the total copy number (CN)
  estimates from a Normal-Uniform mixture.  In particular, we assume the
  CN (suitably transformed) is emitted from a Normal distribution with
  mean given by \code{cnStates}. As outliers are common in
  high-throughput arrays, we allow for unusual values by adding a
  Uniform component to the mixture model that covers the support of the
  CN.  (The support is determined by whether the CN is on the log scale
  as indicated by the \code{is.log} argument).  To estimate the
  probability that CN is an outlier, we calculate CN - CNsmooth where
  CNsmooth is calculated from a running median with window given by
  argument \code{k}.  We assume that the difference (CN-CNsmooth) is a
  mixture of two Normal distributions -- copy number estimates that are
  not outliers should have a Normal distribution with mean zero and
  standard deviation 'sigma1', whereas outliers follow a Normal
  distribution with mean zero and standard deviation 'sigma2', sigma2 >>
  sigma1.  We estimate the responsibilities for the mixture via EM, and
  use these values as a marker-specific estimate of the outlier
  probability.  The emission probability is given by

  pihat * N(mean copy number state, sd) + (1-pihat) * Unif(MIN, MAX).

}

\value{

  Returns an array of the emission probabilities.  The dimensions of the
  array are [feature index,  sample index, state index].

}

\author{
R. Scharpf
}
\seealso{

  \code{\link{cnEmission-methods}}

 See \code{\link{hmm}} method estimates emission probabilities and
 fits the Viterbi algorithm.

 See \code{\link{gtEmission}} for estimating the emission
 probabilities of diallic genotypes for each of the copy number
 states.

}
\examples{
data(oligoSetExample)
oligoSet <- order(oligoSet)
cn.emit <- cnEmission(oligoSet, k=5,
		      cnStates=log2(c(0.1, 1, 2, 3, 4)),
		      is.log=TRUE)

## Note that we require the features to be ordered by
## chromosome and physical position
cn.states <- log2(c(0.1, 1, 2, 3, 4))
oligoSet2 <- oligoSet[sample(seq_len(nrow(oligoSet)), nrow(oligoSet)), ]
tryCatch(cnEmission(oligoSet2, cnStates=cn.states, is.log=TRUE), error=function(e) NULL)
oligoSet3 <- order(oligoSet)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
\keyword{array}
