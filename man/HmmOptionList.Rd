\name{HmmOptionList}
\alias{HmmOptionList}
\alias{hmm.setup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Constructors for HmmOptionList class}
\description{
Constructors for \code{HmmOptionList} class.  An object of the class
contains most of the options for fitting the HMM.
}
\usage{
HmmOptionList(object, copynumberStates = 0:4, states = paste("state",
1:length(copynumberStates), sep = ""), ICE = FALSE, is.log = FALSE,
scaleSds = TRUE, log.initialPr = log(rep(1/length(states),
length(states))), normalIndex = 3L, prGtHom = c(1/3, 0.99, 0.7, 0.6,
0.6), prGtMis = rep(1/length(states), length(states)), prHetCalledHom =
0.001, prHetCalledHet = 0.995, prHomInNormal = 0.8, prHomInRoh = 0.999,
rohStates = logical(), tau = 1e+08, a2n = 1, n2a = 1, a2a = 1, verbose =
2L, ...)
hmm.setup(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class \code{SnpSet} or one of its derivatives.}
  \item{copynumberStates}{numeric. The mean values for each of the
    hidden states. Ignored if \code{object} is not a
    \code{CopyNumberSet}, \code{oligoSnpSet}, or \code{CNSet}.}
  \item{states}{Character vector indicating the names of the hidden states}
  \item{ICE}{Logical. Whether to use crlmm confidence scores estimates
    for computing the emission probabilities of the genotype calls.
    Only supported for some platforms. See details.}
  \item{is.log}{Logical. Whether the copy number estimates in
    \code{object} are on the log scale.  See details}
  \item{scaleSds}{Not implemented. [[Factor to scale the standard
    deviations used to compute the emission probabilities.]]}
  \item{log.initialPr}{Numeric vector of the same length as the number
    of states.  Specifies the initial state probabilities on the log scale.}
  \item{normalIndex}{Integer.  Indicates which element of the
    \code{states} vector corresponds to the 'normal' state.}
  \item{prGtHom}{numeric vector in [0,1],  Elements correspond to the
    probability of a homozygous genotype for each state. Ignored for
    \code{CopyNumberSet} objects.}
  \item{prGtMis}{numeric vector in [0,1],  Elements correspond to the
    probability of a missing genotype for each state. Ignored for
    \code{CopyNumberSet} objects.}
  \item{prHetCalledHom}{numeric in [0,1]. The probability of a heterozygous genotype incorrectly called homozygous. Ignored if
    object is a \code{CopyNumberSet} objects or if ICE is FALSE.}
  \item{prHetCalledHet}{numeric in [0,1],  The probability of a heterozygous genotype correctly called heterozygous. Ignored if
    object is a \code{CopyNumberSet} or if ICE is FALSE.}
  \item{prHomInNormal}{numeric in [0,1].  The probability of a homozygous genotype in a region with copy
    number two. Ignored if ICE is FALSE or if object is a \code{CopyNumberSet}.}
  \item{prHomInRoh}{numeric in [0,1].  The
    probability of observing a homozygous genotype when only homozygous
    genotypes are expected (e.g., hemizygous deletions, copy-neutral
    regions of homozygosity).  Ignored if ICE is FALSE or if object is a
    \code{CopyNumberSet}.}
  \item{rohStates}{integer vector. Indicates which states only
    homozygous genotypes are expected. Ignored if ICE is FALSE or if
    object is a \code{CopyNumberSet}.}
  \item{tau}{   Numeric. Factor for scaling the distance weighted transition
  probability.  See details.}
  \item{a2n}{Integer. Factor for scaling the probability of transitioning
      from an altered state to a normal state. Default is 1
      (recommended).}

  \item{n2a}{Integer. Factor for scaling the probability of transitioning
	      from the normal state to an altered state.  Default is 1
	      (recommended).}

  \item{a2a}{Integer.  Factor for scaling the probability of transitioning
      from an altered state to a different altered stateDefault is 1
      (recommended).}

  \item{verbose}{Logical.  Whether to display all messages and warnings.}

  \item{\dots}{For \code{hmm.setup}, additional arguments passed to
    \code{HmmOptionList} constructor.  Ignored for
    \code{HmmOptionList}.}
}

\details{

  Setting \code{ICE} to TRUE requires availability of the crlmm
  confidence scores for HapMap samples.  Currently only a few platforms
  are supported. See \code{VanillaICE:::icePlatforms()}.

  The transition probability is computed as
  \code{exp(-2 * d/tau)}, where d is the distance between two loci.
  The default value for \code{tau} \code{10^(8)}, but users can adjust this number to
  control the smoothness of the resulting HMM.

  For \code{object} classes that do not support copy number estimation
  (e.g., \code{SnpSet}), \code{is.log} is ignored. For all other
  classes, \code{is.log} is used to determine the plausible range of
  copy number estimates.  When \code{is.log} is TRUE, we assume the
  range of the copy number estimates is [-10, 2.5].  Otherwise, the
  range of the copy number estimates is assumed to be [0, 10].  Values
  more extreme are truncated.  The emission probabilities from the HMM
  are estimated by fitting a Gaussian-Uniform mixture.  The uniform
  component of the mixture model has a support indicated by the above
  ranges.  The Gaussian component has a mean for each state that is
  provided by the \code{copynumberStates} argument. For
  \code{oligoSnpSet} objects, the standard deviation of the Gaussian can
  be marker- and sample-specific if confidence scores are provided in
  the \code{cnConfidence} assay data element of the \code{oligoSnpSet}
  object.  Specifically, \code{1/cnConfidence} is used as the standard
  deviation.  If the \code{cnConfidence} slot is not specified, the
  median absolute deviation of the copy number estimates across
  autosomal chromosomes is used for the standard deviation of both
  autosomal and sex chromosome markers.

}
\value{
Object of class \code{HmmOptionList}.
}
\author{
Rob Scharpf
}

\seealso{
  \code{\linkS4class{HmmOptionList}}, \code{\linkS4class{oligoSnpSet}},
  \code{\linkS4class{SnpSet}}, \code{\linkS4class{SnpSet}},
  \code{\linkS4class{CopyNumberSet}}
}
\examples{
  if(require("crlmm")){
    data(sample.CNSet, package="crlmm")
    hmm.params <- HmmOptionList(sample.CNSet)
    hmm.params
  }
  data(oligoSetExample, package="oligoClasses")
  hmm.params <- HmmOptionList(oligoSet)
  hmm.params
}
\keyword{list}

