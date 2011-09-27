\name{gtEmission}
\alias{gtEmission}
\title{Calculate emission probabilities for genotypes}
\description{
  Calculate emission probabilities for genotypes diallelic genotypes.
}
\usage{
gtEmission(object, hmmOptions, gt.conf, is.snp, cdfName, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    A object of class \code{SnpSet}, \code{oligoSnpSet}, or \code{matrix}.
  }

  \item{hmmOptions}{
    A \code{HmmOptionList}.
  }

  \item{gt.conf}{Ignored unless object is of class \code{matrix}.}

  \item{is.snp}{Ignored unless object is of class \code{matrix}.}

  \item{cdfName}{Ignored unless the argument \code{ICE} for the
    \code{HmmOptionList} object is TRUE.}

  \item{\dots}{Ignored}
}

\details{

  Currently, there are two main approaches for calculating the emission
  probabilities for diallelic genotypes.  Which approach is implemented
  depends on the value of \code{ICE} in the \code{HmmOptionList} object.

  When \code{ICE} is FALSE (the default), the emission probabilities for
  the diallelic genotypes are estimated from a Bernoulli(p_s) with
  \code{p_s} denoting the probability of a homozygous genotype (`AA' or
  `BB') for each state \code{s}. These probabilities can be specified in
  the constructor for the \code{HmmOptionList} class by the argument
  \code{prGtHom}. The corresponding probability for a heterozygous
  genotype (`AB') is 1-p_s for state \code{s}.

  For many calling algorithms, the genotypes are not called.  If many
  no-calls occur in a sample, this can indicate problems with the DNA
  quality.  No calls can also arise when there is poor separation of the
  genotype clusters or when the A and B allele intensities for a sample
  do not cluster into any of the diallelic genotype clusters.  No calls
  should be indicated by the value \code{NA}. We estimate the emission
  probability of a no-call using a Bernoulli(p2_s) with \code{p2_s}
  denoting the probability of a no call for state \code{s}. The numeric
  vector for \code{p2} can be specified in the constructor for the
  \code{HmmOptionList} class through the argument \code{prGtMis}.

  When \code{ICE} is TRUE, we assume that the genotype calls were
  obtained from the \pkg{crlmm}.  This option is only available for a
  few of the platforms that \pkg{crlmm} supports.  The \pkg{crlmm}
  provides an estimate of the confidence for each diallelic genotype
  call.  For the HapMap dataset, we assessed the distribution of the
  confidence scores when the call was correct versus the distribution
  when the call was incorrect (using HapMap genotypes as the gold
  standard).  Using these distributions, we estimate the probability
  that the true diallelic genotype lies in a region of homozygosity
  (possibly suggesting a hemizygous deletion or a region of homozygosity
  induced by uniparental disomy) and the probability that the region has
  a 'normal' proportion of heterozygous genotypes. The constructor
  \code{HmmOptionList} has an argument \code{rohStates} that indicates
  which of the hidden states we expect homozygosity. See the reference
  below for additional details regarding the esimation of emission
  probabilities using the ICE option.


}
\value{
  An array.  The dimensions are features x samples x states.
}
\author{
R.Scharpf
}
\references{Scharpf, RB et al., 2008, Annals of Applied Statistics}
\note{
Only chromosomes 1-24 supported (23=X and 24 = Y).
}

\seealso{
  \code{\link{cnEmission}}, \code{\link{gtEmission-methods}}
}
\examples{
data(oligoSetExample, package="oligoClasses")
oligoSet <- order(oligoSet)
##
## Here, the probability of a missing genotype is 5 times as likely
## for a homozygous-deletion than any of the other states
##
hmmOpts <- HmmOptionList(oligoSet, is.log=TRUE, ICE=FALSE,
			 copynumberStates=0:4,
			 prGtMis=c(1/2, 0.1, 0.1, 0.1, 0.1),
			 prGtHom=c(0.5, 0.99, 0.7, 0.6, 0.6))
gt.emit <- gtEmission(oligoSet, hmmOpts)

## ICE is TRUE
hmmOpts <- HmmOptionList(oligoSet, is.log=TRUE, ICE=TRUE,
			 copynumberStates=c(0,1,2,2,3,4),
			 rohStates=c(FALSE,TRUE,TRUE,FALSE,FALSE),
			 normalIndex=4,
			 prHetCalledHom=0.001,
			 prHetCalledHet=0.995,
			 prHomInNormal=0.7,
			 prHomInRoh=0.999)
tryCatch(gtEmission(oligoSet, hmmOpts), error=function(e) "Annotation platform not supported. See icePlatforms()")
annotation(oligoSet) %in% icePlatforms()

## for illustration
annotation(oligoSet) <- "genomewidesnp6"
gt.emit2 <- gtEmission(oligoSet, hmmOpts)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
\keyword{array}
