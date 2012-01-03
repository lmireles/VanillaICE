\name{hmm3}
\alias{hmm3}

\title{
  HMM for processed data.
}


\description{

  HMM for preprocessed and normalized data. Expects files containing log
  R ratios and B allele frequencies that have already been
  wave-corrected, etc.

}

\usage{
hmm3(filenames, cdfname, universe=c("hg18", "hg19", ""), lrr.colname = "Log.R.Ratio", baf.colname = "B.Allele", samplesPerProcess = 1L, colClasses, TAUP = 1e+08, medianWindow = 5, cnStates = c(-1.5, -0.5, 0, 0, 0.4, 0.8), prOutlier = 0.001, p.hom = 0.6, ...)
}

\arguments{
  \item{filenames}{
    character string of filenames:  must include the complete path such
  that \code{file.exists(filenames)} evaluates to \code{TRUE}.
}

\item{cdfname}{

  character string identifying annotation package.

}

  \item{universe}{
    character string specifying which genome build. Valid options are
    hg18 and hg19.
}


  \item{lrr.colname}{

    character: used to grep for the log R ratios in the header. E.g.,
    \code{grep(lrr.colname, header) } should return a length 1 vector,
    where header is a vector of the column labels.

  }

  \item{baf.colname}{

    character: used to grep for the B allele frequency in the
    header. E.g., \code{grep(baf.colname, header) } should return a
    length 1 vector, where header is a vector of the column labels.

  }

  \item{samplesPerProcess}{

    Number of files per process when parallel computing is enabled.
    This feature is not currently implemented and the value is ignored.

  }

  \item{colClasses}{
    See \code{read.table}.
}

\item{TAUP}{

    Scalar for computing transition probabilities as a function of
    distance.  Large values make the probability of transitioning
    between two different states smaller.  The formula for computing the
    transition probability is

    \code{exp (-2*diff(position)/TAUP)},

    where position is the physical position on the genome.

  }

  \item{medianWindow}{

    The window size for a running median smoother.  When estimating the
    probability that a log R ratio is an outlier, a running median
    smoother with is fit to the log R ratios.  A two-component Gaussian
    mixture model is fitted to the residuals of the smoother.  The
    mixture probability for the 'high variance' distribution is the
    probability that the observation is an outlier.

  }


  \item{cnStates}{

    A numeric vector.  The log R ratios are asssumed to be emitted from
    a Gaussian distribution with parameters conditional on the latent
    copy number state.  The cnStates indicates the mean of the copy
    number estimates (mean log R ratios) for each hidden state. The
    latent copy number states for which the values of cnStates
    correspond are homozygous deletion (total CN=0), hemizygous deletion
    (total CN=1), normal (total CN=2), region of homozygosity (total
    CN=2) single-copy duplication (total CN=2), double-copy duplication
    (total CN=4).

  }

  \item{prOutlier}{

    The probability that the observed B Allele frequency is an outlier
    (or equivalently, the probability of a mistake in the
    genotyping). Currently, the value is assumed to be the same for all
    markers.

  }

  \item{p.hom}{

    A scalar (probability) in the interval [0, 1].  This parameter
    controls the desired size of called copy neutral regions of
    homozgosity (state '4').  In particular, small values of
    \code{p.hom} (e.g., less than 0.5), has the effect of making the
    emission probability for the normal copy number state more similar
    to state '4' for homozygous regions of the genome.  The result of
    which is that longer tracts of homozygosity are required before the
    HMM will transition to state '4'.  As an extreme example, setting
    \code{p.hom = 0} will likely have the effect that no copy-neutral
    regions of homozygosity are called.

  }


  \item{\dots}{

    Currently ignored.

  }
}

%\details{
%
%
%}


\value{

  A object of class \code{RangedDataHMM}.

}

\author{
R. Scharpf
}

\seealso{
\code{\linkS4class{RangedDataHMM}}
}

\examples{
if(require("foreach")){
	path <- system.file("extdata", package="VanillaICE")
	filename <- list.files(path, pattern="LRRandBAF", full.names=TRUE)
	rd <- hmm3(filename, cdfname="gw6crlmm", universe="hg18")
}
}
\keyword{smooth}
\keyword{models}
