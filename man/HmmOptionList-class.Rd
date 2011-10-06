\name{HmmOptionList-class}
\Rdversion{1.1}
\docType{class}
\alias{HmmOptionList-class}
\alias{hmm,CNSet,HmmOptionList-method}
\alias{hmm,CopyNumberSet,HmmOptionList-method}
\alias{hmm,oligoSnpSet,HmmOptionList-method}
\alias{hmm,SnpSet,HmmOptionList-method}

\title{A container for storing options for the hidden Markov model}
\description{
  A container for storing options for the hidden Markov model.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{HmmOptionList(object)}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
Class \code{"\linkS4class{AssayData}"}, by class "list", distance 2.
Class \code{"\linkS4class{list_or_ffdf}"}, by class "list", distance 2.
Class \code{vectorORfactor}, by class "list", distance 3.
}
\section{Methods}{
  \describe{
    \item{hmm}{\code{signature(object = "CNSet", hmm.params = "HmmOptionList")}: ... }
    \item{hmm}{\code{signature(object = "CopyNumberSet", hmm.params = "HmmOptionList")}: ... }
    \item{hmm}{\code{signature(object = "oligoSnpSet", hmm.params = "HmmOptionList")}: ... }
    \item{hmm}{\code{signature(object = "SnpSet", hmm.params = "HmmOptionList")}: ... }
	 }
}
\author{
R. Scharpf
}

\examples{
showClass("HmmOptionList")
}
\keyword{classes}
