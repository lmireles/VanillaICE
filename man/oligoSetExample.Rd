\name{oligoSet}
\alias{oligoSet}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(oligoSet)}
\format{

}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
	Created from the simulated locusLevelData provided in this package.
}
\examples{
	\dontrun{
		## 'oligoSetExample' created by the following
		data(locusLevelData)
		oligoSet <- new("oligoSnpSet",
				copyNumber=log2(locusLevelData[["copynumber"]]/100),
				call=locusLevelData[["genotypes"]],
				callProbability=locusLevelData[["crlmmConfidence"]],
				annotation=locusLevelData[["platform"]])
		oligoSet <- oligoSet[!is.na(chromosome(oligoSet)), ]
		oligoSet <- oligoSet[chromosome(oligoSet) < 3, ]
	}
	data(oligoSet)
	oligoSet
}
\keyword{datasets}
