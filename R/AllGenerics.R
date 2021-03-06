##setGeneric("calculateEmission", function(object, ...) standardGeneric("calculateEmission"))

##setGeneric("computeEmission", function(object, hmmOptions) standardGeneric("computeEmission"))
##setGeneric("computeHmm", function(object, hmmOptions) standardGeneric("computeHmm"))
##setGeneric("plot", function(object, hmm.params, ...) standardGeneric("plot"))
setGeneric("hmm", function(object, hmm.params, use.baf=FALSE, k=5, ...) standardGeneric("hmm"))
setGeneric("hmm2", function(object, hmm.params, use.baf=FALSE, k=5,  ...) standardGeneric("hmm2"))
##setGeneric("snpsetClass", function(object) standardGeneric("snpsetClass"))
##setGeneric("scaleSds", function(object) standardGeneric("scaleSds"))
##setGeneric("emit", function(object, hmm.params) standardGeneric("emit"))
setGeneric("cnEmission", function(object, stdev, ...) standardGeneric("cnEmission"))
setGeneric("gtEmission", function(object, hmm.params, gt.conf, is.snp, cdfName, ...) standardGeneric("gtEmission"))
setGeneric("bafEmission", function(object, ...) standardGeneric("bafEmission"))

setGeneric("sd", useAsDefault=function(x, na.rm=FALSE) stats::sd(x, na.rm))
setGeneric("xyplot2", function(x, data, range, frame=50e3L, ...) standardGeneric("xyplot2"))
setGeneric("xyplot", useAsDefault=function(x, data, ...) lattice::xyplot(x, data, ...))



