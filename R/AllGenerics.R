##setGeneric("calculateEmission", function(object, ...) standardGeneric("calculateEmission"))

##setGeneric("computeEmission", function(object, hmmOptions) standardGeneric("computeEmission"))
##setGeneric("computeHmm", function(object, hmmOptions) standardGeneric("computeHmm"))
##setGeneric("plot", function(object, hmm.params, ...) standardGeneric("plot"))
setGeneric("hmm", function(object, hmm.params, ...) standardGeneric("hmm"))
setGeneric("hmm2", function(object, hmm.params, ...) standardGeneric("hmm2"))
##setGeneric("snpsetClass", function(object) standardGeneric("snpsetClass"))
##setGeneric("scaleSds", function(object) standardGeneric("scaleSds"))
##setGeneric("emit", function(object, hmm.params) standardGeneric("emit"))
setGeneric("cnEmission", function(object, stdev, k=5, cnStates,
				  is.log, is.snp, verbose=TRUE, ...) standardGeneric("cnEmission"))
setGeneric("gtEmission", function(object, hmm.params, gt.conf, is.snp, cdfName, ...) standardGeneric("gtEmission"))
setGeneric("bafEmission", function(object, hmm.params, is.snp, cdfName, ...) standardGeneric("bafEmission"))

setGeneric("sd", useAsDefault=function(x, na.rm=FALSE) stats::sd(x, na.rm))
setGeneric("xyplot2", function(x, data, range, frame=50e3L, ...) standardGeneric("xyplot2"))
setGeneric("xyplot", useAsDefault=function(x, data, ...) lattice::xyplot(x, data, ...))

setGeneric("lrr<-", function(object,value) standardGeneric("lrr<-"))
setGeneric("baf<-", function(object,value) standardGeneric("baf<-"))
