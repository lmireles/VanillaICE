setClass("HmmOptionList", contains="list")
setMethod("initialize", "HmmOptionList",
	  function(.Object, ...){
		  .Object <- callNextMethod()
	  })


setClass("BeadStudioSet", contains="eSet")##,
##	 	 prototype = prototype(
##		 new("VersionedBiobase",
##		     versions=c(classVersion("BeadStudioSet"), BeadStudioSet="0.0.1"))))
setMethod("initialize", "BeadStudioSet",
	  function(.Object,
		   assayData = assayDataNew(baf = baf, lrr = lrr, ...),
		   phenoData = annotatedDataFrameFrom(assayData, byrow=FALSE),
		   featureData = annotatedDataFrameFrom(assayData, byrow=TRUE),
		   experimentData = new("MIAME"),
		   annotation = character(),
		   protocolData = phenoData[,integer(0)],
		   baf = new("matrix"),
		   lrr = matrix(numeric(),
                		   nrow=nrow(baf),
		                   ncol=ncol(baf),
                    		   dimnames=dimnames(baf)),
		   ...) {
	.Object <- callNextMethod(.Object,
				  assayData = assayData,
				  phenoData = phenoData,
				  featureData = featureData,
				  experimentData = experimentData,
				  annotation = annotation,
				  protocolData = protocolData)
	return(.Object)
})


setValidity("BeadStudioSet", function(object) {
	return(all(is.element(c("lrr","baf"), assayDataElementNames(object))))
})

setMethod("lrr", "BeadStudioSet", function(object)
	  assayDataElement(object, "lrr"))
setReplaceMethod("lrr", c("BeadStudioSet", "ANY"),
		 function(object, value) {
			 assayDataElementReplace(object, "lrr", value)
	 })
setMethod("baf", "BeadStudioSet",
	  function(object) {
		  assayDataElement(object, "baf")
	 })
setReplaceMethod("baf", c("BeadStudioSet", "ANY"),
		 function(object, value) {
			 assayDataElementReplace(object, "BAF", value)
	 })



##setClass("DataFrameCN", representation(row.names="character",
##				       names="character"),
##	 contains="list")
##setMethod("dimnames", "DataFrameCN", function(x) list(x@row.names, names(x)))
##setMethod("dim", "DataFrameCN", function(x) c(length(x@row.names), length(x)))


