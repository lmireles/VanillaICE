setClass("HmmOptionList", contains="list")
setMethod("initialize", "HmmOptionList",
	  function(.Object, ...){
		  .Object <- callNextMethod()
	  })

##setClass("DataFrameCN", representation(row.names="character",
##				       names="character"),
##	 contains="list")
##setMethod("dimnames", "DataFrameCN", function(x) list(x@row.names, names(x)))
##setMethod("dim", "DataFrameCN", function(x) c(length(x@row.names), length(x)))


