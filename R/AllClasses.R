setClass("HmmOptionList", contains="list")
setMethod("initialize", "HmmOptionList",
	  function(.Object, ...){
		  .Object <- callNextMethod()
	  })


