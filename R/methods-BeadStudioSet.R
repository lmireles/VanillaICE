


setMethod("xyplot", signature(x="formula", data="BeadStudioSet"),
	  function(x, data, ...){
		  if("range" %in% names(list(...))){
			  xyplot2(x, data, ...)
		  } else {
			  callNextMethod()
		  }
})
