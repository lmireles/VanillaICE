HmmOptionList <- function(object,
			  copynumberStates,
			  states, ##=paste("state", 1:length(copynumberStates), sep=""),
			  ICE=FALSE,
			  is.log=FALSE,
			  scaleSds=TRUE,
			  log.initialPr=log(rep(1/length(states), length(states))),
			  normalIndex,
			  prGtHom, ##=c(1/3, 0.99, 0.7, 0.6, 0.6), ##only used when ICE=FALSE
			  prGtMis=rep(1/length(states), length(states)), ##not applicable when ICE is TRUE (no NA's from crlmm genotypes)
			  prHetCalledHom=0.001, ## ignored unless ICE is TRUE
			  prHetCalledHet=0.995, ## ignored unless ICE is TRUE
			  prHomInNormal=0.8,    ## ignored unless ICE is TRUE
			  prHomInRoh=0.999, ## ignored unless ICE is TRUE
			  rohStates, ## ignored unless ICE is TRUE
			  tau=1e8,
			  a2n=1,
			  n2a=1,
			  a2a=1,
			  verbose=2L, ...){
	## to support previous API
	if("prGenotypeHomozygous" %in% names(list(...)))
		prGtHom <- list(...)[["prGenotypeHomozygous"]]
	if(missing(copynumberStates)){
		if(is.log){
			copynumberStates <- switch(class(object),
						   CNSet=log2(c(0.2, 1, 2, 2, 3, 4)),
						   oligoSnpSet=log2(c(0.1,1,2,2,3,4)),
						   SnpSet=NA,
						   CopyNumberSet=log2(c(0.1,1,2,3,4)))
		} else {
			copynumberStates <- switch(class(object),
						   CNSet=c(0, 1, 2, 2, 3, 4),
						   oligoSnpSet=c(0,1,2,2,3,4),
						   SnpSet=NA,
						   CopyNumberSet=c(0,1,2,3,4))
		}
	}
	if(!ICE){
		prHetCalledHom <- prHetCalledHet <- prHomInNormal <- prHomInRoh <- rohStates <- NA
		if(missing(prGtHom)){
			prGtHom <- switch(class(object),
					  CNSet=c(2/3, 0.99, 0.7, 0.99, 0.7, 0.7),
					  oligoSnpSet=c(2/3, 0.99, 0.7, 0.99, 0.7, 0.7),
					  SnpSet=c(0.99, 0.7),
					  CopyNumberSet=NA)
		}
	} else{
		prGtMis <- prGtHom <- NA
		rohStates <- switch(class(object),
				    SnpSet=c(TRUE, FALSE),
				    oligoSnpSet=c(FALSE,TRUE,FALSE,TRUE,FALSE,FALSE),
				    CNSet=c(FALSE,TRUE,FALSE,TRUE,FALSE,FALSE))
	}
	if(missing(states)){
		states <- switch(class(object),
				 CNSet=c("homDel", "hemDel", "normal", "ROH", "3copy", "4copy"),
				 oligoSnpSet=c("homDel", "hemDel", "normal", "ROH", "3copy", "4copy"),
				 SnpSet=c("ROH", "normal"),
				 CopyNumberSet=c("homDel", "hemDel", "normal", "3copy", "4copy"))
		normalIndex <- grep("normal", states)
	}
	if(missing(normalIndex)){
		stop("'states' was specified by the user, but 'normalIndex' is missing")
	}
	stopifnot(normalIndex %in% seq_along(states))
	res <- list(snpsetClass=class(object),
		    copynumberStates=copynumberStates,
		    states=states,
		    ICE=ICE,
		    is.log=is.log,
		    scaleSds=scaleSds,
		    log.initialPr=log.initialPr,
		    normalIndex=normalIndex,
		    prGtHom=prGtHom,
		    prGtMis=prGtMis,
		    prHetCalledHom=prHetCalledHom,
		    prHetCalledHet=prHetCalledHet,
		    prHomInNormal=prHomInNormal,
		    prHomInRoh=prHomInRoh,
		    rohStates=rohStates,
		    tau=tau,
		    a2n=a2n,
		    n2a=n2a,
		    a2a=a2a,
		    verbose=verbose)
	hmm.opts <- as(res, "HmmOptionList")
	stopifnot(validObject(hmm.opts))
	return(hmm.opts)
}
setValidity("HmmOptionList", function(object){
	##ice <- ICE(object)
	ice <- object$ICE
	S <- length(object$states)
	if(!ice){
		if(object$snpsetClass != "CopyNumberSet"){
			check <- S == length(object$prGtHom)
			if(!check) return(FALSE)
			check <- S == length(object$prGtMis)
			if(!check) return(FALSE)
		}
		if(object$snpsetClass != "SnpSet"){
			check <- S == length(object$copynumberStates)
			if(!check) return(FALSE)
		}
	} else{
		if(!length(object$prHomInNormal) == 1) return(FALSE)
		if(!length(object$prHomInRoh) == 1) return(FALSE)
		if(!length(object$prHetCalledHet)==1) return(FALSE)
		if(!length(object$prHetCalledHom)==1) return(FALSE)
	}
	check <- S == length(object$log.initialPr)
	if(!check) return(FALSE)
	check <- sum(exp(object$log.initialPr)) == 1
	if(!check) return(FALSE)
	TRUE
})

setMethod("show", signature(object="HmmOptionList"),
	  function(object){
		  print(ls.str(object))
	  })



##setMethod("snpsetClass", "HmmOptionList", function(object) object$snpsetClass)
##setMethod("copynumberStates", "HmmOptionList", function(object) object$copynumberStates)
##setMethod("states", "HmmOptionList", function(object) object$states)
##setMethod("ICE", "HmmOptionList", function(object) object$ICE)
##setMethod("copyNumber", "HmmOptionList", function(object) object$copyNumber)
##setMethod("is.log", "HmmOptionList", function(object) object$is.log)
##setMethod("scaleSds", "HmmOptionList", function(object) object$scaleSds)
##setMethod("normalIndex", "HmmOptionList", function(object) object$normalIndex)
##setMethod("log.initialPr", "HmmOptionList", function(object) object$log.initialPr)
##setMethod("prGtHom", "HmmOptionList", function(object) object$prGtHom)
##setMethod("prGtMis", "HmmOptionList", function(object) object$prGtMis)
##setMethod("prHetCalledHom", "HmmOptionList", function(object) object$prHetCalledHom)
##setMethod("prHetCalledHet", "HmmOptionList", function(object) object$prHetCalledHet)
##setMethod("prHomInNormal", "HmmOptionList", function(object) object$prHomInNormal)
##setMethod("prHomInRoh", "HmmOptionList", function(object) object$prHomInRoh)
##setMethod("rohStates", "HmmOptionList", function(object) object$rohStates)
##setMethod("verbose", "HmmOptionList", function(object) object$verbose)
##setMethod("tau", "HmmOptionList", function(object) object$tau)
##setMethod("a2n", "HmmOptionList", function(object) object$a2n)
##setMethod("n2a", "HmmOptionList", function(object) object$n2a)
##setMethod("a2a", "HmmOptionList", function(object) object$a2a)
##setMethod("markerIndex", "HmmOptionList", function(object) object$marker.index)
##setMethod("sampleIndex", "HmmOptionList", function(object) object$sample.index)



##setMethod("initialize", signature(.Object="HmmOptionList"),
##	  function(.Object,
##		   snpsetClass,
##		   copynumberStates,
##		   states,
##		   ICE,
##		   copyNumber,
##		   is.log,
##		   scaleSds,
##		   log.initialPr,
##		   normalIndex,
##		   prGtHom,
##		   prGtMis,
##		   prHetCalledHom,
##		   prHetCalledHet,
##		   prHomInNormal,
##		   prHomInRoh,
##		   rohStates,
##		   verbose){
##		  ##browser()
##		  ##tmp <- callNextMethod()
##	if(missing(is.log)) stop("Must specify whether the copy number is on the log scale using the <is.log> argument.")
##	if(!snpsetClass %in% c("SnpSet", "CopyNumberSet", "oligoSnpSet")){
##		stop("class must be one of SnpSet, CopyNumberSet, or oligoSet")
##	}
##	if(snpsetClass == "SnpSet") copyNumber <- FALSE
##	stopifnot(is.numeric(normalIndex))
##	stopifnot(length(prGtHom) != length(states))
##	if(ICE) stop("not implemented")
##	.Object <- callNextMethod(snpsetClass=snpsetClass,
##				  copynumberStates=copynumberStates,
##				  states=states,
##				  ICE=ICE,
##				  copyNumber=copyNumber,
##				  is.log=is.log,
##				  scaleSds=scaleSds,
##				  log.initialPr=log.initialPr,
##				  normalIndex=normalIndex,
##				  prGtHom=prGtHom,
##				  prGtMis=prGtMis,
##				  prHetCalledHom=prHetCalledHom,
##				  prHetCalledHet=prHetCalledHet,
##				  prHomInNormal=prHomInNormal,
##				  prHomInRoh=prHomInRoh,
##				  rohStates=rohStates,
##				  verbose=verbose)
##})

