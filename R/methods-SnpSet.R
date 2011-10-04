setMethod("hmm2", signature(object="SnpSet", hmm.params="HmmOptionList"),
	  function(object, hmm.params, use.baf=FALSE, ...){
		  ##log.beta <- gtEmission(object, hmm.params)
		  if(use.baf){
			  if(!"baf" %in% ls(assayData(object))){
				  stop("use.baf is true, but baf not in assayData.  See calculateRBaf.")
			  }
			  log.beta.gt <- bafEmission(object, hmm.params, ...)
		  } else {
			  log.beta.gt <- gtEmission(object, hmm.params, ...)
		  }
		  dimnames(log.beta.gt) <- list(NULL,
						sampleNames(object),
						hmm.params$states)
		  viterbi(object, hmm.params, log.E=log.beta.gt)
})

setMethod("hmm", signature(object="SnpSet", hmm.params="HmmOptionList"),
	  function(object, hmm.params, use.baf=FALSE, ...){
		  naindex <- which(is.na(chromosome(object)) | is.na(position(object)))
		  if(length(naindex) > 0)
			  object <- object[-naindex, ]
		  v2 <- hmm.params$verbose
		  naindex <- which(is.na(chromosome(object)) | is.na(position(object)))
		  ice <- hmm.params$ICE
		  if(ice) checkAnnotation(object)
		  if(length(naindex) > 0){
			  warning("NA's in annotation -- check chromosome or position. Ignoring markers with missing annotation")
			  object <- object[-naindex, ]
		  }
		  missingGT <- any(is.na(calls(object)))
		  if(missingGT){
			  if(v2 > 0) message("Some genotypes are NAs.  The default assumes that prGenotypeMissing is the same for each state -- see hmmOptions")
		  }
		  marker.index.list <- split(seq(length=nrow(object)), chromosome(object))
		  chromosomes <- unique(chromosome(object))
		  ix <- order(chromosome(object), position(object))
		  if(any(diff(ix) < 0)) {
			  object <- object[ix, ]
		  }
		  if(is.null(hmm.params$sample.index)){
			  sample.index <- seq(length=ncol(object))
		  } else sample.index <- hmm.params$sample.index
		  if(v2 > 0) {
			  pb <- txtProgressBar(min=0, max=ncol(object), style=3)
		  }
		  res <- vector("list", ncol(object))
		  for(j in seq_along(sample.index)){
			  jj <- sample.index[j]
			  tmp <- vector("list", length(chromosomes))
			  for(k in seq_along(chromosomes)){
				  CHR <- chromosomes[k]
				  i <- marker.index.list[[k]]
				  obj <- object[i, jj]
				  tmp[[k]] <- hmm2(object=obj, hmm.params=hmm.params, use.baf=use.baf, ...)
			  }
			  if(length(tmp) > 1){
				  rdlist <- RangedDataList(tmp)
				  rd <- stack(rdlist)
				  ix <- match("sample", colnames(rd))
				  if(length(ix) > 0) rd <- rd[, -ix]
				  rm(rdlist)
			  } else rd <- tmp[[1]]
			  res[[j]] <- rd
			  rm(tmp, rd); gc()
			  if(v2 > 0) setTxtProgressBar(pb, j)
		  }
		  if(v2 > 0) close(pb)
		  if(length(res) > 1){
			  ##rdlist <- lapply(res, function(x) as(x, "RangedData"))
			  rdlist <- RangedDataList(res)
			  rd <- stack(rdlist)
			  ix <- match("sample", colnames(rd))
			  if(length(ix) > 0) rd <- rd[, -ix]
			  rm(rdlist)
		  } else rd <- res[[1]]
		  return(rd)
	  })

setMethod("gtEmission", signature(object="SnpSet"),
	  function(object, hmm.params, ...){
		  is.ordered <- checkOrder(object)
		  stopifnot(is.ordered)
		  log.emit <- gtEmission(calls(object), hmm.params,
					 is.snp=isSnp(object),
					 gt.conf=confs(object),
					 cdfName=annotation(object), ...)
		  return(log.emit)
	  })

setMethod("bafEmission", signature(object="SnpSet"),
	  function(object, hmm.params, ...){
		  is.ordered <- checkOrder(object)
		  stopifnot(is.ordered)
		  log.emit <- bafEmission(baf(object), hmm.params,
					  is.snp=isSnp(object),
					  cdfName=annotation(object), ...)
		  return(log.emit)
	  })

setMethod("xyplot", signature(x="formula", data="SnpSet"),
	  function(x, data, ...){
		  if("range" %in% names(list(...))){
			  xyplot2(x, data, ...)
		  } else {
			  callNextMethod()
		  }
})

setMethod("baf", signature(object="SnpSet"), function(object) assayDataElement(object, "baf"))

