setMethod("hmm2", signature(object="CopyNumberSet", hmm.params="HmmOptionList"),
	  function(object, hmm.params, ...){
		  verbose <- hmm.params[["verbose"]] > 0
		  log.beta <- cnEmission(object, hmm.params, verbose=verbose, ...)
		  dimnames(log.beta) <- list(featureNames(object),
					     sampleNames(object),
					     hmm.params$states)
		  res <- viterbi(object, hmm.params, log.E=log.beta)
		  return(res)
})

setMethod("hmm", signature(object="CopyNumberSet", hmm.params="HmmOptionList"),
	  function(object, hmm.params, ...){
		  v2 <- hmm.params$verbose
		  naindex <- which(is.na(chromosome(object)) | is.na(position(object)))
		  if(length(naindex) > 0)
			  object <- object[-naindex, ]
		  if(v2 > 0){
			  if("k" %in% names(list(...))){
				  k <- list(...)[["k"]]
			  } else k <- 3
			  message("Using a running median of ", k, " markers to estimate the outlier probability.")
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
				  tmp[[k]] <- hmm2(object=obj, hmm.params=hmm.params, ...)
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
setMethod("sd", signature(x="CopyNumberSet"),
	  function(x, na.rm=FALSE){
		  getSds(x, na.rm=TRUE)
	   })
