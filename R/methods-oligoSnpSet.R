setMethod("hmm2", signature(object="oligoSnpSet", hmm.params="HmmOptionList"),
	  function(object, hmm.params, use.baf=FALSE, ...){
		  verbose <- hmm.params[["verbose"]] > 0
		  log.beta.cn <- cnEmission(object,
					    cnStates=hmm.params[["copynumberStates"]],
					    is.log=hmm.params[["is.log"]],
					    is.snp=isSnp(object), ...)
		  if(use.baf){
			  if(!"baf" %in% ls(assayData(object))){
				  stop("use.baf is true, but baf not in assayData.  See calculateRBaf.")
			  }
			  log.beta.gt <- bafEmission(object, hmm.params, ...)
		  } else {
			  log.beta.gt <- gtEmission(object, hmm.params, ...)
		  }
		  ##log.emission <- emit(object, hmm.params)
		  log.beta <- log.beta.cn+log.beta.gt
		  dimnames(log.beta) <- NULL
		  dimnames(log.beta) <- list(NULL,
					     sampleNames(object),
					     hmm.params$states)
		  res <- viterbi(object, hmm.params, log.E=log.beta)
		  return(res)
})




setMethod("hmm", signature(object="oligoSnpSet", hmm.params="HmmOptionList"),
	  function(object, hmm.params, use.baf=FALSE, ...){
		  v2 <- hmm.params$verbose
		  marker.index <- validChromosomeIndex(object)
		  object <- object[marker.index, ]
		  ice <- hmm.params$ICE
		  if(ice) checkAnnotation(object)
		  missingGT <- any(is.na(calls(object)))
		  if(missingGT){
			  if(v2 > 0) message("Some genotypes are NAs.  The default assumes that prGenotypeMissing is the same for each state -- see hmmOptions")
		  }
		  if(v2 > 0){
			  if("k" %in% names(list(...))){
				  k <- list(...)[["k"]]
			  } else k <- 3
			  message("Using a running median of ", k, " markers to estimate the outlier probability.")
		  }
		  is.ordered <- checkOrder(object)
		  if(!is.ordered) object <- order(object)
		  marker.index.list <- split(seq(length=nrow(object)), chromosome(object))
		  chromosomes <- unique(chromosome(object))
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



##setMethod("emit", signature(object="oligoSnpSet", hmm.params="HmmOptionList"),
##	  function(object, hmm.params){
##		  ICE <- hmm.params$ICE
##		  ##EMIT.THR <- hmm.params[["EMIT.THR"]]
##		  states <- hmm.params$states
##		  verbose <- hmm.params$verbose
##		  normalIndex <- hmm.params$normalIndex
##		  if(all(is.na(cnConfidence(object)))){
##			  message("cnConfidence missing.  Using MAD")
##			  sds <- robustSds2(copyNumber(object))
##			  cnConfidence(object) <- 1/sds
##		  }
##		  ##log.cn.emission <- calculateEmission.copynumber2(object,
##		  ##                                                 hmm.params)
##		  log.cn.emission <- calculateEmission(object, hmm.params)
##		  if(!ICE){
##			  log.gt.emission <- calculateEmission.genotype(object, hmm.params)
##		  } else {
##			  ##assumed order
##			  ## ROH, normal
##			  stop('need to update')
##			  log.gt.emission <- array(NA, dim=c(nrow(object), ncol(object), length(states)),
##						   dimnames=list(featureNames(object),
##						   sampleNames(object),
##						   states))
##			  tmp <- genotypeEmissionCrlmm(object, hmm.params)
##			  rohStates <- which(hmm.params[["rohStates"]])
##			  notRohState <- which(!hmm.params[["rohStates"]])
##			  for(j in rohStates){
##				  log.gt.emission[, , j] <- tmp[, , "ROH"]
##			  }
##			  for(j in notRohState){
##				  log.gt.emission[, , j] <- tmp[, , "normal"]
##			  }
##		  }
##		  log.emission <- log.gt.emission+log.cn.emission
##		  if(any(is.na(log.emission))){
##			  if(verbose==2) message("Converting missing values in the emission matrix to 0")
##			  log.emission[is.na(log.emission)] <- 0
##		  }
##		  return(log.emission)
##})


setMethod("sd", signature(x="oligoSnpSet"),
	  function(x, na.rm=FALSE){
		  getSds(x, na.rm=TRUE)
	   })

setMethod("xyplot2", signature(x="formula",
			       data="eSet",
			       range="RangedDataCNV"),
	  function(x, data, range, frame=0L, ...){
		  ## for now
		  ##if(nrow(range) > 1) frame <- 0L
		  rm <- findOverlaps(range, data, maxgap=frame) ## RangesMatching
		  mm <- matchMatrix(rm)
		  mm.df <- data.frame(mm)
		  mm.df$featureNames <- featureNames(data)[mm.df$subject]
		  marker.index <- mm.df$subject
		  sample.index <- match(sampleNames(range), sampleNames(data))
		  sample.index <- unique(sample.index)
		  data <- data[marker.index, sample.index]
		  mm.df$subject <- match(mm.df$featureNames, featureNames(data))
		  df <- as(data, "data.frame")
		  ##tmp <- todataframe(data, mm.df)
		  list.x <- as.character(x)
		  i <- grep("|", list.x, fixed=TRUE)
		  if(length(i) > 0){
			  zvar <- list.x[[i]]
			  zvar <- strsplit(zvar, " | ", fixed=T)[[1]][[2]]
			  if(zvar == "range"){
				  tmp <- tryCatch(df$range <- mm.df$query, error=function(e) NULL)
			  }
		  }
		  if("gt" %in% colnames(df)){
			  xyplot(x, df,
				 range=range,
				 gt=df$gt,
				 is.snp=df$is.snp,
				 ...)
		  } else {
			  xyplot(x, df,
				 range=range,
				 is.snp=df$is.snp,
				 ...)
		  }
	  })


setMethod("cnEmission", signature(object="oligoSnpSet"),
	  function(object, stdev, k=5, cnStates, is.log, is.snp, ...){
		  ##fn <- featureNames(object)
		  is.ordered <- checkOrder(object)
		  stopifnot(is.ordered)
		  CN <- copyNumber(object)
		  sds <- sd(object)
		  emit <- cnEmission(object=CN, stdev=sds,
				     k=k, cnStates=cnStates,
				     is.log=is.log, is.snp=is.snp, ...)
		  return(emit)
	  })

setMethod("sd", signature(x="oligoSnpSet"),
	  function(x, na.rm=FALSE){
		  getSds(x, na.rm=TRUE)
	   })
