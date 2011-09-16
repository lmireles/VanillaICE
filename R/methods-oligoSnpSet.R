setMethod("hmm2", signature(object="oligoSnpSet", hmm.params="HmmOptionList"),
	  function(object, hmm.params, ...){
		  verbose <- hmm.params[["verbose"]] > 0
		  log.beta.cn <- cnEmission(object, hmm.params, verbose=verbose, ...)
		  log.beta.gt <- gtEmission(object, hmm.params)
		  ##log.emission <- emit(object, hmm.params)

		  log.beta <- log.beta.cn+log.beta.gt
		  dimnames(log.beta) <- list(featureNames(object),
					     sampleNames(object),
					     hmm.params$states)
		  viterbi(object, hmm.params, log.E=log.beta)
})




setMethod("hmm", signature(object="oligoSnpSet", hmm.params="HmmOptionList"),
	  function(object, hmm.params, ...){
		  v2 <- hmm.params$verbose
		  naindex <- which(is.na(chromosome(object)) | is.na(position(object)))
		  if(length(naindex) > 0)
			  object <- object[-naindex, ]
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



setMethod("emit", signature(object="oligoSnpSet", hmm.params="HmmOptionList"),
	  function(object, hmm.params){
		  ICE <- hmm.params$ICE
		  ##EMIT.THR <- hmm.params[["EMIT.THR"]]
		  states <- hmm.params$states
		  verbose <- hmm.params$verbose
		  normalIndex <- hmm.params$normalIndex
		  if(all(is.na(cnConfidence(object)))){
			  message("cnConfidence missing.  Using MAD")
			  sds <- robustSds2(copyNumber(object))
			  cnConfidence(object) <- 1/sds
		  }
		  ##log.cn.emission <- calculateEmission.copynumber2(object,
		  ##                                                 hmm.params)
		  log.cn.emission <- calculateEmission(object, hmm.params)
		  if(!ICE){
			  log.gt.emission <- calculateEmission.genotype(object, hmm.params)
		  } else {
			  ##assumed order
			  ## ROH, normal
			  stop('need to update')
			  log.gt.emission <- array(NA, dim=c(nrow(object), ncol(object), length(states)),
						   dimnames=list(featureNames(object),
						   sampleNames(object),
						   states))
			  tmp <- genotypeEmissionCrlmm(object, hmm.params)
			  rohStates <- which(hmm.params[["rohStates"]])
			  notRohState <- which(!hmm.params[["rohStates"]])
			  for(j in rohStates){
				  log.gt.emission[, , j] <- tmp[, , "ROH"]
			  }
			  for(j in notRohState){
				  log.gt.emission[, , j] <- tmp[, , "normal"]
			  }
		  }
		  log.emission <- log.gt.emission+log.cn.emission
		  if(any(is.na(log.emission))){
			  if(verbose==2) message("Converting missing values in the emission matrix to 0")
			  log.emission[is.na(log.emission)] <- 0
		  }
		  return(log.emission)
})


setAs("oligoSnpSet", "data.frame",
      function(from, to){
	      browser()
	      cn <- copyNumber(from)
	      gt <- calls(from)
	      ##b <- baf(object)[marker.index, sample.index, ]
	      ##r <- logR(object)[marker.index, sample.index, ]
	      md <- mindist(object)[marker.index, sample.index]
##	      if(is.ff){
##		      close(baf(object))
##		      close(logR(object))
##		      close(mindist(object))
##	      }
##	      id <- matrix(c("father", "mother", "offspring"), nrow(b), ncol(b), byrow=TRUE)
##	      empty <- rep(NA, length(md))
	      ## A trick to add an extra panel for genes and cnv
	      ##df <- rbind(df, list(as.integer(NA), as.numeric(NA), as.numeric(NA), as.factor("genes")))
	      ## The NA's are to create extra panels (when needed for lattice plotting)
##	      id <- c(as.character(id), rep("min dist",length(md)))##, c("genes", "CNV"))
	      cn <- as.numeric(cn)
	      gt <- as.numeric(gt)
##	      b <- c(as.numeric(b), empty)
##	      r <- c(as.numeric(r), md)
##	      x <- rep(position(from),
##	      x <- rep(position(object)[marker.index], 4)/1e6
	      is.snp <- rep(isSnp(object)[marker.index], 4)
	      df <- data.frame(x=x, b=b, r=r, id=id, is.snp=is.snp)
	      df2 <- data.frame(id=c(as.character(df$id), "genes", "CNV"),
				b=c(df$b, NA, NA),
				r=c(df$r, NA, NA),
				x=c(df$x, NA, NA),
				is.snp=c(df$is.snp,NA, NA))
	      df2$id <- factor(df2$id, levels=c("father", "mother", "offspring", "min dist", "genes", "CNV"), ordered=TRUE)
	      return(df2)
      })


setMethod("xyplot", signature(x="formula", data="oligoSnpSet"),
	  function(x, data, ...){
		  stopifnot("range" %in% names(list(...)))
		  df <- as(data, "data.frame")
		  ##cn.df <- todf(oligoset, range, )
		  xyplot(x, df, ...)
})
