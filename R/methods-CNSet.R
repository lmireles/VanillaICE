##setMethod("hmm2", signature(object="CNSet", hmm.params="HmmOptionList"),
##	  function(object, hmm.params, ...){
##		  sample.index <- sampleIndex(hmm.params)
##		  marker.index <- markerIndex(hmm.params)
##		  v <- verbose(hmm.params)
##		  if(is.null(sample.index)){
##			  if(v > 0) message("sample.index is NULL.  Using all samples")
##			  sample.index <- seq(length=ncol(object))
##		  }
##		  sample.index.list <- splitIndicesByLength(sample.index, ocSamples())
##		  if(is.null(marker.index)){
##			  if(v > 0) message("marker.index is NULL.  Using all markers")
##			  marker.index <- seq(length=nrow(object))
##		  }
##		  marker.index.list <- split(marker.index, chromosome(object)[marker.index])
##		  oligoClasses:::open(object$SNR)
##		  snr <- object$SNR[]
##		  if(any(snr < 5)) warning(sum(snr < 5), " samples have SNR < 5.  May want to exclude these samples by specifying sample.index.")
##		  oligoClasses:::open(object$SNR)
##		  rm(snr)
##		  fit <- vector("list", length(marker.index.list)*length(sample.index.list))
##		  k <- 1
##		  oligoClasses:::open(object)
##		  for(i in seq_along(marker.index.list)){
##			  ii <- marker.index.list[[i]]
##			  for(j in seq_along(sample.index.list)){
##				  jj <- sample.index.list[[j]]
##				  cnset <- object[ii, jj]
##				  cnset <- cnset[order(position(cnset)), ]
##				  oligoSet <- as(cnset, "oligoSnpSet")
##				  oligoSet <- centerCopyNumber(oligoSet, ...)
##				  sds <- robustSds2(copyNumber(oligoSet), ...)
##				  cnConfidence(oligoSet) <- 1/sds
##				  rm(sds, cnset); gc()
##				  fit[[k]] <- hmm(oligoSet, hmm.params, ...)
##				  k <- k+1
##			  }
##		  }
##		  oligoClasses:::close(object)
##		  if(length(fit) > 1){
##			  rd <- RangedDataList(fit)
##			  rd2 <- stack(rd)
##			  rm(fit,rd); gc()
##			  ix <- match("sample", colnames(rd2))
##			  if(length(ix) > 0) rd2 <- rd2[, -ix]
##		  } else rd2 <- fit[[1]]
##		  return(rd2)
##	  })
setMethod("hmm", signature(object="CNSet", hmm.params="HmmOptionList"),
	  function(object, hmm.params, ...){
		  if("verbose" %in% names(list(...))){
			  verbose <- list(...)[["verbose"]]
		  } else verbose <- TRUE
		  if("by.chromosome" %in% names(list(...))){
			  by.chromosome <- list(...)[["by.chromosome"]]
		  } else by.chromosome <- TRUE
		  if("k" %in% names(list(...))){
			  k <- list(...)[["k"]]
		  } else k <- 5
		  if("sample.index" %in% names(list(...))){
			  sample.index <- list(...)[["sample.index"]]
		  } else sample.index <- seq(length=ncol(object))
		  batch.index <- split(sample.index, batch(object)[sample.index])
		  hmm.params[["verbose"]] <- 0L
		  chrom <- unique(chromosome(object))
		  NN <- length(chrom) * length(batch.index)
		  autosome.index <- which(chromosome(object) < 23)
		  if(verbose){
			  pb <- txtProgressBar(min=0, max=NN, style=3)
		  }
		  results <- vector("list", NN)
		  m <- 1
		  for(i in seq_along(batch.index)){
			  for(j in seq_along(chrom)){
				  CHR <- chrom[j]
				  is.autosome <- CHR < 23
				  J <- batch.index[[i]]
				  I <- which(chromosome(object) == CHR)
				  cnset.batch <- object[I, J]
				  oligoSet <- as(cnset.batch, "oligoSnpSet")
				  oligoSet <- oligoSet[order(position(oligoSet)), ]
				  rm(cnset.batch)
				  oligoSet <- centerAutosomesAt(oligoSet, at=2)
				  hmmOpts <- HmmOptionList(object=oligoSet, verbose=0L)
				  results[[m]] <- hmm(oligoSet, hmmOpts, k=k)
				  if(verbose) setTxtProgressBar(pb, m)
				  m <- m+1
			  }
		  }
		  if(verbose) close(pb)
		  tmp <- stack(RangedDataList(results))
		  index <- match("sample", colnames(tmp))
		  if(length(index) == 1) tmp <- tmp[, -index]
		  rangedData <- tmp
		  return(rangedData)
	  })

setMethod("xyplot", signature(x="formula", data="CNSet"),
	  function(x, data, ...){
		  stopifnot("range" %in% names(list(...)))
		  rd <- list(...)[["range"]]
		  if(!"frame" %in% names(list(...))){
			  w <- width(rd)
			  frame <- w/0.05  * 1/2
		  } else {
			  frame <- list(...)[["frame"]]
		  }
		  marker.index <- featuresInRange(data, rd)
		  sample.index <- match(sampleNames(data), sampleNames(rd))
		  cnset <- cnset[marker.index, sample.index]
		  oligoset <- as(cnset, "oligoSnpSet")
		  xyplot(x, oligoset, ...)
})
