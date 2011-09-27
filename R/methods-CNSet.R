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
			  message("Estimating standard deviation (s) using the median absolute deviation\n",
				  "of the total copy number across autosomal markers.\n",
				  "Emission probabilities estimated using Uniform-Gaussian mixture.\n",
				  "Gaussian component for sample j is N(mu, s_j).")
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

setMethod("xyplot2", signature(x="formula", data="CNSet", range="RangedDataCNV"),
	  function(x, data, range, frame=50e3L, ...){
		  mm <- findOverlaps(range, data, frame=frame)
		  mm.df <- data.frame(mm)
		  mm.df$featureNames <- featureNames(data)[mm.df$subject]
		  marker.index <- unique(mm.df$subject)
		  ##marker.index <- featuresInRange(data, rd, FRAME=frame)
		  sample.index <- match(sampleNames(range), sampleNames(data))
		  data <- data[marker.index, sample.index]
		  ## now we need to know the indices of
		  ## each range after subsetting
		  ## mm.df$subject <- match(mm.df$featureNames, featureNames(data))
		  ## we assume that each range is from a different sample
		  oligoset <- as(data, "oligoSnpSet")
		  df <- as(oligoset, "data.frame")
		  xyplot(x, df,
			 range=range,
			 gt=df$gt,
			 is.snp=df$is.snp,
			 ...)
	  })


##setMethod("xyplot", signature(x="formula", data="CNSet"),
##	  function(x, data, ...){
##		  xyplot2(x, data, ...)
##})
