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
		  chrom <- unique(chromosome(object, na.rm=TRUE))
		  chrom <- chrom[chrom <= 23]
		  NN <- length(chrom) * length(batch.index)
		  autosome.index <- which(chromosome(object) < 23)
		  if(verbose){
			  message("Estimating standard deviation (s) using the median absolute deviation\n",
				  "of the total copy number across autosomal markers.\n")
##				  "Emission probabilities estimated using Uniform-Gaussian mixture.\n",
##				  "Gaussian component for sample j is N(mu, s_j).")
			  pb <- txtProgressBar(min=0, max=NN, style=3)
		  }
		  results <- vector("list", NN)
		  m <- 1
		  for(i in seq_along(batch.index)){
			  J <- batch.index[[i]]
			  if(by.chromosome){
				  for(j in seq_along(chrom)){
					  CHR <- chrom[j]
					  is.autosome <- CHR < 23
					  I <- which(chromosome(object) == CHR)
					  cnset.batch <- object[I, J]
					  oligoSet <- as(cnset.batch, "oligoSnpSet")
					  is.ordered <- checkOrder(oligoSet)
					  if(!is.ordered)
						  oligoSet <- chromosomePositionOrder(oligoSet)
					  rm(cnset.batch)
					  ##oligoSet <- centerAutosomesAt(oligoSet, at=2)
					  ##hmmOpts <- HmmOptionList(object=oligoSet, verbose=0L)
					  hmm.params$verbose <- 0L
					  results[[m]] <- hmm(oligoSet, hmm.params, k=k)
					  if(verbose) setTxtProgressBar(pb, m)
					  m <- m+1
				  }
			  } else {
				  cnset.batch <- object[, J]
				  oligoSet <- as(cnset.batch, "oligoSnpSet")
				  is.ordered <- checkOrder(oligoSet)
				  if(!is.ordered)
					  oligoSet <- chromosomePositionOrder(oligoSet)
				  rm(cnset.batch)
				  hmm.params$verbose <- 0L
				  results[[m]] <- hmm(oligoSet, hmm.params, k=k)
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


