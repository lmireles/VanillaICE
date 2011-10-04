setMethod("cnEmission", signature(object="matrix", stdev="matrix"),
	  function(object, stdev, k=5, cnStates, is.log, is.snp, ...){
		  stopifnot(length(cnStates) > 1)
		  stopifnot(is.numeric(cnStates))
		  stopifnot(all(dim(object) == dim(stdev)))
		  if(any(colSums(is.na(object)) == nrow(object))){
			  stop("Some samples have all missing values. Exclude these samples before continuing.")
		  }
		  S <- length(cnStates)
		  emission.cn <- array(NA, dim=c(nrow(object), ncol(object), S))
		  if(is.log){
			  MIN.CN <- -10
			  MAX.CN <- 2.5
		  } else {
			  MIN.CN <- 0
			  MAX.CN <- 10
		  }
		  if(any(object < MIN.CN, na.rm=TRUE)) object[object < MIN.CN] <- MIN.CN
		  if(any(object > MAX.CN, na.rm=TRUE)) object[object > MAX.CN] <- MAX.CN
		  for(j in 1:ncol(object)){
			  snp.index <- which(is.snp)
			  if(length(snp.index) > 0){
				  cn <- object[snp.index, j, drop=FALSE]
				  s <- stdev[snp.index, j, drop=FALSE]
				  mu.snp <- updateMu(x=cn, mu=cnStates, sigma=s)
				  I <- which(!is.na(as.numeric(cn)))
				  old.tmp <- tmp <- rep(NA, length(as.numeric(cnStates)))
				  cnvector <- as.numeric(cn)[I]
				  prOutlier <- probabilityOutlier(cnvector, k=k)
				  for(l in seq_along(cnStates)){
					  tmp <- (1-prOutlier) * dnorm(x=cnvector,
								       mean=mu.snp[l],##cnStates[l],
								       sd=as.numeric(s)[I]) +
									       prOutlier * dunif(cnvector, MIN.CN, MAX.CN)
					  emission.cn[snp.index[I], j, l] <- tmp
				  }
			  }  ## length(snp.index) > 0
			  np.index <- which(!is.snp)
			  if(length(np.index) > 0){
				  cn <- object[np.index, j, drop=FALSE]
				  s <- stdev[np.index, j, drop=FALSE]
				  mu.np <- updateMu(x=cn, mu=cnStates, sigma=s)
				  I <- which(!is.na(as.numeric(cn)))
				  old.tmp <- tmp <- rep(NA, length(as.numeric(cnStates)))
				  cnvector <- as.numeric(cn)[I]
				  prOutlier <- probabilityOutlier(cnvector, k=k)
				  for(l in seq_along(cnStates)){
					  tmp <- (1-prOutlier) * dnorm(x=cnvector,
								       mean=mu.np[l],##cnStates[l],
								       sd=as.numeric(s)[I]) +
									       prOutlier * dunif(cnvector, MIN.CN, MAX.CN)
					  emission.cn[np.index[I], j, l] <- tmp
				  }
			  }
		  }
		  return(log(emission.cn))
	  })



setMethod("gtEmission", signature(object="matrix"),
	  function(object, hmmOptions, gt.conf, is.snp, cdfName, ...){
		  ICE <- hmmOptions$ICE
		  states <- hmmOptions$states
		  if(!ICE){
			  p <- hmmOptions$prGtHom
			  prGenotypeMissing <- hmmOptions$prGtMis
			  verbose <- hmmOptions$verbose
			  stopifnot(length(p) == length(states))
			  if(!is.numeric(object)) stop("genotypes must be integers (1=AA, 2=AB, 3=BB) or NA (missing)")
			  emission <- array(object, dim=c(nrow(object), ncol(object), length(states))) ##dimnames=list(featureNames(object), sampleNames(object), states))
			  missingGT <- any(is.na(object))
			  for(s in seq(along=states)){
				  tmp <- object
				  tmp[tmp == 1 | tmp == 3] <- p[s]
				  tmp[tmp == 2] <- 1-p[s]
				  index1 <- is.na(tmp) & !is.snp##!isSnp(object)
				  index2 <- is.na(tmp) & is.snp##isSnp(object)
				  if(missingGT){
					  tmp[index2] <- prGenotypeMissing[s]
					  tmp[index1] <- 1/length(states)
				  }
				  emission[, , s] <- tmp
			  }
			  logemit <- log(emission)
			  ##return(logemit)
		  } else {
			  not.valid <- invalidGtConfidence(gt.conf)
			  if(any(not.valid)){
				  stop("Invalid genotype confidence scores.\n",
				       "\tAll confidence scores must be between 0 and 1")

			  }
			  ##stop('need to update ICE option')
			  logemit <- array(NA, dim=c(nrow(object), ncol(object), length(states)))
			  tmp <- genotypeEmissionCrlmm(object, hmmOptions=hmmOptions, gt.conf=gt.conf, cdfName=cdfName)
			  rohStates <- which(hmmOptions[["rohStates"]])
			  notRohState <- which(!hmmOptions[["rohStates"]])
			  if(length(rohStates) > 0){
				  logemit[, , rohStates] <- tmp[, , "ROH"]
			  }
			  if(length(notRohState) > 0){
				  logemit[, , notRohState] <- tmp[, , "normal"]
			  }
		  }
		  return(logemit)
	  })
