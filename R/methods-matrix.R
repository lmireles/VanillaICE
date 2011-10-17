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
			  ##snp.index <- which(is.snp)
			  cn <- object[, j]
			  snp.index <- which(is.snp & !is.na(cn))
			  cn <- cn[snp.index]
			  if(length(snp.index) > 0){
				  ##cn <- object[snp.index, j, drop=FALSE]
				  s <- stdev[snp.index, j]
				  mu.snp <- updateMu(x=cn, mu=cnStates, sigma=s)
				  ##I <- which(!is.na(as.numeric(cn)))
				  ##if(length(I)==0) next()
				  old.tmp <- tmp <- rep(NA, length(as.numeric(cnStates)))
				  ##cnvector <- as.numeric(cn)##[I]
				  ##prOutlier <- as.numeric(probabilityOutlier(cnvector, k=k))
				  prOutlier <- probabilityOutlier(cn, k=k)
				  for(l in seq_along(cnStates)){
					  e <- (1-prOutlier) * dnorm(x=cn, mean=mu.snp[l], sd=s) + prOutlier * dunif(cn, MIN.CN, MAX.CN)
					  emission.cn[snp.index, j, l] <- e
				  }
			  }
		  }
		  is.np <- !is.snp
		  if(any(is.np)){
			  for(j in 1:ncol(object)){
				  np.index <- which(is.np & !is.na(object)[, j])
				  if(length(np.index) > 0){
					  cn <- object[np.index, j]
					  s <- stdev[np.index, j]
					  mu.np <- updateMu(x=cn, mu=cnStates, sigma=s)
					  ##old.tmp <- tmp <- rep(NA, length(as.numeric(cnStates)))
					  prOutlier <- probabilityOutlier(cn, k=k)
					  for(l in seq_along(cnStates)){
						  tmp <- (1-prOutlier) * dnorm(x=cn, mean=mu.np[l], sd=s) + prOutlier * dunif(cn, MIN.CN, MAX.CN)
						  emission.cn[np.index, j, l] <- tmp
					  }
				  }
			  }
		  }
		  logemit <- log(emission.cn)
		  return(logemit)
	  })



setMethod("gtEmission", signature(object="matrix"),
	  function(object, hmm.params, gt.conf, is.snp, cdfName, ...){
		  ICE <- hmm.params$ICE
		  states <- hmm.params$states
		  if(!ICE){
			  p <- hmm.params$prGtHom
			  prGenotypeMissing <- hmm.params$prGtMis
			  verbose <- hmm.params$verbose
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
			  tmp <- genotypeEmissionCrlmm(object, hmm.params=hmm.params, gt.conf=gt.conf, cdfName=cdfName)
			  rohStates <- which(hmm.params[["rohStates"]])
			  notRohState <- which(!hmm.params[["rohStates"]])
			  if(length(rohStates) > 0){
				  logemit[, , rohStates] <- tmp[, , "ROH"]
			  }
			  if(length(notRohState) > 0){
				  logemit[, , notRohState] <- tmp[, , "normal"]
			  }
		  }
		  return(logemit)
	  })

##---------------------------------------------------------------------------
##
## From MinimumDistance
##
##---------------------------------------------------------------------------
tnorm <- function(x, mean, sd, lower=0, upper=1){
	phi <- function(x, mu, sigma) dnorm(x, mu, sigma)
	## cdf of standard normal
	Phi <- function(x, mu, sigma) pnorm(x, mu, sigma)
	res <- phi(x, mean, sd)/(Phi(upper, mean, sd)-Phi(lower, mean, sd))
	ind <- which(x < lower | x > upper)
	if(any(ind)){
		res[ind] <- 0
	}
	res
}
##mosaicProb <- function(bf, sd.mosaic, sd0, sd.5, sd1){
##	initialP <- 0.99## initial probabilty not mosaic
##	iter <- 1
##	##ri=range.index(object)
####	while(any(is.na(rangeIndex))){
####		rangeIndex <- fillInMissing(rangeIndex)
####		iter <- iter+1
####		if(iter > 5) break()
####	}
##	f1 <- tnorm(bf, 0, sd0)
##	f2 <- tnorm(bf, 1, sd1)
##	f3 <- tnorm(bf, 0.25, sd.mosaic, lower=0.05, upper=0.4)
##	f4 <- tnorm(bf, 0.75, sd.mosaic, lower=0.6, upper=0.95)
##	f5 <- tnorm(bf, 0.5, sd.5)
##	tau <- matrix(initialP, nrow(f1), ncol(f1))
##	##LL <- sapply(split(rangeIndex, rangeIndex), length)
##	p.out <- 1e-10
##	## updating f3,f4 might help
##	for(i in 1:10){
##		T1.num <- tau * ((1-p.out)*(1/3*f1 + 1/3*f2 + 1/3*f5) + p.out)
##		T.den <- tau * ((1-p.out)*(1/3*f1 + 1/3*f2 + 1/3*f5) + p.out) + (1-tau)*((1-p.out)*(0.25*f1 + 0.25*f3 + 0.25*f4 + 0.25*f2) + p.out)
##		T1 <- T1.num/T.den
##		tau.next <- apply(T1, 2, mean, na.rm=TRUE)
####		tau.F <- mean(T1[,1], na.rm=TRUE)
####		tau.M <- mean(T1[,2], na.rm=TRUE)
####		tau.O <- sapply(split(T1[, 3], rangeIndex), mean, na.rm=TRUE)
####		tau.F <- rep(tau.F, LL)
####		tau.M <- rep(tau.M, LL)
####		tau.O <- rep(tau.O, LL)
####		tau.next <- cbind(tau.F, tau.M, tau.O)
##		if(abs(sum(tau.next - tau, na.rm=TRUE)) < 1) break()
##		tau <- tau.next
##		##tau.next <- apply(T1, 2, mean, na.rm=TRUE)
##		##tau <- tau.next
##	}
##	return(tau)
##}

updateSigma <- function(x, is.snp, nUpdates=10, sigma0){
	##x <- x[is.snp & x > 0 & x < 1]
	x <- x[is.snp]
	mu <- c(0, 0.5, 1)
	L <- length(mu)
	pi <- rep(1/L, L)
	gamma <- matrix(NA, length(x), L)
	den <- matrix(NA, length(x), L)
	##num <- vector("list", L-1)
	num <- matrix(NA, length(x), L)
	epsilon <- 0.01
	sd.new <- rep(NA,L)
	sigma <- sigma0
	for(iter in seq_len(nUpdates)){
		if(iter > nUpdates) break()
		den[, 1] <- pi[1]*tnorm(x, mean=mu[1], sigma[1])
		den[, 2] <- pi[2]*tnorm(x, mean=mu[2], sigma[2])
		den[, 3] <- pi[3]*tnorm(x, mean=mu[3], sigma[3])
		D <- rowSums(den, na.rm=TRUE)
		for(i in seq_len(L)){
			num <- den[, i] ##pi[i] * dnorm(x, mu[i], sigma)
			gamma[, i] <- num/D
		}
		rs <- rowSums(gamma, na.rm=TRUE)
		gamma <- gamma/rs
		##
		## update the sds
		##
		##sigma2 <- rep(NA,L)
		for(i in seq_len(L)){
			sd.new[i] <- sqrt(sum(gamma[,i]*(x - mu[i])^2,na.rm=TRUE)/(sum(gamma[,i],na.rm=TRUE)))
		}
		##total.gamma <- apply(gamma, 2, sum, na.rm=TRUE)
		pi.new <- apply(gamma, 2, mean, na.rm=TRUE)
		pi <- pi.new
		##dp <- abs(sum(mu - mu.new)) + abs(sum(pi.new-pi))
		dp <- abs(sum(sigma-sd.new,na.rm=TRUE))
		sigma <- sd.new
		if(dp < epsilon) break()
	}
	sds <- sigma
##	is.homozygote <- gamma[, 1] > 0.9 | gamma[, 3] > 0.9
##	index <- which(!is.homozygote)
##	if(length(index) > 1000){
##		x <- x[index]
##		mu <- c(1/3, 0.5, 2/3)
##		sigma0 <- rep(sigma[2], 3)
##		for(iter in seq_len(nUpdates)){
##			if(iter > nUpdates) break()
##			den[, 1] <- pi[1]*tnorm(x, mean=mu[1], sigma[1])
##			den[, 2] <- pi[2]*tnorm(x, mean=mu[2], sigma[2])
##			den[, 3] <- pi[3]*tnorm(x, mean=mu[3], sigma[3])
##			D <- rowSums(den, na.rm=TRUE)
##			for(i in seq_len(L)){
##				num <- den[, i] ##pi[i] * dnorm(x, mu[i], sigma)
##				gamma[, i] <- num/D
##			}
##			rs <- rowSums(gamma, na.rm=TRUE)
##			gamma <- gamma/rs
##			##
##			## update the sds
##			##
##			##sigma2 <- rep(NA,L)
##			for(i in seq_len(L)){
##				sd.new[i] <- sqrt(sum(gamma[,i]*(x - mu[i])^2,na.rm=TRUE)/(sum(gamma[,i],na.rm=TRUE)))
##			}
##			pi.new <- apply(gamma, 2, mean, na.rm=TRUE)
##			pi <- pi.new
##			##dp <- abs(sum(mu - mu.new)) + abs(sum(pi.new-pi))
##			dp <- abs(sum(sigma-sd.new,na.rm=TRUE))
##			sigma <- sd.new
##			if(dp < epsilon) break()
##		} ##for(iter in seq_len(nUpdates))
##		sds <- c(sds[1], sigma, sds[3])
##	} else {## if(length(index) <= 1000)
##		sds <- c(sds[1], sds[2], sds[2], sds[3])
##	}
	return(sds)
}

setMethod("bafEmission", signature(object="matrix"),
	  function(object, hmm.params, is.snp, cdfName, ...){
		  S <- length(hmm.params[["states"]])
		  emission <- array(NA, dim=c(nrow(object), ncol(object), S))
		  ## mixture of 2 truncated normals and 1 normal.
		  ##  Assume mu is known, but sigma is not.
		  sds <- updateSigma(object, is.snp, sigma0=c(0.02, 0.04, 0.02))
		  sd0 <- sds[1]
		  sd1 <- sds[3]
		  sd.5 <- sds[2]
		  ## one way p.out can help is if the state appears
		  ## normal except for a few outliers. Here, bigger
		  ## values for p.out should encourage more 'normal'
		  ## states'.
		  ##
		  ## The downside is that a few hets appearing in a
		  ## region of homozygosity may be more acceptable,
		  ## even though we might prefer a bigger penalty.
		  ##
		  ## Can we estimate p? ... the genotype confidence
		  ## scores may be helpful
		  p.out <- 1e-8
		  q.out <- 1-p.out
		  i <- which(is.snp)
		  obj <- object[i, , drop=FALSE]
		  TN0 <- tnorm(obj, 0, sd0)
		  TN1 <- tnorm(obj, 1, sd1)
		  TN.3 <- tnorm(obj, 1/3, sd.5)
		  TN.6 <- tnorm(obj, 2/3, sd.5)
		  TN.25 <- tnorm(obj, 0.25, sd.5);
		  TN.5 <- tnorm(obj, 0.5, sd.5);
		  TN.75 <- tnorm(obj, 0.75, sd.5)
		  pr2 <- 0.5*TN0 + 0.5*TN1
		  beta.hemizygous <- q.out*pr2 + p.out
		  emission[i, , 1] <- p.out + q.out*dunif(obj, 0, 1) ## 0
		  emission[i, , 2] <- beta.hemizygous     ## 1
		  emission[i, , 3] <- p.out + q.out*(1/3*TN0 + 1/3*TN.5 + 1/3*TN1) ## 2
		  emission[i, , 4] <- beta.hemizygous  ## 2, ROH
		  emission[i, , 5] <- p.out + q.out*(1/4*TN0 + 1/4*TN.3 + 1/4*TN.6 + 1/4*TN1) ## 3
		  emission[i, , 6] <- p.out + q.out*(1/5*TN0 + 1/5*tnorm(obj, 1/4, sd.5) + 1/5*TN.5 + 1/5*tnorm(obj, 0.75, sd.5) + 1/5*TN1)
		  ## assign 1 as the emission probablily for all nonpolymorphic markers
		  ## (across all states)
		  np.index <- which(!is.snp)
		  emission[np.index, , ] <- 1
		  logemit <- log(emission)
		  return(logemit)
	  })
