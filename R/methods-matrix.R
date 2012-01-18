setMethod("cnEmission", signature(object="matrix"),
	  function(object, stdev, k=5, cnStates, is.log, is.snp,
		   normalIndex, verbose=TRUE, chrom, ...){
		  cnEmissionFromMatrix(object=object, stdev=stdev,
				       k=k,
				       cnStates=cnStates,
				       is.log=is.log,
				       is.snp=is.snp,
				       normalIndex=normalIndex,
				       verbose=verbose,
				       chrom=chrom,
				       ...)
	  })

cnEmissionFromMatrix <- function(object, stdev, k=5, cnStates,
				 is.log, is.snp, normalIndex,
				 verbose=TRUE,
				 chrom){
	stopifnot(length(cnStates) > 1)
	stopifnot(is.numeric(cnStates))
	if(missing(stdev)){
		stdev <- .getSds(object)
		## use robust estimate of sample sd
		##s <- apply(object, 2, mad, na.rm=TRUE)
		##stdev <- matrix(s, nrow(object), ncol(object), byrow=TRUE)
	} else if(is(stdev, "matrix")) stopifnot(ncol(stdev) == ncol(object))
	stopifnot(all(dim(object) == dim(stdev)))
	if(any(colSums(is.na(object)) == nrow(object))){
		stop("Some samples have all missing values. Exclude these samples before continuing.")
	}
	S <- length(cnStates)
	emission.cn <- array(NA, dim=c(nrow(object), ncol(object), S))
	rr <- range(object, na.rm=TRUE, finite=TRUE)
	if(is.log){
		MIN.CN <- pmax(-10, rr[1])
		MAX.CN <- pmin(2.5, rr[2])
	} else {
		MIN.CN <- pmax(0, rr[1])
		MAX.CN <- pmin(10, rr[2])
	}
	object <- pmin(object, MAX.CN)
	object <- pmax(object, MIN.CN)
	for(j in seq_len(ncol(object))){
		##snp.index <- which(is.snp)
		cn <- object[, j]
		snp.index <- which(is.snp & !is.na(cn))
		if(length(snp.index) > 0){
			cn <- cn[snp.index]
			if(is(stdev, "matrix")){
				s <- stdev[snp.index, j]
			} else s <- stdev[snp.index]
			mu.snp <- updateMu(x=cn, mu=cnStates, sigma=s, normalIndex=normalIndex)
			old.tmp <- tmp <- rep(NA, length(as.numeric(cnStates)))
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
				if(is(stdev, "matrix")){
					s <- stdev[np.index, j]
				} else s <- stdev[np.index, j]
				mu.np <- updateMu(x=cn, mu=cnStates, sigma=s, normalIndex=normalIndex)
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
}



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
	any.nan <- any(is.nan(sds))
	any.na <- any(is.na(sds))
	allvalid <- !any.nan & !any.na
	stopifnot(allvalid)
	return(sds)
}

setMethod("bafEmission", signature(object="matrix"),
	  function(object, is.snp, prOutlier=1e-3, p.hom=0.95, ...){
		  bafEmissionFromMatrix(object=object,
					is.snp=is.snp,
					prOutlier=prOutlier,
					p.hom=p.hom,...)
	  })

bafEmissionFromMatrix <- function(object, is.snp, prOutlier=1e-3, p.hom=0.95, ...){
	states <- 1:6
	S <- 6
	if("pb" %in% names(list(...))){
		pb <- list(...)[["pb"]]
		pb <- pb/100
		pb[is.na(pb)] <- 0.5
	} else pb <- rep(0.5, nrow(object))
	emission <- array(NA, dim=c(nrow(object), ncol(object), S))
	## mixture of 2 truncated normals and 1 normal.
	##  Assume mu is known, but sigma is not.
	if(ncol(object) > 1 | sum(is.snp) < 1000){
		sds <- updateSigma(object, is.snp, sigma0=c(0.02, 0.04, 0.02))
	} else {
		tmp <- kmeans(object[is.snp, ], centers=c(0,0.5,1))
		sds <- sapply(split(object[is.snp, ], tmp$cluster), mad, na.rm=TRUE)
	}
	minsd <- min(sds)
	if(minsd > 0.02) {
		tmp <- kmeans(object[is.snp, ], centers=c(0,1/4, 1/3, 0.5,2/3, 3/4, 1))
		sds <- sapply(split(object[is.snp, ], tmp$cluster), mad, na.rm=TRUE)
		prOutlier <- prOutlier*10*minsd/.02
	}
	if("scalesd" %in% names(list(...))){
		scalesd <- list(...)[["scalesd"]]
		sds <- sds*scalesd
	}
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
	p.out <- prOutlier
	q.out <- 1-p.out
	if(is(is.snp, "numeric")) is.snp <- as.logical(is.snp)
	i <- which(is.snp)
	obj <- object[i, , drop=FALSE]
	is.hom <- obj < 0.05 | obj > 0.95
	##obj2 <- obj
	##obj <- 0.0052
	TN0 <- tnorm(obj, 0, sd0)
	TN1 <- tnorm(obj, 1, sd1)
	## bound the truncated normals. outliers should be handled by the uniform component
	##		  TN.25 <- tnorm(obj, 0.25, sd.5, upper=0.45, lower=0.05)
	##		  TN.5 <- tnorm(obj, 0.5, sd.5, upper=0.75, lower=0.25)
	##		  TN.75 <- tnorm(obj, 0.75, sd.5, upper=0.95, lower=0.55)
	TN.3 <- tnorm(obj, 1/3, sd.5, lower=0.05, upper=0.45)
	TN.6 <- tnorm(obj, 2/3, sd.5, lower=0.55, upper=0.95)
	TN.25 <- tnorm(obj, 0.25, sd.5, lower=0.05, upper=0.45)
	TN.5 <- tnorm(obj, 0.5, sd.5, lower=0.05, upper=0.95)
	TN.75 <- tnorm(obj, 0.75, sd.5, lower=0.05, upper=0.45)
	p <- pb[i]
	if(any(is.na(p)))
		p[is.na(p)] <- 0.5
	pr2 <- (1-p)*TN0 + p*TN1
	beta.hemizygous <- q.out*pr2 + p.out
	emission[i, , 1] <- p.out + q.out*dunif(obj, 0, 1) ## 0
	emission[i, , 2] <- beta.hemizygous     ## 1
	emission[i, , 4] <- beta.hemizygous  ## 2, ROH
	## if pb is small and genotype is AA
	##   emission for hemiz. del and normal is close
	## if pb is 1/2 and genotype is AA
	##   emission for hemiz. del roughly 2 times normal
	##   -do we want this to be the case?...
	##   -if we force the probs for hom. be the same,
	##    allowing only the prob to het to differ the state would always be normal
	##  Want the probs for hom to differ only by a small amount
	##    such that only long stretches of roh overcome the transition pr.
	q2 <- (1-p)^2
	pq <- p*(1-p)
	p2 <- p^2
	emit.normal <- p.out + q.out*(q2*TN0 + 2*pq*TN.5 + p2*TN1) ## 2
	q3 <- (1-p)^3
	p2q <- p^2*(1-p)
	pq2 <- p*(1-p)^2
	p3 <- p^3
	emit5 <- p.out + q.out*(q3*TN0 + 3*pq2*TN.3 + 3*p2q*TN.6 + p3*TN1) ## 3
	q4 <- (1-p)^4
	pq3 <- p*(1-p)^3
	p2q2 <- p^2*(1-p)^2
	p3q <- p^3*(1-p)
	p4 <- p^4
	emit6 <- p.out + q.out*(q4*TN0 + 4*pq3*tnorm(obj, 1/4, sd.5) + 6*p2q2*TN.5 + 4*p3q*tnorm(obj, 0.75, sd.5) + p4*TN1)
	## homozygous genotypes are less informative
	for(j in seq_len(ncol(obj))){
		i.hom <- which(is.hom[,j])
		##emit.normal[i.hom, j] <- (1-p.hom)*beta.hemizygous[i.hom, j]+p.hom*emit.normal[i.hom, j]
		emit.normal[i.hom, j] <- (1-p.hom)*beta.hemizygous[i.hom, j]+p.hom*emit.normal[i.hom, j]
		emit5[i.hom, j] <- (1-p.hom)*beta.hemizygous[i.hom, j]+p.hom*emit5[i.hom, j]
		emit6[i.hom, j] <- (1-p.hom)*beta.hemizygous[i.hom, j]+p.hom*emit6[i.hom, j]
	}
	emission[i, , 3] <- emit.normal
	emission[i, , 5] <- emit5
	emission[i, , 6] <- emit6
	emit <- emission[is.snp, 1, ]
	emit <- cbind(round(emit, 2), obj)
	colnames(emit)[ncol(emit)] <- "baf"
	## assign 1 as the emission probablily for all nonpolymorphic markers
	## (across all states)
	np.index <- which(!is.snp)
	emission[np.index, , ] <- 1
	logemit <- log(emission)
	return(logemit)
}
