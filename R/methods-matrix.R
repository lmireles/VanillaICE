setMethod("cnEmission", signature(object="matrix"),
	  function(object, stdev, k=5, cnStates, is.log, is.snp,
		   normalIndex, prOutlierCN=0.01, verbose=TRUE, ...){
		  cnEmissionFromMatrix(object=object, stdev=stdev,
				       k=k,
				       cnStates=cnStates,
				       is.log=is.log,
				       is.snp=is.snp,
				       normalIndex=normalIndex,
				       prOutlier=prOutlierCN,
				       verbose=verbose,
				       ...)
	  })

cnEmissionFromMatrix <- function(object, stdev, k=5, cnStates,
				 is.log, is.snp, normalIndex,
				 prOutlier,
				 center.by.type=TRUE,
				 verbose=TRUE, ...){
	stopifnot(length(cnStates) > 1)
	stopifnot(is.numeric(cnStates))
	if(center.by.type){
		## assume that the polymorphic markers and nonpolymorphic markers have the same median copynumber / log r ratio
		snp.index <- which(is.snp)
		mu.snp <- apply(object[snp.index, , drop=FALSE], 2, median)
		object[snp.index, ] <- sweep(object[snp.index, , drop=FALSE],  2, mu.snp)
		if(any(!is.snp)){
			np.index <- which(!is.snp)
			mu.np <- apply(object[np.index, , drop=FALSE], 2, median)
			object[np.index, ] <- sweep(object[np.index, , drop=FALSE], 2, mu.np)
		}
		object <- object+cnStates[normalIndex]
	}
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
			##mu.snp <- updateMu(x=cn, mu=cnStates, sigma=s, normalIndex=normalIndex)
			mu.snp=cnStates
			##old.tmp <- tmp <- rep(NA, length(as.numeric(cnStates)))
			##
			## If the variance of the altered states is much bigger than the variance of the normal
			## state, prOutlier will be near 1 using the approach in the probabilityOutlier function...
			##
			##prOutlier <- probabilityOutlier(cn, k=k)
			for(l in 2:length(cnStates)){
				e <- (1-prOutlier) * dnorm(x=cn, mean=mu.snp[l], sd=s) + prOutlier * dunif(cn, MIN.CN, MAX.CN)
				emission.cn[snp.index, j, l] <- e
			}
			if(is.log){
				emission.cn[snp.index, j, 1] <- (1-prOutlier)*dunif(cn, MIN.CN, -1) + prOutlier*dunif(cn, MIN.CN, MAX.CN)
			} else emission.cn[snp.index, j, 1] <- (1-prOutlier)*dunif(cn, MIN.CN, 1) + prOutlier*dunif(cn, MIN.CN, MAX.CN)
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
				##mu.np <- updateMu(x=cn, mu=cnStates, sigma=s, normalIndex=normalIndex)
				mu.np <- cnStates
				##old.tmp <- tmp <- rep(NA, length(as.numeric(cnStates)))
				##prOutlier <- probabilityOutlier(cn, k=k)
				##for(l in seq_along(cnStates)){
				for(l in 2:length(cnStates)){
					tmp <- (1-prOutlier) * dnorm(x=cn, mean=mu.np[l], sd=s) + prOutlier * dunif(cn, MIN.CN, MAX.CN)
					emission.cn[np.index, j, l] <- tmp
				}
				if(is.log){
					emission.cn[np.index, j, 1] <- (1-prOutlier)*dunif(x=cn, MIN.CN, -1) + prOutlier*dunif(cn, MIN.CN, MAX.CN)
				} else {
					emission.cn[np.index, j, 1] <- (1-prOutlier)*dunif(x=cn, MIN.CN, 1) + prOutlier*dunif(cn, MIN.CN, MAX.CN)
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
	x <- x[x > 0 & x < 1]
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
		den[, 1] <- pi[1]*tnorm(x, mean=mu[1], sigma[1], lower=0)
		den[, 2] <- pi[2]*tnorm(x, mean=mu[2], sigma[2])
		den[, 3] <- pi[3]*tnorm(x, mean=mu[3], sigma[3], upper=1)
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
	  function(object, is.snp, prOutlierBAF=1e-3, p.hom=0.95, ...){
		  bafEmissionFromMatrix(object=object,
					is.snp=is.snp,
					prOutlier=prOutlierBAF,
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
	##if(ncol(object) > 1 | sum(is.snp) < 1000){
	##  Assume mu is known, but sigma is not.
	##sds <- updateSigma(object, is.snp, sigma0=c(0.02, 0.04, 0.02))
	sds <- apply(object, 2, updateSigma, is.snp=is.snp, sigma0=c(0.02,0.04, 0.02))
##	} else {
##		tmp <- kmeans(object[is.snp, ], centers=c(0,0.5,1))
##		sds <- sapply(split(object[is.snp, ], tmp$cluster), mad, na.rm=TRUE)
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
	TN.5 <- tnorm(obj, 0.5, sd.5)
	TN.3 <- tnorm(obj, 1/3, sd0)
	TN.6 <- tnorm(obj, 2/3, sd0)
	TN.25 <- tnorm(obj, 0.25, sd0)
	TN.3 <- tnorm(obj, 1/3, sd0)
	TN.6 <- tnorm(obj, 2/3, sd0)
	TN.25 <- tnorm(obj, 0.25, sd0)
	TN.75 <- tnorm(obj, 0.75, sd0)
	p <- pb[i]
	if(any(is.na(p)))
		p[is.na(p)] <- 0.5
	pr2 <- (1-p)*TN0 + p*TN1
	beta.hemizygous <- (1-p.out/1000)*pr2 + p.out/1000
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
	## small p.hom makes homozygous genotypes less informative
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
	##emit <- cbind(round(emit, 2), obj)
	colnames(emit)[ncol(emit)] <- "baf"
	## assign 1 as the emission probablily for all nonpolymorphic markers
	## (across all states)
	np.index <- which(!is.snp)
	emission[np.index, , ] <- 1
	logemit <- log(emission)
	return(logemit)
}

simulateSingleDupBaf <- function(b, is.snp, from, to, ...){
	stopifnot(is(b, "numeric"))
	names(b) <- NULL
	b.all <- b
	b <- b[is.snp]
	sds <- VanillaICE:::updateSigma(b, is.snp=rep(TRUE, length(b)), sigma0=c(0.02, 0.04, 0.02))
##	if("pb" %in% names(list(...))){
##		pb <- list(...)[["pb"]]
##		pb <- pb/100
##		pb[is.na(pb)] <- 0.5
##		p <- pb
##	} else p <- rep(0.5, length(b))
	p <- 0.5
	q3 <- (1-p)^3
	p2q <- p^2*(1-p)
	pq2 <- p*(1-p)^2
	p3 <- p^3
	##tmp <- cbind(q3, 3*pq2, 3*p2q, p3)
	##stopifnot(all(rowSums(tmp) == 1))
	index <- seq(from, to, by=1)

	z <- sample(1:4, size=length(index), replace=TRUE, prob=c(q3, 3*pq2, 3*p2q, p3))
	nZ <- table(z)
	d1 <- rtnorm(length(index), 0, sds[1], lower=0, upper=1)
	d2 <- rtnorm(length(index), 1/3, sds[1], lower=0, upper=1)
	d3 <- rtnorm(length(index), 2/3, sds[1], lower=0, upper=1)
	d4 <- rtnorm(length(index), 1, sds[3], lower=0, upper=1)
	simB <- rep(NA, length(index))
	simB[z==1] <- sample(d1, nZ[1])
	simB[z==2] <- sample(d2, nZ[2])
	simB[z==3] <- sample(d3, nZ[3])
	simB[z==4] <- sample(d4, nZ[4])
	b.all[index] <- simB
	##het <- rtnorm(length(index), 0.5, sds[2], lower=0, upper=1)
	return(b.all)
}

simulateSingleDupLrr <- function(r, is.snp, cnStates=c(-1.5, -0.5, 0, 0, 0.4, 0.8),
				 from, to, ...){
	stdev <- .getSds(r)
	if(is(stdev, "matrix")) stopifnot(ncol(stdev) == ncol(object))
	stopifnot(all(dim(r) == dim(stdev)))
	if(any(colSums(is.na(r)) == nrow(r))){
		stop("Some samples have all missing values. Exclude these samples before continuing.")
	}
	S <- length(cnStates)
	for(j in seq_len(ncol(r))){
		##snp.index <- which(is.snp)
		cn <- r[, j]
		snp.index <- which(is.snp & !is.na(cn))
		if(length(snp.index) > 0){
			cn <- cn[snp.index]
			if(is(stdev, "matrix")){
				s <- stdev[snp.index, j]
			} else s <- stdev[snp.index]
			prOutlier <- probabilityOutlier(cn, k=k)
			for(l in seq_along(cnStates)){
				e <- (1-prOutlier) * dnorm(x=cn, mean=mu.snp[l], sd=s) + prOutlier * dunif(cn, MIN.CN, MAX.CN)
			}
		}
	}
	simR <- rep(NA, length(index))
	simR[z==1] <- sample(d1, nZ[1])
	simR[z==2] <- sample(d2, nZ[2])
	simB[z==3] <- sample(d3, nZ[3])
	simB[z==4] <- sample(d4, nZ[4])
	##het <- rtnorm(length(index), 0.5, sds[2], lower=0, upper=1)
	return(r.all)
}


simulateDoubleDupBaf <- function(b, is.snp, from, to, ...){
	stopifnot(is(b, "numeric"))
	names(b) <- NULL
	b.all <- b
	b <- b[is.snp]
	sds <- VanillaICE:::updateSigma(b, is.snp=rep(TRUE, length(b)), sigma0=c(0.02, 0.04, 0.02))
##	if("pb" %in% names(list(...))){
##		pb <- list(...)[["pb"]]
##		pb <- pb/100
##		pb[is.na(pb)] <- 0.5
##		p <- pb
##	} else p <- rep(0.5, length(b))
	p <- 0.5
	q4 <- (1-p)^4
	pq3 <- p*(1-p)^3
	p2q2 <- p^2*(1-p)^2
	p3q <- p^3*(1-p)
	p4 <- p^4
	##tmp <- cbind(q3, 3*pq2, 3*p2q, p3)
	##stopifnot(all(rowSums(tmp) == 1))
	index <- seq(from, to, by=1)
	z <- sample(1:5, size=length(index), replace=TRUE, prob=c(q4, pq3, p2q2, p3q, p4))
	nZ <- table(z)
	d1 <- rtnorm(length(index), 0, sds[1], lower=0, upper=1)
	d2 <- rtnorm(length(index), 1/4, sds[1], lower=0, upper=1)
	d3 <- rtnorm(length(index), 1/2, sds[2], lower=0, upper=1)
	d4 <- rtnorm(length(index), 3/4, sds[1], lower=0, upper=1)
	d5 <- rtnorm(length(index), 1, sds[3], lower=0, upper=1)
	simB <- rep(NA, length(index))
	simB[z==1] <- sample(d1, nZ[1])
	simB[z==2] <- sample(d2, nZ[2])
	simB[z==3] <- sample(d3, nZ[3])
	simB[z==4] <- sample(d4, nZ[4])
	simB[z==5] <- sample(d5, nZ[5])
	b.all[index] <- simB
	##het <- rtnorm(length(index), 0.5, sds[2], lower=0, upper=1)
	return(b.all)
}

simulateSingleDelBaf <- function(b, is.snp, from, to, ...){
	index <- seq(from, to, by=1)
	stopifnot(all(diff(index) > 0))
	stopifnot(length(index) > 1)
	stopifnot(is(b, "numeric"))
	names(b) <- NULL
	b.all <- b
	b <- b[is.snp]
	sds <- VanillaICE:::updateSigma(b, is.snp=rep(TRUE, length(b)), sigma0=c(0.02, 0.04, 0.02))
##	if("pb" %in% names(list(...))){
##		pb <- list(...)[["pb"]]
##		pb <- pb/100
##		pb[is.na(pb)] <- 0.5
##		p <- pb
##	} else p <- rep(0.5, length(b))
	p <- 0.5
	q <- 1-p
	##tmp <- cbind(q3, 3*pq2, 3*p2q, p3)
	##stopifnot(all(rowSums(tmp) == 1))

	z <- sample(1:2, size=length(index), replace=TRUE, prob=c(p, q))
	nZ <- table(z)
	d1 <- rtnorm(length(index), 0, sds[1], lower=0, upper=1)
	d2 <- rtnorm(length(index), 1, sds[3], lower=0, upper=1)
	simB <- rep(NA, length(index))
	simB[z==1] <- sample(d1, nZ[1])
	simB[z==2] <- sample(d2, nZ[2])
	b.all[index] <- simB
	##het <- rtnorm(length(index), 0.5, sds[2], lower=0, upper=1)
	return(b.all)
}
