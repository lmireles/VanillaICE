computeBpiEmission.SnpSuperSet <- function(object, hmmOptions, isBPI){
	states <- hmmOptions[["states"]]
	prGtError <- hmmOptions[["prGtError"]]
	ICE <- hmmOptions[["ICE"]]
	emission <- matrix(NA, nrow(object), ncol=2)
	colnames(emission) <- states
	if(ICE){
		pCrlmm <- confs(object)  ## crlmm confidence score
		## take the minimum confidence score in the trio
		pCrlmm <- apply(pCrlmm, 1, min, na.rm=TRUE)
		## set emission probability to min(crlmmConfidence, 0.999)
		I <- as.integer(pCrlmm < (1 - prGtError[["BPI"]]))
		pCrlmm <- pCrlmm*I + (1 - prGtError[["BPI"]])*(1-I)
		emission[,  "BPI"] <- pCrlmm
		##Pr(mendelian inconsistency | BPI) = 0.001
		emission[, "BPI"] <- 1 - pCrlmm
	} else { ##ignore confidence scores
		##Pr(call is consistent with biparental inheritance | BPI) = 0.999
		emission[isBPI==TRUE,  "BPI"] <-  1-prGtError["BPI"]
		##Pr(mendelian inconsistency | BPI) = 0.001
		emission[isBPI==FALSE, "BPI"] <- prGtError["BPI"] ##Mendelian inconsistancy
	}
	##Pr(call is consistent with biparental inheritance | not BPI) = 0.01
	emission[isBPI==TRUE,  "notBPI"] <- prGtError["notBPI"]   ## biparental inheritance, but true state is not Biparental
	##Pr(mendelian inconsistency | not BPI) = 0.99
	emission[isBPI==FALSE, "notBPI"] <- 1-prGtError["notBPI"] ## Mendelian inconsistancy
	log.emission <- log(emission)
	return(log.emission)
}

isBiparental.SnpSuperSet <- function(object, allowHetParent=FALSE){
	##if(length(object$familyMember) < 3) stop("object$familyMember not the right length")
	father <- 1
	mother <- 2
	offspring <- 3
	F <- calls(object[, father])
	M <- calls(object[, mother])
	O <- calls(object[, offspring])
	object <- cbind(F, M, O)
	colnames(object) <- c("father", "mother", "offspring")
	biparental <- isBiparental.matrix(object, allowHetParent=allowHetParent)
	return(biparental)
}

isBiparental.matrix <- function(object, allowHetParent=TRUE){
	F <- object[, 1]
	M <- object[, 2]
	O <- object[, 3]
	##M/F AA, F/M BB, O AB
	##isHet <- offspringHeterozygous(object)  ##offspring is heterozygous
	biparental <- rep(NA, nrow(object))
	biparental[F==1 & M == 3 & O == 2] <- TRUE#Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
	biparental[F==3 & M == 1 & O == 2] <- TRUE#Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
	##M/F AA, F/M BB, O AA or BB
	biparental[F==1 & M == 3 & (O == 1 | O == 3)] <- FALSE#Pr(O | biparental)=0.001, Pr(O | not biparental) = 1-0.01
	biparental[F==3 & M == 1 & (O == 1 | O == 3)] <- FALSE#Pr(O | biparental)=0.001, Pr(O | not biparental) = 1-0.01
	## M/F AA, F/M BB, O AB
	if(allowHetParent) biparental[F == 1 & M == 2 & O == 2] <- TRUE#Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
	if(allowHetParent) biparental[F == 2 & M == 1 & O == 2] <- TRUE#Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
	## F AB, M AA, O BB is not biparental
	## F AA, M AB, O BB is not biparental
	biparental[F == 2 & M == 1 & O == 3] <- FALSE#Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
	biparental[F == 1 & M == 2 & O == 3] <- FALSE#Pr(O | biparental)=0.001, Pr(O | not biparental) = 1-0.01
	## M AA, F AB, O AB
	if(allowHetParent) biparental[F == 2 & M == 3 & O == 2] <- TRUE#Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
	if(allowHetParent) biparental[F == 3 & M == 2 & O == 2] <- TRUE#Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
	## F=AB, M=BB, O=AA is NOT biparental
	biparental[F == 2 & M == 3 & O == 1] <- FALSE#Pr(O | biparental)=0.001, Pr(O | not biparental) = 1-0.01
	biparental[F == 3 & M == 2 & O == 1] <- FALSE#Pr(O | biparental)=0.001, Pr(O | not biparental) = 1-0.01
	return(biparental)
}


getEmission.nps <- function(object, hmmOptions){
	##****************************************************
	##	                                             *
	##  Emission probabilities for nonpolymorphic probes *
	##	                                             *
	##****************************************************
	batch <- unique(object$batch)
	scaleSds <- hmmOptions[["scaleSds"]]
	cnStates <- hmmOptions[["copynumberStates"]]
	verbose <- hmmOptions[["verbose"]]
	if(verbose) message("Computing emission probabilities for nonpolymorphic loci.")
	if(scaleSds){
		##a <- log2(CA(object))
		##sds.a <- apply(a, 2, mad, na.rm=TRUE)
		##sds.a <- sds.a/median(sds.a)
		sds.a <- robustSds(log2(CA(object)))
		sds.a <- matrix(sds.a, nrow(object), ncol(object), byrow=TRUE)
	} else sds.a <- matrix(0, nrow(object), ncol(object))
	emissionProbs <- array(NA, dim=c(nrow(object),
				   ncol(object), length(cnStates)))
	nuA <- getParam(object, "nuA", batch)
	phiA <- getParam(object, "phiA", batch)
	sig2A <- getParam(object, "sig2A", batch)
	if(any(is.na(sig2A))){
		sig2A[is.na(sig2A)] <- median(sig2A, na.rm=TRUE)
	}
	##tau2A <- getParam(object, "tau2A", batch)
	##Assume that on the log-scale, that the background variance is the same...
	##tau2A <- sig2A
	a <- as.numeric(log2(A(object)))
##	if(any(cnStates > 2)){
##		cnStates[cnStates > 2] <- cnStates[cnStates > 2] * 0.85
##	}
	for(k in seq(along=cnStates)){
		CT <- cnStates[k]
		mus.matrix=matrix(log2(nuA + CT*phiA), nrow(object), ncol(object))
		mus <- as.numeric(matrix(log2(nuA + CT*phiA), nrow(object), ncol(object)))
		sds.matrix <- matrix(sqrt(sig2A), nrow(object), ncol(object))
		##sds.matrix <- sds.matrix + sds.a
		sds.matrix <- sds.matrix*sds.a
		sds <- as.numeric(sds.matrix)
		tmp <- matrix(dnorm(a, mean=mus, sd=sds), nrow(object), ncol(object))
		emissionProbs[, , k] <- log(tmp)
	}
	emissionProbs
}


getEmission.snps <- function(object, hmmOptions){
	batch <- unique(object$batch)
	if(length(batch) > 1) stop("batch variable not unique")
	scaleSds <- hmmOptions[["scaleSds"]]
	cnStates <- hmmOptions[["copynumberStates"]]
	verbose <- hmmOptions[["verbose"]]
	if(verbose) message("Computing emission probabilities for polymorphic loci.")
	if(scaleSds){
		##a <- log2(CA(object) + CB(object))
		##sds.a <- apply(a, 2, mad, na.rm=TRUE)
		##sds.a <- sds.a/median(sds.a)
		##sds.a <- robustSds(log2(CA(object) + CB(object)))
		sds.a <- robustSds2(log2(CA(object)))
		##sds.a[sds.a < 1] <- 1
		sds.a <- matrix(sds.a, nrow(object), ncol(object), byrow=TRUE)
	} else sds.a <- sds.b <- matrix(0, nrow(object), ncol(object))
	emissionProbs <- array(NA, dim=c(nrow(object),
				   ncol(object), length(cnStates)))
	corr <- getParam(object, "corr", batch)
	corrA.BB <- getParam(object, "corrA.BB", batch)
	corrB.AA <- getParam(object, "corrB.AA", batch)
	nuA <- getParam(object, "nuA", batch)
	nuB <- getParam(object, "nuB", batch)
	phiA <- getParam(object, "phiA", batch)
	phiB <- getParam(object, "phiB", batch)
	sig2A <- getParam(object, "sig2A", batch)
	sig2B <- getParam(object, "sig2B", batch)
	tau2A <- getParam(object, "tau2A", batch)
	tau2B <- getParam(object, "tau2B", batch)
	a <- as.numeric(log2(A(object)))
	b <- as.numeric(log2(B(object)))
	for(k in seq(along=cnStates)){
		T <- cnStates[k]
		f.x.y <- matrix(0, sum(nrow(object)), ncol(object))
		for(copyA in 0:T){
			copyB <- T-copyA
			sigmaA <- sqrt(tau2A*(copyA==0) + sig2A*(copyA > 0))
			sigmaB <- sqrt(tau2B*(copyB==0) + sig2B*(copyB > 0))
			if(copyA == 0 & copyB > 0) r <- corrA.BB
			if(copyA > 0 & copyB == 0) r <- corrB.AA
			if(copyA > 0 & copyB > 0) r <- corr
			if(copyA == 0 & copyB == 0) r <- 0
			muA <- log2(nuA+copyA*phiA)
			muB <- log2(nuB+copyB*phiB)

			sigmaA <- matrix(sigmaA, nrow=length(sigmaA), ncol=ncol(object), byrow=FALSE)
			sigmaB <- matrix(sigmaB, nrow=length(sigmaB), ncol=ncol(object), byrow=FALSE)
			## scale the variances by a sample-specific estimate of the variances
			## var(I_A, ijp) = sigma_A_ip * sigma_A_jp
			##sigmaA <- sigmaA+sds.a
			##sigmaB <- sigmaB+sds.a
			sigmaA <- sigmaA*sds.a
			sigmaB <- sigmaB*sds.a
			meanA <- as.numeric(matrix(muA, nrow(object), ncol(object)))
			meanB <- as.numeric(matrix(muB, nrow(object), ncol(object)))
			rho <- as.numeric(matrix(r, nrow(object), ncol(object)))
			sd.A <- as.numeric(matrix(sigmaA, nrow(object), ncol(object)))
			sd.B <- as.numeric(matrix(sigmaB, nrow(object), ncol(object)))
			Q.x.y <- 1/(1-rho^2)*(((a - meanA)/sd.A)^2 + ((b - meanB)/sd.B)^2 - 2*rho*((a - meanA)*(b - meanB))/(sd.A*sd.B))
			## For CN states > 1, assume that any of the possible genotypes are equally likely a priori...just take the sum
			## For instance, for state copy number 2 there are three combinations: AA, AB, BB
			##   -- two of the three combinations should be near zero.
			## TODO: copy-neutral LOH would put near-zero mass on both copyA > 0, copyB > 0
			f.x.y <- f.x.y + matrix(1/(2*pi*sd.A*sd.B*sqrt(1-rho^2))*exp(-0.5*Q.x.y), nrow(object), ncol(object))
		}
		emissionProbs[, , k] <- log(f.x.y)
	}
	emissionProbs
}


