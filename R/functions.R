Sweave2pdf <- function(fname, ...){
	fname <- strsplit(fname, ".Rnw")[[1]][[1]]
	suppressWarnings(Sweave(paste(fname, ".Rnw", sep=""), ...))
	texi2dvi(paste(fname, ".tex", sep=""), pdf=TRUE)
}

rowMAD <- function(x, y, ...){
	notna <- !is.na(x)
	sds <- 1.4826*rowMedians(abs(x-rowMedians(x, ...)), ...)
	return(sds)
}

robustSds <- function(x, takeLog=FALSE, ...){
	if(!is.matrix(x)) stop("x is not a matrix")
	if(takeLog) x <- log2(x)
##	if(ncol(x) > 3){
##		sds1 <- rowMAD(x, ...)
##		sds1 <- matrix(sds1, nrow(x), ncol(x))
##		sds2 <- apply(x, 2, "mad", ...)
##		df <- ncol(x)
##		##sds2 <- sds2/median(sds2, na.rm=TRUE)
##		##sds2 <- sds2/min(sds2, na.rm=TRUE)
##		##sds <- t(t(sds1) * sds2)
##	} else {
	sds <- apply(x, 2, "mad", ...)
	sds <- matrix(sds, nrow(x), ncol(x), byrow=TRUE)
	dimnames(sds) <- dimnames(x)
	return(sds)
}

## many markers have a level of noise that is far greater than the
## noise level of the sample.
##
## -- while the uncertainty estimates for the marker may be poorly
##    estimated, shrinking to the noise level of the sample may have
##    the effect of giving uncertainty estimates for the marker that
##    are much too small
robustSds2 <- function(x, DF.PRIOR=10, nSamples, takeLog=FALSE, ...){
	if(!is.matrix(x)) stop("x is not a matrix")
	if(takeLog) x <- log2(x)

	## the posterior is IG(nu_n/2, nu_n * sigma2_n (theta)/2)
	## nu_n = nu_0 + n
	## sigma2_n(theta) = 1/nu_n * [nu_0 * sigma^2_0 + n*s^2_n(theta)],
	## where s^2_n(theta) = Sum(yi-theta)^2/n
	##
	##
	## 1/sigma^2 ~ G(nu_0/2, nu_0/2 * sigma^2_0)
	## E(sigma^2) = sigma^2_0 * nu_0/2 /(nu_0/2 -1)
	##
	## 1/sigma^2 | ... ~ G(nu_n/2, nu_n * sigma^2_n/2)
	##
	## nu_n = nu_0 + n
	## sigma^2_n = 1/nu_n *[nu_0*sigma^2_0 + (n-1)*s^2 + k_0*n/k_n*(ybar-mu_0)^2]
	##
	## prior
	##
	## data
	##
	## posterior inference
	sigma2.0 <- apply(x, 2, MAD, na.rm=TRUE)
	nu0 <- 100
	nn <- ncol(x)
	nu.n <- nu0+nn
	s2 <- rowMAD(x, na.rm=TRUE)
	k0 <- 1 ## ?
	kn <- k0+nn
	ybar <- rowMedian(x, na.rm=TRUE)
	mu0 <- median(ybar)
	##
	## E[1/sigma2_g | ... ] = sigma^2n/2 * nu_n / (nu_n/2 - 1)




	##
##	sigma.marker <- rowMAD(x, na.rm=TRUE)
##	sigma.sample <- apply(x, 2, MAD, na.rm=TRUE)
##	gammahat <- sigma.marker/median(sigma.marker, na.rm=TRUE)
##	df1 <- nSamples-1
##	df2 <- length(sds1)-1
	sds1 <- rowMAD(x, na.rm=TRUE)
##	sds.marker <- (df1*sds.marker + df2*median(sds.marker,na.rm=TRUE))/(df1+df2)
	sds1 <- matrix(sds1, nrow(x), ncol(x))
	sds2 <- apply(x, 2, "mad", constant=2, na.rm=TRUE)
	df <- ncol(x)
	sds2 <- matrix(sds2, nrow(x), ncol(x), byrow=TRUE)
	sds.star <- (sds1 * df + sds2*DF.PRIOR)/(df+DF.PRIOR)
	dimnames(sds.star) <- dimnames(x)
	return(sds.star)
}

viterbi.wrapper <- function(log.emission,
			    log.initial,
			    transitionPr,
			    arm,
			    S,
			    T,
			    result,
			    delta,
			    normal2altered,
			    altered2normal,
			    altered2altered,
			    normalIndex,
			    pAA){
	tmp <- list(as.matrix(as.double(as.matrix(log.emission))),##beta
		    as.double(as.matrix(log.initial)),##initialP
		    as.matrix(as.double(transitionPr)),
		    as.integer(arm),##arm
		    as.integer(S),##number of states
		    as.integer(T),##number of loci
		    result,##placeholder for results
		    as.matrix(as.double(delta)),##delta
		    normal2altered=normal2altered,##c1
		    altered2normal=altered2normal,##c2
		    altered2altered=altered2altered,##c3
		    as.integer(normalIndex),
		    as.double(pAA))##normalState
	.C("viterbi",
	   log.emission=tmp[[1]],
	   log.initial=tmp[[2]],
	   tau=tmp[[3]],
	   arm=tmp[[4]],
	   S=tmp[[5]],
	   T=tmp[[6]],
	   viterbiSeq=tmp[[7]],
	   delta=tmp[[8]],
	   N2A=tmp[[9]],
	   A2N=tmp[[10]],
	   A2A=tmp[[11]],
	   normalIndex=tmp[[12]],
	   pAA=tmp[[13]])  ##can verify that the transition prob. matrix is correct (for last snp)
}

getChromosomeArm <- function(object){
	chrom <- chromosome(object)
	pos <- position(object)
	if(!is.integer(chrom)) {
		chrom <- chromosome2integer(chrom)
	}
	if(!all(chrom %in% 1:24)){
			warning("Chromosome annotation is currently available for chromosomes 1-22, X and Y. Other chromosomes or NA's present.")
			marker.index <- which(chromosome(object) <= 24 & !is.na(chromosome(object)))
			##message("Please add/modify data(chromosomeAnnotation, package='SNPchip') to accomodate special chromosomes")
			##stop()
	} else{
		marker.index <- seq_along(chrom)
	}
	if(!is.integer(pos)) {
		message("Coerced pos to an integer.")
		pos <- as.integer(pos)
	}
	data(chromosomeAnnotation, package="SNPchip", envir=environment())
	chromosomeAnnotation <- as.matrix(chromosomeAnnotation)
	chrAnn <- chromosomeAnnotation
	uchrom <- unique(SNPchip:::integer2chromosome(chrom))
	uchrom <- uchrom[!is.na(uchrom)]
	chromosomeArm <- vector("list", length(uchrom))
	positionList <- split(pos[marker.index], chrom[marker.index])
	chr.names <- unique(chrom)
	chr.names <- as.character(chr.names[!is.na(chr.names)])
	positionList <- positionList[match(chr.names, names(positionList))]
	for(i in seq(along=chr.names)){
		chromosomeArm[[i]] <- as.integer(ifelse(positionList[[i]] <= chrAnn[uchrom[i], "centromereEnd"], 0, 1))
	}
	chromosomeArm <- unlist(chromosomeArm)
	chrom <- chrom[marker.index]
	chromosomeArm <- cumsum(c(0, diff(chromosomeArm) != 0 | diff(chrom) != 0))
	chromosomeArm <- chromosomeArm+1 ##start at 1
	res <- rep(NA, nrow(object))
	res[marker.index] <- chromosomeArm
	return(res)
}

viterbi <- function(object,
		    hmm.params,
		    log.E, ...){
		    ##verbose=TRUE,
		    ##normal2altered=1,
		    ##altered2normal=1,
		    ##altered2altered=1,
		    ##TAUP=1e8, ...){
	##log.E <- hmm.params[["log.emission"]]
	sns <- colnames(log.E)
	if(is.null(sns)) stop("no dimnames for log.emission")
	##log.initial <- hmm.params[["log.initial"]]
	log.initial <- hmm.params$log.initialPr
	verbose <- hmm.params$verbose
	normal2altered <- hmm.params$n2a
	altered2normal <- hmm.params$a2n
	altered2altered <- hmm.params$a2a
	##
	if(normal2altered <= 0) stop("normal2altered must be > 0")
	if(altered2normal <= 0) stop("altered2normal must be > 0")
	if(altered2altered <= 0) stop("altered2altered must be > 0")
	##
	## arm can contain NA's for invalid chromosomes or NA's
	arm <- getChromosomeArm(object)
	normalIndex <- hmm.params$normalIndex##[["normalIndex"]]
	if(normalIndex < 1 | normalIndex > dim(log.E)[[3]]){
		stop("normalIndex in hmm.params not valid")
	}
	TAUP <- hmm.params$tau
	c1 <- normal2altered
	c2 <- altered2normal
	c3 <- altered2altered
	##TT <- nrow(object)
	index <- which(!is.na(arm))
	arm <- arm[index]
	TT <- length(index)
	log.E <- log.E[index, , , drop=FALSE]
	object <- object[index, ]
	##
	states <- hmm.params$states
	names(log.initial) <- states
	S <- length(states)
	delta <- matrix(as.double(0), nrow=TT, ncol=S)
	rangedData <- list()
	for(j in 1:ncol(log.E)){
		rD <- vector("list", length(unique(arm)))
		for(a in seq_along(unique(arm))){
			missingE <- rowSums(is.na(log.E[, j, ])) > 0
			notFinite <- rowSums(!is.finite(log.E[, j, ])) > 0
			missingE <- missingE | notFinite
			I <- arm == a  & !missingE
			if(sum(I) < 2) next()
			T <- sum(I)
			transitionPr <- exp(-2 * diff(position(object)[I])/TAUP)
			##is the lower bound a function of normal2altered, altered2normal, altered2altered?
			minimum <- 1-1/((S-1)*c1) + 0.01
			transitionPr[transitionPr < minimum] <- minimum
			if(any(transitionPr < 0 | transitionPr > 1)) stop("Transition probabilities not in [0,1].  Order object by chromosome and physical position")
			result <- rep(as.integer(0), T)
			viterbiResults <- viterbi.wrapper(log.emission=log.E[I, j, ],
							  log.initial=log.initial,
							  transitionPr=transitionPr,
							  arm=arm[I],
							  S=S,
							  T=T,
							  result=result,
							  delta=delta[I, ],
							  normal2altered=normal2altered,
							  altered2normal=altered2normal,
							  altered2altered=altered2altered,
							  normalIndex=normalIndex,
							  pAA=rep(0, S^2))
			M <- matrix(viterbiResults[["pAA"]], S, S)
			if(!all(is.finite(M))) stop("Infinite values in transition prob. matrix. Check that rows are ordered by physical position")
			if(!all.equal(rowSums(exp(M)), rep(1, S))){
				warning("Rows of the transition probability matrix do not sum to 1")
			}
			viterbiSequence <- viterbiResults[["viterbiSeq"]]
			rl <- Rle(viterbiSequence)
			starts <- start(rl)
			LLR <- rep(999,  length(starts))
			log.emission <- matrix(viterbiResults[["log.emission"]], T, S)
			##** The NA is to stagger the transition probabilities by 1
			##  -- this way, the same index can be used to multiply the transition and emission probabilities
			p <- c(NA, as.numeric(viterbiResults[["tau"]]))
			lP.N2N <- log(1-((1-p)*(S-1)*c1)) ##probability normal -> normal
			lP.N2A <- log((1-p)*c1) ##probability normal -> altered
			P.A2A <- sapply(1-((1-p)*(c2+(S-2)*c3)), function(x) max(x, 0.01))
			lP.A2A <- log(P.A2A) ## probability altered to same altered state
			lP.A2N <- log((1-p)*c2) ##probability altered -> normal
			lP.A2Astar <- log((1-p)*c3) ## probability altered -> different altered state
			##For each seqment, compute the likelihood ratio

			for(k in seq(along=starts)){
				##if(any(LLR < 0)) browser()
				index <- start(rl)[k]:end(rl)[k]
				thisState <- unique(viterbiSequence[index])
				if(thisState == normalIndex){
					LLR[k] <- 0
					next()
				}
				first.index <- min(index)
				last.index <- max(index)
				## 6 Rules  (1 < t < T)
				## 1.  index=1
				## 2.  index=t
				## 3.  index=t,t+1
				## 4.  index=T
				## 5.  index=1,2
				## 6,  index=T-1, T
				##------
				## 1. index=1
				if(first.index == last.index & last.index==1){
					##note the last term cancels
					logLik.vit <- log.initial[thisState]+log.emission[1, thisState] + lP.A2N[2] + log.emission[2, normalIndex]
					logLik.null <- log.initial[normalIndex]+log.emission[1, normalIndex] + lP.N2N[2] + log.emission[2, normalIndex]
					LLR[k] <- logLik.vit-logLik.null
					next()
				}
				##2 index = t
				if(length(index) == 1 & first.index > 1 & last.index < T){
					##note the last term cancels
					logLik.vit <- sum(lP.N2A[index] + log.emission[index, thisState]) + lP.A2N[last.index+1]+log.emission[last.index+1, normalIndex]
					logLik.null <- sum(lP.N2N[index] + log.emission[index, normalIndex]) + lP.N2N[last.index+1] + log.emission[last.index+1, normalIndex]
					LLR[k] <- logLik.vit-logLik.null
					next()
				}
				##if T=t+1?
				## 3: index = t, ..., t*  , t>1, t* < T, t* != t
				if(first.index != last.index & first.index > 1 & last.index < T){
					index2 <- index[-1]
					logLik.vit <- lP.N2A[first.index] +
					              sum(lP.A2A[index2]) +
					              lP.A2N[last.index+1] +
						      log.emission[first.index, thisState] +
						      sum(log.emission[index2, thisState]) +
					              log.emission[last.index+1, normalIndex]
					logLik.null <-
						sum(lP.N2N[index]) +
						lP.N2N[last.index+1]  +
						sum(log.emission[index, normalIndex]) +
						log.emission[last.index+1, normalIndex]
					LLR[k] <- logLik.vit-logLik.null
					next()
				}
				## 4: index = T
				if(first.index == last.index & last.index == T){
					logLik.vit <- lP.N2A[T] + log.emission[T, thisState]
					logLik.null <- lP.N2N[T] + log.emission[T, normalIndex]
					LLR[k] <- logLik.vit-logLik.null
					next()
				}
				## 5: index = 1, 2, ...
				if(first.index != last.index & first.index == 1 & last.index < T){
					index2 <- index[-1]## t=2, ...., t*
					logLik.vit <- log.initial[thisState] + log.emission[first.index, thisState]  + sum(lP.A2A[index2] + log.emission[index2, thisState]) + lP.A2N[last.index+1] + log.emission[last.index+1, normalIndex]
					logLik.null <- log.initial[normalIndex] + log.emission[first.index, normalIndex] + sum(lP.N2N[index2] + log.emission[index2, normalIndex]) + lP.N2N[last.index+1] + log.emission[last.index+1, normalIndex]
					LLR[k] <- logLik.vit-logLik.null
					next()
				}
				if(first.index != last.index & first.index == 1 & last.index == T){
					index2 <- index[-1]## t=2, ...., t*
					logLik.vit <- log.initial[thisState] + log.emission[first.index, thisState]  + sum(lP.A2A[index2] + log.emission[index2, thisState]) ##+ lP.A2N[last.index+1] + log.emission[last.index+1, normalIndex]
					logLik.null <- log.initial[normalIndex] + log.emission[first.index, normalIndex] + sum(lP.N2N[index2] + log.emission[index2, normalIndex]) ##+ lP.N2N[last.index+1] + log.emission[last.index+1, normalIndex]
					LLR[k] <- logLik.vit-logLik.null
					next()
				}
				## 6: index = t, ...T
				if(first.index != last.index & last.index == T){
					index2 <- index[-1]
					logLik.vit <- lP.N2A[first.index] + log.emission[first.index, thisState] + sum(lP.A2A[index2] + log.emission[index2, thisState])
					logLik.null <- lP.N2N[first.index] + log.emission[first.index, normalIndex] + sum(lP.N2N[index2] + log.emission[index2, normalIndex])
					LLR[k] <- logLik.vit - logLik.null
				}
			}
			start.index <- start(rl)
			end.index <- end(rl)
			pos <- position(object)[I]
			##this is tricky since we've added an index to force a segment for each arm.
			start <- pos[start.index]
			end <- pos[end.index]
			##numMarkers <- unlist(numMarkers)
			numMarkers <- width(rl)
			states <- viterbiSequence[start.index]
			ir <- IRanges(start=start, end=end)
			rD[[a]] <- RangedData(ir,
					      ##space=rep(paste("chr", unique(chromosome(object)[I]), sep=""), length(ir)),
					      chrom=rep(unique(chromosome(object)[I]), length(ir)),
					      ##sampleId=sampleNames(object)[j],
					      sampleId=colnames(log.E)[j],
					      state=states,
					      numMarkers=numMarkers,
					      LLR=LLR)
		}
		notnull <- !sapply(rD, is.null)
		rD <- rD[notnull]
		L <- sapply(rD, nrow)
		if(any(L == 1) & any(L > 1)){
			rangedData[[j]] <- c(do.call(c, rD[L == 1]), do.call(c, rD[L > 1]))
		} else {
			rangedData[[j]] <- do.call(c, rD)
		}
	}
##	L <- sapply(rangedData, nrow)
##	if(any(L == 1) & any(L > 1)){
##		rangedData <- c(do.call(c, rangedData[L == 1]), do.call(c, rangedData[L > 1]))
##	} else {
##		rangedData <- do.call(c, rangedData)
##	}
	sampleId <- unlist(lapply(rangedData, function(x) x$sampleId))
	state <- unlist(lapply(rangedData, function(x) x$state))
	numMarkers <- unlist(lapply(rangedData, function(x) x$numMarkers))
	LLR <- unlist(lapply(rangedData, function(x) x$LLR))
	chr <- unlist(lapply(rangedData, function(x) x$chrom))
	starts <- unlist(lapply(rangedData, start))
	ends <- unlist(lapply(rangedData, end))
	ir <- IRanges(start=starts, end=ends)
	rangedData <- RangedDataHMM(ranges=ir,
				    chromosome=chr,
				    sampleId=sampleId,
				    state=state,
				    coverage=numMarkers,
				    LLR=LLR)
	##rangedData <- as(rangedData, "RangedDataCn")
	return(rangedData)
}






##viterbi <- function(emission,
##		    tau,
##		    arm,
##		    initialStateProbs,
##		    verbose=FALSE,
##		    chromosome,
##		    position,
##		    sampleNames,
##		    locusNames,
##		    normalIndex,
##		    returnLikelihood=FALSE,
##		    normal2altered=1,
##		    altered2normal=1,
##		    altered2altered=1){
##	if(class(emission) != "array") stop("emission probabilities must be an array: snps, samples, states. ")
##	if(missing(sampleNames)) sampleNames <- colnames(emission)
##	if(missing(normalIndex)) stop("Must specify integer for normalIndex")
##	if(!is.numeric(normalIndex)) stop("normalIndex should be numeric")
##	viterbiSequence <- matrix(NA, nrow(emission), ncol(emission), dimnames)
##	S <- dim(emission)[3]
##	T <- nrow(emission)
##	if(missing(initialStateProbs)){
##		initialStateProbs <- log(rep(1/S, S))
##	}
##	if(length(initialStateProbs) != S){
##		stop("initialStateProbs (the initial state probabilities, should be a numeric vector of length S, where S is the number of hidden states")
##	}
##	if(!all(initialStateProbs <= 0)){
##		if(all(initialStateProbs >= 0 & initialStateProbs <= 1)){
##			initialStateProbs <- log(initialStateProbs)
##		} else stop("initial state probabilities should be a probability or a log probability")
##	}
##	if(any(is.na(emission))){
##		if(verbose) message("Converting missing values in the emission matrix to 0")
##		emission[is.na(emission)] <- 0
##	}
##	if(any(is.nan(emission))){
##		message("some of the log emission probabilities are NaN.  Replacing with 0")
##		emission[is.nan(emission)] <- 0
##	}
##	if(any(emission < -50)){
##		message("some of the log emission probabilities are very small -- probable outliers.  Replacing with a small value (-10)")
##		emission[emission < -50] <- -50
##	}
##	if(missing(arm)){
##		message("chromosome arm not specified...HMM is not fit separately to each chromosomal arm")
##		arm <- rep(as.integer(1), T)
##	}
##	if(length(arm) != T) {
##		message("arm not the right length.  assuming all values on same chromosomal arm")
##		arm <- rep(as.integer(1), T)
##	}
##	if(missing(tau)){
##		stop("transition probabilities not specified")
##	}
##	if(length(tau) != T) stop("tau must have length T")
##	## The last index is arbitrary and, by default, is NA. Must replace this by a number-- C can not handle NA's
##	tau[is.na(tau)] <- 0
##	delta <- matrix(as.double(0), nrow=T, ncol=S)
##	browser()
##	rangedData <- list()
##	for(j in 1:ncol(results)){
##		rD <- vector("list", length(unique(arm)))
##		for(a in seq(along=unique(arm))){
##			I <- arm == a
##			T <- sum(I)
##			result <- rep(as.integer(0), T)
##			tmp <- list(as.matrix(as.double(as.matrix(emission[I, j, ]))),##beta
##				    as.double(as.matrix(initialStateProbs)),##initialP
##				    as.matrix(as.double(tau[I])),##tau
##				    as.integer(arm[I]),##arm
##				    as.integer(S),##number of states
##				    as.integer(T),##number of loci
##				    result,##placeholder for results
##				    as.matrix(as.double(delta[I, ])),##delta
##				    normal2altered=normal2altered,##c1
##				    altered2normal=altered2normal,##c2
##				    altered2altered=altered2altered,##c3
##				    as.integer(normalIndex),
##				    as.double(rep(0, S^2)))##normalIndex
##			tmp2 <- .C("viterbi",
##				   emission=tmp[[1]],
##				   initialStateProbs=tmp[[2]],
##				   tau=tmp[[3]],
##				   arm=tmp[[4]],
##				   S=tmp[[5]],
##				   T=tmp[[6]],
##				   viterbiSeq=tmp[[7]],
##				   delta=tmp[[8]],
##				   N2A=tmp[[9]],
##				   A2N=tmp[[10]],
##				   A2A=tmp[[11]],
##				   normalIndex=tmp[[12]],
##				   pAA=tmp[[13]])  ##can verify that the transition prob. matrix is correct (for last snp)
##			##check transition probs.
##			M <- matrix(tmp2[["pAA"]], S, S)
##			if(!all(is.finite(M))) stop("Infinite values in transition prob. matrix")
##			if(!all.equal(rowSums(exp(M)), rep(1, S))){
##				warning("Rows of the transition probability matrix do not sum to 1")
##			}
##			viterbiSequence <- tmp2[["viterbiSeq"]]
##			logInitialP <- initialStateProbs
##			rl <- Rle(viterbiSequence)
##			starts <- start(rl)
##			LLR <- rep(NA,  length(starts))
##			logE <- matrix(tmp2[["emission"]], T, S)
##			p <- as.numeric(tmp2[["tau"]])
##			c1 <- normal2altered
##			c2 <- altered2normal
##			c3 <- altered2altered
##			lP.N2N <- log(1-((1-p)*(S-1)*c1)) ##probability normal -> normal
##			lP.N2A <- log((1-p)*c1) ##probability normal -> altered
##			lP.A2A <- log(1-((1-p)*(c2+(S-2)*c3))) ## probability altered to same altered state
##			lP.A2N <- log((1-p)*c2) ##probability altered -> normal
##			lP.A2Astar <- log((1-p)*c3) ## probability altered -> different altered state
##			##For each seqment, compute the likelihood ratio
##			for(k in seq(along=starts)){
##				index <- start(rl)[k]:end(rl)[k]
##				thisState <- unique(viterbiSequence[index])
##				first.index <- min(index)
##				last.index <- max(index)
##				if(last.index < T){
##					next.index <- last.index+1
##					next.state <- viterbiSequence[next.index]
##					lE.next.state <- logE[next.index, next.state]
##					tp.last <- lP.A2N[last.index] ##not next.index!
##				} else{ ##max(index) = T
##					lE.next.state <- NULL
##					tp.last <- NULL
##				}
##				if(min(index) > 1){
##					tps <- c(lP.N2A[min(index)-1], lP.A2A[index[-length(index)]], tp.last)
##					tps.normal <- lP.N2N[c(index-1, max(index))]
##				} else {
##					tps <- c(logInitialP[thisState], lP.A2A[index[-length(index)]], tp.last)
##					tps.normal <- c(logInitialP[normalIndex],
##							lP.N2N[index[-length(index)]])
##				}
##				lEs <- c(logE[index, thisState], lE.next.state)
##				lE.normal <- c(logE[index, normalIndex], lE.next.state)
##				loglik.Vit <- sum(lEs+tps)
##				loglik.null <- sum(lE.normal+tps.normal)
##				LLR[k] <- loglik.Vit-loglik.null
##			}
##			start.index <- start(rl)
##			end.index <- end(rl)
##			pos <- position[I]
##			##this is tricky since we've added an index to force a segment for each arm.
##			start <- pos[start.index]
##			end <- pos[end.index]
##			##numMarkers <- unlist(numMarkers)
##			numMarkers <- width(rl)
##			states <- viterbiSequence[start.index]
##			ir <- IRanges(start=start, end=end)
##			rD[[a]] <- RangedData(ir,
##					      space=rep(paste("chr", unique(chromosome[I]), sep=""), length(ir)),
##					      sampleId=sampleNames[j],
##					      state=states,
##					      numMarkers=numMarkers,
##					      LLR=LLR)
##		}
##		tmp <- do.call(c, rD[sapply(rD, nrow) == 1])
##		tmp2 <- do.call(c, rD[sapply(rD, nrow) > 1])
##		rD <- c(tmp, tmp2)
##		rangedData[[j]] <- rD
##		##results <- list(stateSequence=results, logLikelihoodRatio=lrdiff)
##	}
##	rangedData <- do.call(c, rangedData)
##	return(rangedData)
##}
trioOptions <- function(opts,
			states=c("BPI", "notBPI"),
			##initialP=c(0.99, 0.01),
			##TAUP=1e7,
			prGtError=c(0.001, 0.01),
			##verbose=FALSE,
			allowHetParent=FALSE,
			##normalIndex=1,
			##normal2altered=1,
			##altered2normal=1,
			useCrlmmConfidence=FALSE){
	names(prGtError) <- states
	opts[["states"]] <- states
	opts$prGtError <- prGtError
	opts$allowHetParent <- allowHetParent
	##useCrlmmConfidence=useCrlmmConfidence)
	return(opts)
}

##setMethod("coerce", c("CNSet", "list"), function(from, to){
##	if(ncol(object) < 3){
##		return("No complete trios")
##	}
##	stopifnot(all(c("familyId", "fatherId", "motherId", "individualId") %in% varLabels(object)))
##
##	for(i in 1:nrow(trios)){
##
##
##	}
##})
centerAutosomesAt <- function(x, at, ...){
	stopifnot(!missing(at))
	marker.index <- which(chromosome(x) <= 23)
	cn <- copyNumber(x)[marker.index, ]
	meds <- apply(cn, 2, median, na.rm=TRUE)
	cn.cen <- sweep(cn, 2, meds) + at
	copyNumber(x)[marker.index, ] <- cn.cen
	return(x)
}



hmm.setup <- function(object, ...){  ## whether the save the emission probabilities
	res <- HmmOptionList(object=object, ...)
	res <- as(res, "HmmOptionList")
	return(res)
}

##Code this in C
viterbiR <- function(emission, initialP, tau, arm){
##	emission: log emission probabilities
##	initialP: log initial state probabilities (vector of length S)
##	tau: transition probabilities (original scale)
##	tau.scale: matrix on original scale
	S <- ncol(emission)
	T <- nrow(emission)  ##matrix of T rows and S columns
	delta <- psi <- matrix(NA, T, S)
	delta[1, ] <- initialP + emission[1, ]
	psi[1, ] <- rep(0, S)
	i <- which(rowSums(is.na(emission)) > 0)
	tau.scale <- 1
	#emission[i, ] <- 0
	for(t in 2:T){
		if(t %% 10000 == 0) cat(".")
		if(arm[t] != arm[t-1]){
			delta[t, ] <- initialP + emission[t, ]
			psi[t, ] <- 0
			next()
		}
##		AA <- matrix(NA, nr=S, nc=S)
##		for(i in 1:S){
		AA <- matrix(tau[t-1], nr=S, nc=S)
		epsilon <- (1-tau[t-1])/(S-1)
		##eNormal <- (1-tauNormal[t-1])/(S-1)
		AA[upper.tri(AA)] <- AA[lower.tri(AA)] <- epsilon
		##AA[normalIndex, ] <- rep(eNormal, S)
		##AA[normalIndex, normalIndex] <- tauNormal[t-1]
		AA <- log(AA*tau.scale)
		for(j in 1:S){
			tmp <- delta[t-1, ] + AA[, j]
			delta[t, j] <- max(tmp) + emission[t, j]
			psi[t, j] <- order(tmp, decreasing=TRUE)[1]
		}
	}
	Pstar <- max(delta[nrow(delta), ])
	qhat <- rep(NA, nrow(delta))
	qhat[T] <- order(delta[T, ], decreasing=TRUE)[1]
	for(t in (T-1):1){
		if(arm[t] != arm[t+1]){
			qhat[t] <- order(delta[t, ], decreasing=TRUE)[1]
		} else {
			qhat[t] <- psi[t+1, qhat[t+1]]
		}
	}
	return(qhat)
}

##setMethod("hmm", "oligoSnpSet", function(object, hmmOptions){
##	viterbi(object, hmmOptions)
##})
##
##setMethod("hmm", "CNSet", function(object, hmmOptions){
##	viterbi(object, hmmOptions)
##})
##
##setMethod("hmm", "CopyNumberSet", function(object, hmmOptions){
##	viterbi(object, hmmOptions)
##})
##setMethod("hmm", signature(object="HmmOptionList", cnSet="oligoSnpSet"),
##hmm <- function(object, hmm.params, ...){
##	if(missing(hmm.params)) stop("missing hmm.params.  See hmm.setup")
##	if(missing(object)) stop("missing object.")
##	viterbi(object, hmm.params, ...)
##}


##hmm <- function(object,
##		states,
##		mu=NULL,
##		probs=NULL,
##		takeLog=FALSE,
##		initialP,
##		returnSegments=TRUE,
##		TAUP=1e8,
##		verbose=FALSE,
##		ice=FALSE,
##		envir,
##		normalIndex){
##	if(missing(envir)) envir <- new.env()
##	if(!all(c("position", "chromosome") %in% fvarLabels(object))){
##		stop("'position' and 'chromosome' must be in fvarLabels(object), or transitionPr must be provided")
##	}
##	object <- object[order(chromosome(object), position(object)), ]
##	envir[["locusset"]] <- object
##	if(missing(states)) stop("must specify states")
##	envir[["states"]] <- states
##	if(missing(initialP)) initialP <- rep(1, length(states))/length(states)
##	envir[["initialP"]] <- initialP
##	envir[["mu"]] <- mu
##	if(is.null(probs)){
##		if(ice){
##			probs <- c(0.05, 0.99, 0.7, 0.999)
##		}
##	}
##	envir[["probs"]] <- probs
##	envir[["takeLog"]] <- takeLog
##	envir[["returnSegments"]] <- returnSegments
##	envir[["TAUP"]] <- TAUP
##	envir[["verbose"]] <- verbose
##	envir[["ice"]] <- ice
##	tau <- transitionProbability(chromosome=chromosome(object),
##				     position=position(object),
##				     TAUP=TAUP,
##				     verbose=verbose)
##	arm <- tau[, "arm"]
##	transitionPr <- tau[, "transitionPr"]
##	envir[["transitionPr"]] <- transitionPr
##	envir[["arm"]] <- arm
##	if(takeLog){
##		copyNumber(object) <- log2(copyNumber(object))
##		mu <- log2(mu)
##	}
##	if(verbose) message("Calculating emission probabilities")
##	calculateEmission(object=object,
##			  mu=mu,
##			  probs=probs,
##			  envir=envir,
##			  states=states,
##			  verbose=verbose,
##			  ice=ice)
##	if(!is.null(envir[["emission.gt"]])){
##		emission <- envir[["emission.cn"]]+envir[["emission.gt"]]
##	} else {
##		emission <- envir[["emission.cn"]]
##	}
##	envir[["emission"]] <- emission
##	##emission <- .GlobalEnv[["emission.cn"]]+.GlobalEnv[["emission.gt"]]
##	viterbiResults <- viterbi(initialStateProbs=log(initialP),
##				  emission=emission,
##				  tau=transitionPr,
##				  arm=arm,
##				  verbose=verbose,
##				  normalIndex=normalIndex)
##	dimnames(viterbiResults) <- list(featureNames(object), sampleNames(object))
##	if(returnSegments){
##		viterbiResults <- breaks(x=viterbiResults,
##					 states=states,
##					 position=position(object),
##					 chromosome=chromosome(object),
##					 verbose=verbose)
##	}
##	return(viterbiResults)
##}

setMethod("update", "environment", function(object, ...){
	if(length(ls(object)) == 0) stop("nothing to update")
	hmm(object=object[["locusset"]],
	    states=object[["states"]],
	    mu=object[["mu"]],
	    probs=object[["probs"]],
	    takeLog=object[["takeLog"]],
	    initialP=object[["initialP"]],
	    returnSegments=object[["returnSegments"]],
	    TAUP=object[["TAUP"]],
	    verbose=object[["verbose"]],
	    ice=object[["ice"]],
	    envir=object)
})

findFatherMother <- function(offspringId, object){
	stopifnot(!missing(offspringId))
	family.id <- pData(object)[sampleNames(object) == offspringId, "familyId"]
	father.id <- pData(object)[sampleNames(object) == offspringId, "fatherId"]
	mother.id <- pData(object)[sampleNames(object) == offspringId, "motherId"]
	father.name <- sampleNames(object)[object$familyId == family.id & object$individualId == father.id]
	mother.name <- sampleNames(object)[object$familyId == family.id & object$individualId == mother.id]
	if(length(father.name) > 1 | length(mother.name) > 1){
		stop("More than 1 father and/or more than 1 mother.  Check annotation in phenoData")
	}
	if(length(father.name) < 1 ){
		father.name <- NA
	}
	if(length(mother.name) < 1){
		mother.name <- NA
	}
	fmo.trio <- c(father.name, mother.name, offspringId)
	names(fmo.trio) <- c("father", "mother", "offspring")
	return(fmo.trio)
}


##plot <- function(df, palette, xlim, show.coverage=TRUE, blackBorder=TRUE, sampleLabels.cex=0.5, labelAllSamples=TRUE,...){
##	stopifnot(length(unique(df$chr))==1)
##	mykey <- simpleKey(c("homo-del", "hemi-del", "amp")[palette %in% df$col], points=FALSE,
##		   rectangles=TRUE, col=palette[palette %in% df$col], space="top")
##	mykey$rectangles[["border"]] <- mykey$rectangles[["col"]] <- palette[palette %in% df$col]
##	##df$method=factor(df$method, order=TRUE)
##	if(blackBorder) border <- rep("black", nrow(df)) else border <- df$col
##	if(labelAllSamples) {
##		labels <- df$id
##		ticks.at <- df$y
##	} else {
##		labels <- FALSE
##		ticks.at <- pretty(df$y)
##	}
##	fig <- xyplot(y~midpoint|method, data=df,
##			    panel=function(x, y, x0, x1, chr.size,
##			    col, border, coverage, chr, show.coverage=TRUE, max.y,
##			    ..., subscripts){
##				    panel.grid(h=-1, v=10)
##				    ##yy <- factor(y, order=TRUE)
##				    panel.xyplot(x, y, ..., subscripts)
##				    ##yyy <- as.integer(yy)
##				    h <- 0.75
##				    lrect(xleft=x0[subscripts],
##					  xright=x1[subscripts],
##					  ybottom=y-h/2,
##					  ytop=y+h/2,
##					  border=border[subscripts],
##					  col=col[subscripts], ...)
##				    if(show.coverage)
##					    ltext(x, y,labels=coverage[subscripts], cex=0.6)
##				    ## plot centromere
##				    chr <- unique(as.integer(as.character(df$chr)))
##				    coords <- chromosomeAnnotation[chr, 1:2]/1e6
##				    lrect(xleft=coords[1],
##					  xright=coords[2],
##					  ybottom=0,
##					  ytop=max.y+h/2,
##					  col="grey",
##					  border="grey")
##			    },
##		      x0=df$x0,
##		      x1=df$x1,
##		      col=df$col,
##		      border=border,
##		      alpha=1,
##		      chr.size=df$chr.size,
##		      scales=list(y=list(labels=labels, at=ticks.at, cex=sampleLabels.cex)),
##		      coverage=df$coverage,
##		      xlab="Mb",
##		      ylab="offspring index",
##		      show.coverage=show.coverage,
##		      key=mykey,
##		      par.strip.text=list(lines=0.7, cex=0.6),
##		      prepanel=prepanel.fxn,
##		      xlim=xlim,
##		      max.y=max(df$y), ...)
####		      axis=function(side, text.cex){
####			      panel.axis(side, text.cex=text.cex)}, ...)
##	return(fig)
##}

cnEmission <- function(object, hmmOptions, k=3, verbose=TRUE){
	cnStates <- hmmOptions$copynumberStates##[["copynumberStates"]]
	verbose <- hmmOptions$verbose
	states <- hmmOptions$states
	is.log <- hmmOptions$is.log
	fn <- featureNames(object)
	S <- length(states)
	CN <- copyNumber(object)
	if(any(colSums(is.na(CN)) == nrow(CN))){
		stop("Some samples have all missing values. Exclude these samples before continuing.")
	}
	sds <- getSds(object)
	emission.cn <- array(NA, dim=c(nrow(object), ncol(object), S))
	if(is.log){
		MIN.CN <- -10
		MAX.CN <- 2.5
	} else {
		MIN.CN <- 0
		MAX.CN <- 10
	}
	if(any(CN < MIN.CN, na.rm=TRUE)) CN[CN < MIN.CN] <- MIN.CN
	if(any(CN > MAX.CN, na.rm=TRUE)) CN[CN > MAX.CN] <- MAX.CN
	for(j in 1:ncol(object)){
		cn <- CN[, j, drop=FALSE]
		sd <- sds[, j, drop=FALSE]
		I <- which(!is.na(as.numeric(cn)))
		old.tmp <- tmp <- rep(NA, length(as.numeric(cnStates)))
		cnvector <- as.numeric(cn)[I]
		prOutlier <- probabilityOutlier(cnvector, k=k)
		for(l in seq_along(cnStates)){
			mu <- cnStates[l]
			tmp <- (1-prOutlier) * dnorm(x=cnvector,
						     mean=cnStates[l],
						     sd=as.numeric(sd)[I]) +
							     prOutlier * dunif(cnvector, MIN.CN, MAX.CN)
			emission.cn[I, j, l] <- tmp
		}
	}
	return(log(emission.cn))
}

gtEmission <- function(object, hmmOptions){
	ICE <- hmmOptions$ICE
	states <- hmmOptions$states
	if(!ICE){
		p <- hmmOptions$prGtHom
		prGenotypeMissing <- hmmOptions$prGtMis
		verbose <- hmmOptions$verbose
		stopifnot(length(p) == length(states))
		if(!is.numeric(calls(object))) stop("genotypes must be integers (1=AA, 2=AB, 3=BB) or NA (missing)")
		GT <- calls(object)
		emission <- array(GT, dim=c(nrow(GT), ncol(GT), length(states)), dimnames=list(featureNames(object), sampleNames(object), states))
		missingGT <- any(is.na(GT))
		for(s in seq(along=states)){
			tmp <- GT
			tmp[tmp == 1 | tmp == 3] <- p[s]
			tmp[tmp == 2] <- 1-p[s]
			index1 <- is.na(tmp) & !isSnp(object)
			index2 <- is.na(tmp) & isSnp(object)
			if(missingGT){
				tmp[index2] <- prGenotypeMissing[s]
				## use uniform for nonpolymorphic
				tmp[index1] <- 1/length(states)
			}
			emission[, , s] <- tmp
		}
		logemit <- log(emission)
		##return(logemit)
	} else {
		##stop('need to update ICE option')
		logemit <- array(NA, dim=c(nrow(object), ncol(object), length(states)),
					 dimnames=list(featureNames(object),
					 sampleNames(object),
					 states))
		tmp <- genotypeEmissionCrlmm(object, hmmOptions)
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
}


probabilityOutlier <- function(cn, k=3, verbose){
	## outlier ~ N(0, sigma1), cn ~ N(0, sigma2), sigma2 << sigma1
	## lik= prod_i=1^N Pr(outlier) N(0, sigma1) + (1-Pr(outlier)) N(0, sigma2)
	rmeds <- runmed(cn, k)
	delta <- cn-rmeds
	mu <- 0
	sigma=c(0.5, 0.2)
	tau <- 0.01
	epsilon <- 2; counter <- 1
	while(epsilon > 0.01){
		gamma <- (tau * dnorm(delta, mu, sigma[1]))/ (tau * dnorm(delta, mu, sigma[1]) + (1-tau)*dnorm(delta, mu, sigma[2]))
		## gamma near 1 is likely an outlier
		## gamma near 0 is likely not an outlier
		tau.next <- mean(gamma)
		epsilon <- abs(tau.next - tau)
		tau <- tau.next
		counter <- counter+1
		if(counter > 10) break()
	}
	## plot(delta, gamma)
	return(gamma)
}




checkAnnotation <- function(object){
	if(!annotation(object) %in% icePlatforms()){
		stop("ICE is TRUE, but hapmap crlmm confidence scores for ", annotation(object), " are not available. Using crlmm confidence scores from HapMap samples assayed on the Affy 6.0 platform.")
	}
}


icePlatforms <- function(){
	c("pd.genomewidesnp.6",
	  "genomewidesnp6",
	  "pd.mapping250k.nsp",
	  "pd.mapping250k.sty",
	  "pd.mapping250k.nsp, pd.mapping250k.sty")
}


genotypeEmissionCrlmm <- function(object, hmmOptions){
	if(annotation(object) == "pd.genomewidesnp.6"){
		annotation <- "genomewidesnp6"
	} else annotation <- annotation(object)
	loader(paste(annotation, "Conf.rda", sep=""), .vanillaIcePkgEnv, "VanillaICE")
	hapmapP <- getVarInEnv("reference")
	pHetCalledHom <- hmmOptions[["prHetCalledHom"]]
	pHetCalledHet <- hmmOptions[["prHetCalledHet"]]
	pHomInNormal <- hmmOptions[["prHomInNormal"]]
	pHomInRoh <- hmmOptions[["prHomInRoh"]]
	if(length(annotation(object)) < 1) stop("must specify annotation")
	GT <- as.integer(calls(object))
	GTconf <- confs(object)
	##data(list=paste(annotation, "Conf", sep=""), package="VanillaICE", envir=environment())
	if(length(pHomInNormal) == nrow(GTconf)){  ##convert to vector
		pHomInNormal <- as.numeric(matrix(pHomInNormal, nrow(GTconf), ncol(GTconf), byrow=FALSE))
	} else pHomInNormal <-  rep(pHomInNormal, length(GT))
	hapmapP[, 2] <- 1-exp(-hapmapP[, 2]/1000)
	##p = 1-exp(-X/1000)
	##1000*log(1-p)=X
	##confidence <- 1-exp(-GTconf/1000)
	i11 <- hapmapP[, 1] == 3  ##called homozygous truth homozygous
	i12 <- hapmapP[, 1] == 4  ##called homozygous truth heterozygous
	i21 <- hapmapP[, 1] == 1  ##called het truth hom
	i22 <- hapmapP[, 1] == 2  ##called het truth het
	f11 <- density(hapmapP[i11, 2], from=0, to=1, n=1e3)
	f12 <- density(hapmapP[i12, 2], from=0, to=1, n=1e3)
	f21 <- density(hapmapP[i21, 2], from=0, to=1, n=1e3)
	f22 <- density(hapmapP[i22, 2], from=0, to=1, n=1e3)
	##-------------------------------------------------------------------------
	##distribution of observed call probabilities when the call is homozygous
	##-------------------------------------------------------------------------
	##-------------------------------------------------------------------------
	##P(phat | LOH, gte)
	##-------------------------------------------------------------------------
	##GT <- as.integer(genotypes)
	##confidence <- as.numeric(confidence)
	confidence <- as.numeric(GTconf)
	pTruthIsNormal <- pTruthIsRoh <- rep(NA, length(GT))
	confidence[confidence==0] <- 0.01 ##Otherwise, NA's result
	hom <- which(GT == 1 | GT == 3)
	observedPcalledHom <- cut(confidence[hom], breaks=f11$x, labels=FALSE)
	pTruthIsRoh[hom] <- f11$y[observedPcalledHom]
	het <- which(GT == 2)
	observedPcalledHet <- cut(confidence[het], breaks=f11$x, labels=FALSE)
	pTruthIsRoh[het] <- f21$y[observedPcalledHet]
	##-------------------------------------------------------------------------
	##Calculate P(phat | Normal, HOM)
	##-------------------------------------------------------------------------
	chet1 <- f22$y[cut(confidence[het], breaks=f22$x, labels=FALSE)]
	chet2 <- f21$y[cut(confidence[het], breaks=f21$x, labels=FALSE)]
	##term5[1]=P(true genotype is HET | genotype call is AB, state is normal)
	pTruthIsNormal[het] <- chet1*pHetCalledHet + chet2*(1-pHetCalledHet)
	##chom1=called homozygous truth heterozygous
	chom1 <- f12$y[cut(confidence[hom], breaks=f12$x, labels=FALSE)]
	##chom2=called homozygous truth homozygous
	chom2 <- f11$y[cut(confidence[hom], breaks=f11$x, labels=FALSE)]
	##chom4 <- 0.9999    ##P(HOM|CHOM)
	##probability that the true state is HOM when genotype call is homozygous
	##pHetCalledHom = P(true genotype is HET | calls is AA or BB, state is normal)
	pTruthIsNormal[hom] <- chom1*pHetCalledHom + chom2*(1-pHetCalledHom)
	fNormal <- fLoh <- rep(NA, length(GT))
	fNormal[hom] <- pHomInNormal[hom] * pTruthIsNormal[hom]
	fNormal[het] <- (1-pHomInNormal[het]) * pTruthIsNormal[het]
	fLoh[hom] <- pHomInRoh * pTruthIsRoh[hom]
	fLoh[het] <- (1-pHomInRoh) * pTruthIsRoh[het]
	f <- array(NA, dim=c(nrow(object), ncol(object), 2), dimnames=list(featureNames(object),
							     sampleNames(object),
							     c("normal", "ROH")))
	f[, , "normal"] <- matrix(fNormal, nrow(object), ncol(object))
	f[, , "ROH"] <- matrix(fLoh, nrow(object), ncol(object))
	f[f  == 0] <- min(f[f > 0], na.rm=TRUE)
	f <- log(f)
	return(f)
}

xypanel <- function(x, y,
		    gt,
		    is.snp,
		    range,
		    col.hom="grey60",
		    fill.hom="lightblue",
		    col.het="grey60" ,
		    fill.het="salmon",
		    col.np="grey20",
		    fill.np="green3",
		    ..., subscripts){
	panel.grid(v=0, h=4, "grey", lty=2)
	panel.xyplot(x, y, ...)
	is.snp <- is.snp[subscripts]
	gt <- gt[subscripts]
	hets.index <- which(gt == 2)
	hom.index <- which(gt == 1 | gt == 3)
	if(length(hom.index) > 0)
		lpoints(x[hom.index], y[hom.index], col=col.hom, fill=fill.hom,...)
	if(any(!is.snp))
		lpoints(x[!is.snp], y[!is.snp], col=col.np, fill=fill.np, ...)
	if(length(hets.index) > 0)
		lpoints(x[hets.index], y[hets.index], col=col.het, fill=fill.het, ...)
	j <- panel.number()
	lrect(xleft=start(range)[j]/1e6, xright=end(range)[j]/1e6,
	      ybottom=-10, ytop=10, ...)
}
