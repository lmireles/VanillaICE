##Sweave2pdf <- function(fname, ...){
##	fname <- strsplit(fname, ".Rnw")[[1]][[1]]
##	suppressWarnings(Sweave(paste(fname, ".Rnw", sep=""), ...))
##	texi2dvi(paste(fname, ".tex", sep=""), pdf=TRUE)
##}

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
##robustSds2 <- function(x, DF.PRIOR=10, nSamples, takeLog=FALSE, ...){
##	if(!is.matrix(x)) stop("x is not a matrix")
##	if(takeLog) x <- log2(x)
##
##	## the posterior is IG(nu_n/2, nu_n * sigma2_n (theta)/2)
##	## nu_n = nu_0 + n
##	## sigma2_n(theta) = 1/nu_n * [nu_0 * sigma^2_0 + n*s^2_n(theta)],
##	## where s^2_n(theta) = Sum(yi-theta)^2/n
##	##
##	##
##	## 1/sigma^2 ~ G(nu_0/2, nu_0/2 * sigma^2_0)
##	## E(sigma^2) = sigma^2_0 * nu_0/2 /(nu_0/2 -1)
##	##
##	## 1/sigma^2 | ... ~ G(nu_n/2, nu_n * sigma^2_n/2)
##	##
##	## nu_n = nu_0 + n
##	## sigma^2_n = 1/nu_n *[nu_0*sigma^2_0 + (n-1)*s^2 + k_0*n/k_n*(ybar-mu_0)^2]
##	##
##	## prior
##	##
##	## data
##	##
##	## posterior inference
##	sigma2.0 <- apply(x, 2, MAD, na.rm=TRUE)
##	nu0 <- 100
##	nn <- ncol(x)
##	nu.n <- nu0+nn
##	s2 <- rowMAD(x, na.rm=TRUE)
##	k0 <- 1 ## ?
##	kn <- k0+nn
##	ybar <- rowMedian(x, na.rm=TRUE)
##	mu0 <- median(ybar)
##	##
##	## E[1/sigma2_g | ... ] = sigma^2n/2 * nu_n / (nu_n/2 - 1)
##
##
##
##
##	##
####	sigma.marker <- rowMAD(x, na.rm=TRUE)
####	sigma.sample <- apply(x, 2, MAD, na.rm=TRUE)
####	gammahat <- sigma.marker/median(sigma.marker, na.rm=TRUE)
####	df1 <- nSamples-1
####	df2 <- length(sds1)-1
##	sds1 <- rowMAD(x, na.rm=TRUE)
####	sds.marker <- (df1*sds.marker + df2*median(sds.marker,na.rm=TRUE))/(df1+df2)
##	sds1 <- matrix(sds1, nrow(x), ncol(x))
##	sds2 <- apply(x, 2, "mad", constant=2, na.rm=TRUE)
##	df <- ncol(x)
##	sds2 <- matrix(sds2, nrow(x), ncol(x), byrow=TRUE)
##	sds.star <- (sds1 * df + sds2*DF.PRIOR)/(df+DF.PRIOR)
##	dimnames(sds.star) <- dimnames(x)
##	return(sds.star)
##}

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

computeLoglikForRange <- function(from, to,
				  viterbiSequence, normalIndex,
				  log.initial,
				  log.emission,
				  lP.A2N,
				  lP.N2A,
				  lP.N2N,
				  lP.A2A
				  ){
	index <- seq(from, to)
	thisState <- unique(viterbiSequence[index])
	if(thisState == normalIndex) return(0)
	first.index <- min(index)
	last.index <- max(index)
	T <- length(viterbiSequence)
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
		LLR <- logLik.vit-logLik.null
		return(LLR)
	}
	##2 index = t
	if(length(index) == 1 & first.index > 1 & last.index < T){
		##note the last term cancels
		logLik.vit <- sum(lP.N2A[index] + log.emission[index, thisState]) + lP.A2N[last.index+1]+log.emission[last.index+1, normalIndex]
		logLik.null <- sum(lP.N2N[index] + log.emission[index, normalIndex]) + lP.N2N[last.index+1] + log.emission[last.index+1, normalIndex]
		LLR <- logLik.vit-logLik.null
		return(LLR)
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
		LLR <- logLik.vit-logLik.null
		return(LLR)
	}
	## 4: index = T
	if(first.index == last.index & last.index == T){
		logLik.vit <- lP.N2A[T] + log.emission[T, thisState]
		logLik.null <- lP.N2N[T] + log.emission[T, normalIndex]
		LLR <- logLik.vit-logLik.null
		return(LLR)
	}
	## 5: index = 1, 2, ...
	if(first.index != last.index & first.index == 1 & last.index < T){
		index2 <- index[-1]## t=2, ...., t*
		logLik.vit <- log.initial[thisState] + log.emission[first.index, thisState]  + sum(lP.A2A[index2] + log.emission[index2, thisState]) + lP.A2N[last.index+1] + log.emission[last.index+1, normalIndex]
		logLik.null <- log.initial[normalIndex] + log.emission[first.index, normalIndex] + sum(lP.N2N[index2] + log.emission[index2, normalIndex]) + lP.N2N[last.index+1] + log.emission[last.index+1, normalIndex]
		LLR <- logLik.vit-logLik.null
		return(LLR)
	}
	if(first.index != last.index & first.index == 1 & last.index == T){
		index2 <- index[-1]## t=2, ...., t*
		logLik.vit <- log.initial[thisState] + log.emission[first.index, thisState]  + sum(lP.A2A[index2] + log.emission[index2, thisState]) ##+ lP.A2N[last.index+1] + log.emission[last.index+1, normalIndex]
		logLik.null <- log.initial[normalIndex] + log.emission[first.index, normalIndex] + sum(lP.N2N[index2] + log.emission[index2, normalIndex]) ##+ lP.N2N[last.index+1] + log.emission[last.index+1, normalIndex]
		LLR <- logLik.vit-logLik.null
		return(LLR)
	}
	## 6: index = t, ...T
	if(first.index != last.index & last.index == T){
		index2 <- index[-1]
		logLik.vit <- lP.N2A[first.index] + log.emission[first.index, thisState] + sum(lP.A2A[index2] + log.emission[index2, thisState])
		logLik.null <- lP.N2N[first.index] + log.emission[first.index, normalIndex] + sum(lP.N2N[index2] + log.emission[index2, normalIndex])
		LLR <- logLik.vit - logLik.null
		return(LLR)
	}
	stop("none of conditions in likRatio function satisfied")
}

computeLoglik <- function(viterbiResults,
			  log.initial,
			  log.emission=log.emission,
			  states,
			  normalIndex,
			  nNotMissing,
			  c1, c2, c3,
			  physical.position,
			  CHR,
			  sample.name){
	viterbiSequence <- viterbiResults[["viterbiSeq"]]
	S <- length(states)
	rl <- Rle(viterbiSequence)
	starts <- start(rl)
	LLR <- rep(NA,  length(starts))
	log.emission <- matrix(viterbiResults[["log.emission"]], nNotMissing, S)
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
	for(k in seq_along(starts)){
		LLR[k] <- computeLoglikForRange(from=start(rl)[k],
						to=end(rl)[k],
						viterbiSequence=viterbiSequence,
						normalIndex=normalIndex,
						log.initial=log.initial,
						log.emission=log.emission,
						lP.N2N=lP.N2N,
						lP.N2A=lP.N2A,
						lP.A2N=lP.A2N,
						lP.A2A=lP.A2A)
	}
	start.index <- start(rl)
	end.index <- end(rl)
	##pos <- position(object)[I]
	##this is tricky since we've added an index to force a segment for each arm.
	start <- physical.position[start.index]
	end <- physical.position[end.index]
	##numMarkers <- unlist(numMarkers)
	numMarkers <- width(rl)
	states <- viterbiSequence[start.index]
	ir <- IRanges(start=start, end=end)
	rangedData <- RangedData(ranges=ir,
				 chromosome=rep(CHR, length(LLR)),
				 sampleId=rep(sample.name, length(LLR)),
				 state=states,
				 coverage=numMarkers,
				 LLR=LLR)
	return(rangedData)
}

viterbi <- function(object,
		    hmm.params,
		    log.E, ...){
	sns <- colnames(log.E)
	if(is.null(sns)) stop("no dimnames for log.emission")
	log.initial <- hmm.params$log.initialPr
	verbose <- hmm.params$verbose
	normal2altered <- hmm.params$n2a
	altered2normal <- hmm.params$a2n
	altered2altered <- hmm.params$a2a
	if(normal2altered <= 0) stop("normal2altered must be > 0")
	if(altered2normal <= 0) stop("altered2normal must be > 0")
	if(altered2altered <= 0) stop("altered2altered must be > 0")
##	##
##	## arm can contain NA's for invalid chromosomes or NA's
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
	rangedData <- vector("list", ncol(log.E))
	for(j in 1:ncol(log.E)){
		rD <- vector("list", length(unique(arm)))
		missingE <- rowSums(is.na(log.E[, j, ])) > 0
		notFinite <- rowSums(!is.finite(log.E[, j, ])) > 0
		missingE <- missingE | notFinite
		for(a in seq_along(unique(arm))){
			logicalNotMissing <- I <- arm == a  & !missingE
			if(sum(logicalNotMissing) < 2) next()
			nNotMissing <- T <- sum(logicalNotMissing)
			physical.position <- position(object)[logicalNotMissing]
			transitionPr <- exp(-2 * diff(physical.position)/TAUP)
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
			rD[[a]] <- computeLoglik(viterbiResults,
						 log.initial=log.initial,
						 log.emission=log.emission,
						 states=states,
						 normalIndex=normalIndex,
						 nNotMissing=nNotMissing,
						 c1=c1, c2=c2, c3=c3,
						 physical.position=physical.position,
						 CHR=unique(chromosome(object)[logicalNotMissing]),
						 sample.name=sns[j])
		}
		notnull <- !sapply(rD, is.null)
		rD <- rD[notnull]
		if(length(rD)==1) {
			rangedData[[j]] <- rD[[1]]
		} else{
			rangedData[[j]] <- stack(RangedDataList(rD))
		}
	}
	if(length(rangedData) == 1) {
		rangedData <- rangedData[[1]]
	} else {
		rangedData <- stack(RangedDataList(rangedData))
	}
	return(rangedData)
}






centerAutosomesAt <- function(x, at, ...){
	stopifnot(!missing(at))
	marker.index <- which(chromosome(x) <= 23)
	cn <- copyNumber(x)[marker.index, , drop=FALSE]
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


invalidCnConfidence <- function(x){
	is.na(x) | x <= 0 | is.nan(x) | is.infinite(x)
}

invalidGtConfidence <- function(x){
	is.na(x) | x < 0 | x > 1 | is.nan(x) | is.infinite(x)
}

getSds <- function(object, na.rm=TRUE){
	cn.conf <- cnConfidence(object)
	stopifnot(all(chromosome(object) <= 24))
	notvalid <- invalidCnConfidence(cn.conf)
	CN <- copyNumber(object)
	if(any(notvalid)){
		##if(verbose) message("cnConfidence missing.  Using MAD")
		marker.index <- which(chromosome(object) < 23)
		if(length(marker.index) == 0){
			sds <- matrix(NA, nrow(cn.conf), ncol(cn.conf))
			## sex chromosomes
			marker.index.list <- split(seq_len(nrow(object)), chromosome(object))
			for(i in seq_along(marker.index.list)){
				marker.index <- marker.index.list[[1]]
				tmp <- apply(CN[marker.index, , drop=FALSE], 2, mad, na.rm=TRUE)
				sds[marker.index, ] <- matrix(tmp, length(marker.index), ncol(CN), byrow=TRUE)
			}
		}  else {  ## autosomes present
			s <- apply(CN[marker.index, , drop=FALSE], 2, mad, na.rm=TRUE)
			sds <- matrix(s, nrow(CN), ncol(CN), byrow=TRUE)
		}
		valid <- !notvalid
		if(any(valid)){
			sds[valid] <- 1/cn.conf[valid]
		}
	} else { ## all valid confidence scores
		sds <- 1/cn.conf
	}
	notvalid <- invalidCnConfidence(sds)
	stopifnot(any(!notvalid))
	return(sds)
}

validChromosomeIndex <- function(object){
	index <- which(chromosome(object) <= 24 & !is.na(chromosome(object)) & !is.na(position(object)))
	if(length(index) < 1){
		stop("Chromosome must be 1-24 and physical position \n
                      can not be missing.  See chromosome() and position().\n
                      See integer2chromosome(1:24) for integer codes used by\n
                      VanillaICE.")
	}
	return(index)
}







probabilityOutlier <- function(cn, sigma=c(0.5, 0.2), k=3, verbose=TRUE){
	## outlier ~ N(0, sigma1), cn ~ N(0, sigma2), sigma2 << sigma1
	## lik= prod_i=1^N Pr(outlier) N(0, sigma1) + (1-Pr(outlier)) N(0, sigma2)
	if(length(cn) > k){
		rmeds <- runmed(cn, k)
		delta <- cn-rmeds
	} else delta <- cn-median(cn, na.rm=TRUE)
	mu <- 0
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


genotypeEmissionCrlmm <- function(object, hmm.params, gt.conf, cdfName){
	stopifnot(is(object, "matrix"))
	GT <- as.integer(object)
	##rm(object); gc()
	if(cdfName == "pd.genomewidesnp.6"){
		annotation <- "genomewidesnp6"
	} else annotation <- cdfName
	loader(paste(annotation, "Conf.rda", sep=""), .vanillaIcePkgEnv, "VanillaICE")
	hapmapP <- getVarInEnv("reference")
	pHetCalledHom <- hmm.params[["prHetCalledHom"]]
	pHetCalledHet <- hmm.params[["prHetCalledHet"]]
	pHomInNormal <- hmm.params[["prHomInNormal"]]
	pHomInRoh <- hmm.params[["prHomInRoh"]]
	if(length(cdfName) < 1) stop("must specify annotation")
	##data(list=paste(annotation, "Conf", sep=""), package="VanillaICE", envir=environment())
	if(length(pHomInNormal) == nrow(gt.conf)){  ##convert to vector
		pHomInNormal <- as.numeric(matrix(pHomInNormal, nrow(gt.conf), ncol(gt.conf), byrow=FALSE))
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
	##confidence <- as.numeric(gt.conf)
	confidence <- gt.conf
	rm(gt.conf); gc()
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
	f <- array(NA, dim=c(nrow(object), ncol(object), 2)) ##dimnames=list(featureNames(object),
							##     sampleNames(object),
							  ##   c("normal", "ROH")))
	f[, , 1] <- matrix(fNormal, length(object), ncol(object))
	f[, , 2] <- matrix(fLoh, length(object), ncol(object))
	dimnames(f)[[3]] <- c("normal", "ROH")
	f[f  == 0] <- min(f[f > 0], na.rm=TRUE)
	f <- log(f)
	return(f)
}

xypanel <- function(x, y,
		    gt,
		    is.snp,
		    range,
		    col.hom="grey20",
		    fill.hom="lightblue",
		    col.het="grey20" ,
		    fill.het="salmon",
		    col.np="grey20",
		    fill.np="grey60",
		    show.state=TRUE,
		    ..., subscripts){
	panel.grid(v=0, h=4, "grey", lty=2)
	panel.xyplot(x, y, ...)
	is.snp <- is.snp[subscripts]
	if(!missing(gt)){
		gt <- gt[subscripts]
		hets.index <- which(gt == 2)
		hom.index <- which(gt == 1 | gt == 3)
		if(all(!c("col", "fill") %in% names(list(...)))){
			if(any(!is.snp))
				lpoints(x[!is.snp], y[!is.snp], col=col.np,
					fill=fill.np, ...)
			if(length(hom.index) > 0)
				lpoints(x[hom.index], y[hom.index], col=col.hom,
					fill=fill.hom, ...)

			if(length(hets.index) > 0)
				lpoints(x[hets.index], y[hets.index],
					col=col.het,
					fill=fill.het, ...)
		}
	} else {
		lpoints(x[!is.snp], y[!is.snp], col=col.np,
			fill=fill.np, ...)
		## use whatever col.hom to color SNPs
		lpoints(x[is.snp], y[is.snp], col=col.hom,
			fill=fill.hom, ...)
	}
	j <- panel.number()
	st <- start(range)[j]/1e6
	lrect(xleft=st, xright=end(range)[j]/1e6,
	      ybottom=-10, ytop=10, ...)
	if(show.state){
		## left justify the label to the start of the range
		y.max <- current.panel.limits()$ylim[2]
		ltext(st, y.max, labels=paste("state", state(range)[j]),
		      adj=c(0,1))
	}
}


##
## AD-HOC
updateMu <- function(x, mu, sigma, is.snp, nUpdates=10){
	if(nUpdates==0) return(mu)
	## assume CN is a vector.  Fit EM independently for each
	## sample
	## - might want to do separately for snps, nps
	##
	## compute the responsibilities
	##
	sigma <- sigma[1]
	##pi <- exp(log.initialPr)
	dup.index <- which(duplicated(mu))
	S <- length(mu)
	if(length(dup.index) > 0){
		mu <- unique(mu)
		L <- length(mu)
	} else L <- S
	## fix normal copy number
	mu[3] <- median(x, na.rm=TRUE)
	pi <- rep(1/L, L)
	## fix the sd.  Update the means via em.
	##gamma <- vector("list", L)
	gamma <- matrix(NA, length(x), L)
	den <- matrix(NA, length(x), L)
	##num <- vector("list", L-1)
	num <- matrix(NA, length(x), L)
	epsilon <- 0.01
	for(iter in seq_len(nUpdates)){
		if(iter > nUpdates) break()
		for(i in seq_len(L)){
			den[, i] <- pi[i]*dnorm(x, mean=mu[i], sd=sigma)
		}
		D <- rowSums(den, na.rm=TRUE)
		for(i in seq_len(L)){
			num <- den[, i] ##pi[i] * dnorm(x, mu[i], sigma)
			gamma[, i] <- num/D
		}
		rs <- rowSums(gamma, na.rm=TRUE)
		gamma <- gamma/rs
		total.gamma <- apply(gamma, 2, sum, na.rm=TRUE)
		##
		## update the means with contraints
		##
		mu.new <- rep(NA, length(mu))
		mu.new[3] <- mu[3]
		I <- c(1,2, 4, 5)
		for(i in I){
			if(sum(gamma[, i],na.rm=TRUE) < 0.0001) {
				mu.new[i] <- mu[i]
				next()
			}
			tmp <- sum(gamma[, i] * x, na.rm=TRUE)/total.gamma[i]
			if(i > 1 & i < L){
				## mu[i-1]+sigma < mu[i] < mu[i+1] - sigma
				i1 <- tmp < (mu[i-1] + 1.5*sigma)
				i2 <- tmp > (mu[i+1] - 1.5*sigma)
				if(i1 | i2){
					if(i1)## & tmp < (mu[i+1] -sigma)){
						mu.new[i] <- mu[i-1]+1.5*sigma
					if(i2)
						mu.new[i] <- mu[i+1]-1.5*sigma
				} else mu.new[i] <- tmp
			}
			if(i == 1){
				mu.new[1] <- ifelse(tmp < (mu[2] - 1.5*sigma), tmp, mu[2]-1.5*sigma)
			}
			if(i == L){
				mu.new[L] <- ifelse(tmp > mu[L-1] + 1.5*sigma, tmp, mu[L-1]+1.5*sigma)
			}
		}
		pi.new <- apply(gamma, 2, mean, na.rm=TRUE)
		pi <- pi.new
		##dp <- abs(sum(mu - mu.new)) + abs(sum(pi.new-pi))
		dp <- abs(sum(mu-mu.new))
		mu <- mu.new
		if(dp < epsilon) break()
	}
	if(length(dup.index) > 0){
		tmp <- rep(NA, S)
		tmp[-dup.index] <- mu
		tmp[dup.index] <- tmp[dup.index-1]
		mu <- tmp
	}
	return(mu)
}

arrangeSideBySide <- function(object1, object2){
	grid.newpage()
	lvp <- viewport(x=0,
			y=0.05,
			width=unit(0.50, "npc"),
			height=unit(0.95, "npc"), just=c("left", "bottom"),
			name="lvp")
	pushViewport(lvp)
	nfigs1 <- length(object1$x.limit)
	nfigs2 <- length(object2$x.limit)
	stopifnot(length(nfigs1) == length(nfigs2))
	pushViewport(dataViewport(xscale=c(0,1), yscale=c(0.05,1), clip="on"))
	object1$layout <- c(1, nfigs1)
	print(object1, newpage=FALSE, prefix="plot1", more=TRUE)
	upViewport(0)
	lvp2 <- viewport(x=0.5,
			 y=0.05,
			 width=unit(0.50, "npc"),
			 height=unit(0.95, "npc"), just=c("left", "bottom"),
			 name="lvp2")
	pushViewport(lvp2)
	pushViewport(dataViewport(xscale=c(0,1), yscale=c(0.05,1), clip="on"))
	object2$layout <- c(1, nfigs1)
	print(object2, newpage=FALSE, prefix="plot2", more=TRUE)
}
