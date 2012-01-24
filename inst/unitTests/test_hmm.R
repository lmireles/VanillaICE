test_hmm_oligoSnpSet <- function(){
	## Results are sensitive to the seed, and differs between the unit test and my
	states <- as.integer(c(3, 4, 3, 5, 3, 2, 3, 2, 3, 2, 3))
	nmarkers <- as.integer(c(996, 102, 902, 50, 2467, 102, 1898,
				 99, 900, 20, 160))
	statepath <- rep(states, nmarkers)
	if(FALSE){
		oligoset <- VanillaICE:::artificialData(states, nmarkers)
		save(oligoset, file="../extdata/oligosetForUnitTest.rda")
	} else {
		path <- system.file("extdata", package="VanillaICE")
		load(file.path(path, "oligosetForUnitTest.rda"))
	}
	hmmOpts <- hmm.setup(oligoset, is.log=FALSE)
	##trace(VanillaICE:::cnEmissionFromMatrix, browser)
	##trace(hmm, browser, signature=c("oligoSnpSet", "HmmOptionList"))
	##untrace(hmm, signature=c("oligoSnpSet", "HmmOptionList"))
	##trace(VanillaICE:::hmm2, signature=c("oligoSnpSet", "HmmOptionList"))
	##trace(VanillaICE:::cnEmissionFromMatrix, browser)
	##trace(VanillaICE:::viterbi, browser)
	hmmOpts[["copynumberStates"]] <- c(0, 1, 2, 2, 2.5, 3)
	trace(VanillaICE:::viterbi.wrapper, browser)
	fit <- hmm(oligoset, hmmOpts, k=3)

	object <- oligoset
	hmm.params <- hmmOpts
	mus <- c(0,1,2,2,2.5,3)
	cn <- copyNumber(object)
	sds <- VanillaICE:::.getSds(cn)
	is.log <- FALSE
	rr <- range(cn, na.rm=TRUE, finite=TRUE)
	if(is.log){
		MIN.CN <- pmax(-10, rr[1])
		MAX.CN <- pmin(2.5, rr[2])
	} else {
		MIN.CN <- pmax(0, rr[1])
		MAX.CN <- pmin(10, rr[2])
	}
	emitr <- VanillaICE:::cnEmissionFromMatrix(cn,
						   cnStates=mus,
						   is.snp=isSnp(object),
						   normalIndex=3,
						   is.log=FALSE,
						   prOutlier=0.01)
	b <- baf(object)
	emitb <- VanillaICE:::gtEmission(object, hmm.params)
	log.beta=emitr+emitb
	d <- diff(position(object))
	TAUP <- hmm.params[["tau"]]
	e <- 0.5
	pis <- rep(0, 6)
	pis[3] <- e
	pis[-3] <- (1-e)/5
	log.initial <- log(pis)
	tp <- exp(-2*d)/TAUP
	c1 <- 1
	S <- 6L
	minimum <- 1-1/((S-1)*c1) + 0.01
	tp <- pmax(tp, minimum)
	T <- nrow(object)
	arm <- rep(1L, T)
	vit <- rep(0L, T)
	delta <- matrix(0.0, T, S)
	res <- .C("viterbi",
		  log.emission=as.matrix(as.numeric(log.beta[, 1, ])),
		  log.initial=log.initial,
		  transitionPr=as.matrix(tp),
		  arm=arm,
		  S=S,
		  T=T,
		  result=vit,
		  delta=as.matrix(as.numeric(delta)),
		  normal2altered=1,
		  altered2normal=1,
		  altered2altered=1,
		  normalIndex=3,
		  pAA=rep(0, S^2))
	alpha <- exp(matrix(res[["delta"]], T, S))
	alpha <- alpha/rowSums(alpha, na.rm=T)
	## lik = -log sum(ct) ,where c is the scaling variable
	## rescale: has to be within the viterbi algorithm.
	## use same ct for backward variable
##	lalpha <- matrix(res[["delta"]], T, S)
##	ct <- 1/sum(lalpha[1, ])
##	for(t = 2:T){
##		ct[t] <- 1/sum(alpha[t, ])
##		alpha[t, ] <- prod(ct[1:t])*alpha[t, ]
##	}
##
##	alpha <- alpha/rowSums(alpha) ## now a probability
	g_b <- g_r <- matrix(NA, T, S)## gamma
##	p.out <- prOutlierBAF
##	q.out <- 1-p.out
##	if("pb" %in% names(list(...))){
##		pb <- list(...)[["pb"]]
##		pb <- pb/100
##		pb[is.na(pb)] <- 0.5
##	} else pb <- rep(0.5, nrow(object))
	beta.cn_outlier <- matrix(dunif(cn, 0, 5), T, S)
	p.r <- rep(0.01, S)
	beta.cn <- exp(emitr[, 1, ])
	g_r <- alpha*(beta.cn*(1-p.r))/((1-p.r)*beta.cn + p.r*(beta.cn_outlier))

	p.b <- rep(1e-3, S)
	beta.b_outlier <- matrix(dunif(b, 0, 1), T, S)
	beta.b <- exp(emitb[,1, ])
	g_b <- alpha*(beta.b*(1-p.b))/((1-p.b)*beta.b + p.b*(beta.b_outlier))

	## updates
	## pmix is a 2 x 6 matrix
	tmp <- apply(g_r, 2, sum, na.rm=TRUE)
	tmp2 <- apply(1-g_r, 2, sum, na.rm=TRUE)
	p.r <- tmp2/(tmp+tmp2)

	tmp <- apply(g_b, 2, sum, na.rm=TRUE)
	tmp2 <- apply(1-g_b, 2, sum, na.rm=TRUE)
	p.b <- tmp2/(tmp+tmp2)

	mubar <- rep(NA, 6)
	for(s in 1:6){
		mubar[s] <- sum(g_r[, s]*cn, na.rm=TRUE)/sum(g_r[, s], na.rm=TRUE)  ## mean for each t.
	}
	if(any(diff(mubar) < 0)){ ## must be ordinal
		index <- which(diff(mubar) < 0)
		mubar[index] <- mubar[index+1]
	}
	sdsbar <- rep(NA, 6)
	for(s in 1:6){
		sdsbar[s] <- sqrt(sum(g_r[, s]*(cn - mubar[s])^2, na.rm=TRUE)/sum(g_r[, s], na.rm=TRUE))
	}
	## run viterbi
	## recomput alpha
	## reestimate emission probs. with updated mus and updated sds
	## compute objective function
	## P(O|lambda is sum of alphas)





	## recompute emission probabilities.
	## rerun viterbi.
	o <- findOverlaps(featureData(oligoset), fit[5, ])
	copyNumber(oligoset)[queryHits(o), ]

	T <- 10
	S <- 5
	x <- matrix(0, T, S)
	#x[T, ] <- 1:S
	x[T, ] <- 1:5

	xx <- matrix(x, ncol=1, byrow=TRUE)
	y <- matrix(0, T, S)
	y[1, ] <- 1:5
	yy <- matrix(y, ncol=1,byrow=TRUE)


	checkTrue(identical(state(fit), states))
	checkEquals(coverage2(fit), nmarkers, tolerance=0.01)
	##checkIdentical(coverage2(fit), as.integer(c(996, 102, 902, 54, 2463, 102, 1898, 99, 901, 19, 160)))

	hmmOpts <- hmm.setup(oligoset, is.log=FALSE)
	## fit is sensitive to choice of prOutlier and p.hom
	fit2 <- hmm(oligoset, hmmOpts, use.baf=TRUE, prOutlier=1e-3, p.hom=0.95, k=3)
	checkIdentical(state(fit2), states)
	checkTrue(isTRUE(all.equal(coverage2(fit2), nmarkers, tolerance=0.01)))

	cnemit <- VanillaICE:::cnEmission(object=copyNumber(oligoset),
					  cnStates=c(0, 1, 2, 2, 3, 4),
					  k=3,
					  is.log=FALSE,
					  is.snp=isSnp(oligoset),
					  normalIndex=3)
	bafemit <- VanillaICE:::bafEmission(object=baf(oligoset),
					    is.snp=rep(TRUE, nrow(oligoset)),
					    prOutlier=1e-3,
					    p.hom=0.95)
	LE <- cnemit+bafemit
	log.initial <- log(rep(1/6, 6))
	arm <- VanillaICE:::.getArm(chromosome(oligoset), position(oligoset))
	rdl <- VanillaICE:::viterbi3(arm=arm,
				     pos=position(oligoset),
				     chrom=chromosome(oligoset),
				     LE=LE[, 1, ],
				     log.initial=log.initial,
				     states=1:6,
				     id=sampleNames(oligoset),
				     TAUP=1e8)
	rd <- VanillaICE:::stackRangedData(rdl)
	checkIdentical(state(rd), states)
	checkTrue(isTRUE(all.equal(coverage2(rd), nmarkers, tolerance=0.01)))
}

internalTest_hmm <- function(){
}



