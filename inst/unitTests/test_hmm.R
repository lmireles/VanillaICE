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
	alpha <- alpha/rowSums(alpha) ## now a probability
	g <- array(NA, dim=c(T, 2, S))## gamma
	for(s in seq_len(S)){
		beta1cn <- prOutlier*dunif(cn, MIN.CN, MAX.CN)
		## need also the BAF emission
		den <- beta1cn+(1-prOutlier)*dnorm(cn, mus[s], sds)
		d1 <- beta1cn/den
		d2 <- 1-d1
		g[, 1, s] <- alpha[, s] * d1
		g[, 2, s] <- alpha[, s] * d2
	}
	## updates
	## pmix is a 2 x 6 matrix
	g_1 <- apply(g[, 1, ], 2, sum, na.rm=TRUE)
	g_2 <- apply(g[, 2, ], 2, sum, na.rm=TRUE)
	g_pout <- g_1/(g_1+g_2)
	g_qout <- 1-g_pout
	## betat <- cnEmission(prOutlier, O)  ## T x 6
	## gamma <-                          ## T x 2 x 6.  Second dimension is prOutlier
	## cjk  <-                           ## 2 x 6       prOutlier for each state
	mubar <- matrix(NA, 2, 6)
	mubar1.pout <- (g[, 1, 1] * cn)/sum(g[,1,1], na.rm=TRUE)
	mubar2.qout <- (g[, 2, 1] * cn)/sum(g[, 2, 1], na.rm=TRUE)



	mu_pout <- apply(g[, 1, ]*cnm, 2, sum,na.rm=TRUE)/g_qout
	mu_qout <- apply(g[, 2, ]*cnm, 2, sum, na.rm=TRUE)/g_pout

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



