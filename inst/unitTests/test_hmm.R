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
	fit <- hmm(oligoset, hmmOpts, k=3)
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



