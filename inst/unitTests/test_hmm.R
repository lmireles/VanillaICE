test_hmm_oligoSnpSet <- function(){
	set.seed(1)
	states <- as.integer(c(3, 4, 3, 5, 3, 2, 3, 2, 3, 2, 3))
	nmarkers <- as.integer(c(996, 102, 902, 50, 2467, 102, 1898,
				 99, 900, 10, 160))
	oligoset <- VanillaICE:::artificialData(states, nmarkers)
	hmmOpts <- hmm.setup(oligoset, is.log=FALSE)
	fit <- hmm(oligoset, hmmOpts, k=3)
	checkIdentical(state(fit), states)
	## check that the breakpoints are close
	checkTrue(all.equal(coverage2(fit), nmarkers, tolerance=0.01))

	hmmOpts <- hmm.setup(oligoset, is.log=FALSE)
	## fit is sensitive to choice of prOutlier and p.hom
	fit2 <- hmm(oligoset, hmmOpts, use.baf=TRUE, prOutlier=1e-3, p.hom=0.95, k=3)
	checkIdentical(state(fit2), states)
	checkTrue(all.equal(coverage2(fit2), nmarkers, tolerance=0.01))

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
	checkTrue(all.equal(coverage2(rd), nmarkers, tolerance=0.01))
}



