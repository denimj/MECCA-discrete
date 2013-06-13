simBranchingDiscrete <-
function(n, age, lambda, mu, q12, q21, startState=1, seed=1) {

  ## The 'ran_ecu' function that generates the random numbers has a default seed that is supposed to initialize the generator if no seed is provided.
	qMatrix <- rbind(c(-q12, q12), c(q21, -q21))
	
	wt <- getBdWaitTimes(n=n, age=age, lambda=lambda, mu=mu)
        ## Waiting times is the calculation of the branch lengths between all the nodes that will be created for the clade tip (when it get resolved). It is important to note that the first value of wainting times is the age of the clade, that is, the branch length of the clade tip (unresolved clade tip).
	ss <- noTreeSimFast(wt, qMatrix, startState, seed)
	#res <- list(ns1=ss$nState1, ns2=ss$nState2, seed=ss$seed)
        res <- list(ns1=ss$nState1, ns2=ss$nState2, seed=ss$seed)
	res
	}
