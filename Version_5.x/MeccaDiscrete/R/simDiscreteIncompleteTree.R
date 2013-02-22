simDiscreteIncompleteTree <-
function(phy, tipD, q12, q21, rootState, r, eps, seed=1) {
	mm<-match(phy$tip.label, names(tipD))

	pp <- pruneBL(phy, tipD)
	prunedPhy <- pp[[2]]
	tipLen <- pp[[1]]
	
	qMatrix <- rbind(c(-q12, q12), c(q21, -q21))
	ss <- sim.char(prunedPhy, model.matrix=list(qMatrix), model="discrete", root.state= rootState)[,1,1]

	res <- matrix(nrow=length(ss), ncol=2)
	rownames(res) <- names(ss)
	
	for(i in 1:length(ss)) {
		tl <- tipLen[names(ss)[i]]
		ntaxa <- tipD[names(ss)[i]]
		if(ntaxa==1) {
			## Identify the resolved terminals (only one taxa).
			res[i,] <- 0
			res[i,ss[i]] <- 1
		} else {
			## Identify the unresolved terminals.
			bd <- getBD(r, eps)
			sb <- simBranchingDiscrete(n=ntaxa, age=tl, lambda=bd$b, mu=bd$d, q12=q12, q21=q21, startState=ss[i]) ## I took of the 'seed' argument that was here. But see 'simBranchingDiscrete.R'.
			res[i,] <- c(sb$ns1, sb$ns2)
			seed <- sb$seed
			}
		
		}
		
	return(list(res=res, seed=seed))
}
