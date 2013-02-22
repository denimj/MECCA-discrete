simDiscreteIncompleteTree <-
function(phy, tipD, q12, q21, rootState, r, eps, seed=1) {
	mm <- match(phy$tip.label, names(tipD))

	pp <- pruneBL(phy, tipD)
	prunedPhy <- pp[[2]] ## Phy whose tip clades terminal branch length are equall to 0.
	tipLen <- pp[[1]] ## Vector with the values of the branch length of the clade tips.
	
	qMatrix <- rbind(c(-q12, q12), c(q21, -q21))
	ss <- sim.char(prunedPhy, model.matrix=list(qMatrix), model="discrete", root.state= rootState)[,1,1] ## geiger - Note that each clade tip will receive one value of the simulated char. Therefore, they are treated as normal tips.

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
                  ## I took out the 'getBD' function and put the calculations instead.
                        d <- eps*r/(1-eps)
                        b <- r+d
                        sb <- simBranchingDiscrete(n=ntaxa, age=tl, lambda=b, mu=d, q12=q12, q21=q21, startState=ss[i], seed) ## I took of the 'seed' argument that was here. But see 'simBranchingDiscrete.R'. Note that 'ss[i]' is equal to the startState at the node of a prunned clade. This function is simulating a bd tree with the number of tips equall to ntaxa and its parameters (age, lambda and mu). Then is simulating the evolution of a qualitative character (1 & 2) using the q12 and q21 rates. However, I would expect this to generate a vector of char values with length->ntaxa. But this is not needed, since Stadler T developed a new approach that do not involve simulating the trees and, therefore, is faster than simulating trees and then sampling them.
			res[i,] <- c(sb$ns1, sb$ns2)
                        }
		
		}
        seed <- sb$seed
        colnames(res) <- c("nState1","nState2")
		
	return(list(res=res, seed=seed))
}
