simpleABC <-
function(phy, tipStates, tol=10, ntries=10000, eps=0) {
	mm<-match(phy$tip.label, rownames(tipStates))

	
	
	posterior<-numeric()
	for(i in 1:ntries) {
	# draw parameters from the prior
	q12<-runif(1, min=0.001, max=2)
	q21<-runif(1, min=0.001, max=2)
	
	tipD<-rowSums(tipStates)
	
	# simulate dataset
	simData<-simDiscreteIncompleteTree(phy, tipD, q12, q21, eps)
	x<-log(as.numeric(simData[mm,1])+1)
	y<-log(as.numeric(tipStates[,1])+1)
	
	dd<-dist(rbind(x, y))
	
	if(dd<tol) {
		posterior<-rbind(posterior, c(q12, q21))
		cat(".")
		} else cat("-")
		
	} 
	posterior
}
