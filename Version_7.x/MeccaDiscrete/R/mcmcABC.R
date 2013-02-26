mcmcABC <-
function(phy, tipStates, tol=10, ngen=10000, propRange=2, rtPropRange=0.2, eps=0) {
	mm<-match(phy$tip.label, rownames(tipStates))
	posterior<-matrix(nrow=ngen, ncol=3)
	
	qmin<-log(0.01)
	qmax<-log(2)
	q12<-runif(1, min= qmin, max= qmax)
	q21<-runif(1, min=qmin, max= qmax)
	p1root<-runif(1, min=0, max=1)
	posterior[1,]<-c(exp(q12), exp(q21), p1root)
	tipD<-rowSums(tipStates)

	plot(1, ylim=c(exp(qmin), exp(qmax)), xlim=c(1,ngen), log="y")

	for(i in 2:ngen) {
	# draw parameters from the prior
		nq12<-q12+runif(1, min=-propRange, max= propRange)
		nq21<-q21+runif(1, min=-propRange, max= propRange)
	
	while(1) {
		np1root<-p1root+runif(1, min=-rtPropRange, max= rtPropRange)
		if(np1root>0 & np1root<1) break

		}
	
	rrt<-runif(1)
	if(runif(1)<np1root) rootState=1 else rootState=2
	# simulate dataset
	simData<-simDiscreteIncompleteTree(phy, tipD, exp(nq12), exp(nq21), rootState, eps)
	
	x<-log(as.numeric(simData[,1])+1)
	y<-log(as.numeric(tipStates[mm,1])+1)
	
	dd<-dist(rbind(x, y))
	cat(dd,"\n")
	if(dd<tol) {
		q12<-nq12; q21<-nq21;p1root=np1root;
			cat(".")
		} else cat("-")
	
	posterior[i,]<-c(exp(q12), exp(q21), p1root)
	points(y=posterior[i,], x=rep(i, 3), col=1:3, pch=19, cex=0.5)
		
	} 
	posterior
}
