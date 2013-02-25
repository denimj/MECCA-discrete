noTreeSimFast <-
function(times, qMatrix, startState, seed) {
	q12<-qMatrix[1,2]
	q21<-qMatrix[2,1]
	#nState1<-0; nState2<-0
        nState1<-0; nState2<-0; vezes <- 0; loops <- 0;
        ntimes <- length(times)
        times <- round(times, digits=2)
	
	#out<-.C("noTreeSim", times=as.double(times), q12=as.double(q12), q21=as.double(q21), startState=as.integer(startState), seed=as.integer(seed), nState1=as.integer(nState1), nState2=as.integer(nState2),PACKAGE="MeccaDiscrete")
        out<-.C("noTreeSim", ntimes=as.integer(ntimes), times=as.double(times), q12=as.double(q12), q21=as.double(q21), startState=as.integer(startState), seed=as.integer(seed), nState1=as.integer(nState1), nState2=as.integer(nState2),vezes=as.integer(vezes),loops=as.integer(loops),PACKAGE="MeccaDiscrete")
	return(out)
}
