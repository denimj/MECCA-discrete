getQ <-
function(q, n, model)
{
	if(model=="ER") Q=evenQ(n)*q
	if(model=="SYM") {
		if(length(q)!=n*(n-1)/2) stop("You must supply the correct number of rate categories.")
		Q<-diag(n)
		xx=1
		for(i in 2:n) {
			for(j in 1:(i-1)) {
				Q[i,j]<-Q[j,i]<-q[xx]
				xx<-xx+1
			}
		}	
		for(i in 1:n) diag(Q)[i]<- -sum(Q[i,-i])
	}
	
	if(model=="ARD") {
		if(length(q)!=n*(n-1)) stop("You must supply the correct number of rate categories.")
		Q<-diag(n)
		xx=1
		for(i in 1:n) {
			for(j in (1:n)[-i]) {
				Q[i,j]<-q[xx]
				xx<-xx+1
			}
		}	
		for(i in 1:n) diag(Q)[i]<- -sum(Q[i,-i])
	}
	
	return(Q)
}
