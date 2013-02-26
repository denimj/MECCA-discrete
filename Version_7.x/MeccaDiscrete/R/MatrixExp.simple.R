MatrixExp.simple <-
function(Q)
{
	n<-nrow(Q)
	res<-matrix(0, nrow=n, ncol=n)
	q<-Q[1,2]
	for(i in 1:n)
		res[i, i]<-1/n+(n-1)/n*exp(-n*q)
	res[res==0]<-1/n-1/n*exp(-n*q)
	return(res)
}
