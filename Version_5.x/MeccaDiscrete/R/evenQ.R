evenQ <-
function(n)

{

	q<--diag(n)

	q[q==0]<-1/(n-1)

	return(q)

}
