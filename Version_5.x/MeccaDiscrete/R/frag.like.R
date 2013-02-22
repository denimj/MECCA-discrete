frag.like <-
function(tip.like, bl, q)

{

	nb.states<-ncol(tip.like)

	r<-rep(0, nb.states)

	d<-length(bl)

	p<-list(d)

	for(i in 1:d)

		p[[i]]<-MatrixExp.eig(q*bl[i])

	for(i in 1:nb.states)

		for(j in 1:d) 

			r[i]<-r[i]+logspace_sum(log(p[[j]][i,])+tip.like[j,])

	return(r)

}
