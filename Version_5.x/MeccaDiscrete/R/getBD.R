getBD <-
function(r, eps) {
	d<-eps*r/(1-eps)
	b<-r+d
	return(list(b=b, d=d))
	}
