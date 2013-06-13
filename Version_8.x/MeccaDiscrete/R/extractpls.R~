extractpls <-
function(p, plsdat, Ncomp) {
	x <- numeric(Ncomp);
	for(i in 1:Ncomp) {	
		x[i] <- sum(p * plsdat[, i])
		}
	return(x)
	}
