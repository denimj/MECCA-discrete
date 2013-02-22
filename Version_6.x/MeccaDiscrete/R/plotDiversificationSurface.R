plotDiversificationSurface <-
function (phy, nPoints=10, rInterval=c(0.00001, 0.3), eInterval=NULL, logTransform=T) 
{
    z <- splitEdgeMatrixGeiger(phy, phy$Nnode)
	rootnode = (length(phy$tip.label) + 1)
 	int <- z[z[, 2] > rootnode, ]
    term <- z[z[, 2] < rootnode, ]
    nint <- nrow(int)
    nterm <- nrow(term)
    betaF <- function(r, eps, t1) {
    	if(r<0 | eps<0 | eps>=1) return(-Inf)
        xf <- (exp(r * t1) - 1)/(exp(r * t1) - eps)
        xf
    }

    Lfunc_comb <- function(r, eps) {
    	if(r<0 | eps<0 | eps>=1) return(-Inf)

        (sum(log(1 - betaF(r, eps, term[1:nterm, 4]))) + sum((term[1:nterm, 
            5] - 1) * log(betaF(r, eps, term[1:nterm, 4]))) + nint * 
            log(r) - r * sum(int[1:nint, 4]) - sum(log(1 - (eps * 
            exp(-r * int[1:nint, 3])))))
    }
	
	if(logTransform) {
			x<-exp(seq(from=log(rInterval[1]), to=log(rInterval[2]), length.out=nPoints))
	} else x<-seq(from=rInterval[1], to=rInterval[2], length.out=nPoints)
	
	if(is.null(eInterval)) {
		lnl<-numeric(nPoints)
		for(i in 1:nPoints)
			lnl[i]<-Lfunc_comb(r=x[i], eps=0)
		plot(x, lnl, type="l")	
	} else {
		if(logTransform) {
			y<-exp(seq(from=log(eInterval[1]), to=log(eInterval[2]), length.out=nPoints))
		} else y<-seq(from=eInterval[1], to=eInterval[2], length.out=nPoints)
		lnl<-matrix(nrow=nPoints, ncol=nPoints)
		for(i in 1:nPoints)
			for(j in 1:nPoints)
				lnl[i,j]<-Lfunc_comb(r=x[i], eps=y[j])
		mm<-max(lnl)		
		contour(lnl, levels=c(mm-1, mm-2, mm-3, mm-4, mm-5, mm-10, mm-100), labels=c(1, 2, 3, 4, 5, 10, 100), axes=F, xlab="b-d", ylab="d/b")
		tics<-floor(c(1, nPoints/4, nPoints/2, nPoints*3/4, nPoints))
		axis(1, at=c(0, 0.25, 0.5, 0.75, 1), labels=round(x[tics], 3))
		axis(2, at=c(0, 0.25, 0.5, 0.75, 1), labels=round(y[tics], 3))
	}
}
