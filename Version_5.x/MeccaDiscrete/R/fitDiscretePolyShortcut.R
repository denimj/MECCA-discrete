fitDiscretePolyShortcut <-
function(phyFast, model=c("ER", "SYM", "ARD"), data.names=NULL, plotlnl=F, qLimits=c(0.0000001, 1000), pLimits=c(0.00001, 10))
{
	
	model<-match.arg(model)

	f<-function(x) {
				likelihoodDiscretePolyShortcut(phyFast, exp(x), model)
			}
	nep=0; pLow=-10; pHigh=log(1); pStart=NULL;		
		
				
						
			nRateCats<-getRateCats(factor(1:2), model)
			
			outTries<-list()
			totalbl<-sum(phyFast$edge.length)
			minQ=log(qLimits[1])
			maxQ=log(qLimits[2])
			ntries<-20
			ltry<-numeric(ntries)
			ltry[]<-NA
			lsol<-matrix(nrow= ntries, ncol=nRateCats+nep)
			sp<-numeric(nRateCats)
			qTries<-exp(-7:2)
			
				lower=rep(minQ, nRateCats)
				upper=rep(maxQ, nRateCats)
	
			
			cat("Finding the maximum likelihood solution\n")
			cat("[0        50      100]\n")
			cat("[")
			
			for(j in 1:10) {
				sp<-c(pStart, log(rep(qTries[j], nRateCats)))
				te<-try(outTries[[j]]<-optim(f, par=sp, method="L",  lower=lower, upper=upper), silent=T)
				
				if(class(te)!="try-error") {
					ltry[j]<-outTries[[j]]$value
					lsol[j,]<-exp(outTries[[j]]$par)
				}
				cat(".")

			}
			for(j in 1:10) {
				sp<-c(pStart, runif(nRateCats, minQ, maxQ))
				te<-try(outTries[[10+j]]<-optim(f, par=sp, method="L",  lower=lower, upper=upper), silent=T)
				if(class(te)!="try-error") {
					ltry[10+j]<-outTries[[10+j]]$value
					lsol[10+j,]<-exp(outTries[[10+j]]$par)
				}
				cat(".")

			}
			
			cat("]\n")
			
			ok<-!is.na(ltry)
			
			if(sum(ok)==0) stop("ERROR: No solution found. Does your tree contain zero-length tip branches?")
			
			ltd<-ltry-min(ltry[ok])
			b<-min(which(ltry==min(ltry[ok])))

			gc<-which(ltd<0.1)
			us<-lsol[gc,1]
			usc<-sum((us-min(us))>0.1)			
			b<-min(which(ltry==min(ltry[ok])))
			out<-outTries[[b[1]]]	
			if(usc>1) {out$message="Warning: likelihood surface is flat."}
			
			if(out$convergence!=0) {out$message="Warning: may not have converged to a proper solution."}

			if(out$convergence==0) {out$message="R thinks that this is the right answer."}

			if(model=="ER") {
				res<-list(lnl=-out$value, q=-getQ(exp(out$par), nRateCats, model)[1,1], message=out$message)
			} else {
				res<-list(lnl=-out$value, q=getQ(exp(out$par), nRateCats, model), message=out$message)

			} 
				

		
	
	return(res)

}
