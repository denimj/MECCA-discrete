fitDiscrete <-
function(phy, data, model=c("ER", "SYM", "ARD"), treeTransform=c("none", "lambda", "kappa", "delta", "linearChange", "exponentialChange", "twoRate"), data.names=NULL, plotlnl=F, qLimits=c(0.0001, 1000), pLimits=c(0.00001, 10))
{
	model<-match.arg(model)
	
	treeTransform<-match.arg(treeTransform)
	if(treeTransform=="twoRate" & plotlnl==T) {
				cat("Plotting surfaces not supported for twoRate tree transformation\n")
				plotlnl=F
	}
	
	if(!is.ultrametric(phy)) {
		cat("Warning: some tree transformations in GEIGER might not be sensible for nonultrametric trees.\n")
		}

	if(hasZeroLengthTips(phy)) {
		cat("Warning: your tree has some zero-length tip branches. If the desendent species have different trait values, the likelihood will be zero and this approach will not work.\n")
		}

	td<-treedata(phy, data, data.names, sort=T)

	

	res<-list()

	for(i in 1:ncol(td$data))
	{

	
		if(treeTransform=="none") {
			f<-function(x) {
				likelihoodDiscrete(td$phy, td$data[,i], exp(x), model)
			}
			nep=0; pLow=-10; pHigh=log(1); pStart=NULL;		
		}	
		if(treeTransform=="lambda") {
			f<-function(x) {
				likelihoodDiscrete(td$phy, td$data[,i], exp(x[-1]), lambda=exp(x[1]), model=model)
			}	
			nep=1; pLow=-10; pHigh=log(1); pStart=0.1;
		}
		if(treeTransform=="delta") {
			f<-function(x) {
				likelihoodDiscrete(td$phy, td$data[,i], exp(x[-1]), delta=exp(x[1]), model=model)
				}
			nep=1; pLow=-10; pHigh=log(10);pStart=0.1;
		}
		if(treeTransform=="kappa") {
			f<-function(x) {
				likelihoodDiscrete(td$phy, td$data[,i], exp(x[-1]), kappa=exp(x[1]), model=model)
				}
			nep=1; pLow=-10; pHigh=log(1);pStart=0.1;
		}
		if(treeTransform=="linearChange") {
			f<-function(x) {
				likelihoodDiscrete(td$phy, td$data[,i], exp(x[-1]), endRate=exp(x[1]), linear=T, model=model)
				}
			nep=1; pLow=-10; pHigh=log(10);pStart=0.1;
		}
		if(treeTransform=="exponentialChange") {
			f<-function(x) {
				likelihoodDiscrete(td$phy, td$data[,i], exp(x[-1]), endRate=exp(x[1]), model=model)
				}
			nep=1; pLow=-10; pHigh=log(10);pStart=0.1;
		}
		if(treeTransform=="twoRate") {
			f<-function(x) {
				likelihoodDiscrete(td$phy, td$data[,i], exp(x[-(1:2)]), breakPoint=x[1], endRate=exp(x[2]), model=model)
				}
			mv<-max(branching.times(td$phy))	
			nep=2; pLow=c(mv/1000, -10); pHigh=c(mv, 10);pStart=c(0.1, 0.1);
		}
				
						
			nRateCats<-getRateCats(td$data[,i], model)
			
			outTries<-list()
			totalbl<-sum(td$phy$edge.length)
			minQ=log(0.01/totalbl)
			maxQ=log(1000/totalbl)
			ntries<-20
			ltry<-numeric(ntries)
			ltry[]<-NA
			lsol<-matrix(nrow= ntries, ncol=nRateCats+nep)
			sp<-numeric(nRateCats)
			qTries<-exp(-7:2)
			
			if(nep==0) {
				lower=rep(minQ, nRateCats)
				upper=rep(maxQ, nRateCats)
			} else {
				lower=c(pLow, rep(minQ, nRateCats))
				upper=c(pHigh, rep(maxQ, nRateCats))
			}
			
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

			if(treeTransform=="none" & model=="ER") {
				res[[i]]<-list(lnl=-out$value, q=-getQ(exp(out$par), nRateCats, model)[1,1], message=out$message)
			} else if(treeTransform=="none" & model!="ER") {
				res[[i]]<-list(lnl=-out$value, q=getQ(exp(out$par), nRateCats, model), message=out$message)

			} else if(treeTransform=="twoRate") {
				res[[i]]<-list(lnl=-out$value, q=getQ(exp(out$par[-(1:2)]), nRateCats, model), breakpoint=out$par[1], endRate=exp(out$par[2]), message=out$message)
			} else 	res[[i]]<-list(lnl=-out$value, q=getQ(exp(out$par[-1]), nRateCats, model), treeParam=exp(out$par[1]), message=out$message)

				
			if(!is.null(colnames(td$data))) names(res)[i]<-colnames(td$data)[i] else names(res)[i]<-paste("Trait", i, sep="")
		
		if((model=="SYM" | model=="ARD") & nRateCats>2) {
				cat("Plotting surfaces currently not supported for SYM and ARD models unless your character has only two states.\n")
				plotLnlSurf=F
		}
		
		if(plotlnl) {
			cat("Calculating surface\n")
			if(qLimits[1]<=0) {
				cat("Q must be positive, resetting lower plotting limit to 0.00000001")
				qLimits[1]=0.00000001
			}
			if(treeTransform=="none") {
				if(model=="ER") {
					qx<-exp(seq(log(qLimits[1]), log(qLimits[2]), length=50))
					lnl<-numeric(50)
					for(j in 1:50)
						lnl[j]<- -f(log(qx[j]))
					
					lnlDiff<- -out$value-lnl
					plot(qx, lnl, log="x", type="l", xlab="Rate (q)", ylab="lnL")
				} else {
					qx<-exp(seq(log(qLimits[1]), log(qLimits[2]), length=20))
					qy<-exp(seq(log(qLimits[1]), log(qLimits[2]), length=20))
					lnl<-matrix(nrow=20, ncol=20)
					for(j in 1:20)
						for(k in 1:20)
							lnl[j,k]<- -f(log(c(qx[j], qy[k])))
					
					lnlDiff<- -out$value-lnl
					contour(qx, qy, lnl, xlab="Forward Rate", ylab="Backward Rate")
				}	
			} else {
				px<-exp(seq(log(pLimits[1]), log(pLimits[2]), length=20))
				qy<-exp(seq(log(qLimits[1]), log(qLimits[2]), length=20))
				lnl<-matrix(nrow=20, ncol=20)
				for(j in 1:20)
					for(k in 1:20)
						try(lnl[j,k]<- -f(log(c(px[j], qy[k]))))
				
				lnlDiff<- -out$value-lnl
				contour(px, qy, lnlDiff, levels=c(1, 2, 3, 4, 8, 10, 20, 30, 40, 50, 100, 500, 1000), xlab="Tree Transform parameter estimate", ylab="Rate (q)", main="lnL Surface")		
						
				
			}
			
		}
		
	}
	return(res)

}
