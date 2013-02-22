getDiversificationModel <-
function (z, estimateExtinction=T, startR=0.05, startE=0.5) 
{
    isInt <- is.na(z[,6])
    isTerm <- !isInt
    int <- z[isInt,]
    term <- z[isTerm,]
    
    if(length(dim(int))==0) {int <- rbind(int)}
    if(length(dim(term))==0) {term <- rbind(term)}

    nint <- sum(isInt)
    nterm <- sum(isTerm)
    
    lalphaF <- function(b, d, aa, tt) {
    	num <- log(d)+log(exp((b-d)*tt)-1)
    	den <- log(b*(exp((b-d)*tt))-d)     	
    	# res<-aa*(num-den) # wrong according to dr
	res <- (num-den) # correction by dr - has worked through the code
    	res
    	}	
    
    lbetaF <- function(b, d, aa, tt) {
    	res <- lalphaF(b, d, aa, tt) + log(b) - log(d)
    	res
    }	
    
    lbetaPB <- function(b, aa, tt) {
    	res <- log((exp(b*tt)-1)/(exp(b*tt)))
    }


    Lfunc_comb <- function(r, eps) {
    	if(r<0 | eps<=0 | eps>=1) return(-Inf)

    	#bd <- getBD(r, eps) ## No need to use this function:
    	#b <- bd$b
    	#d <- bd$d
        d <- eps*r/(1-eps)
	b <- r+d
    	    	
    	if(nint==0) {lint=0} else {
    		branchLengths<-int[,3]-int[,4]
    		lint<-nint * 
            log(r) - r * sum(branchLengths) - sum(log(1 - (eps * 
            exp(-r * int[, 3]))))
        }
    	if(nterm==0) {lterm=0} else {
    		a<-term[,"startR"]
    		n<-term[,"endR"]
    		timeInterval<-term[,"startTime"]-term[,"endTime"]
    		endj<-pmin(a, n)
    		sum<-0
    		lnl<-numeric(nterm)
    		for(i in 1:nterm) {
    			lnxx<-numeric(length=endj[i])
    			for(j in 1:endj[i]) { #using likelihoods from Foote et al. 1999, a correction of Raup
    				logAlpha<-lalphaF(b, d, a[i], timeInterval[i])
    				logBeta<-lbetaF(b, d, a[i], timeInterval[i])
    				if(logBeta>0) logBeta=0
    				s1<-lchoose(a[i],j)+lchoose(n[i]-1,j-1)
    				
    				if(logAlpha==-Inf) s2<-0 else s2<-(a[i]-j)*logAlpha
    				s3<-log(((1-exp(logAlpha))*(1-exp(logBeta)))^j)
    				s4<-(n[i]-1)*logBeta
    				s5<-log(1-exp(logAlpha)) # Conditioning on survival to the present
       				lnxx[j]<-s1+s2+s3+s4-s5

    			}
    			lnl[i]<-logspace_sum(lnxx)
    		}
    		lterm<-sum(lnl)
    	}
    	#return(lterm)

    	return(lint+lterm)
	}
	
	Lfunc_comb_pb <- function(r) {
    	if(r<0) return(-Inf)

    	b<-r
    	d<-0
    	eps<-0
    	
    	if(nint==0) {lint=0} else {
    		branchLengths<-int[,3]-int[,4]
    		lint<-nint * 
            log(r) - r * sum(branchLengths) - sum(log(1 - (eps * 
            exp(-r * int[, 3]))))
        }
    	if(nterm==0) {lterm=0} else {
    		a<-term[,"startR"]
    		n<-term[,"endR"]
    		timeInterval<-term[,"startTime"]-term[,"endTime"]
    		endj<-pmin(a, n)
    		sum<-0
    		lnl<-numeric(nterm)
    		for(i in 1:nterm) {
    			lnxx<-numeric(length=endj[i])
    			for(j in 1:endj[i]) { #using likelihoods from Foote et al. 1999, a correction of Raup
    				logAlpha<- -Inf
    				logBeta<-lbetaPB(b, a[i], timeInterval[i])
    				s1<-lchoose(a[i],j)+lchoose(n[i]-1,j-1)
    				
    				if(logAlpha==-Inf) s2<-0 else s2<-(a[i]-j)*logAlpha
    				s3<-log(((1-exp(logAlpha))*(1-exp(logBeta)))^j)
    				s4<-(n[i]-1)*logBeta
    				s5<-log(1-exp(logAlpha)) # Conditioning on survival to the present
       				lnxx[j]<-s1+s2+s3+s4-s5

    			}
    			lnl[i]<-logspace_sum(lnxx)
    		}
    		lterm<-sum(lnl)
    	}
    	
    	return(lint+lterm)
	}

	
    res <- list()
    

    
    if (estimateExtinction == TRUE) {
     	foo <- function(x) 
     		- Lfunc_comb(r=exp(x[1]), eps=exp(x[2]))
		sp<-log(c(startR, startE))
		o<-optim(foo, par=sp, method="N")
		res$LH <- -o$value
    	res$par <- exp(o$par)

    } else {
    	foo <- function(x) 
     		- Lfunc_comb_pb(r=exp(x[1]))
		sp<-log(startR)
		o<-optimize(foo, interval=c(-100, 1))
		res$LH <- -o$objective
    	res$par <- exp(o$minimum)
    	}

	    



    return(res)
}
