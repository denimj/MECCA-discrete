fitDiversification <-
function (phy, richness, estimateExtinction=T, ...) 
{
	
	if(dim(richness)[2]==2) {
		nr<-dim(richness)[1]
		richness<-data.frame(richness[,1], rep(0, nr), richness[,2])
	}
    z <- splitEdgeMatrixGeiger(phy, NA, richness, cutAtStem=T)
    r1 <- getDiversificationModel(z, estimateExtinction, ...)
    if(estimateExtinction) np<-2 else np<-1    
    res <- list()
    res$LH <- r1$LH
    res$aic <- (-2 * r1$LH) + 2*np
    n<-nrow(z)
    res$np<-np
    res$k<-n
    res$aicc<-res$aic+2*np*(np+1)/(n-np-1)
    if(estimateExtinction) eps<-r1$par[2] else eps=0
    
    res$r <- r1$par[1]
    res$eps <- eps
    return(res)
}
