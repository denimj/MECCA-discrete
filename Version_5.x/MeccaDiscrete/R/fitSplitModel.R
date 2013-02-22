fitSplitModel <-
function (phy, estimateExtinction=T,     startR=0.05, startE=0.5) 
{
    phy$node.label <- NULL
    root <- max(phy$edge) - phy$Nnode + 1
    node.list <- 1:max(phy$edge)
    node.list <- node.list[node.list != root]
    x <- branching.times(phy)
    res <- list()
    for (i in 1:length(node.list)) {
        r1 <- NULL
        r2 <- NULL
        z1 <- NULL
        z2 <- NULL
        z <- NULL
        z <- splitEdgeMatrixGeiger(phy, node.list[i])
        z1 <- z[z[, 6] == 1, ]
        z2 <- z[z[, 6] == 2, ]
        
        r1 <- getDiversificationModel(z1, estimateExtinction, startR, startE)
        r2 <- getDiversificationModel(z2, estimateExtinction, startR, startE)

   	 	if(estimateExtinction) np<-4 else np<-2    
    
   
        

        res$node[i] <- node.list[i]
        res$LH[i] <- r1$LH + r2$LH
        res$aic[i] <- (-2 * res$LH[i]) + 2*np*2
        
        if(estimateExtinction) {
        	eps1<-r1$par[2]
        	eps2<-r2$par[2] 
 
        } else {
        	eps1=0
            eps2=0
        }
        
        
        res$r1[i] <- r1$par[1]
        res$e1[i] <- eps1
        res$LH1[i] <- r1$LH
        
        res$r2[i] <- r2$par[1]
        res$e2[i] <- eps2
        res$LH2[i] <- r2$LH

    }
    res <- as.data.frame(res)
    return(res)
}
