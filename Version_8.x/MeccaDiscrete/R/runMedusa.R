runMedusa <-
function (phy, richness, estimateExtinction=T, modelLimit=20, cutAtStem=T, startR=0.05, startE=0.5, ...) 
{

	if(dim(richness)[2]==2) {
		nr<-dim(richness)[1]
		richness<-data.frame(richness[,1], rep(0, nr), richness[,2])
	}
    phy$node.label <- NULL
    
    #reset model limit if tree has fewer than 20 internal branches
	N<-length(phy$tip.label)
	if(modelLimit>(2*N-2)) modelLimit=2*N-2
   
  #holds parameter estimates  
    allRes<-matrix(nrow=modelLimit+1, ncol=6)
    
    root <- max(phy$edge) - phy$Nnode + 1
    node.list<-match(richness[,1], phy$tip.label)
    
    
    xx<-max(node.list)
    intNodes<-(xx+1):max(phy$edge)
    cs<-numeric(length=length(node.list)+length(intNodes))

	for(i in 1:length(node.list))
		cs[i]<-sum(node.list[1:i]==node.list[i])


    node.list<-c(node.list, intNodes)
    cs[node.list %in% intNodes]<-1

    
    

#First pass estimates a birth-death model that integrates phylogenetic and taxonomic likelihoods described in Rabosky et al., 2007      	  
	  
	res <- list()

	baseModel<-fitDiversification(phy, richness, estimateExtinction=T)

      allRes[1,]<-c(0, 0, baseModel$LH, baseModel$np, baseModel$aic, baseModel$aicc)
      
	#Now split the tree on all branches. 

      for (i in 1:length(node.list)) {
        r1 <- NULL
        r2 <- NULL
        z1 <- NULL
        z2 <- NULL
        z <- NULL
        z <- splitEdgeMatrixGeiger(phy, node.list[i], richness, cutAtStem, n2=cs[i])
        z1 <- z[z[, 7] == 1, ]
        z2 <- z[z[, 7] == 2, ]
        
        if(length(dim(z1))==0) z1<-rbind(z1)
        if(length(dim(z2))==0) z2<-rbind(z2)

	#Fit BD model to the two subclades of the tree        

        r1 <- getDiversificationModel(z1,  estimateExtinction, startR, startE)
        r2 <- getDiversificationModel(z2, estimateExtinction, startR, startE)

   	 	if(estimateExtinction) np<-5 else np<-3   
    
   
        k<-2*nrow(richness)-1

        res$node[i] <- node.list[i]
        res$LH[i] <- r1$LH + r2$LH
        res$aic[i] <- (-2 * res$LH[i]) + 2*np
        res$aicc[i]<-res$aic[i]+2*np*(np+1)/(k-np-1)

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
        
        res$np[i]<-np
        res$cs[i]<-cs[i]


      }

      bestModel<-which(res$LH==max(res$LH))[1]
      allRes[2,]<-c(res$node[bestModel], res$cs[bestModel], res$LH[bestModel], res$np[bestModel], res$aic[bestModel], res$aicc[bestModel])
      z<-splitEdgeMatrixGeiger(phy, node.list[bestModel], richness, cutAtStem, n2=cs[bestModel])
      
    
    for(j in 2:modelLimit) {
    	
      res <- list()
      for (i in 1:length(node.list)) {
 
        zNew <- resplitEdgeMatrixGeiger(z, phy, node.list[i], cutAtStem, n2=cs[i])# resplit edge matrix adds a new split to an already-split tree
        LH=0
        np=0
        for(k in 1:max(zNew[,7])){
        	zPart <- zNew[zNew[, 7] == k, ]
        	if(length(dim(zPart))==0) zPart<-rbind(zPart)
        	if(nrow(zPart)!=0) {
        		rPart <- getDiversificationModel(zPart, estimateExtinction, startR, startE)

        		LH<-LH+rPart$LH
        		if(estimateExtinction) np<-np+3 else np<-np+2
        	}
        	
		}
		np<-np-1 # otherwise you charge one extra for the background
    
   
             k<-2*nrow(richness)-1


        res$node[i] <- node.list[i]
        res$LH[i] <- LH
        res$aic[i] <- (-2 * res$LH[i]) + 2*np
        res$aicc[i]<-res$aic[i]+2*np*(np+1)/(k-np-1)

        res$np[i] <- np
        res$cs[i]<-cs[i]

      }

      bestModel<-which(res$LH==max(res$LH))[1]
      allRes[j+1,]<-c(res$node[bestModel], res$cs[bestModel], res$LH[bestModel], res$np[bestModel], res$aic[bestModel], res$aicc[bestModel])
      
      z <- resplitEdgeMatrixGeiger(z, phy, node.list[bestModel], cutAtStem, n2=cs[bestModel])
      
    }
    
    return(allRes)
}
