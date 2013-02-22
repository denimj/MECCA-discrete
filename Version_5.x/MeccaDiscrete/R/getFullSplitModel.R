getFullSplitModel <-
function(phy, estimateExtinction=T, breakList, richness, cutAtStem=T, csList){

	res<-list()
	phy$node.label <- NULL
	
	root <- max(phy$edge) - phy$Nnode + 1
    node.list<-match(richness[,1], phy$tip.label)
    
    
    xx<-max(node.list)
    intNodes<-(xx+1):max(phy$edge)

    x <- branching.times(phy)

	# First break


    z <- splitEdgeMatrixGeiger(phy, breakList[1], richness, cutAtStem=cutAtStem, n2=csList[1])
	if(length(breakList)>1)
	  for(i in 2:length(breakList)) {
	  	kk<-which(node.list==breakList[i])

		z <- resplitEdgeMatrixGeiger(z, phy, breakList[i], cutAtStem=cutAtStem, n2=csList[i])

		}
	
	 for(k in 1:max(z[,7])){
        	zPart <- z[z[, 7] == k, ]
        	if(length(dim(zPart))==0) zPart<-rbind(zPart)

        	if(nrow(zPart)!=0) {
        		rPart <- getDiversificationModel(zPart, estimateExtinction)
        		res[[k]]<-rPart
        		int<-is.na(zPart[,6])
        		if(sum(int)==0)  {
        			nn<-phy$tip.label[zPart[,2]]
        			
        	        res[[k]]$taxa<-nn

        		} else {
        			tt<-zPart[!int, 2]
        			res[[k]]$taxa<-phy$tip.label[tt]

        			}
        
        	}
        	
		}
		list(res=res, z=z)
}
