splitEdgeMatrixGeiger <-
function (phy, node, richness, cutAtStem=T, n2=1) 
{
    bt <- branching.times(phy)
    rootnode <- length(phy$tip.label) + 1
    
    #m1<-match(richness[,1], phy$tip.label)
    #m2<-match(phy$edge[,1], names(bt))

    interior<-phy$edge[,2] %in% phy$edge[,1]
 	tips<-!interior 
 	
 	z<-matrix(ncol=7, nrow=sum(interior)+nrow(richness))
    colnames(z)<-c("anc", "dec", "startTime", "endTime", "startR", "endR", "Group")
  
 	for(i in 1:sum(interior)) {
 		anc<-phy$edge[interior,1][i]
 		dec<-phy$edge[interior,2][i]
 		startTime<-bt[names(bt)==anc]
 		endTime<-startTime-phy$edge.length[interior][i]
 		sr<-1
 		er<-NA
 		group<-NA
 		z[i,]<-c(anc, dec, startTime, endTime, sr, er, group)
 		} 
  
  	m1<-match(richness[,1], phy$tip.label)

  
	for(i in 1:nrow(richness)) {
		wt<-phy$tip.label[m1][i]
		we<-which(phy$edge[,2]==m1[i])
		anc<-phy$edge[we,1]
 		dec<-phy$edge[we,2]
 		startTime<-bt[names(bt)==anc]
 		endTime<-richness[i,2]
 		sr <- 1
 		er<-richness[i,3]
 		group<-NA
 		z[i+sum(interior),]<-c(anc, dec, startTime, endTime, sr, er, group)
		}
	 
	for(i in 1:max(z[,2])) {
		nn<-which(z[,2]==i)
		o<-order(z[nn,4], decreasing=T)
		if(length(nn)>1) for(j in 2:length(nn)) {
			z[nn,][o[j],3]<-z[nn,][o[j-1],4]
			z[nn,][o[j],5]<-z[nn,][o[j-1],6]
			}
		}
		
	dec2<-numeric(length(z[,2]))
	for(i in 1:max(z[,2])) {
		ss<-z[,2]==i
		nn<-sum(ss)
		oo<-order(z[ss,4], decreasing=T)
		dec2[ss]<-oo
		}
	
		
	if(is.na(node)) {
		z[,7]<-1
	} else {
		setZ<-function(z, n, n2, atTipFlag=T) {
			if(n %in% z[,2]) {
				if(atTipFlag) {
       				ok<-z[,2]==n & dec2 >= n2
       			} else {
       				ok<-z[,2]==n
       			}
       			z[ok,7]<-2
				ns<-node.sons(phy, n)
				if(length(ns)!=0) for(k in 1:length(ns))
					z=setZ(z, ns[k], n2, atTipFlag=F)
			} else if(n %in% z[,1]) z[,7]<-2
			z
		}
		
		z<-setZ(z, node, n2)
        
    	z[is.na(z[,7]),7]<-1
    }            
    
    if(!cutAtStem) {
    	row<-which(z[, 1] == node)
    	if(length(row)>0)
			z[row,7]<-1
    	}
       
    return(z)
}
