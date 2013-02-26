pruneBL <-
function(phy, tipDiversity)
{
	#m<-match(names(tipDiversity), phy$tip.label)
	tipBranches<-numeric(length=length(tipDiversity))
	for(i in 1:length(tipDiversity)) {
		if(tipDiversity[i]>1) {
			nt<-which(phy$tip.label==names(tipDiversity)[i])
			nn<-which(phy$edge[,2]==nt)
			tipBranches[i]<-phy$edge.length[nn]
			phy$edge.length[nn]<-0
		}
		
	}
	names(tipBranches)<-names(tipDiversity)

	return(list(tipBranches, phy))
}
