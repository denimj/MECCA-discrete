pruneZeroStem <-
function(phy, tipDiversity)
{
	m<-match(names(tipDiversity), phy$tip.label)
	for(i in 1:length(m)) {
		if(tipDiversity[m][i]>1) {
			nt<-which(phy$tip.label==names(tipDiversity)[m][i])
			nn<-which(phy$edge[,2]==nt)

			phy$edge.length[nn]<-0
		}
		
	}
	return(phy)
}
