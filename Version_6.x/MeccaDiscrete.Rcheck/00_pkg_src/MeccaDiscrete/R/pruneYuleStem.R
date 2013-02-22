pruneYuleStem <-
function(phy, tipDiversity)
{
	m<-match(names(tipDiversity), phy$tip.label)
	for(i in 1:length(m)) {
		if(tipDiversity[m][i]>1) {
			nt<-which(phy$tip.label==names(tipDiversity)[m][i])
			nn<-which(phy$edge[,2]==nt)
			t<-phy$edge.length[nn]
			n<-tipDiversity[m][i]
			stem<-t/sum(1/(1:n))
			phy$edge.length[nn]<-stem
		}
		
	}
	return(phy)
}
