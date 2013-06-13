hasZeroLengthTips <-
function(phy)
{
	ntips<-length(phy$tip.label)
	tips<-phy$edge[,2]<=ntips
	nn<-phy$edge.length[tips]==0
	if(sum(nn)>0) return(T)
	return(F)
	}
