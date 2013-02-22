summaryMedusa <-
function(phy, richness, out, cutoff=4, plotTree=T, useCorrection=F, cutAtStem=T, estimateExtinction=T, ...) {
	
	if(dim(richness)[2]==2) {
		nr<-dim(richness)[1]
		richness<-data.frame(richness[,1], rep(0, nr), richness[,2])
	}
	breaks<-numeric()
	csList<-numeric()
	i=2
	ml<-dim(out)[1]-1
	if(useCorrection) col=6 else col=5
	while(1) {
		if(i>ml) break;
		if((out[i-1,col]-out[i,col])<cutoff) break;
		breaks[i-1]<-out[i,1]
		csList[i-1]<-out[i,2]
		i<-i+1
		}
	rr<-getFullSplitModel(phy, estimateExtinction= estimateExtinction, breakList=breaks, richness=richness, cutAtStem=cutAtStem, csList=csList)
	
	if(plotTree) {
		mm<-match(phy$edge[,2], rr$z[,2])
		ec<-rr$z[mm,7]
		plot(phy, edge.color=ec, ...)
		}
	
	
	rr$res	
	}
