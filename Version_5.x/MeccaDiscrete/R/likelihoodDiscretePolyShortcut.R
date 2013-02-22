likelihoodDiscretePolyShortcut <-
function(phyFast, q, model="ER")
{
	
	Q<-getQ(q, n=2, model=model)
	
	tmp <- as.numeric(phyFast$edge)
	nb.tip <- max(tmp) #number of tips
	nb.node <- -min(tmp) #number of internal nodes
	nb.states <- 2 #fast version restricted to 2-state characters
	l <- matrix(-Inf, nrow=nrow(phyFast$edge), ncol=nb.states) #makes matrix that will store likelihoods of each character state at each node
	root <- numeric(nb.states) #makes vector that will store root likelihoods
	
	new2old.phylo(phyFast)->phyFast
	
	for(i in 1:nrow(phyFast$edge)) #for each edge
		if(as.numeric(phyFast$edge[i,2])>0) {
			whichTip<-as.numeric(phyFast$edge[i,2])
			div<-phyFast$tipDiversity[whichTip]
			trait1<-phyFast$trait1[whichTip]
			ll<-phyFast$tipLength[whichTip]
			
			trait2<-div-trait1
			
			
			#cat(phyFast$tip.label[whichTip], phyFast$tipDiversity[whichTip], trait1, trait2, "\n")
			# for tips with state 1
			tip.like<-matrix(c(0,-Inf), nrow=1)
			lt1<-trait1*frag.like(tip.like, ll, Q)
		
			# for tips with state 2
			tip.like<-matrix(c(-Inf,0), nrow=1)
			lt2<-trait2*frag.like(tip.like, ll, Q)


			l[i,]<-lt1+lt2
		}
		
		# fix for zero-length tips, which return NaNs instead of 0 in some cases
		
		l[is.na(l)]<- 0
		
			#if the edge is connected to a terminal taxon, you set the likelihood of the tip value equal to 1 and all others equal to zero.
		times <- branching.times(old2new.phylo(phyFast)) #get node to tip distances
		-1*(1:(max(as.numeric(names(times)))-min(as.numeric(names(times)))+1))->names(times)
		times = max(times) - times #convert into root to node tips

	while(1) {
		
		if(sum(as.numeric(phyFast$edge[,2])>0)==2) break 
	
		#obtain ancestors of current tips
		x <- match(phyFast$edge[,1], phyFast$edge[,2])[as.numeric(phyFast$edge[,2])>0] #finds nodes connected to terminal taxa
		#find last node with two tip descendent
		a <- max(x[duplicated(x)])
		t <- which(phyFast$edge[,1]== phyFast$edge[a,2])
		bl <- phyFast$edge.length[t]
		age = times[which(names(times)== phyFast$edge[a,2])]
		l[a,] <- frag.like(l[t,], bl, Q)
		#next line effectively prunes out the tips just used
		phyFast$edge[a,2]<-1
		phyFast$edge[t,2]<-0

	
	}
	t <- which(as.numeric(phyFast$edge[,2])>0)
	bl <- phyFast$edge.length[t]
	root <- frag.like(l[t,], bl, Q)
	neglnl=-logspace_sum(root)+log(nb.states)
	return(neglnl)
}
