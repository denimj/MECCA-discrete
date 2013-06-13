collapseTreeForMECCA <-
function(phy, trait=NULL, timeFraction=0.5, endClades=0) {

### Question! What does 'endClades' arg do?
  
  	rich <- rep(1, length(phy$tip.label)) #This is the number of tips
  	names(rich) <- phy$tip.label
	rootDepth  <- max(branching.times(phy)) #Gives the distance from the tips (any) to the root. Should not work in unrooted trees. Therefore, the phy argument needs a rooted tree.
	threshold <- rootDepth*timeFraction # gives the threshold for collapsing - younger than this gets collapsed - This case collapsing all clades at 'halfway' of the phylogeny.

	if(timeFraction==0 & endClades==0)
		return(NULL)

	while(1) { ## it's going to iterate through until there are no nodes left that satisfy the threshold requirement ## Notice that the condition of the 'while' is always satisfied.
		bt <- branching.times(phy)

		if(timeFraction != 0) {
			ct <- which(bt < threshold) # which branches events are too young ## Remember: bt are computed from tip to root.
			if(length(ct) == 0) break; ## ct is a counter and the cycle need to take elements from it.
			nodeOrder <- order(bt[ct], decreasing=T) ## order all the branching events that are to be collapsed ## The index 'ct' only gives the number of elements of 'bt' to be ordered. Given that you ordered 'bt', the lenght of 'ct' is all values that you need. Notice that when you order a vector with a vector with less elements, the resulting vector will have the lenght of the shorter one. This should work the same as 'order(bt[1:lenght(ct)], decreasing = T)' but it may work faster, for you are using less functions. The other option would be use 'which', however, 'which' is a function that look to all elements of the vector and, therefore, is not fast.
			collapseNode <- names(bt[ct][nodeOrder])[1] ## get the Node number of the first node to be collapsed ## 'bt[ct]' gives all the nodes that need to be collapsed. The index of 'nodeOrder' put the elements of the vector ordered without using the function again (save time). The names of the 'bt' vector are the number (labels) of the nodes.
		} else if(endClades != 0) {
			if(length(bt) <= (endClades-1)) break;
			nodeOrder <- order(bt) ##Notice that we need the same object here but we are in a different if statement, therefore we need to create it again. But now we have ALL the nodes of phy instead of just the ones we want to collapse. That is not a problem, for we have the 'ct' above (line 16).
			collapseNode <- names(bt[nodeOrder])[1]

		}
		ns <- tips(phy, as.numeric(collapseNode)) # The previous version was using the function 'node.leaves' however this will be deprecated of the geiger package.
		whichTip <- which(phy$tip.label == ns[1]) ## pull out the first tip from the clade to be collapsed. ## Gives the index of the first tip in the subset selected by 'node.leaves'.
		whichCollapse <- which(phy$tip.label %in% ns[-1]) ## gets all the others for removal ## The '%in%' means set membership. Therefore, which tips of phy are in the set formed by ns[-1], given that ns[-1] is the set of tips to be removed, minus the first element.
		rich[whichTip] <- rich[whichTip] + sum(rich[whichCollapse]) # computes the total richess for the dropped clade and puts it in the rich vector at the position for the retained (first) tip ##  Why not use 'length(ns)'?
		names(rich)[whichTip] <- paste(ns, collapse=".") ## renames richness ##The tip that will be retained will have the name of all tips that where collapsed. These is a good way to keep track of watch have been collapsed.

		rich <- rich[-whichCollapse] # removes the collpased tipes from richness ## Also do the job of setting the counter.
		
		phy$tip.label[whichTip] <- paste(ns, collapse=".") ## Here we are also changing the tip labels of the phy.
		for(j in 2:length(ns))
			phy <- drop.tip(phy, ns[j])
		
	}

	strsplit(phy$tip.label, "\\.") -> terminals; ## This really cool function will split the vector in different objects using the pattern given. Notice that 'strsplit' will create a list everytime that it split one element of the 'phy$tip.label' character string.
	
	
	if(is.null(trait)) {
		
		species <- numeric(length(terminals)) ## This is not the number of original tips. This is the number of tips after the last prunning cycle. Notice that 'terminals' is a list and not a vector. The function 'numeric' will create a vector of '0' with the given length. 
		
		for (i in 1:length(terminals) ) {
			species[i] <- length(terminals[[i]]) ## This set in 'species' the number of species of each tip, after the prunning cycle (i.e., the number of original tips).
		}
		
		phy$tip.label <- paste("meccaTree", seq(1:length(phy$tip.label)), sep = "_");	
		richness <- data.frame("clade" = phy$tip.label, "richness" = species) ##Creating a table with the number of species of each tip.

		
	} else {
	cladeMean <- numeric(length(terminals));
	cladeVar <- numeric(length(terminals));
	species <- numeric(length(terminals))
	
	
	for (k in 1:length(terminals) ) {
		cladeMean[k] <- mean(trait[terminals[[k]]]);
		cladeVar[k] <- var(trait[terminals[[k]]]);
		species[k] <- length(terminals[[k]])
		}
		
	phy$tip.label <- paste("meccaTree", seq(1:length(phy$tip.label)), sep = "_");	
	cladeVar[is.na(cladeVar)] <- 0;

	richness <- data.frame("clade" = phy$tip.label, "cladeMean" = cladeMean, "cladeVariance" = cladeVar, "richness" = species)
	}

	return(list(phy=phy, richness=richness))
}
