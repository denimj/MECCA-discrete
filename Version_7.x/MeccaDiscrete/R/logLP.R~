logLP <-
function(phy, b, d) {
# function that computes the phylogenetic likelihood for a backbone tree following Rabosky et al 2007

#Arguments
	# phy <- a phylogenetic tree
	# b <- the proposed birth rate
	# d <- the proposed death rate 	
	
	r <- b - d # computes the diversification rate
	a <- d / b # computes the extinction fraction (also called relative extinction by Rabosky)
	N <- phy$Nnode - 1 # computes the number of internal branches
		
	Xi <- cbind(phy$edge, phy$edge.length)
	Xi <- Xi[-which(Xi[ ,2] <= length(phy$tip.label)), ] ## Only the internal branches. This takes off all the terminal branches.
	Xi <- cbind(Xi, numeric(length(Xi[ ,3]))) ## Add one empty column to Xi.
	bt <- branching.times(phy) ## A vector with the distance from each node to the tip. names(bt) equall to the node numbers in the first column of phy$edge
	
	for (i in 1: length(Xi[ ,4])) {
		Xi[i, 4] <- bt[which(names(bt)==Xi[i, 1])]
		}
		
	part3 <- numeric(length(Xi[ ,4]))
	
	for (p in 1:length(part3)) {
		part3[p] <- log(1 - a * exp(-r*Xi[p, 4]))
		}
		
	logLp <- (N * log(r) - r * sum(Xi[ ,3]) - sum(part3))
	
	return(logLp)
	}
