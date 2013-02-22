logLP <-
function(phy, b, d) {
# function that computes the phylogenetic likelihood for a backbone tree following Rabosky et al 2007

#Arguments
	# phy <- a phylogenetic tree
	# b <- the proposed birth rate
	# d <- the proposed death rate 	
	
	r <- b - d; # computes the diversification rate
	a <- d / b; # computes the extinction fraction
	N <- phy$Nnode - 1; # computes the number of internal branches
		
	Xi <- cbind(phy$edge, phy$edge.length);
	Xi <- Xi[-which(Xi[ ,2] <= length(phy$tip.label)), ];
	Xi <- cbind(Xi, numeric(length(Xi[ ,3])));
	bt <- branching.times(phy);
	for (i in 1: length(Xi[ ,4])) {
		Xi[i, 4] <- bt[which(names(bt)==Xi[i, 1])]; 
		}
		

	
	part3 <- numeric(length(Xi[ ,4]));
	
	for (p in 1:length(part3)) {
		
		part3[p]<- log(1 - a * exp(-r*Xi[p, 4]));
		
		}
		
	logLp <- (N * log(r) - r * sum(Xi[ ,3]) - sum(part3));
	
	return(logLp);
	
	}
