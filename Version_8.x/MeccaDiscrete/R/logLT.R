logLT <-
function(phy, richness, b, d) {
# function that computes the taxonomic likelihood for terminals with species richness in a backbone tree following Rabosky et al 2007

#Arguments
	# phy <- a phylogenetic tree
	# richness <- a vector of species richness values
	# b <- the proposed birth rate
	# d <- the proposed death rate 	

	
	r <- b-d;
	a <- d/b;
	t <- cladeAge(phy);
	n <- richness[match(phy$tip.label, names(richness))];
	beta <- (exp(r * t) - 1) / (exp(r * t) - a);
		
	logLT <- sum(log(1-beta)) + sum((n-1)*log(beta));
	return(logLT);
	
	}
