cladeAge <-
function(phy) {
# a function to pull tip lineage ages from a phylogenetic tree
# Arguments: a phylogenetic tree
# Values: a named vector of ages for terminal lineages		
			
	tips <- which(phy$edge[ ,2] <= length(phy$tip.label))
	cladeAges <- phy$edge.length[tips]
	tipNumbers <- phy$edge[tips, 2]
	names(cladeAges) <- phy$tip.label[tipNumbers]
	m <- match(phy$tip.label, names(cladeAges)) 
	cladeAges <- cladeAges[m]
	return(cladeAges)
}
