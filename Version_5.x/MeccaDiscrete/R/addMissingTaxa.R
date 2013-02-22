addMissingTaxa <-
function(phy, tipDiversity, tipTraits, nsims=1000, method=c("resolvedFixed", "resolvedVariableOrdered",  "unresolvedFixed"), bRate=NA, dRate=0) 
{
	method<-match.arg(method)
	phyFast<-NA
	if(method=="resolvedFixed") {
		
		###########################################################
		# This is the "resolved" method of Davies and Harmon      #
		# for when traits are fixed within clades                 #
		###########################################################
		
		if(is.na(bRate)) stop("Must supply at least a birth rate (bRate) with resolvedFixed method.")
		
		# find number of original tips
		origtips<-length(phy$tip.label)
		
		# identify which tree edges are tip clades
		tt<-phy$edge[,2]<=origtips
		
		# save old tips and their names
		oldtip<-phy$edge.length[tt]
		names(oldtip)<-phy$tip.label[phy$edge[tt,2]]
	
		# cut off all tip branches leaving behind zero-length stems
		phy<-pruneZeroStem(phy, tipDiversity)
		
		# sort tipDiversity so that smaller clades get added first
		# this speeds up the code a tiny bit because the tree grows slowly at first
		# matters the most when there are a couple really large unresolved clades in the tree
		tipDiversity<-sort(tipDiversity)

		tipPhenotypes<-numeric()

		# loop through each tip one-by-one
		for(i in 1:length(tipDiversity)) {
			
		  # if diversity was one, the branch was not pruned off by pruneZeroStem so nothing needs to be done except add to the phenotypes
		  
		  if(tipDiversity[i]==1) {
		  	nt<-which(phy$tip.label==names(tipDiversity)[i]) 
			nn<-which(phy$edge[,2]==nt)
			b<-tipDiversity[i]
			t<-oldtip[names(tipDiversity)[i]]
			# assign the phenotypes
			# get the right row
			rr<-which(rownames(tipTraits)==names(tipDiversity)[i])
			
			# find the most common trait
			# if there is a tie pick the first one
			trt<-which(tipTraits[rr,]==max(tipTraits[rr,]))[1]
			
			names(trt)<-names(tipDiversity)[i]
			tipPhenotypes<-c(tipPhenotypes, trt)

			
		  	}	
		  if(tipDiversity[i]>1) {
		  	
		  	# locate the node in the tree structure (which tip nt, which edge nn), its diversity and age
			nt<-which(phy$tip.label==names(tipDiversity)[i]) 
			nn<-which(phy$edge[,2]==nt)
			b<-tipDiversity[i]
			t<-oldtip[names(tipDiversity)[i]]
			
			# simulate new tree for this missing clade
			newtree<-sim.bd.taxa.age(n=b, numbsim=1, lambda=bRate, mu=dRate, frac=1, age=t, mrca=F)[[1]]
			
			# calculate and add root edge
			rootedge<-t-max(branching.times(newtree))
 			phy$edge.length[nn]<-rootedge			
			# add unique tip labels
			newtree$tip.label=paste(names(tipDiversity)[i],".", b:1,  sep="")
			
			# paste on to the tree at the proper place
			phy<-bind.tree(phy, newtree, where=nt)
			
			# assign the phenotypes
			# get the right row
			rr<-which(rownames(tipTraits)==names(tipDiversity)[i])
			
			# find the most common trait
			# if there is a tie pick the first one
			trt<-which(tipTraits[rr,]==max(tipTraits[rr,]))[1]
			
			# make new phenotypes and stick on to results
			newPhen<-rep(trt, tipDiversity[i])
			names(newPhen)<-newtree$tip.label
			tipPhenotypes<-c(tipPhenotypes, newPhen)
			
			}
		}
	} else if(method=="resolvedVariableOrdered") {
		
		#########################################################################
		# This is the "resolved" method of Davies and Harmon      				#
		# for when traits are assigned to taxa in clades 		  				#
		# by ladderizing the nodes and assigning traits	in that order		  	#
		#########################################################################
		
		if(is.na(bRate)) stop("Must supply at least a birth rate (bRate) with resolvedVariableOrdered method.")
		
		# find number of original tips
		origtips<-length(phy$tip.label)
		
		# identify which tree edges are tip clades
		tt<-phy$edge[,2]<=origtips
		
		# save old tips and their names
		oldtip<-phy$edge.length[tt]
		names(oldtip)<-phy$tip.label[phy$edge[tt,2]]
	
		# cut off all tip branches leaving behind zero-length stems
		phy<-pruneZeroStem(phy, tipDiversity)
		
		# sort tipDiversity so that smaller clades get added first
		# this speeds up the code a tiny bit because the tree grows slowly at first
		# matters the most when there are a couple really large unresolved clades in the tree
		tipDiversity<-sort(tipDiversity)

		# need to add checker to see if tipDiversity and tipTraits are compatable
		# rowSums(tipTraits) == tipDiversity
		# accounting for mismatch in names

		tipPhenotypes<-numeric()

		nStates<-ncol(tipTraits)


		# loop through each tip one-by-one
		for(i in 1:length(tipDiversity)) {
			
		  # if diversity was one, the branch was not pruned off by pruneZeroStem so nothing needs to be done	
		  if(tipDiversity[i]==1) {
		  	
		  	# assign the phenotype
			# get the right row
			rr<-which(rownames(tipTraits)==names(tipDiversity)[i])
			
			# find the most common trait
			trt<-which(tipTraits[rr,]==max(tipTraits[rr,]))[1]
			
			# make new phenotypes and stick on to results
			newPhen<-trt
			names(newPhen)<-names(tipDiversity)[i]
			tipPhenotypes<-c(tipPhenotypes, newPhen)
			
		  	} else if(tipDiversity[i]>1) {
		  	
		  	# locate the node in the tree structure (which tip nt, which edge nn), its diversity and age
			nt<-which(phy$tip.label==names(tipDiversity)[i]) 
			nn<-which(phy$edge[,2]==nt)
			b<-tipDiversity[i]
			t<-oldtip[names(tipDiversity)[i]]
			
			# simulate new tree for this missing clade
			newtree<-sim.bd.taxa.age(n=b, numbsim=1, lambda=bRate, mu=dRate, frac=1, age=t, mrca=F)[[1]]
			
			# calculate and add root edge
			rootedge<-t-max(branching.times(newtree))
 			phy$edge.length[nn]<-rootedge			
			# add unique tip labels
			newtree$tip.label=paste(names(tipDiversity)[i],".", b:1,  sep="")
			
			ladTree<-read.tree(text=write.tree(newtree))
			ladTree <-ladderize(ladTree)
			
			# paste on to the tree at the proper place
			phy<-bind.tree(phy, ladTree, where=nt)
			
			# assign the phenotypes
			# get the right row
			rr<-which(rownames(tipTraits)==names(tipDiversity)[i])
			
			# pull out the state freq and order them
			traitFreq<-sort(tipTraits[rr,], decreasing=T)
			
			# create a vector of tipTraits
			newPhen<-rep(1:nStates, times= traitFreq)
			
			# get tip names in "ladderized" order
			ntips<-length(ladTree$tip.label)
			tipNums<-ladTree$edge[,2][ladTree$edge[,2]<=ntips]
			names(newPhen)<-ladTree$tip.label[tipNums]
			
			#randomly assign to tips
			newPhen<-sample(newPhen)
			names(newPhen)<-newtree$tip.label
			
			# paste together
			tipPhenotypes<-c(tipPhenotypes, newPhen)
			
			}
		}
	} else if(method=="unresolvedFixed") {
		
		###########################################################
		# This is the "unresolved" method of Davies and Harmon    #
		# for when traits are fixed within clades                 #
		###########################################################
		
		
		# find number of original tips
		origtips<-length(phy$tip.label)
		
		# identify which tree edges are tip clades
		tt<-phy$edge[,2]<=origtips
		
		# save old tips and their names
		oldtip<-phy$edge.length[tt]
		names(oldtip)<-phy$tip.label[phy$edge[tt,2]]
	
		# cut off all tip branches leaving behind zero-length stems
		phy<-pruneYuleStem(phy, tipDiversity)
		
		# phyFast will be a special tree that can be used with the rapid fitDiscreteShortcut
		phyFast<-phy
		mm<-match(phyFast$tip.label, names(tipDiversity))
		phyFast$tipDiversity<-tipDiversity[mm]
		phyFast$tipLength<-numeric(length=length(phyFast$tip.label))
		names(phyFast$tipLength)<-phyFast$tip.label
		
		phyFast$trait1<-numeric(length=length(phyFast$tip.label))
		
		# sort tipDiversity so that smaller clades get added first
		# this speeds up the code a tiny bit because the tree grows slowly at first
		# matters the most when there are a couple really large unresolved clades in the tree
		tipDiversity<-sort(tipDiversity)

		tipPhenotypes<-numeric()

		# loop through each tip one-by-one
		for(i in 1:length(tipDiversity)) {
			
		  # if diversity was one, the branch was not pruned off by pruneZeroStem 
		# still need to set tipLength of phyFast to real branch length
		# locate the node in the tree structure (which tip nt, which edge nn), its diversity and age

			nt<-which(phy$tip.label==names(tipDiversity)[i]) 
			nn<-which(phy$edge[,2]==nt)
			b<-tipDiversity[i]
			t<-oldtip[names(tipDiversity)[i]]
		
		
			nnt<-which(phyFast$tip.label==names(tipDiversity)[i]) 

		  if(tipDiversity[i]==1) {	

			phyFast$tipLength[nnt]<-0 # to prevent adding extra length when calculating likelihoods
			rr<-which(rownames(tipTraits)==names(tipDiversity)[i])
			trt<-which(tipTraits[rr,]==max(tipTraits[rr,]))[1]
			if(trt==1) phyFast$trait1[nnt]<-1
			
			names(trt)<-names(tipDiversity)[i]
			tipPhenotypes<-c(tipPhenotypes, trt)

			
			} else if(tipDiversity[i]>1) {
		  
			stemLength<-phy$edge.length[nn]
			nn<-tipDiversity[i]
			totHist<-nn*t/sum(1/(1:nn))
			newEdge<-(totHist-stemLength)/nn
			
			nnt<-which(phyFast$tip.label==names(tipDiversity)[i]) 
			phyFast$tipLength[nnt]<-newEdge

			
			newtree<-list()
			newtree$edge<-matrix(nrow=nn, ncol=2)
			newtree$edge[,1]<-nn+1
			newtree$edge[,2]<-1:nn
			newtree$edge.length<-numeric(length=nn)
			newtree$edge.length[]<-newEdge
			newtree$Nnode=1
			newtree$tip.label=paste(names(tipDiversity)[i],".", nn:1,  sep="")
			class(newtree)<-"phylo"
		
			phy<-bind.tree(phy, newtree, where=nt)
			
						
			# assign the phenotypes
			# get the right row
			rr<-which(rownames(tipTraits)==names(tipDiversity)[i])
			
			# find the most common trait
			# if there is a tie pick the first one
			trt<-which(tipTraits[rr,]==max(tipTraits[rr,]))[1]
			
			# make new phenotypes and stick on to results
			newPhen<-rep(trt, tipDiversity[i])
			names(newPhen)<-newtree$tip.label
			tipPhenotypes<-c(tipPhenotypes, newPhen)
			
			if(trt==1) phyFast$trait1[nnt]<-tipDiversity[i]
			
			}
		}
	
		}
	list(phy=phy, tipPhenotypes= tipPhenotypes, phyFast=phyFast)
	
	}
