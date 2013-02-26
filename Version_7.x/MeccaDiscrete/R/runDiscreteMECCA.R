runDiscreteMECCA <-
function(phy, richness, startState, Ngens = 10000, printFreq = 100, qPrior=c(0.001, 1), scale = 1, divSampleFreq = 0,  summarizeSummaries = FALSE, outputName ="discreteMecca", seed=1) {  ## we need the priors in here!
  
## To run MECCA, the following arguments are required
	# phy - an incompletely resolved phylogenetic tree
	# richness - a named vector of species richness values for the tips in the phylogenetic tree
	# cal - the output of calibrateMECCA ## What is the calibrateMECCA? Try to search for it in the Slater paper!
	# Ngens - the number of generations that MECCA will be run for

	## house keeping block - check if everything is in the correct order, get rid of node labels and generate a vector of clade age values ##

        wn <- getOption("warn")
        options(warn=-1)
        nc <- geiger:::.treedata(phy, richness) ## Before this was using name.check() but this function is not available anymore.
        na <- sum(is.na(richness))
        if(!nc == "OK") { stop("names in tree and data do not match") }
        else {
           if(!na == 0) { stop("richness contains missing values") }
        }
        options(warn=wn)
  
	phy$node.label <- NULL

	m <- match(phy$tip.label, names(richness))
	richness <- richness[m]
	       
	cladeAges <- cladeAge(phy)
	
	dcrit <- startState$dcrit
	obs <- startState$obsTraits
	plsobs <- startState$plsObserved
	pls <- startState$plsLoadings
	ncomp <- ncol(pls)
	
	## end of name checking and data ordering ##
        
	## we generate a matrix to hold our sampled parameters and their associated summary statistics | will hold the pls transformed summary data
	
	distSim <- file(paste(outputName, "distSimFile.txt", sep = "_"), "w")
	
	cat("SigmaSq", "RootState", paste("pls", 1:length(plsobs), sep = ""), file = distSim, sep = "\t")
	cat("\n", file=distSim)
	
	## will hold the untransformed summary data
	charSim <- file(paste(outputName, "charSimFile.txt", sep = "_"), "w")
	cat("q12", "q21", names(obs), file = charSim, sep = "\t")
	cat("\n", file= charSim)
	
	bdSim <- file(paste(outputName, "bdSimFile.txt", sep = "_"), "w")
	writeLines(paste("Lambda", "Mu",  "lkl"), con = bdSim)

	## ABCtoolkit requires a matrix of observed data in the same order as the posterior sample. MECCA will provide this for you as yourMECCAoutputName$observed

	## pls transformed observed
	names(plsobs) <- c(paste("pls", 1:length(plsobs), sep = ""))
	
	write.table(t(plsobs), paste(outputName, "distObsFile.txt", sep = "_"), row.names = F, sep = "\t", quote = F)
			
	write.table(t(obs), paste(outputName, "ObsFile.txt", sep = "_"), row.names = F, sep = "\t", quote = F)
	
	# TRACK PROPOSALS ETC
	#proposalMat <- matrix(data = NA, nrow = Ngens, ncol = 3 + length(obs));
	#colnames(proposalMat) <- c("rate", "root", "distance", names(obs));
	## end of matrix creation ##

	nAcceptBD = 0
	nAcceptQ = 0
	## set up blank trace plot
	# plot(0, xlim=c(0, Ngens), ylim=c(0,0.1), cex=0.001,main="Trace Plot",xlab="generation", ylab="params")
	## we initialize MECCA with the mean of the within-tolerance parameters sampled during calibration
 
	thetaQ12 <- startState$startingQ12
	thetaQ21 <- startState$startingQ21
	thetaB <- startState$startingB
	thetaD <- startState$startingD

	Lk <- logLP(phy, thetaB, thetaD) + logLT(phy, richness, thetaB, thetaD)
		
	## and now we can run MECCA ##
	#cat("Generation", "\t", "Accepted", "\n\n");

	for(ngen in 1:Ngens) {
 
 	while(1) {
 		## use a sliding window proposal range with width derived from the calibration step ##
		thetaBprop<-getProposal(thetaB, startState$widthB * scale, min = 0, max = Inf)
		thetaDprop<-getProposal(thetaD, startState$widthD * scale, min = 0, max = Inf)
			
		if(thetaBprop > thetaDprop) ## we will only proceed with the proposed values if they are greater than 0 and b > d
		break; 
	}
		
	## update q's 
	thetaQ12prop <- getProposal(thetaQ12, startState$widthQ12 * scale, min = qPrior[1], max = qPrior[2])
	thetaQ21prop <- getProposal(thetaQ21, startState$widthQ21 * scale, min = qPrior[1], max = qPrior[2])
		
	propLk <- logLP(phy, thetaBprop, thetaDprop) + logLT(phy, richness, thetaBprop, thetaDprop) # compute the joint taxonomic/phylogenetic likelihood of the tree/richness data based on propB and propD

	if(propLk == - Inf) {
		LKratio <- 0
	} else {
		LKratio <- exp(propLk - Lk)
	}
		
	#### At this point we need to either accept or reject the proposed B and D values - we accept if LKratio >/= 1 and accept with probability p if not. If we accept, we simulate trees under these values and simulate trait evolution on them for abc - if we don't accept, we simulate data and trees using the current states and store the likelihoods and summaries and distances. We might still reject the proposals for sigmasq and root state if traits aren't close though 
	if(LKratio >= 1 || runif(1) < LKratio) {			
		thetaB <- thetaBprop; thetaD <- thetaDprop; Lk <- propLk
		nAcceptBD <- nAcceptBD + 1
	} 
	
	writeLines(paste(thetaB, thetaD, Lk), con = bdSim)
	
	##### now that we have updated our diversification parameters or rejected the proposed move we can move on to simulating trees and traits ####	
	
	p1root <- thetaQ21prop/(thetaQ12prop + thetaQ21prop)
	
	if(runif(1)<p1root) {rootState <- 1} else {rootState <- 2}
	simData <- simDiscreteIncompleteTree(phy, richness, q12= thetaQ12prop, q21= thetaQ21prop, rootState=rootState, r= thetaB-thetaD, eps=thetaD/thetaB, seed) ## DCS: This line was using 'mtree' before, this object does not exist.
		 
	sims <- simData$res[,1]/rowSums(simData$res)

 	 plsSim <- extractpls(sims, pls, ncomp)
 	 
 	 dSim <- abs(dist(rbind(plsobs, plsSim))[1])
 	 
 	 #proposalMat[ngen, ] <- c(thetaSprop, thetaMprop, dSim, sims);
   	
   		if (dSim <= dcrit) {## if it's odd assess q's
   		
				thetaQ12 <- thetaQ12prop ## accept the proposals
				thetaQ21 <- thetaQ21prop ## accept the proposals

				nAcceptQ = nAcceptQ + 1
				cat(thetaQ12, thetaQ21, plsSim, file = distSim, sep = "\t")
				cat("\n", file = distSim)
				cat(thetaQ12, thetaQ21, sims, file = charSim, sep = "\t")
				cat("\n", file= charSim)
		} else { # if dSim > dcrit
				
				p1root<-thetaQ21/(thetaQ12 + thetaQ21)
				
				if(runif(1)<p1root) {rootState<-1} else {rootState<-2}
	
				Sim <- simDiscreteIncompleteTree(phy, richness, q12= thetaQ12, q21= thetaQ21, rootState=rootState, r= thetaB-thetaD, eps=thetaD/thetaB) ##DCS: Line was using 'mtree'. This object does not exist.
				
				sims <- Sim$res[,1]/rowSums(Sim$res)
				
 	 			plsSim <- extractpls(sims, pls, ncomp)
 	 			
				cat(thetaQ12, thetaQ21, plsSim, file = distSim, sep = "\t")
				cat("\n", file = distSim)
				cat(thetaQ12, thetaQ21, sims, file = charSim, sep = "\t")
				cat("\n", file= charSim)
			
		} 

	if(ngen %% printFreq == 0) {
			cat("Generation", ngen, "/", round(nAcceptBD/ngen, 2), round(nAcceptQ/(ngen/2), 2), "\n\n");
			}
	} ## END MCMC
			
	close(charSim); close(bdSim); close(distSim)
}
