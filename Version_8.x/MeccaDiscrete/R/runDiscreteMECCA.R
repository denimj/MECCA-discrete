runDiscreteMECCA <-
function(phy, richness, startState, Ngens = 10000, printFreq = 100, qPrior=c(0.001, 1), scale = 1, scale.dcrit = 1, divSampleFreq = 0,  summarizeSummaries = FALSE, outputName ="discreteMecca", seed=1) {  ## we need the priors in here!

## To run MECCA, the following arguments are required
	# phy - an incompletely resolved phylogenetic tree
	# richness - a named vector of species richness values for the tips in the phylogenetic tree
	# cal - the output of calibrateMECCA
	# Ngens - the number of generations that MECCA will be run for
        # scale.dcrit - one scale to modify the critical distance between the data and the proposal in order to accept or reject new values of q12 and q21. Values <1 will make less proposals to be accepted and values >1 will make more proposals to be accepted.

	## house keeping block - check if everything is in the correct order, get rid of node labels and generate a vector of clade age values ##

        wn <- getOption("warn")
        options(warn=-1)
        nc <- geiger:::.treedata(phy, richness)
        na <- sum(is.na(richness))
        if(!nc == "OK") {
        	stop("names in tree and data do not match")
        	} else {
           	if(!na == 0) { stop("richness contains missing values") }
        }
        options(warn=wn)
  
	phy$node.label <- NULL
	m <- match(phy$tip.label, names(richness))
	richness <- richness[m]
	       
	cladeAges <- cladeAge(phy) ## Vector of the branch lengths of the tips.
	
	dcrit <- startState$dcrit
        dcrit <- dcrit * scale.dcrit
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
	writeLines(paste("Lambda", "Mu",  "lkl", sep = "\t"), con = bdSim)

	## ABCtoolkit requires a matrix of observed data in the same order as the posterior sample. MECCA will provide this for you as yourMECCAoutputName$observed --> DCS: yourMECCAoutputName$charSimFile.txt

	## pls transformed observed. Tables from things generated previously in StartMecca.
	names(plsobs) <- c(paste("pls", 1:length(plsobs), sep = ""))
	write.table(t(plsobs), paste(outputName, "distObsFile.txt", sep = "_"), row.names = F, sep = "\t", quote = F)
	write.table(t(obs), paste(outputName, "ObsFile.txt", sep = "_"), row.names = F, sep = "\t", quote = F)
	
	# Keep track of the proposal fit, acceptance rates and other info:
        proposalMatBD <- file(paste(outputName, "track_prop_BD.txt", sep = "_"), "w")
	cat("Birth_root", "Birth_prop","Death_root","Death_prop","LKratio","accepted", file = proposalMatBD, sep = "\t")
        cat("\n", file = proposalMatBD)
        
        proposalMatQij <- file(paste(outputName, "track_prop_Qij.txt", sep = "_"), "w")
	cat("Q12_root", "Q12_prop","Q21_root","Q21_prop","PLS_1","accepted", file = proposalMatQij, sep = "\t")
        cat("\n", file = proposalMatQij)

        ############### end of matrix creation ######################

	nAcceptBD <- numeric(0)
	nAcceptQ <- numeric(0)
        acceptedBD <- logical(0)
        acceptedQij <- logical(0)
	## we initialize MECCA with the mean of the within-tolerance parameters sampled during calibration
 
	thetaQ12 <- startState$startingQ12
	thetaQ21 <- startState$startingQ21
	thetaB <- startState$startingB
	thetaD <- startState$startingD

	Lk <- logLP(phy, thetaB, thetaD) + logLT(phy, richness, thetaB, thetaD) ## This will be the prior likelihood of the tree? Using the prior values from startmecca?
		
	## and now we can run MECCA ##

	for(ngen in 1:Ngens) {

          ## Reset the values for the acceptance flags.
          acceptedBD <- FALSE
          acceptedQij <- FALSE
 
 		while(1) {
## use a sliding window proposal range with width derived from the calibration step ##
## Note that the sliding window is wider than the sd from the prior distribution.
			thetaBprop <- getProposal(thetaB, startState$widthB * scale, min = 0, max = Inf)
			thetaDprop <- getProposal(thetaD, startState$widthD * scale, min = 0, max = Inf)
			if(thetaBprop > thetaDprop)
			break
		}
		
	## update q's -> Apply sliding window.
		thetaQ12prop <- getProposal(thetaQ12, startState$widthQ12 * scale, min = qPrior[1], max = qPrior[2])
		thetaQ21prop <- getProposal(thetaQ21, startState$widthQ21 * scale, min = qPrior[1], max = qPrior[2])
		
		propLk <- logLP(phy, thetaBprop, thetaDprop) + logLT(phy, richness, thetaBprop, thetaDprop) # compute the joint taxonomic/phylogenetic likelihood of the tree/richness data based on propB and propD

		if(propLk == - Inf) {
			LKratio <- 0
			} else {
			LKratio <- exp(propLk - Lk)
		}
		
## At this point we need to either accept or reject the proposed B and D values - we accept if LKratio >/= 1 and accept with probability p if not. [DSC This probability p is important for the acceptance algorithm not to be too greedy.] If we accept, we simulate trees under these values and simulate trait evolution on them for abc - if we don't accept, we simulate data and trees using the current states and store the likelihoods and summaries and distances [DSC Remember that the frequency is important in the posterior distribution, therefore the value of LKratio that was not beated by the prop need to be used to simulate - we need to simulate again because of the stocasticy of the process]. We might still reject the proposals for sigmasq (sigma square) and root state if traits aren't close though. [DSC The traits will be compared using the matrix distances and the dcrit value.]
 
		if(LKratio >= 1 || runif(1) < LKratio) {			
			thetaB <- thetaBprop
			thetaD <- thetaDprop
			Lk <- propLk
			nAcceptBD <- nAcceptBD + 1
                        acceptedBD <- TRUE
		}
	
		writeLines(paste(thetaB, thetaD, Lk, sep = "\t"), con = bdSim)
	
## now that we have updated our diversification parameters or rejected the proposed move we can move on to simulating trees and traits.	
## DSC -> We need to simulate tree and traits and to reject or accept the q12 and q21 proposals before moving to the next generation.

		p1root <- thetaQ21prop/(thetaQ12prop + thetaQ21prop)
	
		if(runif(1) < p1root) {rootState <- 1} else {rootState <- 2}
		
		simData <- simDiscreteIncompleteTree(phy, richness, q12= thetaQ12prop, q21= thetaQ21prop, rootState=rootState, r= thetaB-thetaD, eps=thetaD/thetaB, seed)
		sims <- simData$res[,1]/rowSums(simData$res)
 	 	plsSim <- extractpls(sims, pls, ncomp)
 	 	dSim <- abs(dist(rbind(plsobs, plsSim))[1])

   		if (dSim <= dcrit) {  ## if it's odd assess q's
   			thetaQ12 <- thetaQ12prop  ## accept the proposals
			thetaQ21 <- thetaQ21prop  ## accept the proposals
			nAcceptQ <- nAcceptQ + 1
                        acceptedQij <- TRUE
			cat(thetaQ12, thetaQ21, plsSim, file = distSim, sep = "\t")
			cat("\n", file = distSim)
			cat(thetaQ12, thetaQ21, sims, file = charSim, sep = "\t")
			cat("\n", file= charSim)
			} else {  # if dSim > dcrit
				p1root <- thetaQ21/(thetaQ12 + thetaQ21)
				## reject the proposal and simulate using the other values.
				if(runif(1) < p1root) {rootState <- 1} else {rootState <- 2}
				Sim <- simDiscreteIncompleteTree(phy, richness, q12= thetaQ12, q21= thetaQ21, rootState=rootState, r= thetaB-thetaD, eps=thetaD/thetaB)
				sims <- Sim$res[,1]/rowSums(Sim$res)
 	 			plsSim <- extractpls(sims, pls, ncomp)
				cat(thetaQ12, thetaQ21, plsSim, file = distSim, sep = "\t")
				cat("\n", file = distSim)
				cat(thetaQ12, thetaQ21, sims, file = charSim, sep = "\t")
				cat("\n", file= charSim)
			}
          
          ## Writing the values of the proposals, fit and etc.
          ## Everytime that the proposal is accepted, the values of thetaX and thetaXprop will be equal.
          ## Would be simple to check for this. However, it is better to flag when acceptance occurr.
          cat(thetaB,thetaBprop,thetaD,thetaDprop,LKratio,acceptedBD, file = proposalMatBD, sep = "\t")
          cat("\n", file = proposalMatBD)
          cat(thetaQ12,thetaQ12prop,thetaQ21,thetaQ21prop,dSim,acceptedQij, file = proposalMatQij, sep = "\t")
          cat("\n", file = proposalMatQij)
                
		if(ngen %% printFreq == 0) {
			cat("Generation", ngen, "\n\n")
		}
	} ## END MCMC
			
	close(charSim); close(bdSim); close(distSim); close(proposalMatBD); close(proposalMatQij)
}
