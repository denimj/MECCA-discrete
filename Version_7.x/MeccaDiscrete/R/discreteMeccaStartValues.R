discreteMeccaStartValues <-
function(calibrationOutput, phy, tipStates, tolerance = 0.01, plsComponents) {

	
	m <- match(phy$tip.label, rownames(tipStates));
	tipStates <- tipStates[m,];
	
	obs<-tipStates[,1]/rowSums(tipStates)
		
	stdobs <- obs; # observed stats for standardization and pls
	
	bdcal <- as.matrix(calibrationOutput$calResBD);
	bmcal <- as.matrix(calibrationOutput$calResChar);	
	
	Ncomp <- ncol(plsComponents); # gives the number of PLS components
	Nsims <- nrow(bdcal)

	bdretained <- bdcal[sort(bdcal[,3]), ][seq(1, Nsims * tolerance), ]; # retained diversification parameters 
		
	obsPLS <- numeric(Ncomp);
	
	#print(length(stdobs));
	#print(length(plsComponents[, 1]));
	
	for(i in 1:length(obsPLS)) {		
		obsPLS[i] <- sum(stdobs * plsComponents[, i]);
	} # creates the pls asjusted oberved data
	print("done");
	
	plsSim <- matrix(NA, nrow = Nsims, ncol = Ncomp);
	params <- bmcal[ , 1:2];
	stats <- bmcal[ , - (1:2)];
	
	for(i in 1:Nsims) {
	
	plsSim[i, ] <- extractpls(stats[i,], plsComponents, Ncomp); ## this needs adjusting for different numbers of pls components
	
	} ## now this holds the pls adjusted distances 
	
	dmat <- numeric(Nsims);

	for (i in 1:Nsims) {
	
		dmat[i] <- dist(rbind(obsPLS, plsSim[i, ]));
		} # computes the euclidean distances for the simulations
	
	pdmat <- cbind(params, dmat);
	bmretained <- pdmat[order(pdmat[,3]),][seq(1, Nsims * tolerance) , ];
	
	## and now we compute the starting values (mean) and proposal width (sd) for the MCMC ##
	
	
	mean(bmretained[,1]) -> startingQ12; ## takes the mean of the prior values of sigma that resulted in variances within our distance 
	sd(bmretained[,1]) -> widthQ12;	## takes the sd of the prior values of sigma that resulted in variances within our distance
	mean(bmretained[,2]) -> startingQ21; ## takes the mean of the prior values of sigma that resulted in variances within our distance 
	sd(bmretained[,2]) -> widthQ21;	## takes the sd of the prior values of sigma that resulted in variances within our distance
	
	mean(bdretained[,1]) -> startingB;	## takes the mean of the prior values of lambda that resulted in variances within our distance
	sd(bdretained[,1]) -> widthB;	## takes the sd of the prior values of lambda that resulted in variances within our distance
	mean(bdretained[,2]) -> startingD; ## mean of prior death rates giving good taxa
	sd(bdretained[,2]) -> widthD;  ## sd of prior death rates giving good taxa



		return(list("startingQ12" = startingQ12, "widthQ12"=widthQ12, "startingQ21"= startingQ21, "widthQ21"=widthQ21,"startingB" = startingB, "widthB" = widthB, "startingD" = startingD, "widthD" = widthD, "dcrit" = max(bmretained[,3]), "obsTraits" = obs, "plsObserved" = obsPLS, "plsLoadings" = plsComponents));
		
	}
