getDivMLE <-
function(phy, richness) {
	require(laser)
# a function to compute the maximum likelihood estimators of birth and death from a richness tree. uses the function fitNDR_1rate from Laser (Rabosky 2006)
	print("getting starting diversification parameters");

# Arguments
	# phy <- a phylogenetic tree
	# richness <- a vector of species richness values for the tips in phy
	
# Values

	# lik <- the log Likelihood of the MLE
	# div <- the  MLE of diversification rate
	# lambda <- the  MLE of birth rate  
	# mu <- the  MLE of extinction rate
            
	phy2 <- getTipdata(richness, phy);
	ext <- seq(0, 0.99, 0.005); # gives extinction fractions from 0 to 99%
	
	lik<-numeric(length(ext));
	div<-numeric(length(ext));
	lam<-numeric(length(ext));
	
	for(i in 1:length(ext)) {
	
		fitNDR_1rate(phy2, eps = ext[i], rbounds = c(0.0001, .5), combined = TRUE) -> fit;
		lik[i] <- fit$LH;
		div[i] <- fit$r;
		lam[i] <- fit$lambda;
	
		}
names(lik) <- ext;
names(div) <- ext;
names(lam) <- ext;


return(list("lik" = max(lik), "div" = div[which(lik==max(lik))],"lamda" = lam[which(lik==max(lik))], "mu" = ext[which(lik==max(lik))] * lam[which(lik==max(lik))]
));	
}
