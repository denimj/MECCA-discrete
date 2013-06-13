acceptance <-
function(current, proposed, distribution, params) {
	
	if(distribution == "raw" | distribution == "log") {
		return("accept");
	}
	
	if(distribution == "lognorm") {
		k <- length(current);
		priorden <- numeric(k)
		for(i in 1:k) {
		priorden[i] <- dnorm(proposed[i], mean = params[1], sd = params[2]) / dnorm(current[i], mean = params[1], sd = params[2]);
		}
		priorden <- prod(priorden);
		p <- runif(1);
		
		if(priorden >= p) {
			
			return("accept")
			} else {
				
				return("reject")
				}
		
		}	
	
	}
