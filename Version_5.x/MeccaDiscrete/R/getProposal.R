getProposal <-
function(currentState, psi, min, max) {

# This function returns a proposal using a sliding window	
# Arguments;
	
	# current state <- the current state of the parameter;
	# psi <- a tuning parameter - here, some proportion of the standard deviation from the calibration step;
	# min <- the lower bound for the prior distribution;
	# max <- the upper bound for the prior distribution;
	
# Value;	
	# prop <- the proposal; 
	
	prop <- currentState + ((runif(1) - 0.5) * psi);
	
	if (prop < min) { prop <- min + (min - prop) } # if lower than lower bound, bounce back;
	if (prop > max) { prop <- max - (prop - max) } # if higher than upper bound, bounce back;

	return(prop);
	
	}
