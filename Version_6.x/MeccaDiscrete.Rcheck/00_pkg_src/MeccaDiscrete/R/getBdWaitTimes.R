getBdWaitTimes <-
function (n, age, lambda, mu, rho=1) 
{
    #edge <- c(-1, -2)
    #leaves <- c(-2)
    #timecreation <- c(0, 0)
    #time <- 0
    #maxspecies <- -2
    #edge.length <- c(0)
    #index = 2
  ## Rho is the sampling probability. The probability that a given tip (individual) are sampled. Then, values below 1 equal to incomplete samplings.

  ## Waiting times are the intervals between the root of a tree to each of the specitaion events in time. They are measured as time slices of the complete phylogeny. 
  
  ## The two lines below seem not to make sense. Then I did:
    #specevents <- vector() ##Since the vector is empty, the value of 'age' will be equall to 'specevents' directly.
    #specevents <- c(specevents, age)
  specevents <- age

  # The calculations below were based on the work from Stadler T. You would expect the need to simulate all the trees and then sample. However, this matematical procedure avoids the need to do this and make the process faster. 
    for (j in 1:(n - 1)) {
        r <- runif(1, 0, 1)
        tor <- age
        lamb1 <- rho * lambda
        mu1 <- mu - lambda * (1 - rho)
        if (lambda > mu) {
            temp <- 1/(lamb1 - mu1) * log((lamb1 - mu1 * exp((-lamb1 + 
                mu1) * tor) - mu1 * (1 - exp((-lamb1 + mu1) * 
                tor)) * r)/(lamb1 - mu1 * exp((-lamb1 + mu1) * 
                tor) - lamb1 * (1 - exp((-lamb1 + mu1) * tor)) * 
                r))
        }
        else {
            temp <- -((tor * r)/(-1 - lambda * rho * tor + lambda * 
                rho * tor * r))
        }
        specevents <- c(specevents, temp) ## This is appending the values of temp in the specevents vector. The first value of the vector is 'age' parameter.
    }
    specevents <- sort(specevents, decreasing = TRUE)
	specevents    
}
