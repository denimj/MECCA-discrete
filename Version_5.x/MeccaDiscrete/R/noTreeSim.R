noTreeSim <-
function(times, qMatrix, startState, seed) {
	currState <- numeric(length=dim(qMatrix)[1])
	currState[startState] <- 1
	
	for(i in 2:length(times)) {
		# pick waiting time
		wt <- times[i-1]-times[i]
				
		# simulate character change on all branches up to the wait time
		pp <- MatrixExp(qMatrix * wt)
		
		#how many change?
		ch12 <- rbinom(n=1, size=currState[1], prob=pp[1,2])
		ch21 <- rbinom(n=1, size=currState[2], prob=pp[2,1])

		# add these to the totals
		currState[1] <- currState[1]-ch12+ch21
		currState[2] <- currState[2]-ch21+ch12
		
		# which type is affected by the event?
		prT1 <- currState[1]/sum(currState)
		if(runif(1)<prT1) subj <- 1 else subj <- 2
				
		# speciation
		currState[subj] <- currState[subj]+1

		}
	
	# last interval
	wt <- times[i]
				
	# simulate character change on all branches up to the wait time
	pp <- MatrixExp(qMatrix * wt)
		
	#how many change?
	ch12 <- rbinom(n=1, size=currState[1], prob=pp[1,2])
	ch21 <- rbinom(n=1, size=currState[2], prob=pp[2,1])

	# add these to the totals
	currState[1] <- currState[1]-ch12+ch21
	currState[2] <- currState[2]-ch21+ch12
	# All the results in currState have length() = 1.
	currState[3] <- seed
	
	return(currState)
	}
