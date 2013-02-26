inspectCalOut <-
function(bdcal, bmcal, nhot = 0) {
	
	k <- 4 + nhot
	#quartz() ##DCS: Do not need this line, par will do the job.
	par(mfrow=c(ceiling(k/2),2))
	L <- nrow(bdcal)
	for(i in 1:2) { hist(bdcal[,i], main =NULL,  breaks = 100, xlab = colnames(bdcal)[i]); }
	for(i in 1:(2+nhot)) { hist(bmcal[,i], main =NULL,  breaks = 100, xlab = colnames(bmcal)[i]); }
	X11() #DCS: This will work in Linux now... Good...
	plot(seq(1:L), bdcal[1:L, 3], type ="l", xlab = "generation", ylab = "logLk");
	
	}
