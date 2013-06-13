getRateCats <-
function(data, model)
{
	if(model=="ER") return(1)
	n<-nlevels(factor(data))
	if(model=="SYM") return(n*(n-1)/2)
	if(model=="ARD") return(n*(n-1))
	
}
