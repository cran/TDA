bottleneck <-
function(Diag1, Diag2, dimension){
	
	if (class(Diag1)!="diagram" && class(Diag1)!="matrix" && !is.data.frame(Diag1) ) 
		stop("Diag1 should be a diagram or a matrix")	
	if (class(Diag2)!="diagram" && class(Diag2)!="matrix" && !is.data.frame(Diag1) )  
		stop("Diag2 should be a diagram or a matrix")	
	if (!is.vector(dimension) || length(dimension)!=1) stop("dimension should be a integer")	

	Diag1=Diag1[which(Diag1[,1]==dimension),2:3]
	Diag2=Diag2[which(Diag2[,1]==dimension),2:3]
	
	if (length(Diag1)==0) Diag1=c(dimension,0,0)
	if (length(Diag2)==0) Diag2=c(dimension,0,0)
	
	if (class(Diag1)!="matrix") Diag1=t(Diag1) #in case there is only 1 point
	if (class(Diag2)!="matrix") Diag2=t(Diag2) #in case there is only 1 point		

	out=1
	out1=.C("bottleneck", as.double(t(Diag1)), as.integer(nrow(Diag1)), as.double(t(Diag2)), as.integer(nrow(Diag2)), as.double(out), dup=TRUE, package="TDA")[[5]]
	
	return(out1)
}
