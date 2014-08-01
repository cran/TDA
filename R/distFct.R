distFct<-
function(X, Grid){
	
	if (!is.numeric(X) && !is.data.frame(X)) stop("X should be a matrix of coordinates")
	if (!is.numeric(Grid) && !is.data.frame(Grid)) stop("Grid should be a matrix of coordinates")
		
	X=as.matrix(X) 
	Grid=as.matrix(Grid)
	distances=knnx.dist(X, Grid, k=1, algorithm=c("kd_tree"))	
	dOut=as.vector(distances)
	return(dOut)
}
