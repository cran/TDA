dtm<-
function(X, Grid, m0){

	if (!is.numeric(X) && !is.data.frame(X)) stop("X should be a matrix of coordinates")
	if (!is.numeric(Grid) && !is.data.frame(Grid)) stop("Grid should be a matrix of coordinates")
	if (!is.vector(m0) || length(m0)!=1) stop("m0 should be a number between 0 and 1")
	
    X=as.matrix(X) 
	Grid=as.matrix(Grid)
	k0=ceiling(m0*dim(X)[1])
	distances=knnx.dist(X, Grid, k=k0, algorithm=c("kd_tree"))	
	d2=apply(distances^2, 1, sum)
	dOut=sqrt(d2/k0)
	return(dOut)
}
