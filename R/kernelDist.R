kernelDist<-
function(X, Grid, h){

	if (!is.numeric(X) && !is.data.frame(X)) stop("X should be a matrix of coordinates")
	if (!is.numeric(Grid) && !is.data.frame(Grid)) stop("Grid should be a matrix of coordinates")
	if (!is.vector(h) || length(h)!=1) stop("h should be a positive number")
	
    X=as.matrix(X) 
	Grid=as.matrix(Grid)
	out=rep(0, nrow(Grid))
	outR= .C("kdeDist", as.double(X), as.integer(nrow(X)), as.integer(ncol(X)), as.double(Grid), as.integer(nrow(Grid)), as.double(h), "out"=as.double(out),dup=FALSE, package="persistence")$out
	return(outR)
}
