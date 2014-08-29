kde <-
function(X,Grid,h, printProgress=FALSE){

	if (!is.numeric(X) && !is.data.frame(X)) stop("X should be a matrix of coordinates")
	if (!is.numeric(Grid) && !is.data.frame(Grid)) stop("Grid should be a matrix of coordinates")
	if (!is.vector(h) || length(h)!=1) stop("h should be a positive number")
	if (!is.logical(printProgress)) stop("printProgress should be a logical variable")
	    
    X=as.matrix(X) 
	Grid=as.matrix(Grid)
	out=rep(0, nrow(Grid))
	outR= .C("kde", as.double(X), as.integer(nrow(X)), as.integer(ncol(X)), as.double(Grid), as.integer(nrow(Grid)), as.double(h), as.integer(printProgress),
	"out"=as.double(out), dup=FALSE, package="TDA")$out
	return(outR)
}
