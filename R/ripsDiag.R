ripsDiag <-
function(X, maxdimension, maxscale, dist="euclidean", printStatus=FALSE){

	if (!is.numeric(X) && !is.data.frame(X)) stop("X should be a matrix of coordinates")
	if (!is.vector(maxdimension) || length(maxdimension)!=1) stop("maxdimension should be an integer")
	if (!is.vector(maxscale) || length(maxscale)!=1) stop("maxscale should be a number")
	if (!is.logical(printStatus)) stop("printStatus should be logical")

	if (dist=="euclidean"){
		#in 32bit architectures ripsL2Diag doesn't work
		#Diag=ripsL2Diag(X,maxdimension,maxscale, printStatus)		
		distX=as.matrix(dist(X))
		Diag=ripsArbitDiag(distX,maxdimension,maxscale, printStatus)		
	} else if (dist=="arbitrary"){
		Diag=ripsArbitDiag(X,maxdimension,maxscale, printStatus)
	} else stop ("dist should be 'euclidean' or 'arbitrary'")

	return(Diag)
}
