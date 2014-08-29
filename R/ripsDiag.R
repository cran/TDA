ripsDiag <-
function(X, maxdimension, maxscale, dist="euclidean", library="GUDHI", printProgress=FALSE){
	
	if (library=="dionysus" || library=="DIONYSUS") library="Dionysus"
	if (library=="gudhi" || library=="Gudhi") library="GUDHI"
	if (library!="Dionysus" && library!="GUDHI") stop("library should be either 'Dionysus' or 'GUDHI'")   
	
	
	if (!is.numeric(X) && !is.data.frame(X)) stop("X should be a matrix of coordinates")
	if (!is.vector(maxdimension) || length(maxdimension)!=1) stop("maxdimension should be an integer")
	if (!is.vector(maxscale) || length(maxscale)!=1) stop("maxscale should be a number")
	if (!is.logical(printProgress)) stop("printProgress should be logical")

	if (dist=="euclidean"){
		#in 32bit architectures ripsL2Diag doesn't work
		#Diag=ripsL2Diag(X,maxdimension,maxscale, library,printProgress)		
		if (library=="GUDHI")
		{
			Diag=ripsL2Diag(X,maxdimension,maxscale, library,printProgress)			
		} else{
			distX=as.matrix(dist(X))
			Diag=ripsArbitDiag(distX,maxdimension,maxscale, printProgress)		
		}
	} else if (dist=="arbitrary"){
		Diag=ripsArbitDiag(X,maxdimension,maxscale, printProgress)
	} else stop ("dist should be 'euclidean' or 'arbitrary'")

	return(Diag)
}
