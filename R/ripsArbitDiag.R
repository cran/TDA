ripsArbitDiag <-
function(distX,maxdimension, maxscale, printProgress=FALSE){
	
	if (!is.numeric(distX) && !is.matrix(distX) && !is.data.frame(distX)) stop("distX should be a matrix of distances")
	if (!is.vector(maxdimension) || length(maxdimension)!=1) stop("maxdimension should be an integer")
	if (!is.vector(maxscale) || length(maxscale)!=1) stop("maxscale should be a number")
	if (!is.logical(printProgress)) stop("printProgress should be logical")

		
	write.table(distX,"inputTDA.txt", row.names=F, col.names=F, sep=" " )	
	out1=.C("ripsArbit", as.integer(maxdimension+1), as.double(maxscale), as.integer(printProgress), dup=FALSE, package="TDA")
	Diag=as.matrix(read.table("outputTDA.txt", sep=""))
	
	N=dim(Diag)[1]
	remove=NULL  # we remove points with lifetime=0
	for (i in 1:N){
		if (Diag[i,2]==Diag[i,3]) remove=c(remove,i)
	}
	#remove points with lifetime=0
	if (!is.null(remove)) Diag=Diag[-remove,]  
	#Diag[which(Diag==Inf)]=maxscale	

	if (class(Diag)!="matrix") Diag=t(Diag) #in the case there is only 1 point
	colnames(Diag)=c("dimension","Birth", "Death")
	class(Diag)="diagram"
	attributes(Diag)$maxdimension=maxdimension
	attributes(Diag)$scale=c(0, maxscale)
	attributes(Diag)$call=match.call()
	Diag[1,3]=maxscale

	return(list("diagram"=Diag))
}
