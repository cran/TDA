ripsL2Diag <-
function(X,maxdimension, maxscale, printStatus=FALSE){

	if (!is.numeric(X) && !is.data.frame(X)) stop("X should be a matrix of coordinates") 
	if (!is.vector(maxdimension) || length(maxdimension)!=1) stop("maxdimension should be an integer")
	if (!is.vector(maxscale) || length(maxscale)!=1) stop("maxscale should be a number")
	if (!is.logical(printStatus)) stop("printStatus should be logical")


	write.table(X,"inputDionysus.txt", row.names=F, col.names=F, sep=" " )	
	out1=.C("rips", as.integer(maxdimension+1), as.double(maxscale), as.integer(printStatus) ,dup=FALSE, package="persistence")
	Diag=as.matrix(read.table("outputDionysus.txt", sep=""))
	
	N=dim(Diag)[1]
	remove=NULL  # we remove points with lifetime=0
	for (i in 1:N){
		if (Diag[i,2]==Diag[i,3]) remove=c(remove,i)
	}
	#remove points with lifetime=0
	if (!is.null(remove)) Diag=Diag[-remove,]  
	#Diag[which(Diag==Inf)]=maxscale	
	colnames(Diag)=c("dimension","Birth", "Death")
	
	if (class(Diag)!="matrix") Diag=t(Diag) #in the case there is only 1 point
	
	class(Diag)="diagram"
	attributes(Diag)$maxdimension=maxdimension
	attributes(Diag)$scale=c(0, maxscale)
	attributes(Diag)$call=match.call()
	Diag[1,3]=maxscale
		
	return(Diag)
}
