ripsL2Diag <-
function(X,maxdimension, maxscale, library="GUDHI", printProgress=FALSE){

	if (library=="dionysus" || library=="DIONYSUS") library="Dionysus"
	if (library=="gudhi" || library=="Gudhi") library="GUDHI"
	if (library!="Dionysus" && library!="GUDHI") stop("library should be either 'Dionysus' or 'GUDHI'")   

	if (!is.numeric(X) && !is.data.frame(X)) stop("X should be a matrix of coordinates") 
	if (!is.vector(maxdimension) || length(maxdimension)!=1) stop("maxdimension should be an integer")
	if (!is.vector(maxscale) || length(maxscale)!=1) stop("maxscale should be a number")
	if (!is.logical(printProgress)) stop("printProgress should be logical")


	if (library=="Dionysus")
	{
		write.table(X,"inputTDA.txt", row.names=F, col.names=F, sep=" " )	
		out1=.C("rips", as.integer(maxdimension+1), as.double(maxscale), as.integer(printProgress) ,dup=FALSE, package="TDA")
		Diag=as.matrix(read.table("outputTDA.txt", sep=""))		
	} else{
		max_num_pairs=5000  # to be added as an option
		points=as.vector(t(X))
		dim=ncol(X)
		num_points=nrow(X)
		diagram=matrix(-1, ncol=3, nrow=max_num_pairs)
		
		out=.C("rips_persistence_diagram_GUDHI", as.double(points),
	                                       as.integer(dim),
	                                       as.integer(num_points),
	                                       as.double(maxscale),
	                                       as.integer(maxdimension+1), #GUDHI needs +1
	                                       as.double(diagram),
	                                       as.integer(printProgress),
	                                       dup=FALSE, package="TDA")
		Diag=as.matrix(read.table("outputTDA.txt", sep=""))
	}
	
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
			
		return(list("diagram"=Diag))
	
}
