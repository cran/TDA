maxPersistence<-
function(FUN, parameters, X, Xlim, Ylim=NA, Zlim=NA, by, sublevel=TRUE, B=30, alpha=0.05, parallel=FALSE, printProgress=FALSE){

	if (!is.function(FUN)) stop("FUN should be a function")	
	if (!is.vector(parameters) || !is.numeric(parameters)) stop("parameters should be a numeric vector")
	if (!is.numeric(X) && !is.data.frame(X)) stop("X should be a matrix of coordinates")
	if (!is.vector(Xlim) || length(Xlim)!=2) stop("Xlim should be vector of length 2")
	if (!is.na(Ylim) && (!is.vector(Ylim) || length(Ylim)!=2) ) stop("Ylim should be vector of length 2")
	if (!is.na(Zlim) && (!is.vector(Zlim) || length(Zlim)!=2) ) stop("Zlim should be vector of length 2")
	if (!is.vector(by) || length(by)!=1) stop("by should be a positive number")
	if (!is.logical(sublevel)) stop("sublevel should be logical")
	if (!is.vector(B) || length(B)!=1) stop("B should be an integer")
	if (!is.vector(alpha) || length(alpha)!=1) stop("alpha should be a number between 0 and 1")
     if (!is.logical(parallel)) stop("parallel should be logical")
     if (!is.logical(printProgress)) stop("printProgress should be logical")


	# in case there is only 1 point
	if (is.vector(X)) X=t(X)

	Grid=gridBy(Xlim,Ylim,Zlim,by)$grid
	
	Kseq=length(parameters)
	eps=numeric(Kseq)
	numberSignificant=numeric(Kseq)
	significantPers=numeric(Kseq)
	
	Pers=list()
	
	if (printProgress) cat("Progress: ")
	for (i in 1:Kseq){
		
		Diag= gridDiag(X, FUN, Xlim, Ylim, Zlim, by=by, sublevel=sublevel, printStatus=F, diagLimit=NULL, parameters[i])
		Diag[1,3]=Diag[1,2] #remove first component with infinite persistence
		Pers[[i]]=cbind(Diag[,1], Diag[,3]-Diag[,2])
		colnames(Pers[[i]])=c("dimension", "Persistence")
			
		if (B>0){
		eps[i] = bootstrapBand(X, FUN, Grid, B=B, alpha=alpha, parallel=parallel, printStatus=F, parameters[i])$width
		} else eps[i]=0
		numberSignificant[i]=sum( Pers[[i]][,2]> (2*eps[i]) )
		significantPers[i]= sum(pmax(0, Pers[[i]][,2]-(2*eps[i])))

		if (printProgress) cat(round(i/Kseq,2), " ")
	}
	if (printProgress) cat("\n")

	# Two criterions
	Param1=parameters[which(numberSignificant==max(numberSignificant))]
	Param2=parameters[which(significantPers==max(significantPers))]
	
	out=list("parameters"=parameters, "sigNumber"=numberSignificant, "sigPersistence"=significantPers, "bands"=eps, "Persistence"=Pers)
	
	class(out)="maxPersistence"
	attributes(out)$call=match.call()

	return(out)
}
