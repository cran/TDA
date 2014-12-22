maxPersistence<-
function(FUN, parameters, X, lim, by, maxdimension=length(lim)/2-1, sublevel=TRUE, library="Dionysus", B=30, alpha=0.05, bandFUN="bootstrapBand", distance="bottleneck", dimension=1, p=1, parallel=FALSE, printProgress=FALSE){

	if (!is.function(FUN)) stop("FUN should be a function")	
	if (!is.character(bandFUN)) stop("bandFUN should be a string: either 'bootstrapBand' or 'bootstrapDiagram'")	
	if (!is.vector(parameters) || !is.numeric(parameters)) stop("parameters should be a numeric vector")
	if (!is.numeric(X) && !is.data.frame(X)) stop("X should be a matrix of coordinates")
	if (2*ncol(X)!=length(lim)) stop("dimension of X does not match with lim")	
	if (!is.vector(by) || length(by)!=1) stop("by should be a positive number")
	if (!is.logical(sublevel)) stop("sublevel should be logical")
	if (!is.vector(B) || length(B)!=1) stop("B should be an integer")
	if (!is.vector(alpha) || length(alpha)!=1) stop("alpha should be a number between 0 and 1")
     if (!is.logical(parallel)) stop("parallel should be logical")
     if (!is.logical(printProgress)) stop("printProgress should be logical")


	# in case there is only 1 point
	if (is.vector(X)) X=t(X)

	Grid=gridBy(lim, by=by)$grid
	
	Kseq=length(parameters)
	eps=numeric(Kseq)
	numberSignificant=numeric(Kseq)
	significantPers=numeric(Kseq)
	
	Pers=list()
	
	if (printProgress)
	{
		cat("0   10   20   30   40   50   60   70   80   90   100\n")
		cat("|----|----|----|----|----|----|----|----|----|----|\n")
		cat("*")		
	}
	percentageFloor=0

	for (i in 1:Kseq){
		
		Diag= gridDiag(X=X, FUN=FUN, lim=lim, by=by, maxdimension=maxdimension, sublevel=sublevel, library=library, location=FALSE, printProgress=FALSE, diagLimit=NULL, parameters[i])$diagram
		Diag[1,3]=Diag[1,2] #remove first component with infinite persistence
		Pers[[i]]=cbind(Diag[,1], Diag[,3]-Diag[,2])
		colnames(Pers[[i]])=c("dimension", "Persistence")
			
		if (B>0){
			if (bandFUN=="bootstrapDiagram"){
				eps[i] = bootstrapDiagram(X=X, FUN=FUN, lim=lim, by=by, sublevel=sublevel, library=library, B=B, alpha=alpha, distance=distance, dimension=dimension, p=p, printProgress=FALSE, parameters[i])
				selectDimension=which(Pers[[i]][,1]==dimension)
				numberSignificant[i]=sum( Pers[[i]][selectDimension,2]> (2*eps[i]) )
				significantPers[i]= sum(pmax(0, Pers[[i]][selectDimension,2]-(2*eps[i])))
			}else {
				eps[i] = bootstrapBand(X, FUN, Grid, B=B, alpha=alpha, parallel=parallel, printProgress=F, parameters[i])$width
				numberSignificant[i]=sum( Pers[[i]][,2]> (2*eps[i]) )
				significantPers[i]= sum(pmax(0, Pers[[i]][,2]-(2*eps[i])))
			}
		} else eps[i]=0
		

		if (printProgress && floor((100*i/Kseq-percentageFloor)/2)>0)
		{
			for (j in 1:(floor((100*i/Kseq-percentageFloor)/2)))
			{
				cat("*")
				percentageFloor=percentageFloor+2
			}

		}
	}
	if (printProgress)
	{ 
		cat("\n")
	}


	# Two criterions
	Param1=parameters[which(numberSignificant==max(numberSignificant))]
	Param2=parameters[which(significantPers==max(significantPers))]
	
	out=list("parameters"=parameters, "sigNumber"=numberSignificant, "sigPersistence"=significantPers, "bands"=eps, "Persistence"=Pers, "bandFUN"=bandFUN, "dimension"=dimension)
	
	class(out)="maxPersistence"
	attributes(out)$call=match.call()

	return(out)
}
