bootstrapDiagram <- 
function(X, FUN, lim, by, sublevel=TRUE, library="Dionysus", B=30, alpha=0.05, distance="bottleneck", dimension=1, p=1, printProgress=FALSE, ...){
     
	if (!is.numeric(X) && !is.data.frame(X)) stop("X should be a matrix of coordinates")
	if (class(FUN)!="function") stop("FUN should be function")
	if (2*ncol(X)!=length(lim)) stop("dimension of X does not match with lim")
	if (!is.vector(by) || length(by)!=1) stop("by should be a positive number")
	if (!is.logical(sublevel)) stop("sublevel should be logical")
	if (!is.numeric(B) || length(B)!=1 || B<1) stop("B should be a positive number")
	if (!is.numeric(alpha)) stop("alpha should be a number")
	if (!is.numeric(dimension) || length(dimension)!=1) stop("dimension should be a positive integer")
	if (!is.vector(p) || length(p)!=1 || p < 1) stop("p should be a positive integer")	
     # if (!is.logical(parallel)) stop("parallel should be logical")
     if (!is.logical(printProgress)) stop("printProgress should be logical")


	if (distance!="wasserstein" && distance!="bottleneck") stop("distance should be a string: either 'bottleneck' or 'wasserstein'")

     X=as.matrix(X)
     n = nrow(X)

	maxdimension=dimension
     parallel=FALSE
     
     Diag=gridDiag(X=X, FUN=FUN, lim=lim, by=by, maxdimension=maxdimension, sublevel=sublevel, library=library, printProgress=FALSE, diagLimit=NULL, ...)$diagram

     if (parallel) {boostLapply=mclapply
     	} else boostLapply=lapply

     if (printProgress) cat("Bootstrap: ")
     
     if (distance=="wasserstein") {
	     width=boostLapply(1:B, FUN=function(i){
	          I = sample(1:n,replace=TRUE,size=n)
	          Y = as.matrix(X[I,])
	          Diag1 = gridDiag(X=Y, FUN=FUN, lim=lim, by=by, maxdimension=maxdimension, sublevel=sublevel, library=library, printProgress=FALSE, diagLimit=NULL, ...)$diagram
	          	width1 = wasserstein(Diag,Diag1, p=p,dimension=dimension)
	          if (printProgress) cat(i," ")
	     	  return(width1)
	     })
     } else{
     	width=boostLapply(1:B, FUN=function(i){
          I = sample(1:n,replace=TRUE,size=n)
          Y = as.matrix(X[I,])
          Diag1 = gridDiag(X=Y, FUN=FUN, lim=lim, by=by, maxdimension=maxdimension, sublevel=sublevel, library=library, printProgress=FALSE, diagLimit=NULL, ...)$diagram
          width1 = bottleneck(Diag,Diag1, dimension=dimension)
          if (printProgress) cat(i," ")
     	  return(width1)
     	})
     }
     	     		
     if (printProgress) cat("\n")
     width=unlist(width)
	 width = quantile(width,1-alpha)
     
     return(width)
}