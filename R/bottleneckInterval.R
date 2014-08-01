bottleneckInterval <- 
function(X, FUN, Xlim, Ylim=NA, Zlim=NA, by=(Xlim[2]-Xlim[1])/20, sublevel=TRUE, B=30, alpha=0.05, dimension=1, printStatus=FALSE, ...){
     
     if (!is.numeric(X) && !is.data.frame(X)) stop("X should be a matrix of coordinates")
     if (class(FUN)!="function") stop("FUN should be function")

	if (!is.vector(Xlim) || length(Xlim)!=2) stop("Xlim should be vector of length 2")
	if (!is.na(Ylim) && (!is.vector(Ylim) || length(Ylim)!=2) ) stop("Ylim should be vector of length 2")
	if (!is.na(Zlim) && (!is.vector(Zlim) || length(Zlim)!=2) ) stop("Zlim should be vector of length 2")
	if (!is.vector(by) || length(by)!=1) stop("by should be a positive number")
	if (!is.logical(sublevel)) stop("sublevel should be logical")
     
     if (!is.numeric(B) || length(B)!=1 || B<1) stop("B should be a positive number")
     if (!is.numeric(alpha)) stop("alpha should be a number")

     if (!is.numeric(dimension) || length(dimension)!=1) stop("dimension should be a positive integer")
     # if (!is.logical(parallel)) stop("parallel should be logical")
     if (!is.logical(printStatus)) stop("printStatus should be logical")


     X=as.matrix(X)
     n = nrow(X)
     
     parallel=FALSE
     
     Diag=gridDiag(X, FUN, Xlim=Xlim, Ylim=Ylim, by=by, sublevel=sublevel, printStatus=FALSE, ...)

     if (parallel) {boostLapply=mclapply
     	} else boostLapply=lapply

     if (printStatus) cat("Bootstrap: ")
     width=boostLapply(1:B, FUN=function(i){
          I = sample(1:n,replace=TRUE,size=n)
          Y = as.matrix(X[I,])
          Diag1 = gridDiag(Y, FUN, Xlim=Xlim, Ylim=Ylim, by=by, sublevel=sublevel, printStatus=FALSE, ...)
          width1 = bottleneck(Diag,Diag1, dimension=dimension)
          if (printStatus) cat(i," ")
     	  return(width1)
     	}     	
     	)
     if (printStatus) cat("\n")
     width=unlist(width)
	 width = quantile(width,1-alpha)
     
     return(width)
}