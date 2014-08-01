bootstrapBand <- 
function(X, FUN, Grid, B=30, alpha=0.05, parallel=FALSE, printStatus=FALSE, ...){
     
     if (!is.numeric(X) && !is.data.frame(X)) stop("X should be a matrix of coordinates")
     if (class(FUN)!="function") stop("FUN should be function")
     if (!is.numeric(Grid) && !is.data.frame(Grid)) stop("Grid should be a matrix of coordinates")
     if (!is.numeric(B) || length(B)!=1 || B<1) stop("B should be a positive number")
     if (!is.numeric(alpha)) stop("alpha should be a number")
     if (!is.logical(parallel)) stop("parallel should be logical")
     if (!is.logical(printStatus)) stop("printStatus should be logical")

     ff = FUN(X,Grid, ...)      
     X=as.matrix(X)
     n = nrow(X)
     m = nrow(Grid)
     width = rep(0,B)
     if (parallel) {boostLapply=mclapply
     	} else boostLapply=lapply


     if (printStatus) cat("Bootstrap: ")
     width=boostLapply(1:B, FUN=function(i){
          I = sample(1:n,replace=TRUE,size=n)
          Y = as.matrix(X[I,])
          bootF = FUN(Y,Grid, ...)
          width1 = max(abs(ff-bootF))
          if (printStatus) cat(i," ")
     	  return(width1)
     	}     	
     	)
     if (printStatus) cat("\n")
     width=unlist(width)
	 width = quantile(width,1-alpha)
     
     UPband=ff+width
     LOWband=ff-width
     #LOWband[which(LOWband<0)]=0   #set negative values of lower band =0
     Band=cbind(LOWband, UPband)
     
     out=list("width"=width, "fun"=ff, "band"=Band)
          
     return(out)
}