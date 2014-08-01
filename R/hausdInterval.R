hausdInterval <- 
function(X, m, B=30, alpha=0.05, parallel=FALSE, printStatus=FALSE){
     
     if (!is.numeric(X) && !is.data.frame(X)) stop("X should be a matrix of coordinates")
     if (!is.numeric(m)) stop("m should be a number")
     if (!is.numeric(B)) stop("B should be a number")
     if (!is.numeric(alpha)) stop("alpha should be a number")
     if (!is.logical(parallel)) stop("parallel should be logical")
     if (!is.logical(printStatus)) stop("printStatus should be logical")

     X=as.matrix(X)
     n = nrow(X)
     width = rep(0,B)
     if (parallel) {boostLapply=mclapply
     	} else boostLapply=lapply


     if (printStatus) cat("Bootstrap: ")
     width=boostLapply(1:B, FUN=function(i){
          I = sample(1:n,replace=FALSE,size=m)
          Y = as.matrix(X[I,])
          LL = max(knnx.dist(Y, X, k=1, algorithm="kd_tree"))
          if (printStatus) cat(i," ")
     	  return(LL)
     	}
     	)
     if (printStatus) cat("\n")
     width=unlist(width)
	 width = 2*quantile(width,1-alpha)
     
     out=width
          
     return(out)
}