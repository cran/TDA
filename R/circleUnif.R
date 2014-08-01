circleUnif<-
function(n,r=1){

     if (!is.vector(n) || length(n)!=1) stop("n should be a integer")	
     if (!is.vector(r) || length(r)!=1) stop("r should be a number")	

     ########################
     th = runif(n,0,2*pi)
     x1 = r*cos(th)
     x2  = r*sin(th)
     X = cbind(x1,x2)
     return(X)
}
