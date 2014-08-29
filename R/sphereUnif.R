sphereUnif<-
function(n,d,r=1){

     if (!is.vector(n) || length(n)!=1) stop("n should be a integer")	
     if (!is.vector(d) || length(d)!=1) stop("d should be a integer")	
     if (!is.vector(r) || length(r)!=1) stop("r should be a number")
     
     ########################
     X = array(rnorm(n*(d+1)),dim=c(n,d+1))
     for (i in 1:n)
     {
       while(norm(X[i,],"2")==0)
       {
         X[i,]=rnorm(d+1)
       }
       X[i,]=r*X[i,]/norm(X[i,],"2")
     }
     return(X)     
}
