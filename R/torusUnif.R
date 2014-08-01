torusUnif<-
function(n,a,c){

	if (!is.vector(n) || length(n)!=1 || !is.numeric(n)) stop("n should be a integer")	
	if (!is.vector(a) || length(a)!=1) stop("a should be a number")	
	if (!is.vector(c) || length(c)!=1) stop("c should be a number")	

	theta=NULL
	while (length(theta)<n){
		xvec=runif(1,0,2*pi)
		yvec=runif(1,0,1/pi)
		fx=(1+(a/c)*cos(xvec))/(2*pi)
		if (yvec<fx) theta=c(theta, xvec)
	}

	phi=runif(n,0,2*pi)
	x=(c+a*cos(theta))*cos(phi)
	y=(c+a*cos(theta))*sin(phi)
	z=a*sin(theta)
	
	out=cbind(x,y,z)
	return(out)
}
