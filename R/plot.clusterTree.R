plot.clusterTree <-
function(x, type="lambda",color=NULL, add=FALSE, ...){
	
	if (class(x)!="clusterTree") stop("Tree should be an object of class clusterTree")
	if (!is.null(color) && (!is.vector(color) || length(color)!=1) ) stop("color should have length 1")
     if (!is.logical(add)) stop("add should be logical")
	
	id=x$id
	base=x$Xbase
	
	if (type=="lambda"){
		bottom=x$Ybottom
		top=x$Ytop
	} else if (type=="kappa"){
		bottom=x$Kbottom
		top=x$Ktop
	} else stop("type should be 'lambda' or 'kappa'")

	sons=x$sons
	
	Ylim=max(top[which(!is.na(top))])
	
	if (!add)
	plot(c(0,0),c(0,Ylim), type="n", xlim=c(0,1), ylab=type, xlab="", ...)
	
	#vertical lines
	if (is.null(color)) color=id
	segments(base,bottom,base,top, col=color, lwd=3)

	## now the horizontal lines
	if (length(sons)>0){	
		for (i in 1:length(sons)){
			if (!is.null(sons[[i]]) && length(sons[[i]]>0) && !is.na(sons[[i]])){
				x=c(min(base[sons[[i]]]) , max(base[sons[[i]]]) )
				y=bottom[sons[[i]][1]]
				segments(x[1],y,x[2],y, lwd=3, col=ifelse(is.null(color),1,color))
			}	
		}
	}
}

