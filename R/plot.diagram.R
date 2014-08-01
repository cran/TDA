plot.diagram <-
function(x, diagLim=NULL, dimension=NULL, col=NULL, rotated=FALSE, barcode=FALSE, band=NULL, add=FALSE,...){

    if (class(x)!="diagram" && class(x)!="matrix" && !is.data.frame(x) )  
    	stop("x should be a diagram, or a P by 3 matrix")
	if (!is.null(diagLim) && (!is.vector(diagLim) || length(diagLim)!=2) ) stop("diagLim should be a vector of length 2")
	if (!is.null(dimension) && (!is.vector(dimension) || length(dimension)!=1) ) stop("dimension should be an integer")
     if (!is.logical(rotated)) stop("rotated should be logical")
     if (!is.logical(barcode)) stop("barcode should be logical")
	if (!is.null(band) && (!is.vector(band) || length(band)!=1) ) stop("band should be a  number")
     if (!is.logical(add)) stop("add should be logical")

	################################################
	if (is.null(diagLim) && class(x)=="diagram"){
		diagLim=attributes(x)$scale
	} else if (is.null(diagLim)) {
		diagLim=c(min(x[,2:3]), max(x[,2:3]))
	}
	
	sublevel=FALSE
	if (colnames(x)[2]=="Birth"){
		sublevel=TRUE	
	}
	
	if (!is.null(dimension)) x=x[which(x[,1]==dimension),]
	
	symb=x[,1]
	for (i in 1:length(symb)){
		if (symb[i]==0) symb[i]=16
		else if (symb[i]==1) symb[i]=2
		else if (symb[i]==2) symb[i]=5
		else if (symb[i]==5) symb[i]=1
		}

	if (is.null(col)){
		col=x[,1]+1						# betti0 black, betti1 red
		for (i in 1:length(symb)){
			if (symb[i]==5) col[i]=4		# betti2 blue
			if (symb[i]==3) col[i]=3		# betti3 green
			}
	}

	### barcode plot
	if (barcode){
		if (length(col)==1) col=rep(col,nrow(x))
		## first we sort the bars
		maxD=max(x[,1])
		minD=min(x[,1])
		if (maxD>0){
			sortedDiag=x
			sortedCol=col
			posD=which(x[,1]==minD)
			lD=0
			for (dd in (minD):maxD){
				oldlD=lD
				posD=which(x[,1]==dd)
				if (length(posD)!=0){
					lD=oldlD+length(posD)				
					sortedDiag[(oldlD+1):(lD),]=x[posD,]
					sortedCol[(oldlD+1):(lD)]=col[posD]
				}
			}
			x=sortedDiag
			col=sortedCol
		}	
			
		## now we plot the bars
		left=x[,2]
		right=x[,3]
		n = length(left)

		Bmax = max(right)
		Bmin = min(left)
		plot(c(Bmin,Bmax),c(1,n+1), type="n", xlab="", ylab="", 
		xlim=c(Bmin,Bmax),ylim=c(0,n+1),
		xaxt="n", yaxt="n", ...)
		axis(1)
		
		lwid=rep(2,n)
		ltype=rep(1,n)
		if (!is.null(band)){
			for(i in 1:n){
				if ((x[i,3]-x[i,2])<=band) {
					ltype[i]=3
					lwid[i]=1.5
				}
			}
		}
		
		segments(left,1:n,right,1:n, lwd=lwid, lty=ltype, col=col)
		
			
	} else{  ### diagram plot

		if (rotated==TRUE){

			if (add==FALSE) plot(0, 0,type="n", axes=F, xlim=diagLim, ylim=diagLim, xlab=" ", ylab=" ", ...)
			if (!is.null(band)) 	polygon(c(0,diagLim[2]+1, diagLim[2]+1,0),c(0,0,band,band),col="pink", lwd=0.01, border="white")

			points((x[,2]+x[,3])/2, (x[,3]-x[,2])/2 ,col=col,pch=symb,lwd=2,cex=1)
		} else{

			if (add==FALSE) plot(0, 0,type="n", axes=F, xlim=diagLim, ylim=diagLim, xlab=" ", ylab=" ", ...)
			if (!is.null(band)) 	polygon(c(0,diagLim[2]+1,diagLim[2]+1,0),c(0,diagLim[2]+1, diagLim[2]+1+band,band),col="pink", lwd=0.01, border="white")


			points(x[,2],x[,3],pch=symb,lwd=2, cex=1, col=col)
			abline(0,1)
		}
		
		if (add==FALSE){
			axis(1)
			axis(2)
			if (sublevel) {
				if (!rotated){
					title(main="", xlab="Birth", ylab="Death")
				}else title(main="", ylab="(Death-Birth)/2", xlab="(Death+Birth)/2")
			} 
			if (!sublevel) {
				if (!rotated){
					title(main="", xlab="Death", ylab="Birth")
				}else title(main="", ylab="(Birth-Death)/2", xlab="(Death+Birth)/2")
			} 
		}
	}
}
