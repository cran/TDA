gridDiag <-
function(X, FUN, Xlim, Ylim=NA, Zlim=NA, by=(Xlim[2]-Xlim[1])/20, sublevel=TRUE, printStatus=FALSE, diagLimit=NULL, ...){

	if (!is.function(FUN)) stop("FUN should be a function")	
	if (!is.numeric(X) && !is.data.frame(X)) stop("X should be a matrix of coordinates")
	if (!is.vector(Xlim) || length(Xlim)!=2) stop("Xlim should be vector of length 2")
	if (!is.na(Ylim) && (!is.vector(Ylim) || length(Ylim)!=2) ) stop("Ylim should be vector of length 2")
	if (!is.na(Zlim) && (!is.vector(Zlim) || length(Zlim)!=2) ) stop("Zlim should be vector of length 2")
	if (!is.vector(by) || length(by)!=1) stop("by should be a positive number")
	if (!is.logical(sublevel)) stop("sublevel should be logical")
	if (!is.logical(printStatus)) stop("printStatus should be logical")
	if (!is.null(diagLimit) && (!is.vector(diagLimit) || length(diagLimit)!=1) ) stop("diagLimit should be a positive number")	

	# in case there is only 1 point
	if (is.vector(X)) X=t(X)
	
	if (ncol(X)>3) stop("gridDiag currently works in dimension 1,2, or 3")

	Grid=gridBy(Xlim, Ylim, Zlim, by=by)
	p=FUN(X,Grid$grid,...)
		
	dim=Grid$dim
	
	#convert p into grid format
	Dim1=dim[1]
	Dim2=dim[2]
	Dim3=dim[3]
	ppp=array(p,dim)   
	ppp=apply(ppp,2,c)
	if (Dim2==1 & Dim3==1) {
		ppp=t(ppp)
		ppp=rbind(c(dim, rep(NA,Dim1-3)), ppp)
	} else{
		ppp=rbind(c(dim,rep(NA,Dim2-3)),ppp)  #add dimension of the grid
	}
	gridValues=ppp
	gridValues[1,]=ppp[1,]	
	if (sublevel==FALSE) gridValues[-1,]=max(p)-ppp[-1,]	
	
	#write input.txt and read output.txt
	write.table(gridValues,"inputDionysus.txt", row.names=F, col.names=F, sep=" " )
	computeGrid=.C("grid", as.integer(printStatus),dup=FALSE, package="persistence")
	out=read.table("outputDionysus.txt", sep="\n")
	
	##convert output.txt in matrix format
	vecOut=as.vector(out$V1)
	whichDimens=c((1:length(vecOut))[-grep(" ",vecOut)], length(vecOut)+1)
	dim=NULL
	for (i in 1:(length(whichDimens)-1)){
		dim=c(dim, rep(i-1, whichDimens[i+1]-whichDimens[i]-1))
		}
	out2=data.frame(dim,life=out[grep(" ",vecOut),1])
	life2=matrix(NA, ncol=2, nrow=length(dim))
	for (i in 1:length(dim)){
		life2[i,]=as.numeric(unlist(strsplit(as.character(out2[i,2]), " "))	)
		}
	
	Diag=cbind(dim,life2)
	if (sublevel==FALSE) {
		colnames(Diag)=c("dim", "Death", "Birth")
		Diag[,2:3]=max(p)-Diag[,2:3]
		Diag[1,3]= ifelse(is.null(diagLimit), 0, diagLimit) 
	} else {
		colnames(Diag)=c("dim", "Birth", "Death")
		Diag[1,3]=ifelse(is.null(diagLimit), max(p), diagLimit) 
	}
	if (sublevel==FALSE) Diag[,2:3]=Diag[,3:2]

	if (class(Diag)!="matrix") Diag=t(Diag) #in the case there is only 1 point
	class(Diag)="diagram"
	attributes(Diag)$maxdimension=max(Diag[,1])
	nonInf=which(Diag[,2]!=Inf & Diag[,3]!=Inf)
	attributes(Diag)$scale=c(min(Diag[nonInf,2:3]), max(Diag[nonInf,2:3]))
	attributes(Diag)$call=match.call()
	return(Diag)
}
