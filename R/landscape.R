landscape<-
function(Diag, dimension=1, KK=1, tseq=seq( min(Diag[,2:3]), max(Diag[,2:3]) , length=500) )
{
    
    if (class(Diag)!="diagram" && class(Diag)!="matrix" && !is.data.frame(Diag) )  
    	stop("Diag should be a diagram, or a P by 3 matrix")
	if (!is.vector(dimension) || length(dimension)!=1) stop("dimension should be an integer")
	if (!is.vector(KK) || length(KK)!=1) stop("KK should be an integer")
	if (!is.vector(tseq) || !is.numeric(tseq)) stop("tseq should be a numeric vector")
    
    isNA=length(which(Diag[, 1] == dimension))
    if (isNA==0) return(rep(0, length(tseq))) #in case there are no features with dimension "dimension"
    	
	Diag=Diag[which(Diag[,1]==dimension),]
	if (class(Diag)!="matrix") Diag=t(Diag) #in the case there is only 1 point
	
	Npoints=nrow(Diag)

    fab = matrix(NA, nrow = length(tseq), ncol = Npoints)
    lambda = numeric()
    for (j in 1:Npoints) {    
        fab[,j]=sapply(1:length(tseq), FUN=function(i){
        	max(min(tseq[i] - Diag[j, 2], Diag[j,3] - tseq[i]), 0)        	    
        })        
    }
    lambda=sapply(1:length(tseq),  FUN=function(i){
    	sort(fab[i, ], decreasing = TRUE)[KK]  	
    })
    return(lambda)
}