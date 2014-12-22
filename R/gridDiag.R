gridDiag <-
function(X, FUN, lim, by, maxdimension=length(lim)/2-1, sublevel=TRUE, library="Dionysus", location=FALSE, printProgress=FALSE, diagLimit=NULL, ...){


  if (library=="dionysus" || library=="DIONYSUS") library="Dionysus"
  if (library=="phat" || library=="Phat") library="PHAT"
  if (library!="Dionysus" && library!="PHAT") stop("library should be either 'Dionysus' or 'PHAT'")   
  
	if (!is.function(FUN)) stop("FUN should be a function")	
	if (!is.numeric(X) && !is.data.frame(X)) stop("X should be a matrix of coordinates")
  if (!is.numeric(lim) || (length(lim) %% 2 != 0)) stop("lim should be either a matrix or a vector of even elements")
	if (!is.numeric(by) || !is.vector(by) || (length(by)!=1 && length(by)!=length(lim)/2) || !all(by>0) ) stop("by should be either a positive number or a positive vector of length equals dimension of grid")
	if (!is.vector(maxdimension) || length(maxdimension)!=1 || maxdimension<0) stop("maxdimension should be a nonnegative number")  
	if (!is.logical(sublevel)) stop("sublevel should be logical")
	if (!is.logical(printProgress)) stop("printProgress should be logical")
	if (!is.null(diagLimit) && (!is.vector(diagLimit) || length(diagLimit)!=1) ) stop("diagLimit should be a positive number")	

	# in case there is only 1 point
	if (is.vector(X)) X=t(X)
	
	if (maxdimension >= length(lim)/2)
	{
	  maxdimension = length(lim)/2-1
	}

	Grid=gridBy(lim, by=by)
	p=FUN(X,Grid$grid,...)
		
	dim=Grid$dim

  gridValues=p
  if (sublevel==FALSE) gridValues=-p
	
	#write input.txt and read output.txt
  if (ncol(X)<=3)
  {
  	computeGrid=.C("grid",FUNvaluesInput=as.double(gridValues),gridDimensionInput=as.integer(length(dim)),gridNumberInput=as.integer(dim),maxdimensionInput=as.integer(maxdimension),decompositionInput=as.character("5tetrahedra"),libraryInput=as.character(library),locationInput=as.integer(location),printInput=as.integer(printProgress),
                   dup=FALSE, package="TDA")
  }
  else
  {
    computeGrid=.C("grid",FUNvaluesInput=as.double(gridValues),gridDimensionInput=as.integer(length(dim)),gridNumberInput=as.integer(dim),maxdimensionInput=as.integer(maxdimension),decompositionInput=as.character("barycenter"),libraryInput=as.character(library),locationInput=as.integer(location),printInput=as.integer(printProgress),
                   dup=FALSE, package="TDA")
  }
	
  out=read.table("outputTDA.txt", sep="\n")

  ##convert outputTDA.txt in matrix format
 	vecOut=as.vector(out$V1)
  if (location==FALSE)
  {
 	  whichDimens=c((1:length(vecOut))[!grepl(" ",vecOut)], length(vecOut)+1)
  } else
  {
    whichDimens=(1:length(vecOut))[!grepl("[ C]",vecOut)]
  }
 	dim=NULL
   if (length(whichDimens)>1)
   {
   	for (i in 1:(length(whichDimens)-1)){
   		dim=c(dim, rep(as.numeric(vecOut[whichDimens[i]]), whichDimens[i+1]-whichDimens[i]-1))
   		}
   }
  if (location==FALSE)
  {
 	  life=out[grep(" ",vecOut),1]
  } else
  {
    life=(out[grep(" ",vecOut),1])[1:length(dim)]
  }
  life2=matrix(as.numeric(unlist(strsplit(as.character(life)," "))),ncol=2, nrow=length(dim), byrow=TRUE)
	Diag=cbind(dim,life2)
  if (location==TRUE)
  {
    loc=out[(grep("L",vecOut)+1):(grep("C",vecOut)-1),1]
    loc2=matrix(as.numeric(unlist(strsplit(as.character(loc)," "))),ncol=2, nrow=length(loc), byrow=TRUE)
    BirthLocation=Grid$grid[loc2[,1],]
    DeathLocation=Grid$grid[loc2[,2],]
    
    if (library=="Dionysus")
    {
      cycle=c("",as.character(out[(grep("C",vecOut)+1):nrow(out),1]))
      CycleLocation = lapply(strsplit(as.character(cycle)," "),function(c){Grid$grid[as.numeric(c),]})
    }
  }
	  
  if (nrow(Diag)>0) {
  	if (sublevel==FALSE) {
  		colnames(Diag)=c("dim", "Death", "Birth")
  		Diag[,2:3]=-Diag[,2:3]
  		Diag[1,3]= ifelse(is.null(diagLimit), min(p), diagLimit)
  	} else {
  		colnames(Diag)=c("dim", "Birth", "Death")
  		Diag[1,3]=ifelse(is.null(diagLimit), max(p), diagLimit) 
  	}
  	if (sublevel==FALSE) Diag[,2:3]=Diag[,3:2]
  }

	if (class(Diag)!="matrix") Diag=t(Diag) #in the case there is only 1 point
	class(Diag)="diagram"
	attributes(Diag)$maxdimensionension=max(Diag[,1])
	nonInf=which(Diag[,2]!=Inf & Diag[,3]!=Inf)
	attributes(Diag)$scale=c(min(Diag[nonInf,2:3]), max(Diag[nonInf,2:3]))
	attributes(Diag)$call=match.call()
  if (location==FALSE)
  {
    out=list("diagram"=Diag)
  } else if (library=="PHAT")
  {
    out=list("diagram"=Diag,"birthLocation"=BirthLocation,"deathLocation"=DeathLocation)
  } else
  {
    out=list("diagram"=Diag,"birthLocation"=BirthLocation,"deathLocation"=DeathLocation,"cycleLocation"=CycleLocation)
  }
	return(out)
}
