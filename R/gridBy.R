gridBy<-
function(Xlim=c(0,1), Ylim=NA, Zlim=NA, by=(Xlim[2]-Xlim[1])/10){
	if (is.numeric(Xlim) && length(Xlim)==2 && is.na(Ylim)){
		grid=matrix(seq(Xlim[1], Xlim[2], by=by), ncol=1)	
		dim=c(length(grid),1,1)
	} else if (is.numeric(Xlim) && 
				length(Xlim)==2 && 
				is.numeric(Ylim) && 
				length(Ylim)==2 && 
				is.na(Zlim)){
					xgrid=seq(Xlim[1], Xlim[2], by=by)
					ygrid=seq(Ylim[1], Ylim[2], by=by)
					dim=c(length(xgrid), length(ygrid),1)
					grid=as.matrix(expand.grid(xgrid, ygrid)[1:2])
			} else if (is.numeric(Xlim) && 
						length(Xlim)==2 && 
						is.numeric(Ylim) && 
						length(Ylim)==2 && 
						is.numeric(Zlim) && 
						length(Zlim)==2) {
							xgrid=seq(Xlim[1], Xlim[2], by=by)
							ygrid=seq(Ylim[1], Ylim[2], by=by)
							zgrid=seq(Zlim[1], Zlim[2], by=by)
							k1=length(xgrid)
							k2=length(ygrid)
							k3=length(zgrid)
							thirdCol=rep(zgrid, each=k1*k2)
							dim=c(k1,k2,k3)
							grid=as.matrix(cbind(expand.grid(xgrid,ygrid)[1:2], thirdCol))
					}
	colnames(grid)=NULL
	out=list("dim"=dim, "Xlim"=Xlim, "Ylim"=Ylim, "Zlim"=Zlim, "grid"=grid)
	return(out)	
}