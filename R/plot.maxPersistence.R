plot.maxPersistence<-
function(x,...){
	
	parameter=x$parameters
	Pers=x$Persistence
	eps=x$bands
	Kseq=length(parameter)
	maxPers=0
	for (k in 1:Kseq){
		maxPers=max(maxPers, Pers[[k]][,2])	
	}


	plot(parameter, rep(maxPers, Kseq), type="n", ylim=c(0,1.18*maxPers), ylab="Persistence", axes=F, ... )
	axis(1)
	axis(2)
	
	for (i in 1:Kseq){	
		
		symb=Pers[[i]][,1]
		for (j in 1:length(symb)){
			if (symb[j]==0) symb[j]=16
			else if (symb[j]==1) symb[j]=2
			else if (symb[j]==2) symb[j]=5
		}

		col=Pers[[i]][,1]+1						# betti0 black, betti1 red
		for (j in 1:length(symb)){
			if (col[j]==3) col[j]=4		# betti2 blue
			if (symb[j]==3) col[j]=3		# betti3 green
			}

		points(rep(parameter[i], nrow(Pers[[i]])), Pers[[i]][,2], col=col, pch=symb, lwd=2)	
		points(parameter[i], 1.2*maxPers, col=1, pch=16, lwd=2)	
	}
	
	eps=pmin(eps, rep(1.1*maxPers/2, length(eps)) )
	polygon(c(parameter, parameter[Kseq:1]), c(2*eps, rep(0,Kseq)), col=alpha(2,0.1))
	abline(h=1.18*maxPers, lty=2)
	
	# legend(parameter[floor(0.6*Kseq)], maxPers, c("components", "loops", "voids"), pch=c(16,2,5), pt.lwd=2, col=c(1,2,4))
		
}
