findKtree <-
function(bb, parent, sons, compBranch, n){
	Ktop=numeric(bb)
	Kbottom=numeric(bb)
	for (i in 1:bb){
		if (is.na(parent[i])){
			Kbottom[i]=NA
			Ktop[i]=NA
		} else if (parent[i]==0){ 
				Kbottom[i]=0	
				if (i <= length(sons)){
					Ksons=sons[[i]]
					if (!is.null(Ksons) && !is.na(Ksons)){
					Ktop[i]=(length(compBranch[[i]])-length(unlist(compBranch[Ksons])) )/n
					} else Ktop[i]=length(compBranch[[i]])/n
				} else Ktop[i]=length(compBranch[[i]])/n
			} else{
				Kbottom[i]=Ktop[parent[i]]
				if (i <= length(sons)){
					Ksons=sons[[i]]
					if (!is.null(Ksons) && length(Ksons)!=0 && !is.na(Ksons)){
					Ktop[i]=Kbottom[i]+(length(compBranch[[i]])-length(unlist(compBranch[Ksons])) )/n
					} else Ktop[i]=Kbottom[i]+length(compBranch[[i]])/n
				} else Ktop[i]=Kbottom[i]+length(compBranch[[i]])/n
			}
		}

	return(list("Ktop"=Ktop,"Kbottom"=Kbottom))
}
