plotRule <-
function(Tree){

	if (length(Tree$id)>1){
	goWhile=TRUE
	
		while(goWhile){
	
			# max density of each cluster	
			Tops=numeric(length(Tree$id))
			for (j in Tree$id){
				Tops[j]=max(Tree$density[Tree$DataPoints[[j]]])
			}
		
			uniqueParents=unique(Tree$parent)
			uniqueParentsNo0=setdiff(uniqueParents, 0)
			
			for (i in uniqueParentsNo0){
				
				bros=Tree$sons[[i]]
				newID=bros[order(Tops[bros], decreasing=TRUE)]  # TODO what if same height?
						
				if (sum(bros!=newID)!=0){
				
					NewTree=Tree
					
					#update new IDs
					NewTree$Ytop[bros]=Tree$Ytop[newID]	
					NewTree$Ktop[bros]=Tree$Ktop[newID]	
					
					for (j in 1:length(bros)){
						NewTree$parent[which(Tree$parent==bros[j])]=newID[j]
					}
						
					for (j in 1:length(newID)){
						if (!is.null(Tree$sons[bros[j]][[1]])){
							NewTree$sons[[newID[j]]]= Tree$sons[bros[j]][[1]]
						} else NewTree$sons[[newID[j]]]= NA
					}
					
					for (j in 1:length(newID)){
						NewTree$DataPoints[[newID[j]]]= Tree$DataPoints[[bros[j]]]
					}
					
					## Now we modify Xbase and silos
					for (i in 1:length(NewTree$id)){
						if (NewTree$parent[i]==0){
						Bros=which(NewTree$parent==0)	
						rank=which(Bros==i)
						NewTree$silo[[i]]=siloF(c(0,1), length(Bros), rank)	
						NewTree$Xbase[i]=sum(NewTree$silo[[i]])/2
						} else{
						Bros=which(NewTree$parent==NewTree$parent[i])	
						rank=which(Bros==i)
						NewTree$silo[[i]]=siloF(NewTree$silo[[NewTree$parent[i]]], length(Bros), rank)	
						NewTree$Xbase[i]=sum(NewTree$silo[[i]])/2			
						}
					}
					
					Tree=NewTree	
					break
				} 		
			}
			
			if (i==rev(uniqueParentsNo0)[1]) goWhile=FALSE
					 
		}
	}
	return(Tree)
}
