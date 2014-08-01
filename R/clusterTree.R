clusterTree <-
function(X, k, h=NULL, density="knn", dist="euclidean", d=NULL, Nlambda=100, printStatus=FALSE){
	
	if (!is.numeric(X) && !is.data.frame(X)) stop("X should be a matrix of coordinates")
	if (!is.vector(k) || length(k)!=1) stop("k should be an number")
	if (!is.null(h) && (!is.vector(h) || length(h)!=1)) stop("h should be a real value")
	if (!is.null(Nlambda) && (!is.vector(Nlambda) || length(Nlambda)!=1)) stop("Nlambda should be a number")
	if (!is.logical(printStatus)) stop("printStatus should be logical")
	if (density=="kde" && dist!="euclidean") stop("kde is only possible with dist='euclidean' ")
	if (!is.null(d) && (!is.vector(d) || length(d)!=1)) stop("d should be a number")
				
	if (dist=="euclidean"){ 
		n=dim(X)[1]
		d=dim(X)[2]
	} else if (dist=="arbitrary") { 
		n=dim(X)[1]
		if (is.null(d)) d=1
		distMat=X
	} else {
		stop("Wrong dist")
	}
	
	
	## start adjacency matrix
	adjMat=diag(n)

	## Compute density estimator: knn or kde
	if (density=="knn" && dist=="euclidean"){
		knnInfo=get.knn(X, k=k, algorithm="kd_tree")
		for(i in 1:n){
			adjMat[i, knnInfo$nn.index[i,] ]=1
		}
		r.k= apply(knnInfo$nn.dist, 1, max)		
		v.d=pi^(d/2) /gamma(d/2+1)
		hat.f=k/(n*v.d*r.k^d)	
	} else if (density=="knn" && dist=="arbitrary"){
		orderMat=apply(distMat,2,order)[2:(k+1),]
		for(i in 1:n){
			adjMat[i,orderMat[,i]]=1		
		}		
		r.k=apply( apply(X,2,sort)[2:(k+1),], 2, max)
		v.d=pi^(d/2) /gamma(d/2+1)		
		hat.f=k/(n*v.d*r.k^d)	
	} else if (density=="kde"){	#kde estimate
		knnInfo=get.knn(X, k=k, algorithm="kd_tree")
		for(i in 1:n){
			adjMat[i, knnInfo$nn.index[i,] ]=1
		}
		hat.f=kde(X,X,h)	
	} else stop("Wrong density")
	
	# ordered value of the density
	ord.hat.f=order(hat.f)

	# starting graph	
	G=graph.adjacency(adjMat, mode="undirected")	  ## could be changed to directed
	
	# Lambda grid
	Lambda=hat.f[ord.hat.f]
	if (is.null(Nlambda)) { 
		Nlambda=n
	} else {
		Nlambda=min(n,Nlambda)
		Lambda=seq(min(Lambda),max(Lambda), length=Nlambda)
	}

	exclude=numeric()
	
	## in CLUSTERS we store the clusters found for each level of lambda_j
	CLUSTERS=list()
	if (printStatus) cat("Percentage: ")
	for (j in 1:Nlambda){
		OldExcluded=exclude
		lambda=Lambda[j]
		present=which(hat.f>=lambda)  
		exclude=setdiff(1:n,present)  # points with density less than lambda
		NewExcluded=setdiff(exclude,OldExcluded)   # the new excluded point  
		G[NewExcluded, present]=FALSE     # remove edges of the new excluded point
		clust=clusters(G)
		CLUSTERS[[j]]=list("no"=clust$no,"mem"=clust$mem, "present"=present, "exclude"=exclude)		
		if (printStatus) cat(round(j/Nlambda, 3), " ")
	}
	if (printStatus) cat("\n")
	
	## Now assign ID, Generation and Components to each new cluster
	id=0                   # id assigned to each new cluster
	components=list()      # data points contained in each new cluster
	generation=numeric()   #generation of each new cluster
	for (j in 1:Nlambda){
		presentMembership=unique(CLUSTERS[[j]]$mem[CLUSTERS[[j]]$present])
		for (i in presentMembership){		
			id=id+1
			components[[id]]=which(CLUSTERS[[j]]$mem==i)
			generation[id]=j				
		}
	}
	
	# find the father of each cluster
	father=numeric()
	startF=which(generation==2)[1]
	for (i in startF:length(components)){
			for (j in which(generation==(generation[i]-1)) ){
				if (setequal(intersect(components[[i]], components[[j]]), components[[i]])){
						father[i]=j
						break
					}	
			}
	}
	father[is.na(father)]=0    #for the roots set father = 0
	
	## Convert the clusters into BRANCHES
	bb=0                 # count number of branches
	branch=numeric()     # map each cluster into a branch
	base=numeric()       # x coordinate of the base of each branch
	top=numeric()        # y-top of each branch
	bottom=numeric()     # y-bottom of each branch
	compBranch=list()    # data points corresponding to this branch
	silo=list()          # x coordinates of each branch silo
	rank=numeric()       # rank among brothers
	parent=numeric()     # father of each branch
	sons=list()          # sons of each branch 
	
	# if there is more than 1 root create a fake single root
	if (sum(generation==1)>1){
		bb=bb+1
		silo[[bb]]=c(0,1)
		base[bb]=0.5
		compBranch[[bb]]=1:n
		rank[bb]=1
		parent[bb]=0   # parent of roots is set to be 0
		top[bb]=0
		bottom[bb]=0   # y-bottom of roots is set to be 0
	}
	for (i in 1:length(father)){

		# the first generation is treated separately		
		# multiple roots are children of the fake single root
		if (sum(generation==1)>1 & generation[i]==1){ 
			Bros=which(generation==generation[i])	
			bb=bb+1
			branch[i]=bb
			rank[bb]=sum(generation[1:i]==generation[i] & father[1:i]==father[i])  #same gen, same father.
			silo[[bb]]=siloF(c(0,1), length(Bros), rank[bb])	
			base[bb]=sum(silo[[bb]])/2
			top[bb]=min(hat.f[components[[i]]])
			compBranch[[bb]]=components[[i]]		
			parent[bb]=1   # parent of roots is set to be 1
			bottom[bb]=0   # y-bottom of roots is set to be 0
			## add this branch to the list of sons of its parent
			if (length(sons)<parent[bb]) {sons[[ parent[bb] ]]=bb
				} else {sons[[ parent[bb] ]]=c(sons[[ parent[bb] ]],bb)}		
		} else if (sum(generation==1)==1 & generation[i]==1){ #is there is 1 root
			bb=bb+1
			branch[i]=bb
			silo[[bb]]=c(0,1)
			base[bb]=0.5
			top[bb]=min(hat.f[components[[i]]])
			compBranch[[bb]]=components[[i]]		
			parent[bb]=0   # parent of roots is set to be 0
			bottom[bb]=0   # y-bottom of roots is set to be 0
		} else{	
			Bros=which(generation==generation[i] & father==father[i])	
			## if the cluster has brothers, then there is a split and new branches are created
			if (length(Bros)>1){	
				bb=bb+1
				branch[i]=bb
				parent[bb]=branch[father[i]]
				rank[bb]=sum(generation[1:i]==generation[i] & father[1:i]==father[i])	
				silo[[bb]]=siloF(silo[[parent[bb]]], length(Bros), rank[bb])	
				base[bb]=sum(silo[[bb]])/2
				top[bb]=min(hat.f[components[[i]]]) 
				bottom[bb]=top[parent[bb]]
				compBranch[[bb]]=components[[i]]
				## add this branch to the list of sons of its parent
				if (length(sons)<parent[bb]) {sons[[ parent[bb] ]]=bb
					} else {sons[[ parent[bb] ]]=c(sons[[ parent[bb] ]],bb)}			
			}
			
			## if the cluster does not have brothers, then no new branches are created
			## and this cluster is assigned to an old branch	
			if (length(Bros)==1){  
				for (j in which(generation==(generation[i]-1))){
					if (setequal(intersect(components[[i]], components[[j]]), components[[i]]))
					belongTo=branch[j]
				}
				top[belongTo]=min(hat.f[components[[i]]]) #update top of the branch
				branch[i]=belongTo
			}
		}
	}	
	
	## info for the kappa tree
	kTree=findKtree(bb, parent, sons, compBranch, n)
	Ktop=kTree$Ktop
	Kbottom=kTree$Kbottom
	
	out=list("n"=n, "id"=1:bb, "Xbase"=base, "Ybottom"=bottom, "Ytop"=top, "Kbottom"=Kbottom, "Ktop"=Ktop, "silo"=silo, "sons"=sons, "parent"=parent, "DataPoints"=compBranch, "density"=hat.f)
	class(out)="clusterTree"
	
	out1=plotRule(out) ## relabel branches according to plotting rules.
	
	return(out1)
}
