
#####################################################################
##  functions needed for my unlabelled tree metric
#####################################################################
# require(digest)

# note plotting the tree with show.node.label=TRUE allows you to check that the nodes are labelled the right way. 
# could consider comparing: orig trees, trees with all tips retracted (towards backbone), trees retracted again, etc. 
# could consider weightings but might get to violating triangle ineq

# new with hashing because labels get too big too fast 


#' Create the labels for a tree
#'
#' @param tree A tree in phylo format (or convertible with as(tree,"phylo")
#' @param allChar Whether to use all character precision (default is false)
#' @return A vector of labels for the tree including tips and internal nodes
#' @examples
#' treelabels(rtree(10)
treelabels <- function(tree,allChar=FALSE) { 
     if (class(tree) != "phylo") tree = as(tree,"phylo")
	 if (is.null(tree$tip.label))
        stop("This tree has no tips")
	num.tips=length(tree$tip.label)
	labels=NA + 0*(1:(nrow(tree$edge)+1))
	names(labels)[1:num.tips]=tree$tip.label;
	names(labels)[(num.tips+1): length(labels)]=paste("node",1:(length(labels)-num.tips),sep="")
	labels[1:num.tips]=ifelse(allChar=TRUE,"1",1)     # tips are 0 for now 
	NodeIDS= (num.tips + 1) : (2*num.tips -1)
	while (any(is.na(labels))) { 
		IsReady = NodeIDS[ vapply(NodeIDS,function(x) !any(is.na(labels[tree$edge[which(tree$edge[,1]==x),2]])) & is.na(labels[x])  ,FUN.VALUE=TRUE) ]
		TheseLabels = unlist(sapply(IsReady, function(x) getLabelFromPairs(labels[tree$edge[tree$edge[,1]==x,2]],useAllCharacters=allChar)))
		labels[IsReady]=TheseLabels
		}
		return(labels)
	}

#' Plot a tree showing the labels 
#' @param tree A tree in phylo format
#'  @return No variable is returned
#'  @examples
#'  plotlabels(rtree(10))
plotlabels <- function(tree) {
    nn=length(tree$tip.label)
tree$node.label=treelabels(tree)[(nn+1):(2*nn-1)]
tree$tip.label=1+0*(1:nn)
    plot(tree,show.node.label=TRUE, edge.width=8,cex=2,edge.color="grey")
}


getLabelFromPairs <-function( twolabels,useAllCharacters=FALSE ) {
if (useAllCharacters==TRUE) {
    minMax=charMinMax(twolabels) 
	k = minMax$max
	j = minMax$min
	newLabel=charSum(c(charDiv2(charProd(c(k, charSum(c(k,"-1"))))),charSum(c(j,"1")))) #this is k*(k+1)/2 + j + 1 (include the 1 if 0 is a tip)
      return(newLabel)
} else {
    l1=twolabels[1]; l2=twolabels[2]; 
        if (nchar(l1) < 14 & nchar(l2) < 14) { 
                k = max(as.numeric(l1),as.numeric(l2))
                j = min(as.numeric(l1),as.numeric(l2))
                return( k*(k-1)/2 + j + 1) # NOTE if 1 is a tip and there are no 0s allowed (full binary tree) the correct expression is 1/2 k (k-1) + j + 1. 
                } else { return(digest(sort(c(l1,l2)))) } 
}

}


#this multiplies a char number for a single digit char
charProdSimple <- function(twolabels) {
	twolabels=twolabels[sort(nchar(twolabels),index.return=TRUE)$ix] #sort by length
	lshort=twolabels[1]; llong=twolabels[2];
	N=nchar(llong)
	if (nchar(lshort)!=1) {stop("not single digit")}
	digitsLong=as.numeric(strsplit(llong,split="")[[1]])
	dS=as.numeric(lshort)
	remainder=0
	charProdSimple=""
	for (i in N:1) {
		dL=digitsLong[i]
		dProd=as.character(dL*dS+remainder)
		if (nchar(dProd)==2) {
			dProdSplit=strsplit(dProd,"")[[1]]
			remainder=as.numeric(dProdSplit[1])
			charProdSimple=paste(as.character(dProdSplit[2]),charProdSimple,sep="")
		} else {
			charProdSimple=paste(dProd,charProdSimple,sep="")
			remainder=0
		}
	}
	if (remainder!=0) {charProdSimple=paste(as.character(remainder),charProdSimple,sep="")}
	return(charProdSimple)
}


#' Get distance between two vectors of labels
#' @param x A vector of tree labels (as created by treelabels)
#' @param y A vector of tree labels (as created by treelabels)
#' @return D1 distance: the symmetric set difference between x and y
#' @examples
#' labeldistance(treelabels(rtree(10)), treelabels(rtree(12))
labeldistance <- function(x,y) {
	# for each unique element of x, how many times does it come up in x, and how many in y?
	uni.x=unique(x)
	ynotinx=setdiff(y,x) # things in y not in x 
	dcounts=vapply(c(uni.x,ynotinx), function(k) abs(length(which(y==k))-length(which(x==k))),FUN.VALUE=1)
	return(sum(dcounts))
	# for each unique element of y NOT already counted in x, how many times does it happen in y?
	# distance is the sum of all of those numbers.
	 
}

#' Get distance between two trees
#' @param tree1 A phylogenetic tree, phylo format or convertible 
#' @param tree2 A phylogenetic tree, phylo format or convertible 
#' @return labeldistance between the label sets of the two trees
#' @examples
#' distunlab(rtree(10), rtree(12))
distunlab<-function(tree1,tree2) {
	lab1=treelabels(tree1); lab2=treelabels(tree2);
	return(labeldistance(lab1,lab2))
}

##4 veclabeldistance: instead of the symmetric set difference use the L2 norm of vectors of the #s of each individual unique label. 

# samesize = FALSE: divide all entries of hte label count by n, and add a compensatory espsilon | na -nb| term to the distance. if nottips =TRUE discount tip # entirely. 

#' Get D2 distance between two vectors of labels
#' @param lab1 A vector of tree labels (as created by treelabels)
#' @param lab2 A vector of tree labels (as created by treelabels)
#' @param samesize Logical: are the trees the same size? 
#' @param eps: if the trees are not the same size (samesize=F), 
#' the components in the distance will be divided by the number of tips. 
#' Then eps*(difference in tip number) is added to preserve the metric property. 
#' @param nottips: exclude the number of 1s from the distance (default F)
#' @return D2 distance: the L2 norm comparing the vectors of labels
#' @examples
#' veclabeldistance(treelabels(rtree(10)), treelabels(rtree(12))
veclabeldistance <- function(lab1, lab2, samesize=TRUE,eps=0,nottips=FALSE) {
	Unis=unique(c(lab1,lab2)); 
	if (nottips) {Unis = Unis[-which(Unis==min(Unis))]}
	ntips1=sum(lab1==min(lab1)); ntips2=sum(lab2==min(lab2)); 
	if (samesize==TRUE) {
	components= vapply(Unis, function(x) abs(sum(lab1==x)-sum(lab2==x)), FUN.VALUE=0)
	} 
	if (samesize==FALSE) {
		components= vapply(Unis, function(x) abs(sum(lab1==x)/ntips1-sum(lab2==x)/ntips2), FUN.VALUE=0)
		}
	return(sqrt(sum(components^2))+eps*abs(ntips1-ntips2))
	}
	
	
veclabel<- function(lab,Nmax=50) {
	return(vapply(1:Nmax, function(x) sum(lab==x),FUN.VALUE=1))
} # NOTE this does not include the number of tips. 

## 4 all pairwise distances between a list of trees 

#' All pairwise D1 distances for a list of trees
#' @param trees A list of trees
#' @param listoflabels Alternatively, a list of labels. 
#' If listoflabels is NA (the default), this function
#'  first computes the labels for the set of trees. 
#' If listoflabels is not NA, the labels are compared directly and the trees are ignored.
#' @return An object of class 'dist' containing tree-tree D1 distances
#' @examples
#' dd=multiDistUnlab(rmtree(10,12))
multiDistUnlab <- function( trees,listoflabels=NA) {
	 if (is.na(listoflabels)) {	num_trees <- length(trees)} else {num_trees <- length(listoflabels) }

    if (num_trees < 2) {
        stop("multiDistUnlab expects at least two trees")
    }
   # if (!is.na(trees) & is.null(names(trees))) 
   #     names(trees) <- 1:num_trees
   # else if (length(unique(names(trees))) != num_trees) {
   #     warning("duplicates detected in tree labels - using generic names")
   #     names(trees) <- 1:num_trees
   #  }
#     lab <- names(trees)
	distances <- matrix(0, num_trees, num_trees)
if (is.na(listoflabels)) {   listoflabels <- lapply(trees, treelabels) }

	# listoflabels=lapply(trees, treelabels) 
	
	sapply(1:(num_trees - 1), function(i) {
               print(paste("Computing remaining distances to tree ", i,sep=""))
                sapply((i + 1):num_trees, function(j) {
                  distances[i, j] <<- distances[j, i] <<- labeldistance(listoflabels[[i]],listoflabels[[j]])
                })
            })
    return(as.dist(distances))
	
	}
	
	
#' All pairwise D2 distances for a list of trees
#' @param trees A list of trees
#' @param listoflabels Alternatively, a list of labels. If listoflabels is NA
#'  (the default), this function first computes the labels for the set of trees.
#'  If listoflabels is not NA, the labels are compared directly and the trees are ignored.
#' @param samesize Logical: whether the trees are the same size
#' @param eps: if the trees are not the same size (samesize=F), the components 
#' in the distance will be divided by the number of tips.
#' Then eps*(difference in tip number) is added to preserve the metric property. 
#' @param nottips: exclude the number of 1s from the distance (default F)
#' @return An object of class 'dist' containing tree-tree D2 distances (uses veclabeldistance) 
#' @examples
#' dd=vecmultiDistUnlab(rmtree(10,12))
vecMultiDistUnlab <- function(trees,listoflabels=NA,samesize=TRUE,eps=0,nottips=FALSE) {
if (is.na(listoflabels)) {	num_trees <- length(trees)} else {num_trees <- length(listoflabels) }
if (is.na(listoflabels)) {  print("creating labels"); listoflabels <- lapply(trees, treelabels) }
	distances <- matrix(0, num_trees, num_trees)
	sapply(1:(num_trees - 1), function(i) {
                print(paste("Computing remaining distances to tree ", i,sep=""))
                sapply((i + 1):num_trees, function(j) {
                  distances[i, j] <<- distances[j, i] <<- veclabeldistance(listoflabels[[i]],listoflabels[[j]],samesize=samesize,eps=eps,nottips=nottips)
                })
            })
    return(as.dist(distances))
	
	}
	
	
	
getns <- function(trees) { vapply(trees, function(tree) length(tree$tip.label), FUN.VALUE=1)} 

#' Prune a tree to the desired size, randomly uniformly dropping tips
#' @param tree A tree in phylo format
#' @param size an integer, the size of the final tree
#' @return a tree in phylo format, with the desired number of tips (size)
#' @examples
#' prunetosize(rtree(20),12)
prunetosize <- function(tree,size) {
	nt=length(tree$tip.label)
	if (nt <= size) {return(tree)} else {
		todrop =sample(1:nt, nt-size); 
		return(drop.tip(tree, todrop))
		}
	}
	
	################### WEIGHT THE LOWER LABELS HIGHER THAN THE HIGH ONES: 
	
# use 1/sqrt(label) or tiny number to *weight* contributions to v in vec distance: 	

#' Get weighted D2 distance between two vectors of labels
#' @param lab1 A vector of tree labels (as created by treelabels)
#' @param lab2 A vector of tree labels (as created by treelabels)
#' @param weights A function in the form ifelse(is.numeric(x),VALUE1,VALUE2) 
#' where VALUE1 is the weight applied to labels that have not been hashed 
#' (integers; fewer than 12 digits) and VALUE2 is the value if hashed. 
#' Default is (1/(1+log(x)) and 0.0001
#' @return weighted D2 distance
#' @examples
#' weightlabeldistance(treelabels(rtree(10)), treelabels(rtree(12))
	weightlabeldistance <- function(lab1, lab2, weights=function(x) ifelse(is.numeric(x),1/(1+log(x)),0.0001)) {
	Unis=unique(c(lab1,lab2)); 
		
	components = vapply(Unis, function(x) weights(x)*abs(sum(lab1==x)-sum(lab2==x)),FUN.VALUE=0)
		
	return(sqrt(sum(components^2)))
	}

	
#' All weighted pairwise distances for a list of trees
#' @param trees A list of trees
#' @param listoflabels Alternatively, a list of labels. If listoflabels is NA 
#' (the default), this function first computes the labels for the set of trees.
#' If listoflabels is not NA, the labels are compared directly and the trees are ignored.
#' @param weights A function in the form ifelse(is.numeric(x),VALUE1,VALUE2) 
#' where VALUE1 is the weight applied to labels that have not been hashed (fewer
#' than 12 digits) and VALUE2 is the value if hashed. Default is (1/(1+log(x)) and 1e-8
#' @return An object of class 'dist' containing tree-tree distances
#' @examples
#' dd=weightmultiDistUnlab(rmtree(10,12))
weightMultiDistUnlab <- function(trees,listoflabels=NA,weights=function(x) ifelse(is.numeric(x),1/(1+log(x)),1e-8)) {
if (is.na(listoflabels)) {	num_trees <- length(trees)} else {num_trees <- length(listoflabels) }
if (is.na(listoflabels)) {   listoflabels <- lapply(trees, treelabels) }
	distances <- matrix(0, num_trees, num_trees)
	sapply(1:(num_trees - 1), function(i) {
                print(paste("Computing remaining distances to tree ", i,sep=""))
                sapply((i + 1):num_trees, function(j) {
                  distances[i, j] <<- distances[j, i] <<- weightlabeldistance(listoflabels[[i]],listoflabels[[j]],weights)
                })
            })
    return(as.dist(distances))
	
	}
	
	
	
########################### DEPRECATED DO NOT USE, FOR REFERENCE ONLY 
## 1: get internal node labels by the tree isomorphism system
# treelabelsOrig <- function(tree,tips.zero=TRUE) {
#	n=tree$Nnode+1 # number of tips in the tree#
#	IntNodes=(n+1):(2*n-1)
#	Pointers=t(vapply(IntNodes, function(x) tree$edge[tree$edge[,1]==x,2],FUN.VALUE=c(0,0)))	# rownames(Pointers)=IntNodes
#	pidm=matrix(NA,nrow=nrow(Pointers),ncol=ncol(Pointers))
#	pidm[Pointers<=n]=0 # tips are a 1 , or a 0? 
#	nodelabs=0*(1:(n-1)) # initialise 
##	isReady=which(!is.na(pidm[,1]) & !is.na(pidm[,2])) # correspond to rows of pidm#
#	for (Node in isReady) {
#		DecTags=pidm[Node,]; k=max(DecTags); j=min(DecTags);
#		NewLabel=k*(k+1)/2 + j +1 ; # HAD +1 here  for system where 0 is a tip 
#		pidm[Pointers==Node+n]=NewLabel; # or  as.integer(names(Pointers[,1])[Node])
#		nodelabs[Node]=NewLabel;
#	}
#}
#k=max(pidm[1,]);j=min(pidm[1,]); RL=k*(k+1)/2+j +1 ; # root label. HAD +1 here if 0 is a tip
#nodelabs[1]=RL
#idvec=c(pidm[,1],pidm[,2],RL)
#return(list(labels=idvec,nodelabs=nodelabs))
# }



