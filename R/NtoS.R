##################################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@****************************************************************************@@#
#@@**       0                                                        0       **@@#
#@@**      0|0               %--%--%--%--%--%--%--%--%              0|0      **@@#
#@@**     0 | 0              |                       |             0 | 0     **@@#
#@@**    0  |  0             %     A L G E B R A     %            0  |  0    **@@#
#@@**   0   |   0            |                       |           0   |   0   **@@#
#@@**  0----+----0           %          O F          %          0----+----0  **@@#
#@@**   0   |   0            |                       |           0   |   0   **@@#
#@@**    0  |  0             %       T R E E S       %            0  |  0    **@@#
#@@**     0 | 0              |                       |             0 | 0     **@@#
#@@**      0|0               %--%--%--%--%--%--%--%--%              0|0      **@@#
#@@**       0                                                        0       **@@#
#@@****************************************************************************@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
##################################################################################


##################################################################################
##################################################################################
####     %--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%     ####
####     |                     F U N C T I O N S                        |     ####
####     %--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%     ####
##################################################################################
##################################################################################

#Takes a matrix (or an array) with two columns (or two elements). 
#Each row is a pair of subtree indices. Remember the first element is larger than the second
#Returns the index of the tree originated joining the subtrees
T.from.subTs=function(T1T2) {
	
	if (class(T1T2)=="numeric") {
		T1=T1T2[1]
		T2=T1T2[2]
	} else {
	T1=T1T2[,1]
	T2=T1T2[,2]
	}
	
	T.from.subTs=T1*(T1-1)/2+T2+1
	return(T.from.subTs)
}


#This is a bijection : N --> Z. It is the inverse of the one below
N.to.Z=function(n) {
	if (min(n)>=0) {
		N.to.Z=n/2*((n+1)%%2) - (n+1)/2*(n%%2)
		return(N.to.Z)
	}
	else {stop(sprintf("not a natural or 0", n))}
}


#This is a bijection Z --> N. It is the inverse of the one above
Z.to.N=function(z) {
	Z.to.N=2*z*(z>=0) -(2*z+1)*(z<0)
	return(Z.to.N)
}


#This is a bijection N --> Q+. It is the inverse of the one below.
factors.and.exponents=function(N) {
	N.factors=as.numeric(factorize(N))
	N.unique.factors=unique(N.factors)
	N.factors.exponents=c()
	for (i in 1:length(N.unique.factors)) {
		N.factors.exponents[i]=sum(N.factors==N.unique.factors[i])
	}
	return(list(factors=N.unique.factors, exponents=N.factors.exponents))
}
N.to.Qp=function(n) {
	fAndE=factors.and.exponents(n+1)
	n1.unique.factors=fAndE$factors
	n1.exponents=fAndE$exponents
	newExponents=N.to.Z(n1.exponents)
	N.to.Qp=prod(n1.unique.factors^newExponents)
	num=prod(n1.unique.factors^(newExponents*(newExponents>0)))
	den=prod(n1.unique.factors^(-newExponents*(newExponents<0)))
	return(list(q=N.to.Qp,num=num,den=den))
}


#This is a bijection Q+ --> N. It is the inverse of the one above.
Qp.to.N=function(q, num=NA, den=NA) {
	if (is.na(num) & is.na(den)) {
		fracString=attr(fractions(q, cycle=1000, max.denominator=999999999999),"fracs")
		num=as.numeric(strsplit(fracString, "/")[[1]][1])
		den=as.numeric(strsplit(fracString, "/")[[1]][2])
	}

	if (is.na(den)) {
		den=1
	}
	
	G=GCD(num,den)
	if (G>1) {
		num=num/G
		den=den/G
	} 
	
	numFandE=factors.and.exponents(num)
	denFandE=factors.and.exponents(den)
	
	Qp.to.N=prod(numFandE$factors^Z.to.N(numFandE$exponents))*prod(denFandE$factors^Z.to.N(-denFandE$exponents))-1		

	return(Qp.to.N)	
}


#This is a bijection N --> Q. 
N.to.Q=function(n) {
	if (n==0) {
		N.to.Q=0
		return(N.to.Q)
	} else {
		N.to.Q=(-1)^n * N.to.Qp(floor((n-1)/2))
		return(N.to.Q)
	}
}


#Takes an array of tree indices and returns the indices of the first subtrees
subTs.from.T=function(T) {
	L=length(T)
	
	T1=ceiling((-1+sqrt(1+8*(T-1)))/2)
	T2=T-T1*(T1-1)/2-1

	subTs.from.T=matrix(c(T1,T2),nrow=L,ncol=2)
	return(subTs.from.T)
}


#Takes an index n. Returs a matrix with two columns
#The first column is a node number, the second column is its ancestor
#the root is the node number who has 0 as ancestor
nthDetails=function(n, anc=0, counter=1) {

	if (n==1) {
		node=counter
		ancestor=anc

		nthDetails=matrix(c(node, ancestor), ncol=2, nrow=1)
	}
	else if (n==2) {
		node=c(counter, counter+1, counter+2)
		ancestor=c(anc, counter, counter)

		nthDetails=matrix(c(node, ancestor), ncol=2, nrow=3)
	}
	else {
		t1t2=subTs.from.T(n)
		t1t2anc=counter
		
		counterT1=counter+1
		t1Tree=nthDetails(t1t2[1], t1t2anc, counterT1)
		
		counterT2=counter+1+nrow(t1Tree)
		t2Tree=nthDetails(t1t2[2], t1t2anc, counterT2)

		ancestor=anc
		node=counter
		currentNode=matrix(c(node, ancestor), nrow=1, ncol=2)

		t1row=nrow(t1Tree)
		t2row=nrow(t2Tree)
		L=t1row+t2row+1
		nthDetails=matrix(nrow=L,ncol=2)
		nthDetails[1,]=currentNode
		nthDetails[2 : (t1row+1),]=t1Tree
		nthDetails[(t1row+2) : L,]=t2Tree	
	}

	return(nthDetails)
}


#Takes a matrix and an array representing element permutation
#Returns the matris whose elements have been permutated following the array
#example: if p=c(1,2,3) then 1 is substituted by 2, 2 by 3 and 3 by 1.
swap=function(p,M) {
	M1=M
	L=length(p)
		
	for (i in 1:(L-1)) {
		M1[M==p[i]]=p[i+1]
	}
	
	M1[M==p[L]]=p[1]

	swap=M1
	return(swap)
}


#Takes a tree index. It returns a similar matrix to nthTree() but with different labelling
#It first create the edge matrix with nthDetails(), then it corrects the labelling
#to be conform to the phylo4 class.
nthPhyloMat=function(n) {
	treeMat=nthDetails(n)
	nodes=treeMat[,1]
	ancestor=treeMat[,2]

	L=length(nodes)
	currentValues=ancestor[duplicated(ancestor)]
	numberInternalNodes=length(currentValues)
	newValues=L:(L-numberInternalNodes+1)
	minimum=L-numberInternalNodes+1

	toAvoid=currentValues[which(currentValues>=minimum)]
	if (length(toAvoid)==0) {toAvoid=1}
	for (i in 1:length(toAvoid)) {
		newValues=newValues[newValues!=toAvoid[i]]
	}

	valuesToChange=currentValues[which(currentValues<minimum)]
	
	nthPhyloMat=matrix(c(ancestor, nodes), nrow=L, ncol=2)
	
	for (i in 1:length(valuesToChange)) {
		nthPhyloMat=swap(c(valuesToChange[i],newValues[i]), nthPhyloMat)
	}

	root=nthPhyloMat[,2][which(nthPhyloMat[,1]==0)]
	nthPhyloMat=swap(c(root,minimum), nthPhyloMat)

#	internalNodes=(L-numberInternalNodes+1):L
#	root=phyloMat[,2][phyloMat[,1]==0]
#	phyloMat=swap(c(root,internalNodes[1]), phyloMat)
	
	S=sort(nthPhyloMat[,2],index.return=TRUE)
	nthPhyloMat[,2]=S$x
	nthPhyloMat[,1]=nthPhyloMat[,1][S$ix]

	return(nthPhyloMat)	
}


#Takes an edge matrix of a tree. It returns the number of generations and theindividuals in each generation
generations <- function(phyloMat) {

	L=length(phyloMat[,1])
	mothers=phyloMat[,1][duplicated(phyloMat[,1])]
	numberInternalNodes=length(mothers)
	
	nodesGeneration=list()
	nodesGeneration[[1]]=1:(L-numberInternalNodes)
	generations=0
	s=1
	while (s>0) {
		generations=generations+1
		dupMothers=phyloMat[,1][nodesGeneration[[generations]]]
		mothers=dupMothers[!duplicated(dupMothers)]
		nodesGeneration[[generations+1]]=mothers[which(mothers!=0)]
		s=length(nodesGeneration[[generations+1]])
	}
	
	nodesGeneration[[generations+1]]=NULL
	return(list(number=generations, nodes=nodesGeneration))
}


#Takes a tree index. It returns the correspondinf tree in class phylo4
#NOTE: no good for plotting, the internal labelling, although it is correct,
#confuses treePlot {phylobase}. For plotting use nthPhylo() instead
nthPhylo4=function(n) {	
	phyloMat=nthPhyloMat(n)	

	nthPhylo4=phylo4(phyloMat, edge.length=rep(1,length(phyloMat[,1])))
	return(nthPhylo4)	
}

#Takes a tree index. It returns the corresponding tree in phylo format

#' Obtain the tree for integer n
#' @param n An integer
#' @return phylo object tree corresponding to n
#' @examples
#' nthPhylo(10)
#' @import phylobase 
#' @import gmp 
#' @import MASS
#' @import numbers
#' @import igraph 
nthPhylo <- function(n) {
edges=nthPhyloMat(n)
root=edges[edges[,1]==0,2]
edgesNoRoot=edges[-root,]
nthPhylo=makelabphylotree(edgesNoRoot,rep(1,nrow(edgesNoRoot)),Root=root)
return(nthPhylo)
}


# this uses igraph (functions graph, graph.dfs) and ape (rtree)
makelabphylotree <- function(Edges, Lengths, Root,FLAGS=NULL) {
G <- graph(edges=t(Edges));
orderT <- graph.dfs(G,Root)$order; 
 newLengths=0*Lengths; 

oE<-order(Edges[,2]); Edges <- Edges[oE, ];  Lengths <-Lengths[oE];

newEdges <- matrix(NA, nrow(Edges), ncol(Edges)) 
ooT<-order(orderT[-1]);
newEdges[ooT,] <- Edges
newLengths[ooT]= Lengths
if (!is.null(FLAGS)) {
    FLAGS <- FLAGS[oE];
    newFlags=0*FLAGS;
    newFlags[ooT]=FLAGS
}


Nnode <- (length(orderT)-1)/2; 
Ntips <- Nnode+1; 
pt <- rtree(Ntips); 
pt$Nnode <- Nnode; 
pt$edge <-  newEdges;
pt$edge.length=newLengths;

tiplabels <- paste("t",1:Ntips,sep="") # tip numbers themselves

if (!is.null(FLAGS)) {
    tipind=which(newEdges[,2]<=Ntips)
alllabels <- rep("b",nrow(newEdges))

   alllabels[tipind] <- paste("t",1:Ntips,sep="") #  the one listed first in edges has the first tip label
    # and so on. NEED TO NOW RE-ORDER THE FLAGS (keep track of tip meta-data; not needed in treetop)
tiplabels <- alllabels[tipind]
tipinternalnums <- newEdges[tipind,2]
myflags=0*(1:Ntips)
 myflags[newEdges[tipind,2]]=newFlags[tipind]
    tiplabels[myflags==1]=paste("ps_",tiplabels[myflags==1],sep="")
    
}
pt$tip.label=tiplabels
return(pt)
}
