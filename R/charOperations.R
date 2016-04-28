#load("//icnas4.cc.ic.ac.uk/gp111/downloads/usaAndTropFig3trees.Rdata")

treelabels <- function(tree,tips.zero=TRUE) { 
     if (class(tree) != "phylo") tree = as(tree,"phylo")
	 if (is.null(tree$tip.label))
        stop("This tree has no tips")
	num.tips=length(tree$tip.label)
	labels=NA + 0*(1:(nrow(tree$edge)+1))
	names(labels)[1:num.tips]=tree$tip.label;
	names(labels)[(num.tips+1): length(labels)]=paste("node",1:(length(labels)-num.tips),sep="")
	labels[1:num.tips]="0"     # tips are 0 for now 
	NodeIDS= (num.tips + 1) : (2*num.tips -1)
	while (any(is.na(labels))) { 
		IsReady = NodeIDS[ vapply(NodeIDS,function(x) !any(is.na(labels[tree$edge[which(tree$edge[,1]==x),2]])) & is.na(labels[x])  ,FUN.VALUE=TRUE) ]
		TheseLabels = unlist(sapply(IsReady, function(x) getLabelFromPairs(labels[tree$edge[tree$edge[,1]==x,2]])))
		labels[IsReady]=TheseLabels
		#print(sum(is.na(labels)))
		}
		return(labels)
	}



getLabelFromPairs <-function( twolabels ) {
	minMax=charMinMax(twolabels) 
	k = minMax$max
	j = minMax$min
	newLabel=charSum(c(charDiv2(charProd(c(k, charSum(c(k,"1"))))),charSum(c(j,"1")))) #this is k*(k+1)/2 + j + 1 (include the 1 if 0 is a tip)
	return(newLabel)
}

charSum <- function(labels) {
	nLabels=length(labels)
	
	labels=labels[sort(nchar(labels),index.return=TRUE)$ix] #sort by length

	N=nchar(labels[nLabels])
	
	digits=matrix(ncol=nLabels,nrow=N)
	for (i in 1:nLabels) {
		digitsTemp=as.numeric(strsplit(labels[i],split="")[[1]])
		digits[,i]=c(rep(0,N-length(digitsTemp)),digitsTemp)
	}

	sumTemp=digits%*%rep(1,nLabels)
	
	sum=digitise(sumTemp)

	charSum=paste(as.character(sum),collapse="")
	
	return(charSum)
}

digitise <- function(array) {
	l=length(array)
	
	tens=floor(array/10)
	units=array-tens*10
	
	digitisedTemp=c(tens,0)+c(0,units)
	
	while(any(digitisedTemp>=10)) {
		tens=floor(digitisedTemp/10)
		units=digitisedTemp-tens*10
		
		if (tens[1]>=1) {
			digitisedTemp=c(tens,0)+c(0,units)
		} else {
			digitisedTemp=c(tens[2:(l+1)],0)+units
		}
	}
	
	while (digitisedTemp[1]==0 & length(digitisedTemp)>1) {
		digitisedTemp=digitisedTemp[2:length(digitisedTemp)]
	}

	digits=paste(as.character(digitisedTemp),collapse="")

	return(digits)
}

#this multiplies two char numbers
charProd <- function(twolabels) {
	nLabels=2				#length(twolabels)
	
	twolabels=twolabels[sort(nchar(twolabels),index.return=TRUE)$ix] #sort by length

	N=nchar(twolabels[2])		#nchar(twolabels[nLabels])
	n=nchar(twolabels[1])
	
	digitsShort=cbind(as.numeric(strsplit(twolabels[1],split="")[[1]]))
	digitsLong=rbind(as.numeric(strsplit(twolabels[2],split="")[[1]]))
	
	if (log10(N)+log10(n)>8) {
		prodTemp=rep(0,N+n-1)
		for (s in 1:(N+n-1)) {
			indShort=max(1,s-N+1):min(n,s)
			indLong=rep(s+1,length(indShort))-indShort
			prodTemp[s]=sum(digitsShort[indShort]*digitsLong[indLong])
		}
	} else {
		prodMat=digitsShort%*%digitsLong
		d=row(prodMat)+col(prodMat)
		prodMatDiags=split(prodMat,d)
	
		prodTemp=c()
		for (i in 1:(N+n-1)) {
			prodTemp[i]=sum(prodMatDiags[[i]])
		}
	}
	
	prodDigits=digitise(prodTemp)
	
	charProd=paste(as.character(prodDigits),collapse="")
	
	return(charProd)
}


charDiv2 <- function(number) {
	N=nchar(number)
	digits=as.numeric(strsplit(number,split="")[[1]])
	if (digits[N]%%2 !=0 ){stop("number not even")}

	ratio=floor(digits/2)
	reminder=digits%%2

	reminderToAdd=c(0,reminder[1:(N-1)]*5)
	
	divTemp=ratio+reminderToAdd
	
	charDiv2=digitise(divTemp)

	return(charDiv2)
}

charMinMax <- function(twolabels) {
	twolabels=twolabels[sort(nchar(twolabels),index.return=TRUE)$ix] #sort by length
	l1=twolabels[1]; l2=twolabels[2];
	N=nchar(l2)
	templ1=paste(paste(rep("0",N-nchar(l1)), collapse=""),l1,sep="")
	digits1=as.numeric(strsplit(templ1,split="")[[1]])
	digits2=as.numeric(strsplit(l2,split="")[[1]])
	diff=digits1-digits2
	test=which(diff != 0)
	if (length(test)==0) {return(list(max=l2,min=l1))}
	else if (diff[test][1]>0) {return(list(max=l1,min=l2))}
	else {return(list(max=l2,min=l1))}
}
