#load("//icnas4.cc.ic.ac.uk/gp111/downloads/usaAndTropFig3trees.Rdata")


charSum=function(labels) {
	nLabels=length(labels)
	
	for (i in 1:nLabels) {
		if (substring(labels[i],1,1)!="+"&substring(labels[i],1,1)!="-") {
			labels[i]=paste("+",labels[i],sep="")
		}
	}
	
	labels=labels[sort(nchar(labels),index.return=TRUE)$ix] #sort by length

	N=nchar(labels[nLabels])-1
	
	digits=matrix(ncol=nLabels,nrow=N)
	
	for (i in 1:nLabels) {
		charDigits=strsplit(labels[i],split="")[[1]]
		sign=(charDigits[1]=="+")-(charDigits[1]=="-")
		digitsTemp=as.numeric(charDigits[2:length(charDigits)])
		digits[,i]=sign*c(rep(0,N-length(digitsTemp)),digitsTemp)
	}

	sumTemp=digits%*%rep(1,nLabels)
	
	sum=digitise(sumTemp)

	charSum=paste(as.character(sum),collapse="")
	
	return(charSum)
}

digitise=function(array) {
	sign=sign(array[1])
	
	l=length(array)
	
	tens=sign(array)*floor(abs(array)/10)
	units=array-tens*10
	
	digitisedTemp=c(tens,0)+c(0,units)
	
	while(any(abs(digitisedTemp)>=10)) {
		tens=(digitisedTemp>0)*floor(abs(digitisedTemp)/10)
		units=digitisedTemp-tens*10
		
		if (tens[1]>=1) {
			digitisedTemp=c(tens,0)+c(0,units)
		} else {
			digitisedTemp=c(tens[2:(l+1)],0)+units
		}
	}
	
	digitisedTemp=discardFirstZeros(digitisedTemp)

	while (sum(abs(digitisedTemp[sign(digitisedTemp)!=sign]))>0) {
		toChange=which(sign(digitisedTemp)!=sign)
		check=abs(sign(digitisedTemp[which(sign(digitisedTemp)!=sign)]))
		digitisedTemp[toChange]=(10*sign+digitisedTemp[toChange])*check
		digitisedTemp[toChange-1]=digitisedTemp[toChange-1]-1*check*sign
	}
	
	digitisedTemp=discardFirstZeros(digitisedTemp)

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

#
#charDiv2 <- function(number) {
#	N=nchar(number)
#	digits=as.numeric(strsplit(number,split="")[[1]])
#	if (digits[N]%%2 !=0 ){stop("number not even")}
#
#	ratio=floor(digits/2)
#	reminder=digits%%2
#
#	reminderToAdd=c(0,reminder[1:(N-1)]*5)
#	
#	divTemp=ratio+reminderToAdd
#	
#	charDiv2=digitise(divTemp)
#
#	return(charDiv2)
#}
#


#This divides by 2 an even char number
charDiv2<- function(number) {
	N=nchar(number)
	digits=as.numeric(strsplit(number,split="")[[1]])
	if (digits[N]%%2 !=0 ){stop("number not even")}
	charDiv2=""
	remainder=0
	for (i in 1:N) {
		toBeDivided=10*remainder+digits[i]
		charDiv2=paste(charDiv2,as.character(floor(toBeDivided/2)),sep="")
		remainder=toBeDivided%%2
	}
	nCharDiv2=nchar(charDiv2)
	digitsCharDiv2=strsplit(charDiv2,"")[[1]]
	if (nCharDiv2>1 & digitsCharDiv2[1]=="0") {charDiv2=paste(digitsCharDiv2[2:nCharDiv2],collapse="")}
	return(charDiv2)
}

#This finds min and max between char numbers
charMinMax<- function(twolabels) {
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
