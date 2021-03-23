plotCutDistribution <-
function(assignedSeqs,REmapList, chr="chr1",xlim, title="RE cut frequency distribution", xlab="Chromosome Location (bp)",ylab="Frequency", round=TRUE, n.sequence)
{
	if(missing(assignedSeqs)){
		stop("assignedSeqs is required, which is generated from assignSeq2REsite!")
	}
	
	chr = sub("chr", "", chr, ignore.case =TRUE)
	temp = assignedSeqs$passed.filter		
	temp = temp[as.character(temp$Chr) ==chr,]
	temp1 = rowsum(temp$Weight, group=temp$REid)
	temp2 = cbind(rownames(temp1), unlist(temp1))
	if (round == TRUE){
		temp2[,2] = round(as.numeric(temp2[,2])-0.01)
	}
	colnames(temp2) = c("REid", "total.weight")
	temp1 = merge(temp2, temp[,2:7], by="REid")
	if (missing(REmapList)){
		message("REmapList is not specified, so the map distribution will not draw along with the assigned sequence distribution!")
	}else{	
		ref = REmapList[[chr]]
		
		total = 0
		for(i in length(REmapList)){
			total = total + length(start(REmapList[[names(REmapList)[i]]]))
		}
		
		if (missing(n.sequence)){
			expectedCount =1
		}else{
			expectedCount = n.sequence/total
			expectedCount = max(expectedCount, 1)
			#expectedCount = min(expectedCount, max(as.numeric(temp1$total.weight)))
		}		
	}
	if (missing(xlim) || length(xlim) !=2){
		plot(as.numeric(temp1$REstart), as.numeric(temp1$total.weight), main=title, xlab=xlab, ylab=ylab, type="h")
		if (!missing(REmapList))
		{
			plot(as.numeric(temp1$REstart), as.numeric(temp1$total.weight), main=title, xlab=xlab, ylab=ylab, type="h", ylim=c(0, max(expectedCount, max(as.numeric(temp1$total.weight)))))
			points(start(ref), rep(expectedCount, length(start(ref))), col="red", pch=25,cex=1)
		}else{
			plot(as.numeric(temp1$REstart), as.numeric(temp1$total.weight), main=title, xlab=xlab, ylab=ylab, type="h")
		}
	}else{
		if (!missing(REmapList)){
			plot(as.numeric(temp1$REstart), as.numeric(temp1$total.weight), main=title, xlab=xlab, ylab=ylab, type="h", ylim=c(0, max(expectedCount, max(as.numeric(temp1$total.weight)))))
			points(start(ref), rep(expectedCount, length(start(ref))),  col="red", xlim=xlim, pch=25, cex=1)
		}else{
			plot(as.numeric(temp1$REstart), as.numeric(temp1$total.weight), main=title, xlab=xlab, ylab=ylab, type="h", xlim=xlim)
		}
	}			
}

