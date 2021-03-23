distanceHistSeq2RE <-
function(assignedSeqs,longestDist=1000,title="histogram of distance to assigned RE site", xlab="Distance to assigned RE site",ylab="Frequency",ylim="")
{
		temp = assignedSeqs$passed.filter
		r = hist(temp$Distance[abs(temp$Distance)<longestDist],main=title, xlab=xlab,ylab=ylab)
		d = unique(temp$Distance)
		temp1=cbind(d[order(d)],as.vector(rowsum(temp$Weight, group=temp$Distance)))
		breaks = r$breaks
		for (i in 2:length(breaks))
		{
			r$counts[i-1] = sum(temp1[temp1[,1]<=breaks[i] & temp1[,1] >breaks[i-1],2])
		}
		if (missing(ylim) || length(ylim) !=2)
		{
			plot(r,main=title, xlab=xlab, ylab=ylab, ylim=c(0,max(r$counts)+5))
		}
		else
		{
			plot(r,main=title, xlab=xlab, ylab=ylab, ylim=ylim)
		}
}

