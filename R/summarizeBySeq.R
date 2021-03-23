summarizeBySeq <-
function(assignedSeqs, by=c("Weight", "SEQid"))
{
	if (missing(assignedSeqs)) {
        	stop("Missing requirement argument assignedSeqs!")
	}
	by = match.arg(by)
	temp2 = assignedSeqs$passed.filter
	if (by=="Weight" || by=="W")
	{
			temp1 = rowsum(temp2$Weight, group=temp2$SEQid)
			temp = cbind(rownames(temp1), unlist(temp1[,1]))
			colnames(temp) = c("SEQid", "total.weight")
		}
	else
	{
		temp = as.data.frame(table(temp2$SEQid))
		colnames(temp) = c("SEQid", "total.count")
	}
	temp
}

