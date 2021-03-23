summarizeByRE <-
function(assignedSeqs, by=c("Weight", "REid"),sampleName="",round=TRUE)
{
	if (missing(assignedSeqs)) {
        	stop("Missing requirement argument assignedSeqs!")
	}
	by = match.arg(by)
	temp1 = assignedSeqs$passed.filter
	if (by=="Weight" || by =="W")
	{
			temp = rowsum(temp1$Weight, group=temp1$REid)
			temp = cbind(rownames(temp), unlist(temp[,1]))
			if (round == TRUE)
			{
				temp[,2] = round(as.numeric(temp[,2])-0.001)
			}
			if (sampleName == "")
			{
				colnames(temp) = c("REid", "total.weight")
			}
			else
			{
				colnames(temp) = c("REid", sampleName)
			}
			temp
	}
	else
	{
		temp = as.data.frame(table(temp1$REid))
		if (sampleName == "")
			{
				colnames(temp) = c("REid", "total.count")
			}
			else
			{
				colnames(temp) = c("REid", sampleName)
			}
		temp
	}
}

