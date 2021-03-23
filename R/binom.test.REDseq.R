binom.test.REDseq <-
function(REsummary,col.count =2, multiAdj = TRUE, multiAdjMethod="BH", prior.p=0.000001)
{
	REsummary[,col.count] = round(as.numeric(REsummary[,col.count])-0.01)
	total = sum(as.numeric(REsummary[,col.count]))
	count = unique(as.numeric(REsummary[,col.count]))
	pvalue = do.call(rbind, lapply(count, function(x) {
		temp = binom.test(x,total, p=prior.p)
		c(x,temp$p.value, temp$est)
	}))
	group = colnames(REsummary)[col.count]
	colnames(pvalue) = c(paste(group,"count", sep="."), "p.value", "cut.frequency")	
	colnames(REsummary)[col.count] = paste(group,"count", sep=".")
	REsummary = merge(REsummary, pvalue)

	if (multiAdj)
	{
		  procs = c(multiAdjMethod)
        res <- mt.rawp2adjp(REsummary[,dim(REsummary)[2]-1], procs)
        adjp = unique(res$adjp)
        colnames(adjp)[1] = "p.value"
        colnames(adjp)[2] = paste(multiAdjMethod,"adjusted.p.value", sep = ".")
        merge(REsummary, adjp, all.x = TRUE)
		}
	else
	{
			REsummary
	}
}

