compareREDseq <-
function(REsummary,col.count1=2, col.count2=3, multiAdj = TRUE, multiAdjMethod="BH", maxP=1, minCount=1)
{
	if(missing(REsummary))
	{
		stop("REsummary is required!")
	}
	REsummary[,col.count1] = round(as.numeric(REsummary[,col.count1])-0.01)
	REsummary[,col.count2] = round(as.numeric(REsummary[,col.count2])-0.01)
	total1 = sum(as.numeric(REsummary[,col.count1]))
	total2 = sum(as.numeric(REsummary[,col.count2]))
	group1 = colnames(REsummary)[col.count1]
	group2 = colnames(REsummary)[col.count2]
	REsummary = REsummary[REsummary[,col.count1] >= minCount | REsummary[,col.count2] >= minCount,]
	REsummary = cbind(REsummary, c(total1), c(total2))
	colnames(REsummary)[dim(REsummary)[2]] =paste(group2,"total",sep=".")
	colnames(REsummary)[dim(REsummary)[2]-1] =paste(group1,"total", sep=".")
	unique.pair = unique(cbind(as.numeric(REsummary[,col.count1]), as.numeric(REsummary[,col.count2])))
	colnames(unique.pair) = c(colnames(REsummary)[col.count1],colnames(REsummary)[col.count2])
	pvalue =do.call(rbind,lapply(1:dim(unique.pair)[1], function(i) {
		temp = matrix(c(unique.pair[i,2], unique.pair[i,1], total2 - unique.pair[i,2], total1-unique.pair[i,1]),
      	 		nrow = 2,  dimnames =  list(c(group1, group2),  c("cut", "not cut")))
		r1 = fisher.test(temp)
		c(unique.pair[i,],r1$p.value, r1$est)
	}))
	colnames(pvalue) = c(paste(group1, "count", sep="."),paste(group2, "count",sep= "."),"p.value", "odds.ratio")
	colnames(REsummary)[col.count1] = paste(group1, "count", sep=".")
	colnames(REsummary)[col.count2] = paste(group2, "count", sep=".")
	
	REsummary = merge(REsummary, pvalue)	
	if (multiAdj)
	{
				procs = c(multiAdjMethod)
				ptemp = REsummary[,dim(REsummary)[2]-1]
        			res <- mt.rawp2adjp(ptemp, procs)
        			adjp = unique(res$adjp)
        			colnames(adjp)[1] = "p.value"
        			colnames(adjp)[2] = paste(multiAdjMethod, "adjusted.p.value", sep = ".")
        			temp = merge(REsummary, adjp, all.x = TRUE)
				temp[is.na(temp[, dim(temp)[2]]),dim(temp)[2]] = 1
				temp[temp[,dim(temp)[2]] <=maxP,]
	}
	else
	{
				REsummary[REsummary[,dim(REsummary)[2]-1] <=maxP,]
	}
}

