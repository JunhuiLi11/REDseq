assignSeq2REsite <- function (input.RD.list, REmap.RD.list, cut.end = 'both',cut.offset = 1, seq.length = 36, allowed.offset = 5, min.FragmentLength = 60, max.FragmentLength = 300, partitionMultipleRE = c("unique", "average", "estimate","best", "random"), chr="all")
{
	cat(date(), "Validating input ...\n")
	partitionMultipleRE = match.arg(partitionMultipleRE)
	if (missing(input.RD.list)) {
		stop("Missing required argument input.RD.list!")
	}
	if (missing(REmap.RD.list)) {
		stop("Missing requirement argument REmap.RD.list!")
	}
	if(cut.end %in% c('3','5','both')){
		if(cut.end == 'both'){
			cutend=c('3','5')
		}else{
			cutend=cut.end
		}
	}else{
		stop("cut.end should be one of '3', '5' and 'both'")
	}
	
	allChr.map <- unique(names(REmap.RD.list))
	if(length(chr)==1 && chr=="all"){
		allChr.input <- unique(names(input.RD.list))
	}else{
		allChr.input <- chr
	}
	
	cat(date(), "Prepare map data ...\n")
	common.chr <- intersect(allChr.map,allChr.input)
	REmap.RD.list <- REmap.RD.list[common.chr]
	r1.keep.all <- NULL
	r.passed.uRE.all <- NULL
	r.notpassed.all <- NULL
	r.passed.mRE.all <- NULL
	for(ce in cutend){
		r1 = lapply(allChr.input, function(chr) {
			if (chr %in% allChr.map) {
				cat(date(), "Align to chromosome ", chr, " ...\n")
				thisInput = input.RD.list[[chr]]
				if (class(thisInput) != "GRanges") {
					stop("No valid input.RD.list passed in. It needs to be GRanges object")
				}
				thisMap = REmap.RD.list[[chr]]
				if (class(thisMap) != "GRanges") {
					stop("No valid REmap.RD.list passed in. It needs to be GRanges")
				}
				this.strand = as.character(strand(thisInput))
				
				if(ce=='3'){
					Input.IR = IRanges(start = end(thisInput), width = 1, names = names(thisInput),strand = this.strand)
					names(Input.IR) <- paste0(names(Input.IR),"_3")
				}else if(ce=='5'){
					Input.IR = IRanges(start = start(thisInput), width = 1, names = names(thisInput),strand = this.strand)
					names(Input.IR) <- paste0(names(Input.IR),"_5")
				}
				
				Map.IR = IRanges(start = start(thisMap), end = end(thisMap), names = names(thisMap))
				
				rm(thisInput)
				rm(thisMap)
				
				#nearestRE = Map.IR[selectHits(nearest(Input.IR, Map.IR, select = "all"),select="last")]
				nearestRE = Map.IR[nearest(Input.IR, Map.IR)]
				
				Distance <- rep(NA,length=length(Input.IR))
				plus.loc <- which(GenomicRanges::mcols(Input.IR)$strand == "+")
				
				Distance[plus.loc] <- start(Input.IR[plus.loc]) - start(nearestRE[plus.loc])
				minus.loc <- which(GenomicRanges::mcols(Input.IR)$strand == "-")
				Distance[minus.loc] <- start(Input.IR[minus.loc]) - end(nearestRE[minus.loc])
				
				data.frame(SEQid = names(Input.IR), REid = names(nearestRE),
									 Chr = c(chr), strand = this.strand, SEQstart = start(Input.IR),
									 SEQend = end(Input.IR), REstart = start(nearestRE),
									 REend = end(nearestRE), Distance = Distance,Weight = c(1))
			}
		})
		
		names(r1) <- common.chr
		cat(date(), "Finished 1st round of aligning! Start the 2nd round of aligning ...\n")
		newMap <- list()
		newInput <- list()
		r1.keep.list <- NULL
		#chr="chr1"
		for(chr in names(r1)){
			sub.r1 <- r1[[chr]]
			r1.keep = sub.r1[(sub.r1$Distance <= (cut.offset + allowed.offset) &
													as.character(sub.r1$strand) %in% c("1", "+") & sub.r1$Distance >= 0) |
												 (sub.r1$Distance >= (cut.offset - seq.length - allowed.offset) &
														as.character(sub.r1$strand) %in% c("-1", "-") & sub.r1$Distance <= 0),]
			r1.keep.list[chr] <- list(r1.keep)
			newMap[chr] = list(unique(data.frame(REid = r1.keep$REid, Chr = r1.keep$Chr,
																					 REstart = r1.keep$REstart, REend = r1.keep$REend)))
			newInput[chr] = list(sub.r1[!sub.r1$SEQid %in% r1.keep$SEQid, ])
		}
		rm(r1)
		allChr.map = unique(names(newMap))
		allChr.input = unique(names(newInput))
		common.chr <- intersect(allChr.map,allChr.input)
		
		#r2 = future_lapply(allChr.input, function(chr) {
		r2 = lapply(allChr.input, function(chr) {
			if (chr %in% allChr.map) {
				cat(date(), "Align to chromosome ", chr, " ... \n")
				thisInput = newInput[[chr]]
				thisMap = newMap[[chr]]
				
				Input.IR = IRanges(start = thisInput$SEQstart, end = thisInput$SEQend, names = thisInput$SEQid,strand=thisInput$strand)
				id.strand = cbind(as.character(thisInput$SEQid), as.character(thisInput$strand))
				colnames(id.strand) = c("SEQid", "strand")
				Map.IR = IRanges(start = thisMap$REstart, end = thisMap$REend, names = thisMap$REid)
				matches = findOverlaps(Input.IR, Map.IR, maxgap = max.FragmentLength, select = "all")
				matchmatrix = as.matrix(matches)
				qname = Input.IR[matchmatrix[, 1]]
				tname = Map.IR[matchmatrix[, 2]]
				loc.plus <- which(GenomicRanges::mcols(qname)$strand == "+")
				loc.minus <- which(GenomicRanges::mcols(qname)$strand == "-")
				Distance <- rep(NA,length=length(qname))
				Distance[loc.plus] <- start(qname[loc.plus]) - start(tname[loc.plus])
				Distance[loc.minus] <- start(qname[loc.minus]) - end(tname[loc.minus])
				SEQid = cbind(names(qname))
				colnames(SEQid) = "SEQid"
				this.strand = merge(id.strand, SEQid, by = "SEQid")
				data.frame(SEQid = names(qname), REid = names(tname),
									 Chr = c(chr), strand = this.strand[, 2], SEQstart = start(qname),
									 SEQend = end(qname), REstart = start(tname),
									 REend = end(tname), Distance = Distance, Weight = c(1))
			}
		})
		names(r2) <- common.chr
		cat(date(), "Start filtering ... \n")
		
		r1.keep <- NULL
		r.notpassed <- NULL
		r.passed.mRE <- NULL
		r.passed.uRE <- NULL
		#chr="chrX"
		for(chr in names(r2)){
			sub.r2 <- r2[[chr]]
			cat(date(), "Filter chromosome ", chr, " ... \n")
			r.passed = sub.r2[(sub.r2$Distance <= 0 & abs(sub.r2$Distance) <= max.FragmentLength &
													 abs(sub.r2$Distance) >= min.FragmentLength) | (sub.r2$Distance >= 0 &
																																						abs(sub.r2$Distance) <= (max.FragmentLength - seq.length) &
																																						abs(sub.r2$Distance) >= (min.FragmentLength - seq.length)), ]
			
			if(nrow(r.passed)==0){
				r.notpassed.sub = sub.r2
				r.notpassed <- rbind(r.notpassed,r.notpassed.sub)
			}else{
				r.notpassed.sub = sub.r2[!sub.r2$SEQid %in% r.passed$SEQid, ]
				seqID.count = as.data.frame(table(r.passed$SEQid))
				colnames(seqID.count) = c("SEQid", "count")
				seqID.mRE = as.character(seqID.count[seqID.count$count > 1,]$SEQid)
				r.passed.mRE.sub = r.passed[as.character(r.passed$SEQid) %in% seqID.mRE, ]
				r.passed.uRE.sub = r.passed[!as.character(r.passed$SEQid) %in% seqID.mRE, ]
				r.notpassed <- rbind(r.notpassed,r.notpassed.sub)
				r.passed.mRE <- rbind(r.passed.mRE,r.passed.mRE.sub)
				r.passed.uRE <- rbind(r.passed.uRE,r.passed.uRE.sub)
				rm(seqID.mRE)
			}
			r1.keep <- rbind(r1.keep,r1.keep.list[[chr]])
		}
		r1.keep.all <- rbind(r1.keep.all,r1.keep)
		r.passed.uRE.all <- rbind(r.passed.uRE.all,r.passed.uRE)
		r.notpassed.all <- rbind(r.notpassed.all,r.notpassed)
		r.passed.mRE.all <- rbind(r.passed.mRE.all,r.passed.mRE)
	}
	
	
	r.unique <- rbind(r1.keep.all, r.passed.uRE.all)
	if (partitionMultipleRE == "unique") {
		list(passed.filter = r.unique, notpassed.filter = rbind(r.notpassed.all,r.passed.mRE.all))
	}else {
		cat(date(), "Partitioning reads over RE sites within ", max.FragmentLength, "...\n")
		if (partitionMultipleRE == "average" || partitionMultipleRE == "random") {
			nMapped = rowsum(r.passed.mRE.all$Weight, group = r.passed.mRE.all$SEQid)
			nMapped = cbind(rownames(nMapped), nMapped[, 1])
			colnames(nMapped) = c("SEQid", "nMapped")
			r.passed.mRE.all.1 = merge(r.passed.mRE.all, nMapped, by = "SEQid")
			if (partitionMultipleRE == "average") {
				r.passed.mRE.all.1$Weight = 1/as.numeric(as.character(r.passed.mRE.all.1$nMapped))
				list(passed.filter = rbind(r.unique, r.passed.mRE.all.1[,1:10]), notpassed.filter = r.notpassed.all, mREwithDetail = r.passed.mRE.all.1)
			}else {
				r.passed.mRE.all.1$Weight = 0
				i = 1
				while (i < dim(r.passed.mRE.all.1)[1]) {
					j = sample(seq(1, as.numeric(as.character(r.passed.mRE.all.1$nMapped[i]))),1)
					r.passed.mRE.all.1$Weight[i + j - 1] = 1
					i = i + as.numeric(as.character(r.passed.mRE.all.1$nMapped[i]))
				}
				r.passed.mRE.all.1 = r.passed.mRE.all.1[r.passed.mRE.all.1$Weight == 1, ]
				list(passed.filter = rbind(r.unique, r.passed.mRE.all.1[,1:10]), notpassed.filter = r.notpassed.all, mREwithDetail = r.passed.mRE.all.1)
			}
		}else {
			cat(date(), "get count for each RE ...\n")
			reID.count = as.data.frame(table(r1.keep.all$REid))
			colnames(reID.count) = c("REid", "count")
			r.passed.mRE.all.1 = merge(as.data.frame(r.passed.mRE.all),as.data.frame(reID.count), by = "REid")
			temp = rowsum(r.passed.mRE.all.1$count, group = r.passed.mRE.all.1$SEQid)
			temp = cbind(rownames(temp), temp[, 1])
			colnames(temp) = c("SEQid", "total.count")
			r.passed.mRE.all = merge(r.passed.mRE.all.1, temp, by = "SEQid")
			r.passed.mRE.all$Weight = r.passed.mRE.all$count/as.numeric(as.character(r.passed.mRE.all$total.count))
			if (partitionMultipleRE == "estimate") {
				list(passed.filter = rbind(r.unique, r.passed.mRE.all[,1:10]), notpassed.filter = r.notpassed.all, mREwithDetail = r.passed.mRE.all)
			}else {
				nMapped = as.data.frame(table(r.passed.mRE.all$SEQid))
				colnames(nMapped) = c("SEQid", "nMapped")
				r.passed.mRE.all = merge(r.passed.mRE.all, nMapped, by = "SEQid")
				r.passed.mRE.all$Weight = 0
				i = 1
				while (i < dim(r.passed.mRE.all)[1]) {
					r.passed.mRE.all$Weight[r.passed.mRE.all$count == max(r.passed.mRE.all$count[i:(i + as.numeric(as.character(r.passed.mRE.all$nMapped[i])) - 1)]) & r.passed.mRE.all$SEQid == r.passed.mRE.all$SEQid[i]] = 1
					i = i + as.numeric(as.character(r.passed.mRE.all$nMapped[i]))
				}
				r.passed.mRE.all = r.passed.mRE.all[r.passed.mRE.all$Weight ==1, ]
				list(passed.filter = rbind(r.unique, r.passed.mRE.all[,1:10]), notpassed.filter = r.notpassed.all, mREwithDetail = r.passed.mRE.all)
			}
		}
	}
}
