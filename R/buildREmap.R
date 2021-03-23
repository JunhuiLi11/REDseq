buildREmap <- function (REpatternFilePath, format = "fasta", BSgenomeName, outfile, chr="all") 
{
	if (missing(REpatternFilePath)) {
		stop("missing required parameter REpatternFilePath!")
	}
	if (!file.exists(REpatternFilePath)) {
		stop("REpatternFilePath specified as ", REpatternFilePath, " does not exsist!")
	}
	if (format != "fasta" && format != "fastq") {
		stop("format needs to be either fasta or fastq!")
	}
	if (missing(BSgenomeName) || class(BSgenomeName) != "BSgenome") {
		stop("BSgenomeName is required as BSgenome object!")
	}
	if (file.exists(outfile)) {
		outfile = paste("REDseq_buildREmap", outfile, sep = "_")
	#stop("outfile specified as ", outfile, " already exists! Please rename the outfile!")
	}
	dict = readDNAStringSet(REpatternFilePath, format, use.names = TRUE)
	
	seqnames <- seqnames(BSgenomeName)
	if(length(chr) == 1 && chr=="all"){
		chr <- seqnames
	}else{
		if(!all(chr %in% seqnames)){
			stop("chr should be pattern with the chromosomes of given species or equals to 'all', please check it!")
		}
	}
	REmap <- list()
	for(seqname in chr){
		cat(date(), seqname, " start ...\n")
		searchPattern(dict, BSgenomeName = BSgenomeName, outfile = outfile, chr=seqname)
		outfile <- paste0(outfile,"_",seqname)
		if(file.exists(outfile)){
			GRangesMap <- toGRanges(outfile, format = "BED", header = TRUE, sep = "\t")
			### BED file is 0 based. However, the outfile is 1 based
			GenomicRanges::start(GRangesMap) <- GenomicRanges::start(GRangesMap) - 1
			REmap[seqname] <- list(GRangesMap)
			file.remove(outfile)
		}
	}
	return(REmap)
}
