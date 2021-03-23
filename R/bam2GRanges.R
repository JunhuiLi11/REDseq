bam2GRanges <- function (bamfile, bamindex = bamfile, subChr = 'all', pairedEndReads = FALSE, removeDuplicateReads = FALSE, minMapq = 10, maxFragmentWidth = 1000, blacklist = NULL, what = "mapq")
{
	bamGR <- list()
	chrom.lengths <- GenomeInfoDb::seqlengths(Rsamtools::BamFile(bamfile))
	chrom.in.data <- names(chrom.lengths)
	
	if(length(subChr)==1 && subChr=="all"){
		chrom.name <- chrom.in.data
	}else{
		if(!all(subChr %in% chrom.in.data)){
			stop("chr should be pattern with the chromosomes of given species or equals to 'all', please check it!")
		}
		chrom.name <- subChr
	}
	for(chr in chrom.name){
		bam_chr <- chromstaR::readBamFileAsGRanges(bamfile, bamindex = bamfile, chromosomes = chr, pairedEndReads = pairedEndReads,remove.duplicate.reads = removeDuplicateReads, min.mapq = minMapq, max.fragment.width = maxFragmentWidth, blacklist = blacklist, what = what)
		names(bam_chr) <- paste0(seqnames(bam_chr),"_",c(1:length(bam_chr)))
		bamGR[as.character(chr)] <- list(bam_chr)
	}
	return(bamGR)
}