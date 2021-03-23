writeHits <-
function(seqname, matches, strand, file="", append=FALSE)
 {
 	if(missing(seqname))
	{
		stop("seqname is required as character!")
	}
	if (missing(matches) || class(matches) != "XStringViews")
	{
		stop("matches is required as XStringViews object!")
	}
	if (missing(strand))
	{
		stop("strand is required as + or - !")
	}
	if (file.exists(file) && !append)
 			warning("existing file ", file, " will be overwritten with 'append=FALSE'")
 		if (!file.exists(file) && append)
 			warning("new file ", file, " will have no header with 'append=TRUE'")
 		hits <- data.frame(chrom=rep.int(seqname, length(matches)),
 		chromStart=start(matches),
 		chromEnd=end(matches),
 		name=paste(names(matches),paste(seqname,1:length(matches),sep="."), sep="."),
		score= rep(1,length(matches)),
		strand=rep.int(strand, length(matches)),
 		check.names=FALSE)
 		write.table(hits, file=file, append=append, quote=FALSE, sep="\t",
 row.names=FALSE, col.names=!append)
 }

