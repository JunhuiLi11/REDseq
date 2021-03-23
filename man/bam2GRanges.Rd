\name{bam2GRanges}
\alias{bam2GRanges}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Build a list of GRanges from bam file
}
\description{
Build a list of GRanges from bam file
}
\usage{bam2GRanges(bamfile,bamindex=bamfile,subChr='all',
 pairedEndReads=FALSE,removeDuplicateReads=FALSE,minMapq=10,
 maxFragmentWidth=1000,blacklist=NULL,what="mapq")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bamfile}{
	A sorted BAM file.
}
  \item{bamindex}{
	BAM index file. Can be specified without the .bai ending. If the index file does not exist it will be created and a warning is issued.
}
  \item{subChr}{
	Chromosome to search pattern(s), default is 'all'
}
  \item{pairedEndReads}{
	Set to TRUE if you have paired-end reads in your BAM files
}
  \item{removeDuplicateReads}{
	Whether or not to remove duplicate reads
}
  \item{minMapq}{
	Minimum mapping quality when importing from BAM files, default is 10
}
  \item{maxFragmentWidth}{
	Maximum allowed fragment length, default is 1000
}
  \item{blacklist}{
	A GRanges-class or a bed(.gz) file with blacklisted regions. Reads falling into those regions will be discarded
}
  \item{what}{
	A character vector of fields that are returned. Uses the Rsamtools::scanBamWhat function. See Rsamtools::ScanBamParam to see what is available.
}
}
\value{
Output a list of GRanges
}
\author{
Junhui Li
}

\examples{
	library(REDseq)
	bamfile = system.file("extdata", "example.bam", package="REDseq")
	GRlistBam <- bam2GRanges(bamfile, bamindex = bamfile, subChr = 'chr12', pairedEndReads = FALSE, removeDuplicateReads = FALSE, minMapq = 10, maxFragmentWidth = 1000, blacklist = NULL, what = "mapq")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
