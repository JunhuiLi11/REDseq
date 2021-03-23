\name{buildREmap}
\alias{buildREmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Build a list of genome wide cut site map for a Restriction Enzyme (RE)
}
\description{
Build a genome-wide cut map for a Restriction Enzyme (RE)
}
\usage{
buildREmap(REpatternFilePath, format = "fasta", BSgenomeName, outfile, chr="all")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{REpatternFilePath}{
	File path storing the recognition pattern of a RE
}
  \item{format}{
	format of the pattern file, either "fasta" (the default) or "fastq
}
  \item{BSgenomeName}{
	BSgenome object, please refer to available.genomes in BSgenome package for details
}
  \item{outfile}{
	temporary output file for writing the matched chromosome location to
}
  \item{chr}{
	chromosome to search pattern(s), default is 'all'
}
}
\value{
Output REmap as a list of GRanges
}
\author{
Lihua Julie Zhu
}

\examples{
	library(REDseq)
	REpatternFilePath = system.file("extdata", "examplePattern.fa", package="REDseq")
	library(BSgenome.Celegans.UCSC.ce2)
	REmap.list <- buildREmap( REpatternFilePath, BSgenomeName=Celegans, outfile=tempfile(),chr='all')
	REmap.list[["chrI"]]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
