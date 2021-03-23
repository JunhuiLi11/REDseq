\name{writeHits}
\alias{writeHits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
write the hits of pattern search to a bed file
}
\description{
write the hits of pattern search to a bed file, internal function used by searchPattern
}
\usage{
writeHits(seqname, matches, strand, file = "", append = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seqname}{
	Chromosome name
}
  \item{matches}{
	XStringViews object storing matched chromosome locations
}
  \item{strand}{
	strand of the match
}
  \item{file}{
	file path where the hits is written to
}
  \item{append}{
	TRUE if append to existing file, false if start a new file
}
}
\value{
	results are saved in the file specified by outfile
}
\references{
http://bioconductor.org/packages/2.8/bioc/vignettes/BSgenome/inst/doc/GenomeSearching.pdf
}
\author{
Lihua Julie Zhu
}

\seealso{
searchPattern, buildREmap
}
\examples{
	library(REDseq)
	x <- DNAString("AAGCGCGATATG")
	m <- matchPattern("GCGC", x)
	REDseq:::writeHits(seqname="chr1", m, strand="+",file="exampleWriteHits.bed", append=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
