\name{searchPattern}
\alias{searchPattern}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 search pattern(s) in a genome
}
\description{
internal function for searching pattern(s) in a genome used by buildREmap leveraging BSgenome package
}
\usage{
searchPattern(dict0, BSgenomeName, outfile = "",chr='all')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dict0}{
	DNAStringSet object, see examples
}
  \item{BSgenomeName}{
BSgenome object, please refer to available.genomes in BSgenome package for details
}
  \item{outfile}{
	The file path to write the search results as a bed file
}
  \item{chr}{
	chromosome to search pattern(s), default is 'all'
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
buildREmap
}
\examples{
	library(REDseq)
	patternFile = system.file("extdata", "examplePattern.fa", package="REDseq")
	example.dict0 <- readDNAStringSet(patternFile, "fasta")
	library(BSgenome.Celegans.UCSC.ce2)
	REDseq:::searchPattern(dict0 = example.dict0, BSgenomeName=Celegans, outfile=tempfile(),chr='all')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal  }
