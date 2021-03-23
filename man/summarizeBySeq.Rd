\name{summarizeBySeq}
\alias{summarizeBySeq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Output count/weight summary by sequences
}
\description{
Output count/weight summary by sequences with each row representing each sequences
}
\usage{
summarizeBySeq(assignedSeqs, by =c("Weight", "SEQid"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{assignedSeqs}{
	output from assignSeq2REsite	
}
  \item{by}{
Weight if sum up the weight for each sequence, SEQid if sum the occurrence of each sequence
}
}
\value{
	a matrix with SEQid as the first column and total count/weight as the second column
}
\author{
Lihua Julie Zhu
}

\seealso{
summarizeByRE, assignSeq2REsite
}
\examples{
	library(REDseq)
	data(example.assignedREDseq)
	summarizeBySeq(example.assignedREDseq, by="Weight")
	summarizeBySeq(example.assignedREDseq,by="SEQid")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc}
