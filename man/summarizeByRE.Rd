\name{summarizeByRE}
\alias{summarizeByRE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Output count/weight summary by restriction enzyme cut site ID (REid)
}
\description{
Output count/weight summary by REid with each row representing each REid
}
\usage{
summarizeByRE(assignedSeqs, by=c("Weight", "REid"),sampleName="",round=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{assignedSeqs}{
	output from assignSeq2REsite
}
  \item{by}{
Weight if sum up the weight for each REid, REid if sum the occurrence of each REid.
}
  \item{sampleName}{
The name of the sample used as the count column name.
}
  \item{round}{
TRUE: the sum of the weight is rounded up if the fraction part is greater than 0.5. FALSE: as it is.
}
}
\value{
a matrix with REid as the first column and total count/weight as the second column, that can be used for the downstream analysis with DEseq or edgeR.
}
\author{
Lihua Julie Zhu
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
summarizeBySeq, assignSeq2REsite
}
\examples{
	library(REDseq)
	data(example.assignedREDseq)
	summarizeByRE(example.assignedREDseq,by="REid",sampleName="example")
	summarizeByRE(example.assignedREDseq,by="Weight",sampleName="example")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
