\name{plotCutDistribution}
\alias{plotCutDistribution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot cut frequencies of RE sites along a given chromosome
}
\description{
plot cut frequencies of RE sites along a chromosome, which gives a bird-eye view of genome-wide frequent-cut regions and RE inaccessible regions.
}
\usage{
plotCutDistribution(assignedSeqs,REmapList, chr="chr1",xlim, 
title="RE cut frequency distribution",
xlab="Chromosome Location (bp)",ylab="Frequency", 
round=TRUE, n.sequence)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{assignedSeqs}{
result returned from assignSeq2REsite
}
 \item{REmapList}{
REmapList used in assignSeq2REsite and generated from buildREmap
}
  \item{chr}{
chromosome to be plotted
}
  \item{xlim}{
range of x to be plotted
}
  \item{title}{
an overall title for the plot
}
  \item{xlab}{
a title for the x axis
}
  \item{ylab}{
a title for the y axis
}
  \item{round}{
TRUE: the sum of the weight is rounded up if the fraction part is greater than 0.5. FALSE: as it is.
}
 \item{n.sequence}{
total uniquely mapped sequences in the dataset for estimating the expected count for each RE site. If omitted, the expected count for each RE site will be set as 1 as default.
}
}
\author{
Lihua Julie Zhu
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
assignSeq2REsite, distanceHistSeq2RE
}
\examples{
	library(REDseq)
	data(example.assignedREDseq)
	data(example.map.list)
	plotCutDistribution(example.assignedREDseq,example.map.list,
chr="2", xlim =c(3012000, 3020000))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graph }
