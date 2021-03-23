\name{distanceHistSeq2RE}
\alias{distanceHistSeq2RE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the distance distribution from sequence to the associated RE sites
}
\description{
Give an overview of the distance distribution from all assigned sequences to the associated RE sites. If average or estimate is used for assigning sequences to RE sites, the count for histogram drawing will be adjusted with the weight assigned.}
\usage{
distanceHistSeq2RE(assignedSeqs, longestDist = 1000, 
title = "histogram of distance to assigned RE site", 
xlab = "Distance to assigned RE site", ylab = "Frequency", ylim="")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{assignedSeqs}{
result returned from assignSeq2REsite
}
  \item{longestDist}{
longest distance to keep in the plot}
  \item{title}{
an overall title for the plot 
}
  \item{xlab}{
a title for the x axis
}
  \item{ylab}{
a title for the y axis
}
\item{ylim}{
range of y to be plotted
}
}
\author{
Lihua Julie Zhu
}


\seealso{
assignSeq2REsite, distanceHistSeq2RE
}
\examples{
library(REDseq)
data(example.assignedREDseq)
distanceHistSeq2RE(example.assignedREDseq,ylim=c(0,20))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graph }
