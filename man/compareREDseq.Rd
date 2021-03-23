\name{compareREDseq}
\alias{compareREDseq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compare two RED Sequencing Dataset
}
\description{
For early stage experiment without replicates, compareREDseq outputs differentially cut RE sites between two experimental conditions using Fisher's Exact Test.
}
\usage{
compareREDseq(REsummary, col.count1 = 2, col.count2 = 3, multiAdj = TRUE,
 multiAdjMethod = "BH", maxP = 1, minCount = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{REsummary}{
A matrix with a RE id column, 2 count/weight column, see examples
}
  \item{col.count1}{
The column where the total count/weight for the 1st experimental condition is
}
  \item{col.count2}{
The column where the total count/weight for the 2nd experimental condition is
}
  \item{multiAdj}{
Whether apply multiple hypothesis testing adjustment, TURE or FALSE
}
  \item{multiAdjMethod}{
Multiple testing procedures, for details, see mt.rawp2adjp in multtest package
}
  \item{maxP}{The maximum p-value to be considered to be significant}
  \item{minCount}{For a RE site to be included, the tag count from at least one of the experimental condictions >= minimumCount}
}
\value{
\item{p.value}{the p-value of the test}
 \item{*.count}{weight/count from the input column col.count1 and col.count2}
\item{*.total}{total weight/count from input column col.count1 and col.count2}
 \item{REid}{the id of the restriction enzyme from the input}
 \item{odds.ratio}{an estimate of the odds ratio for 2nd experimental condition vs. 1st experimental condition}
 \item{*.adjusted.p.value}{applicable if multiAdj=TRUE, adjusted p.value using the method * specified in multiAdjMethod}
}
\author{
Lihua Julie Zhu
}


\seealso{
binom.test.REDseq
}
\examples{
library(REDseq)
x= cbind(c("RE1", "RE2", "RE3", "RE4"), c(10,1,100, 0),c(5,5,50, 40))
colnames(x) = c("REid", "control", "treated")
compareREDseq(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Statistics }
