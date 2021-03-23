\name{binom.test.REDseq}
\alias{binom.test.REDseq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Binomial test for REDseq dataset
}
\description{
 For any early stage experiment with one experimental condition and one biological replicate, binom.test.REDseq computes p-value for each RE site in the genome.
}
\usage{
binom.test.REDseq(REsummary, col.count = 2, multiAdj = TRUE, 
multiAdjMethod = "BH", prior.p = 0.000001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{REsummary}{
A matrix  returned from summarizeByRE with a RE id column, a count/weight column. See examples
}
  \item{col.count}{
The column where the total count/weight is
}
  \item{multiAdj}{
Whether apply multiple hypothesis testing adjustment, TURE or FALSE
}
  \item{multiAdjMethod}{
Multiple testing procedures, for details, see mt.rawp2adjp in multtest package
}
  \item{prior.p}{
It is the probability of assigning a mapped sequence tag to a given RE site. Assuming each RE site gets cut equally, then the prior.p = 1/number of total RE sites in the genome.
}
}
\value{
\item{p.value}{p-value of the test}
 \item{*.count}{weight/count from the input REsummary}
 \item{REid}{the id of the restriction enzyme from the input REsummary }
 \item{cut.frequency}{cut frequency}
 \item{*.adjusted.p.value}{applicable if multiAdj=TRUE, adjusted p.value using * method specified in multiAdjMethod}
}
\author{
Lihua Julie Zhu
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
compareREDseq
}
\examples{
	library(REDseq)
	REsummary = cbind(c("RE1", "RE2", "RE3"), c(10,1,100))
	colnames(REsummary) = c("REid", "control")
	binom.test.REDseq(REsummary)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Statistics}
