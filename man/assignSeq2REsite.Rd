\name{assignSeq2REsite}
\alias{assignSeq2REsite}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Assign mapped sequence tags to corresponding restriction enzyme (RE) cut sites
}
\description{
Given the sequence tags aligned to a genome as a list of GRanges, and a list map built using the buildREmap function, assignSeq2REsite first identifies RE sites that have mapped sequence tags around the cut position taking consideration of user-defined offset, sequence length and strand in the aligned sequences. These RE sites are used as seeds for assigning the remaining tags depending on which of five strategies the users select for partitioning sequences associated with multiple RE sites, i.e., unique, average,estimate, best and random.
}
\usage{
assignSeq2REsite(input.RD.list, REmap.RD.list, cut.end = 'both',cut.offset = 1,
seq.length = 36, allowed.offset = 5, min.FragmentLength = 60, max.FragmentLength = 300,
partitionMultipleRE = c("unique", "average", "estimate","best", "random"), chr="all")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.RD.list}{
List of GRanges as mapped sequences: see example below
}
  \item{REmap.RD.list}{
List of GRanges as restriction enzyme (RE) cut site map: see example below
}
  \item{cut.end}{
Single end('3' or '5') or both end('both') need to be analyzed, default is 'both'. Please set min.FragmentLength = 0, max.FragmentLength = 1, seq.length = 1 with both calls.
}
  \item{cut.offset}{
The cut offset from the start of the RE recognition sequence: index is 0 based, i.e.,1 means the RE cuts at position 2. 
}
  \item{seq.length}{
Sequence length: 36 means that the sequence tags are 36-base long. 
}
  \item{allowed.offset}{
Offset allowed to count for imperfect sticky end repair and primer addition.
}
  \item{min.FragmentLength}{
Minimum fragment length of the sequences size-selected for sequencing
}
  \item{max.FragmentLength}{
Maximum fragment length of the sequences size-selected for sequencing
}
  \item{partitionMultipleRE}{
The strategy for partitioning sequences associated with multiple RE sites.
For strategy unique, only sequence tags that are associated with a unique RE site within the distance between min.FragmentLength and max.FragmentLength are kept for downstream analysis. For strategy average, sequence tags are partitioned equally among associated RE sites. For strategy estimate, sequence tags are partitioned among associated RE sites with a weight function, which is determined using the count distribution of the RE seed sites described in the description section above. For strategy best, sequence tags are assigned to the most probable RE sties with the same weight function as that in strategy estimate. For strategy random, the sequence tags are randomly assigned to one of the multiple associated RE sites.
}
  \item{chr}{
Chromosomes to assign mapped sequence tags to corresponding restriction enzayme(RE) cut sites, default is 'all'
}
}
\value{
	\item{passed.filter}{Sequences assigned to RE(s), see the example r.unique$passed.filter}
	\item{notpassed.filter}{Sequences not assigned to any RE, see example r.unique$notpassed.filter}
	\item{mREwithDetail}{Detailed assignment information for sequences associated with multiple RE sites. Only available when partitionMultipleRE is set to average or estimate, see r.estimate$mREwithDetail in the examples 
}}
\references{1. Roberts, R.J., Restriction endonucleases. CRC Crit Rev Biochem, 1976. 4(2): p. 123-64. \cr2.Kessler, C. and V. Manta, Specificity of restriction endonucleases and DNA modification methyltransferases a review (Edition 3). Gene, 1990. 92(1-2): p. 1-248. \cr3. Pingoud, A., J. Alves, and R. Geiger, Restriction enzymes. Methods Mol Biol, 1993. 16: p. 107-200. \cr}
\author{
 Lihua Julie Zhu
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
buildREMap, example.REDseq, example.map, example.assignedREDseq
}
\examples{
	library(REDseq)
	#data(example.REDseq)
	#data(example.map)
	#example.map.list <- list()
	#example.map.list['2'] <- list(example.map)
	#example.REDseq.list <- list()
	#example.REDseq.list['2'] <- list(example.REDseq)
	data(example.REDseq.list)
	data(example.map.list)

	r.unique = assignSeq2REsite(example.REDseq.list, example.map.list, cut.end ='5',
	cut.offset = 1, seq.length = 36, allowed.offset = 5, 
	min.FragmentLength = 60, max.FragmentLength = 300, 
	partitionMultipleRE = "unique")

	r.average = assignSeq2REsite(example.REDseq.list, example.map.list, cut.end ='5',
	cut.offset = 1, seq.length = 36, allowed.offset = 5, min.FragmentLength = 60,
	max.FragmentLength = 300, partitionMultipleRE = "average")

	r.random = assignSeq2REsite(example.REDseq.list, example.map.list, cut.end ='5',
	cut.offset = 1, seq.length = 36, allowed.offset = 5, min.FragmentLength = 60,
	 max.FragmentLength = 300, partitionMultipleRE = "random")

	r.best = assignSeq2REsite(example.REDseq.list, example.map.list, cut.end ='5',
	cut.offset = 1, seq.length = 36, allowed.offset = 5, min.FragmentLength = 60,
	max.FragmentLength = 300, partitionMultipleRE = "best")

	r.estimate = assignSeq2REsite(example.REDseq.list, example.map.list, cut.end ='5',
	cut.offset = 1, seq.length = 36, allowed.offset = 5, min.FragmentLength = 60,
	max.FragmentLength = 300, partitionMultipleRE = "estimate")

	r.estimate$passed.filter
	r.estimate$notpassed.filter
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc}
