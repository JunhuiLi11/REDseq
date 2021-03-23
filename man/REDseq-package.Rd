\name{REDseq-package}
\alias{REDseq-package}
\alias{REDseq}
\docType{package}
\title{
REDseq
}
\description{
REDSeq is a Bioconductor package for building genomic map of restriction enzyme sites REmap, assigning sequencing tags to RE sites using five different strategies, visualizing genome-wide distribution of differentially cut regions with the REmap as reference and the distance distribution of sequence tags to corresponding RE sites, generating count table for identifying statistically significant RE sites using edgeR or DEseq.
}
\details{
\tabular{ll}{
Package: \tab REDseq\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-05-10\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Lihua Julie Zhu

Maintainer:  
Lihua Julie Zhu <julie.zhu@umassmed.edu>
}
\references{
1. Roberts, R.J., Restriction endonucleases. CRC Crit Rev Biochem, 1976. 4(2): p. 123-64.2. Kessler, C. and V. Manta, Specificity of restriction endonucleases and DNA modification methyltransferases a review (Edition 3). Gene, 1990. 92(1-2): p. 1-248. \cr3. Pingoud, A., J. Alves, and R. Geiger, Restriction enzymes. Methods Mol Biol, 1993. 16: p. 107-200. \cr4. Anders, S. and W. Huber, Differential expression analysis for sequence count data. Genome Biol, 2010. 11(10): p. R106. \cr5. Robinson, M.D., D.J. McCarthy, and G.K. Smyth, edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics, 2010. 26(1): p. 139-40. \cr6. Zhu, L.J., et al., ChIPpeakAnno: a Bioconductor package to annotate ChIP-seq and ChIP-chip data. BMC Bioinformatics, 2010. 11: p. 237. \cr
7. Pages, H., BSgenome package. http://bioconductor.org/packages/2.8/bioc/\cr
vignettes/BSgenome/inst/doc/GenomeSearching.pdf \cr
8. Zhu, L.J., et al., REDseq: A Bioconductor package for Analyzing High Throughput Sequencing Data from Restriction Enzyme Digestion. (In preparation) \cr}

\keyword{ package }
\seealso{
buildREmap, assignSeq2REsit, plotCutDistribution, distanceHistSeq2RE, summarizeByRE, summarizeBySeq, compareREseq, binom.test.REDseq, bam2GRanges
}
\examples{
	if(interactive()){
	library(ChIPpeakAnno)
	REpatternFilePath = system.file("extdata", "examplePattern.fa", package="REDseq")
	library(BSgenome.Celegans.UCSC.ce2)
	buildREmap( REpatternFilePath, BSgenomeName=Celegans, outfile=tempfile(),chr='all')
	bamfile = system.file("extdata", "example.bam", package="REDseq")
	GRlistBam <- bam2GRanges(bamfile, bamindex = bamfile, subChr = 'chr12', pairedEndReads = FALSE, removeDuplicateReads = FALSE, minMapq = 10, maxFragmentWidth = 1000, blacklist = NULL, what = "mapq")
	library(REDseq)
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
	data(example.assignedREDseq)
	plotCutDistribution(example.assignedREDseq,example.map.list, 
chr="2", xlim =c(3012000, 3020000))
	distanceHistSeq2RE(example.assignedREDseq,ylim=c(0,20))
	summarizeByRE(example.assignedREDseq,by="Weight",sampleName="example")
	REsummary  =summarizeByRE(example.assignedREDseq,by="Weight")
	binom.test.REDseq(REsummary)
}
}
