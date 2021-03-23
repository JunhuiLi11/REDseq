\name{example.assignedREDseq}
\alias{example.assignedREDseq}
\docType{data}
\title{
an example assigned REDseq dataset
}
\description{
an example assigned REDseq dataset generated from assignSeq2REsite
}
\usage{data(example.assignedREDseq)}
\format{
  The format is:
List of 3\cr
 $ passed.filter   :'data.frame': Sequences that passed the filters:\cr
  ..$ SEQid   :Sequence ID \cr
  ..$ REid    : Restriction Enzyme Site ID \cr
  ..$ Chr     : Chromosome \cr
  ..$ strand  : Strand \cr
  ..$ SEQstart: Sequence Start \cr
  ..$ SEQend  : Sequence End \cr
  ..$ REstart : Restriction Enzyme Site Start \cr
  ..$ REend   : Restriction Enzyme Site End \cr
  ..$ Distance: Distance from SEQstart to REstart \cr
  ..$ Weight  : Weighted count for this REid and this SEQid \cr
 $ notpassed.filter:'data.frame' : Sequences that did not pass the filters \cr
  ..$ SEQid   :Sequence ID \cr
  ..$ REid    : Restriction Enzyme Site ID \cr
  ..$ Chr     : Chromosome \cr
  ..$ strand  : Strand \cr
  ..$ SEQstart: Sequence Start \cr
  ..$ SEQend  : Sequence End \cr
  ..$ REstart : Restriction Enzyme Site Start \cr
  ..$ REend   : Restriction Enzyme Site End \cr
  ..$ Distance: Distance from SEQstart to REstart \cr
  ..$ Weight  : Weighted count for this REid and this SEQid \cr
 $ mREwithDetail   :'data.frame': Detailed information about the sequences that are associated with multiple REid - for debugging: \cr
  ..$ SEQid   :Sequence ID \cr
  ..$ REid    : Restriction Enzyme Site ID \cr
  ..$ Chr     : Chromosome \cr
  ..$ strand  : Strand \cr
  ..$ SEQstart: Sequence Start \cr
  ..$ SEQend  : Sequence End \cr
  ..$ REstart : Restriction Enzyme Site Start \cr
  ..$ REend   : Restriction Enzyme Site End \cr
  ..$ Distance: Distance from SEQstart to REstart \cr
  ..$ Weight  : Weighted count for this REid and this SEQid \cr
  ..$ count   : count of seed for this REid and SEQid \cr
  ..$ total.count: total number of seeds that are associated with this SEQid \cr
  }
\examples{
library(REDseq)
data(example.assignedREDseq)
## maybe str(example.assignedREDseq) ; plot(example.assignedREDseq) ...
}
\keyword{datasets}
