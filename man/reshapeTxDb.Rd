\name{reshapeTxDb}
\Rdversion{1.1}
\alias{reshapeTxDb}

\title{
Reshape a TrascriptDb object with additional information
}
\description{

}
\usage{
reshapeTxDb <- function(txdb,disjoin=TRUE,with.junctions=TRUE,probelen,ignore.strand=TRUE,junction.overlap=5L,
                        exclude.non.std=TRUE,exclude.mir=TRUE,what=c('exon','cds'),mcpar,verbose=FALSE,test.genes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{txdb}{An object of class \code{TranscriptDb}}
  \item{disjoin}{Logical. Create disjoint exons ranges}
  \item{with.junctions}{Logical. Add 'biological junctions'}
  \item{probelen}{Reads length}
  \item{ignore.strand}{Logical. Ignore strand information}
  \item{junction.overlap}{Minimum overlap of a read with either exons for a
    junction read in order to be counted}
  \item{exclude.non.std}{Logical. Shall the 'non standard' chromosomes be removed?}
  \item{exclude.mir}{Logical. Remove microRNAs}
  \item{what}{What kind of exons to take into account. General "exons" or "cds"}
  \item{mcpar}{Parameter object from \code{\link{BiocParallel}} for parallel computation}
  \item{verbose}{Logical. Print information during computation}
  \item{test.genes}{Character vector. Subset of gene names to consider. Basically for debugging.}
  
}

\value{
A GRangesList containing GRanges for each gene
}
\references{
Chen Suo, Stefano Calza, Agus Salim and Yudi Pawitan. Joint estimation of isoform expression and isoform-
specific read distribution using RNA-Seq data across samples - 2011 Submitted.
}
\author{
Stefano Calza <stefano.calza@med.unibs.it>, Chen Suo, Agus Salim and Yudi Pawitan
}

\seealso{
\code{\link{makeTranscriptDb}}, \code{\link{exonsBy}}
}
\examples{
\dontrun{

seqs <- seqnames(seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene))
sel <- rep(FALSE,length(seqs))
names(sel) <- seqs
sel["chr14"] <- TRUE
isActiveSeq(TxDb.Hsapiens.UCSC.hg19.knownGene) <- sel

library(TxDb.Hsapiens.UCSC.hg19.knownGene)

param = SerialParam()

txdb <- reshapeTxDb(TxDb.Hsapiens.UCSC.hg19.knownGene,probelen = 50L,mcpar=param)

}
}

\keyword{RNA-seq}
\keyword{annotation}% __ONLY ONE__ keyword per line
