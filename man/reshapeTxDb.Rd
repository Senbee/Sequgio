\name{reshapeTxDb}
\alias{reshapeTxDb}

\title{
Reshape a TranscriptDb object with additional information
}
\description{
  This function process a TranscriptDb-class object disjoining
  overlapping genomic regions and adding additional
  information needed for design matrices building and fitting.
}
\usage{
reshapeTxDb(txdb,disjoin=TRUE,with.junctions=TRUE,probelen,ignore.strand=TRUE,junction.overlap=5L,exclude.non.std=TRUE,exclude.ids=NULL,include.only.ids=NULL,
what=c('exon','cds'),mcpar,verbose=FALSE,test.genes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{txdb}{An object of class TranscriptDb.}
  \item{disjoin}{Logical. Create disjoint exons ranges}
  \item{with.junctions}{Logical. Add 'biological junctions'}
  \item{probelen}{Reads length}
  \item{ignore.strand}{Logical. Ignore strand information}
  \item{junction.overlap}{Minimum overlap of a read with either exons for a
    junction read in order to be counted}
  \item{exclude.non.std}{Logical. Shall the 'non standard' chromosomes be removed?}
  \item{exclude.ids}{Character. Which features to exclude (Defaults to none)}
    \item{include.only.ids}{Character. Which features to include (Defaults to
    all)}
  \item{what}{What kind of exons to take into account. General "exons" or "cds"}
  \item{mcpar}{An object of (virtual) class "BiocParallelParam"
    specifying parameters for parallel computing using \code{\link{BiocParallel}}}
  \item{verbose}{Logical. Print information during computation}
  \item{test.genes}{Character vector. Subset of gene names to consider. Basically for debugging.}
  
}

\value{
A GRangesList containing GRanges for each gene
}
\references{
  Suo C, Calza S, Salim A, Pawitan Y. Joint estimation of isoform
  expression and isoform-specific read distribution using multisample
  RNA-Seq data.   Bioinformatics. 2014 Feb 15;30(4):506-13
}
\author{
Stefano Calza <stefano.calza@unibs.it>, Chen Suo, Agus Salim and Yudi Pawitan
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
\keyword{annotation}
