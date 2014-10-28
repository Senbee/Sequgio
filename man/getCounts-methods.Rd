\name{getCounts}
\alias{getCounts}
\alias{getCounts-methods}
\alias{getCounts,seqCounts-method}
\title{
Extract read counts
}
\description{
Retrive the number of read counts in each exon or junction.
}
\usage{

\S4method{getCounts}{seqCounts}(Object)

}
\arguments{
  \item{Object}{an object of class seqCounts.}
}

\value{
A list contains the number of reads in exons and junctions. 
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
\code{\link{makeXmatrix}}, \code{\link{fitModels}}
}
\examples{
\dontrun{
library(RNAseqData.HNRNPC.bam.chr14)
data("TxDb")

bflist <- BamFileList(RNAseqData.HNRNPC.bam.chr14_BAMFILES,asMates=TRUE)

seqObj <- setCounts(bflist,txdb,fileName='test')

param <- SerialParam()

doCounts(seqObj,mcpar=param)
bigM <- getCounts(seqObj)
allCounts <- as.matrix(bigM)
}
}
\keyword{RNA-seq}
\keyword{reads}
