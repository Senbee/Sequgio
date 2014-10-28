\name{setCounts}
\alias{setCounts}
\alias{setCounts-methods}
\alias{setCounts,seqCounts-method}
\title{
Create a container for read counts, a bigmatrix instance.
}
\description{
Create a container for read counts, a bigmatrix instance.
}
\usage{

setCounts(bfl,txdb,sample.names=NULL,fileName=NULL,rootpath=getwd())

}
\arguments{
  \item{bfl}{a BamFile or BamFileList object}
  \item{txdb}{a reshaped (\code{\link{reshapeTxDb}}) TranscriptDb
    object}
  \item{sample.names}{Optional. Character vector of sample names,
    i.e. counts matrix column names}
  \item{fileName}{Optional. File name for bigmatrix description file}
  \item{rootpath}{Optional. Path to where to store the description file}
  
}

\value{
An object of class \code{seqCounts}. 
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
