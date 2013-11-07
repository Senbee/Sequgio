\name{getCounts}
\alias{getCounts}
\title{
Extract read counts
}
\description{
Obtain the number of read counts in each exon or junction.
}
\usage{

getCounts(target,txdb,type="BAM",bam.params=NULL,minoverlap=5L,ignore.strand=TRUE,mapq.filter=NULL,use.samtools=FALSE,unique.only=FALSE,regesp=NULL,mcpar,only.proper=FALSE)

}
\arguments{
  \item{target}{a two columns data.frame contains the file names to be
    imported and sample names of the dataset.}
  \item{txdb}{An object produced by the function
    \code{\link{reshapeTxDb}}}
  \item{type}{Currently only BAM files are supported}
  \item{bam.params}{Parameters for BAM files input. See details.}
  \item{minoverlap}{Integer. Minimum overlap for reads to span
    junctions.}
  \item{ignore.strand}{Logical. Ignore strand information}
  \item{mapq.filter}{Specify optional \code{MAPQ} value for filtering}
  \item{use.samtools}{Use samtools program to filter BAM instead of
    Rsamtools functions.}
  \item{unique.only}{Logical. Keep only reads uniquely aligned?}
  \item{regesp}{Character. A regular expression to be applied to
    template name (\code{QNAME} in SAM file specifications)}
  \item{mcpar}{An object of (virtual) class "BiocParallelParam"
    specifying parameters for parallel computing using \code{\link{BiocParallel}}}
  \item{only.proper}{Logical. For pair-end reads: keep only proper
    matched reads?}
}

\value{
A list contains the number of reads in exons and junctions. 
}
\references{
Chen Suo, Stefano Calza, Agus Salim and Yudi Pawitan. Joint estimation of isoform expression and isoform-
specific read distribution using RNA-Seq data across samples - 2013 Submitted.
}
\author{
Stefano Calza <stefano.calza@med.unibs.it>, Chen Suo, Agus Salim and Yudi Pawitan
}

\seealso{
\code{\link{makeXmatrix}}, \code{\link{fitModels}}
}
\examples{
\dontrun{
library(RNAseqData.HNRNPC.bam.chr14)
data("TxDb")

param <- SerialParam()

target <- data.frame(filenames=RNAseqData.HNRNPC.bam.chr14_BAMFILES,
                     index=RNAseqData.HNRNPC.bam.chr14_BAMFILES,
                     samplenames=RNAseqData.HNRNPC.bam.chr14_RUNNAMES, stringsAsFactors=FALSE)

bigM <- getCounts(target,txdb,mcpar=param)
allCounts <- as.matrix(bigM)
}
}
\keyword{RNA-seq}
\keyword{reads}
