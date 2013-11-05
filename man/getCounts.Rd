\name{getCounts}
\Rdversion{1.1}
\alias{getCounts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract read counts
}
\description{
Obtain the number of read counts in each exon or junction.
}
\usage{
getCounts <- function(target, junction.names, ExByGene, minoverlap=25L,
                      type="BAM", bam.params, ignore.strand=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{target}{a two columns data.frame contains the file names to be imported and sample names of the dataset.}
  \item{junction.names}{A vector constains names of all junctions, indexed by their chromosome number and position.}
  \item{TxByGe}{An object of GRangesList.}
  \item{ExByGene}{An object of GRangesList.}
  \item{type}{Character. Type of input files}
  \item{bam.params}{A series of paramters passing into the function. See examples below.}
  \item{...}{}
}

\value{
A list contains the number of reads in exons and junctions. 
}
\references{
Chen Suo, Stefano Calza, Agus Salim and Yudi Pawitan. Joint estimation of isoform expression and isoform-
specific read distribution using RNA-Seq data across samples - 2011 Submitted.
}
\author{
Stefano Calza <stefano.calza@med.unibs.it>, Chen Suo, Agus Salim and Yudi Pawitan
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{makeXmatrix}}, \code{\link{fitModels}}
}
\examples{
\dontrun{

library(RNAseqData.HNRNPC.bam.chr14)

target <- data.frame(filenames=RNAseqData.HNRNPC.bam.chr14_BAMFILES,
                     index=RNAseqData.HNRNPC.bam.chr14_BAMFILES,
                     samplenames=RNAseqData.HNRNPC.bam.chr14_RUNNAMES, stringsAsFactors=FALSE)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{RNA-seq}
\keyword{reads}% __ONLY ONE__ keyword per line
