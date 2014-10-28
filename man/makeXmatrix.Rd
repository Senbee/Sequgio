\name{makeXmatrix}
\Rdversion{1.1}
\alias{makeXmatrix}

\title{
Create design matrix
}
\description{
Extract useful information from an annotation file to create design matrix used in expression estimation.
}
\usage{
makeXmatrix(object,method=c("SE","PE"),mulen,verbose=FALSE,mcpar,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class GRangesList as returned by \code{\link{reshapeTxDb}}}
  \item{method}{Specify what method to use, if a design matrix structure for "single end" ("SE") or "paired-end" ("PE") library. Defualts to "SE"}
  \item{mulen}{Numeric. And estimate of average fragment length}
    \item{verbose}{verbose output. Defaults to no messages.}
  \item{mcpar}{An object of (virtual) class "BiocParallelParam"
    specifying parameters for parallel computing using
    \code{\link{BiocParallel}}}
  \item{...}{Additional arguments for specific methods}
}
\value{
A list containing the design matrix for each unit
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
\code{\link{getCounts}}, \code{\link{fitModels}}
}
\examples{
data("TxDb")
param = SerialParam()
Design <- makeXmatrix(txdb[1:5],method="PE",mulen=200,mcpar=param)
}
\keyword{RNA-seq}
\keyword{annotation}
