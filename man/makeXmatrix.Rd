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
function(object,method=c("SE","PE"),probelen,mulen,mcpar,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class GRangesList as returned by \link{reshapeTxDb}}
  \item{method}{Specify what method to use, if a design matrix structure for "single end" ("SE") or "paired-end" ("PE") library. Defualts to "SE"}
  \item{probelen}{Integer specifying probe length}
  \item{mulen}{Numeric. And estimate of average fragment length}
  \item{mcpar}{Parameters for parallel computing based on \code{\link{BiocParallel}}}
}
\value{
A list containing the design matrix for each unit
}
\references{
Chen Suo, Stefano Calza, Agus Salim and Yudi Pawitan. Joint estimation of isoform expression and isoform-
specific read distribution using RNA-Seq data across samples - 2013 Submitted.
}
\author{
Stefano Calza <stefano.calza@med.unibs.it>, Chen Suo, Agus Salim and Yudi Pawitan
}


\seealso{
\code{\link{getCounts}}, \code{\link{fitModels}}
}
\examples{
\dontrun{

data("TxDb")
param = SerialParam()
Design <- makeXmatrix(txdb,method="PE",probelen=50L,mulen=200,mcpar=param)

}
}
\keyword{RNA-seq}
\keyword{annotation}
