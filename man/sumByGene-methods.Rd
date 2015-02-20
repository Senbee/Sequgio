\name{sumByGene}
\alias{sumByGene}
\alias{sumByGene-methods}
\alias{sumByGene,sumByGene-method}
\title{
Summarize isoform expression by gene
}
\description{
This function would simply summarize (sum) expression values within gene
}
\usage{
sumByGene(object,txdb)
}
\arguments{
  \item{object}{List. Results from fitModels}
  \item{txdb}{Either a GRangesList containing GRanges for each gene plus some
    annotation informations, or a data.frame with at least gene_id and
    tx_name fields.}
}

\value{
  A data.frame with gene level expression
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
\code{\link{reshapeTxDb}}, \code{\link{fitModels}}
}
\examples{

data(Counts)
data(Design)
data(TxDb)

iGenes <- names(Design)

fit1 <- fitModels(iGenes[21],Design,allCounts)
fit2 <- lapply(iGenes[21:22],fitModels,design=Design,counts=allCounts)

subByGene(fit2,txdb)

}

\keyword{RNA-seq}

