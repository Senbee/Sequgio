\name{fitModels}
\alias{fitModels}
\title{
Fit a Poisson model for expression and read distribution estimation
}
\description{
Fit a joint statistical model that iteratively estimate non-
uniform isoform-specific read distribution and gene isoform
expression estimation.
}
\usage{
fitModels(iGene, design, counts, probeLen = 50L, minoverlap=5L,robust = TRUE, use.joint = TRUE, verbose = FALSE, ls.start=FALSE,ridge.lambda = 0, maxit = 50L, std.err = FALSE, Q1 = 0.75, DF = 3L, fix.lambda = exp(20),max.exons = 100, use.trueiLen=FALSE, trueiLen, attr.iLen=FALSE, attr.df=FALSE,useC=FALSE)
}
\arguments{
  \item{iGene}{Charater. Gene IDs}
  \item{design}{A list constains the design matrix for every gene.}
  \item{counts}{A matrix contains the read counts. }
  \item{probeLen}{Integer. Probe length in base pair (bp). Default is
    set to50.}
  \item{minoverlap}{Integer. Minimum overlap on either side for junction reads}
  \item{robust}{Use robust estimation}
  \item{use.joint}{Logical. If TRUE, a joint statistical model is fitted, taking non-uniform read distribution into account. Default to TRUE}
  \item{verbose}{Logical. If TRUE, gene ID is printed after final estimation is done. Default to FALSE}
  \item{ls.start}{Logical. If FALSE, robust regression is performed. Default to TRUE}
  \item{ridge.lambda}{Penalization term}
  \item{maxit}{Integer. Maximum number of iterations in joint estimation of read distribution and expression level. Default is set to 50.}
  \item{std.err}{Logical. If TRUE, standard error of expression level is included in the output. Default to FALSE}
  \item{Q1}{Numeric. The percentile of residuals used as a threshold in assigning weight. Default is set to 0.75.}
  \item{DF}{Numeric. Degress of freedom of the smoothing. Default is set to 3. }
  \item{fix.lambda}{Numeric. A fixed large penalization value used when there are too many exons in a gene. Default is set to exp(20).}
  \item{max.exons}{Numeric. If the number of exons in a gene exceeds this number, a fixed large lambda value is used in smoothing read intensity. Default is set to 100.}
  \item{use.trueiLen}{Logical. If True, estimated true read distribution is employed, thus only expression level needs to be estimated. Default to FALSE}
 \item{trueiLen}{True read distribution}
   \item{attr.iLen}{Logical. If TRUE estimated read intensity is included in the output. Default to FALSE}
  \item{attr.df}{Logical. If TRUE degree of freedom used in smoothing read intensity is included in the output. Default to FALSE}
  \item{useC}{Logical. If TRUE C code would be called. Default to FALSE
    (Experimental and still under testing)}
}

\value{
  A list contains expression leve, estimated read distribution, standard error of expression level and etc. 
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
\code{\link{makeXmatrix}}, \code{\link{getCounts}}
}
\examples{

data(Counts)
data(Design)

iGenes <- names(Design)

fit1 <- fitModels(iGenes[21],Design,allCounts)
fit2 <- lapply(iGenes[21:22],fitModels,design=Design,counts=allCounts)
}

\keyword{RNA-seq}
\keyword{poisson}
