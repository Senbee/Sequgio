\name{Sequgio-package}
\alias{Sequgio-package}
\alias{Sequgio}
\docType{package}
\title{
Gene isoform expression and isoform-specific read distribution estimation based on RNA-seq data
}
\description{
There are two components to this package. These are: i) Construct design matrices used in expression estimation and extract read counts data from raw input filesf. ii) Output the expression levels and optionally the read distribution, standard error of the expression estimates and etc.
}
\details{
\tabular{ll}{
Package: \tab Sequgio\cr
Type: \tab Package\cr
Version: \tab 0.1.22\cr
Date: \tab 2011-12-07\cr
License: \tab GPL version 2.0 or newer\cr
Depend: \tab GenomicFeatures, Biostrings, lattice, Rsamtools, latticeExtra, RColorBrewer, AnnotationDbi, Biobase\cr
LazyLoad: \tab yes\cr
}
Index:\cr
\tabular{ll}{
makeXmatrix \tab create design matrix\cr
getCounts \tab import raw input files and extract read counts\cr
fitModels \tab Fit a Poisson model to estimate the gene isoform expression and read distribution\cr

}

}
\author{
Stefano Calza, Chen Suo, Agus Salim and Yudi Pawitan
Maintainer: Stefano Calza <stefano.calza@med.unibs.it>
}
\references{
Chen Suo, Stefano Calza, Agus Salim and Yudi Pawitan. Joint estimation of isoform expression and isoform-
specific read distribution using RNA-Seq data across samples - 2011 Submitted.
}


