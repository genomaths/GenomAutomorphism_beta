% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allClasses.R
\docType{class}
\name{ListCodonMatrix-class}
\alias{ListCodonMatrix-class}
\alias{ListCodonMatrix}
\alias{valid.ListCodonMatrix}
\title{A Convenient Class to Store Codon Coordinates in given
Genetic Code cube.}
\usage{
ListCodonMatrix(object, cube, group, seq_alias = NULL, names = NULL)

valid.ListCodonMatrix(x)
}
\arguments{
\item{object}{A list of CodonMatrix-class objects}

\item{x}{A 'ListCodonMatrix-class' object}
}
\description{
ListCodonMatrix-class objects are generated by function \link{codon_matrix}.
}
\examples{
## ListCodonMatrix-class objects are generated by function 'codon_matrix'.
## Let's create DNAStringSet-class object
base <- DNAStringSet(x = c( seq1 ='ACGTGATCAAGT',
                            seq2 = 'GTGTGATCCAGT',
                            seq3 = 'TCCTGATCAGGT'))

x1 <- codon_matrix(base)
x1

## Extract the first element
x1[1]
x1$codon.1
x1[[1]]
  
}
\keyword{internal}
