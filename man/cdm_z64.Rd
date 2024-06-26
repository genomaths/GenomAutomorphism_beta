% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_cdm_z64.R
\docType{data}
\name{cdm_z64}
\alias{cdm_z64}
\title{Codon Distance Matrices for the Standard Genetic Code on Z4}
\format{
A list object.
}
\usage{
data("cdm_z64", package = "GenomAutomorphism")
}
\description{
This is a list of 24 codon distance matrices created with function
\link{codon_dist_matrix} in the set of 24 genetic-code cubes on Z4
(using the default weights and assuming the standard genetic code (SGC).
The data set is created to speed up the computation when working with DNA
sequences from superior organisms. Since distance matrices are symmetric,
it is enough to provide the lower matrix. Each matrix is given as
named/labeled vector (see the example).
}
\examples{
## Load the data set
data("cdm_z64", package = "GenomAutomorphism")
cdm_z64

## The lower matrix (given as vector) for cube "TCGA" (picking out the 20 
## first values). Observe that this vector is labeled. Each numerical value 
## corresponds to the distance between the codons specified by the 
## name/label on it. For example, the distance between codons TTT and TCT 
## is: 0.0625.

head(cdm_z64[[ "TCGA" ]], 20)

}
\keyword{datasets}
