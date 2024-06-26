% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allClasses.R, R/show-methods.R
\docType{class}
\name{MatrixSeq-class}
\alias{MatrixSeq-class}
\alias{MatrixSeq}
\alias{show,MatrixSeq-method}
\alias{show-MatrixSeq}
\title{Definition of MatrixSeq-class}
\usage{
MatrixSeq(seqs, matrix, names, aaindex, phychem, accession)

\S4method{show}{MatrixSeq}(object)
}
\arguments{
\item{seqs, matrix, names, aaindex, phychem, accession}{See detail section}

\item{object}{An object from 'MatrixSeq' class}
}
\value{
Given the slot values, it defines a MatrixSeq-class.

A MatrixSeq-class object

Print/show of a MatrixSeq-class object.
}
\description{
This is a very simple flexible class to store DNA and
aminoacid aligned sequences together with their physicochemical properties.
That is, a place where each aminoacid or codon from the sequence is
represented by numerical value from a physicochemical index.
}
\details{
\describe{
\item{\strong{seqs}: }{A string character vector of DNA or aminoacid
sequences.}
\item{\strong{matrix}: }{A numerical matrix or a numerical vector
(in the constructor) carrying the specified aminoacid physicochemical
indices for aminoacid in the DNA or aminoacid sequence(s).}
\item{\strong{names}: }{Alias/names/IDs DNA or aminoacid sequences.}
\item{\strong{aaindex}: }{Aminoacid index database where the
physicochemical index can be found.}
\item{\strong{phychem}: }{Description of the physicochemical index applied
to represent the DNA or aminoacid sequences.}
\item{\strong{accession}: }{Accession number or ID of the applied
physicochemical index in the database.}
}
}
\examples{
aln <- c(S1 = "ATGCGGATTAGA", S2 = "ATGACGATCACA", S3 = "ATGAGATCACAG")
cd <- DNAMultipleAlignment(aln)
r1 <- peptide_phychem_index(unmasked(cd), acc = "EISD840101")
r1

## Extract the second aminoacid sequence
r1[2]

## Using the sequence given name 
r1$S1

## Extract the second aminoacid value from the first sequence
r1[1,2]

## Change the name the second sequence
names(r1) <- c('S1', 'Seq1', 'S1')
r1

## Extract the amino acid sequences
slot(r1, 'seqs')

}
\author{
Robersy Sanchez \url{https://genomaths.com}
}
\keyword{internal}
