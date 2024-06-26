% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_mutscore.R
\name{get_mutscore}
\alias{get_mutscore}
\alias{get_mutscore,character,character-method}
\alias{get_mutscore,BaseSeq,missing-method}
\alias{get_mutscore,DNAStringSet,missing-method}
\alias{get_mutscore,DNAMultipleAlignment,missing-method}
\title{Get Mutation Score from an AAindex or a Mutation/Distance Matrix}
\usage{
get_mutscore(aa1, aa2, ...)

\S4method{get_mutscore}{character,character}(
  aa1,
  aa2,
  acc = NULL,
  aaindex = NULL,
  mutmat = NULL,
  alphabet = c("AA", "DNA"),
  num.cores = 1L,
  tasks = 0L,
  verbose = FALSE,
  ...
)

\S4method{get_mutscore}{BaseSeq,missing}(
  aa1,
  aa2,
  acc = NULL,
  aaindex = NULL,
  mutmat = NULL,
  alphabet = c("AA", "DNA"),
  stat = mean,
  numcores = 1L,
  num.cores = 1L,
  tasks = 0L,
  output = c("dist", "matrix", "vector"),
  na.rm = TRUE,
  verbose = TRUE,
  ...
)

\S4method{get_mutscore}{DNAStringSet,missing}(
  aa1,
  aa2,
  acc = NULL,
  aaindex = NULL,
  mutmat = NULL,
  alphabet = c("AA", "DNA"),
  stat = mean,
  num.cores = 1L,
  tasks = 0L,
  verbose = TRUE,
  output = c("dist", "matrix", "vector"),
  na.rm = TRUE,
  ...
)

\S4method{get_mutscore}{DNAMultipleAlignment,missing}(
  aa1,
  aa2,
  acc = NULL,
  aaindex = NULL,
  mutmat = NULL,
  alphabet = c("AA", "DNA"),
  stat = mean,
  num.cores = 1L,
  tasks = 0L,
  verbose = TRUE,
  output = c("dist", "matrix", "vector"),
  na.rm = TRUE,
  ...
)
}
\arguments{
\item{aa1, aa2}{A simple character representing an amino acids or a
character string of letter from the amino acid alphabet or base-triplets
from the DNA/RNA alphabet. If \strong{\emph{aa1}} is an object from any
of the classes: \link{BaseSeq}, \code{\link[Biostrings]{DNAStringSet}}, or
\code{\link[Biostrings]{DNAMultipleAlignment}}, then argument
\strong{\emph{aa2}} is not required.}

\item{...}{Not in use.}

\item{acc}{Accession id for a specified mutation or contact potential
matrix.}

\item{aaindex}{Database where the requested accession id is locate. The
possible values are:  "aaindex2" or "aaindex3".}

\item{mutmat}{A mutation or any score matrix provided by the user.}

\item{alphabet}{Whether the alphabet is from the 20 amino acid (AA) or
four (DNA)/RNA base alphabet. This would prevent mistakes, i.e.,
the strings "ACG" would be a base-triplet on the DNA alphabet or simply
the amino acid sequence of alanine, cysteine, and glutamic acid.}

\item{num.cores, tasks}{Parameters for parallel computation using package
\code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
use, i.e. at most how many child processes will be run simultaneously (see
\code{\link[BiocParallel]{bplapply}} and the number of tasks per job (only
for Linux OS).}

\item{verbose}{Optional. Only if num.cores > 1. If TRUE, prints the
function log to stdout.}

\item{stat}{Statistic that will be used to summarize the scores of the
DNA sequences provided. Only if \strong{\emph{aa1}} is an object from any of
the classes: \link{BaseSeq}, \code{\link[Biostrings]{DNAStringSet}}, or
\code{\link[Biostrings]{DNAMultipleAlignment}}.}

\item{numcores}{An integer to setup the number of parallel workers via
\code{\link[parallel]{makeCluster}}.}

\item{output}{Optional. Class of the returned object. Only if
\strong{\emph{aa1}} is an object from any of the classes: \link{BaseSeq},
\code{\link[Biostrings]{DNAStringSet}}, or
\code{\link[Biostrings]{DNAMultipleAlignment}}.}

\item{na.rm}{a logical evaluating to TRUE or FALSE indicating whether NA
values should be stripped before the computation proceeds.}
}
\value{
A single numeric score or a numerical vector, or if
\strong{\emph{aa1}} is an object from any of the classes: \link{BaseSeq},
\code{\link[Biostrings]{DNAStringSet}}, or
\code{\link[Biostrings]{DNAMultipleAlignment}}, then depending on the
user selection the returned object will be:
\enumerate{
\item A lower diagonal numerical vector of the sequence pairwise scores.
\item A \code{\link[stats]{dist}}-class object.
\item A whole score matrix.
}
}
\description{
This function is applied to get the mutation or contact
potential scores representing the similarity/distance between amino acids
corresponding to substitution mutations. The scores are retrieved from a
mutation matrix or a statistical protein contact potentials matrix from
\href{https://www.genome.jp/aaindex/}{AAindex} (ver.9.2).

Alternatively, the mutation scores can be estimated based on an user
mutation matrix, for example, see \link{aminoacid_dist} and \link{codon_dist_matrix}.
}
\details{
If a score matrix is provided by the user, then it must be a
symmetric matrix 20x20.
}
\examples{
## A single amino acids substitution mutation
get_mutscore("A", "C", acc = "MIYS930101", aaindex = "aaindex2")

## A tri-peptide mutation
get_mutscore(aa1 = "ACG", aa2 = "ATG", acc = "MIYS930101",
            aaindex = "aaindex2", alphabet = "AA")

## A single base-triple mutation, i.e., a single amino acid substitution 
## mutation
get_mutscore(aa1 = "ACG", aa2 = "CTA", acc = "MIYS930101",
            aaindex = "aaindex2", alphabet = "DNA")

## Peptides can be also written as:
get_mutscore(aa1 = c("A","C","G"), aa2 = c("C","T","A"), 
            acc = "MIYS930101", aaindex = "aaindex2", alphabet = "AA")
}
\seealso{
\link{aa_mutmat}, \link{aaindex2} and \link{aaindex3}.
}
\author{
Robersy Sanchez \url{https://genomaths.com}
}
