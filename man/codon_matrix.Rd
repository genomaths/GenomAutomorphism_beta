% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/codon_matrix.R
\name{codon_matrix}
\alias{codon_matrix}
\alias{codon_matrix,BaseSeqMatrix-method}
\alias{codon_matrix,DNAStringSet-method}
\alias{codon_matrix,DNAMultipleAlignment-method}
\title{Codon Coordinate Matrix}
\usage{
codon_matrix(base, ...)

\S4method{codon_matrix}{BaseSeqMatrix}(base, num.cores = 1L, tasks = 0L, verbose = TRUE, ...)

\S4method{codon_matrix}{DNAStringSet}(
  base,
  cube = c("ACGT", "AGCT", "TCGA", "TGCA", "CATG", "GTAC", "CTAG", "GATC", "ACTG",
    "ATCG", "GTCA", "GCTA", "CAGT", "TAGC", "TGAC", "CGAT", "AGTC", "ATGC", "CGTA",
    "CTGA", "GACT", "GCAT", "TACG", "TCAG"),
  group = c("Z4", "Z5"),
  num.cores = 1L,
  tasks = 0L,
  verbose = TRUE
)

\S4method{codon_matrix}{DNAMultipleAlignment}(
  base,
  cube = c("ACGT", "AGCT", "TCGA", "TGCA", "CATG", "GTAC", "CTAG", "GATC", "ACTG",
    "ATCG", "GTCA", "GCTA", "CAGT", "TAGC", "TGAC", "CGAT", "AGTC", "ATGC", "CGTA",
    "CTGA", "GACT", "GCAT", "TACG", "TCAG"),
  group = c("Z4", "Z5"),
  num.cores = 1L,
  tasks = 0L,
  verbose = TRUE
)
}
\arguments{
\item{base}{A \code{\link[Biostrings]{DNAMultipleAlignment}}, a
\code{\link[Biostrings]{DNAStringSet}}, or a \link{BaseSeqMatrix}.}

\item{...}{Not in use yet.}

\item{num.cores, tasks}{Parameters for parallel computation using package
\code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
use, i.e. at most how many child processes will be run simultaneously (see
\code{\link[BiocParallel]{bplapply}} and the number of tasks per job (only
for Linux OS).}

\item{verbose}{If TRUE, prints the function log to stdout}

\item{cube}{A character string denoting one of the 24 Genetic-code cubes,
as given in references (3-4).}

\item{group}{A character string denoting the group representation for the
given base or codon as shown in reference (3-4).}
}
\value{
A \link{ListCodonMatrix} class object with the codon coordinate on its
metacolumns.
}
\description{
This function build the coordinate matrix for each sequence from an aligned
set of DNA codon sequences.
}
\details{
The purpose of this function is making the codon coordinates from multiple
sequence alignments (MSA) available for further downstream statistical
analyses, like those reported in references (1) and (2).
}
\examples{
## Load the MSA of Primate BRCA1 DNA repair genes
data("brca1_aln")

## Get the DNAStringSet for the first 33 codons and apply 'codon_matrix'
brca1 <- unmasked(brca1_aln)
brca1 <- subseq(brca1, start = 1, end = 33)
codon_matrix(brca1)

## Get back the alignment object and apply 'codon_matrix' gives us the 
## same result.
brca1 <- DNAMultipleAlignment(as.character(brca1))
codon_matrix(brca1)

}
\references{
\enumerate{
\item Lorenzo-Ginori, Juan V., Aníbal Rodríguez-Fuentes, Ricardo Grau
Ábalo, and Robersy Sánchez Rodríguez. "Digital signal processing in the
analysis of genomic sequences." Current Bioinformatics 4, no. 1 (2009):
28-40.
\item Sanchez, Robersy. "Evolutionary analysis of DNA-protein-coding
regions based on a genetic code cube metric." Current Topics in Medicinal
Chemistry 14, no. 3 (2014): 407-417.
\item Robersy Sanchez, Jesus Barreto (2021) Genomic Abelian Finite
Groups.
\href{https://doi.org/10.1101/2021.06.01.446543}{doi: 10.1101/2021.06.01.446543}
\item M. V Jose, E.R. Morgado, R. Sanchez, T. Govezensky, The 24 possible
algebraic representations of the standard genetic code in six or in three
dimensions, Adv. Stud. Biol. 4 (2012) 119-152.\href{https://is.gd/na9eap}{PDF}.
\item R. Sanchez. Symmetric Group of the Genetic-Code Cubes. Effect of the
Genetic-Code Architecture on the Evolutionary Process MATCH Commun. Math.
Comput. Chem. 79 (2018) 527-560.
}

\enumerate{
\item 
\item 
}
}
\seealso{
\link{codon_coord}, \link{base_coord} and \link{base2int}.
}
\author{
Robersy Sanchez \url{https://genomaths.com}
}
