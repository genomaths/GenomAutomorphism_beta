% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base2int.R
\name{base2int}
\alias{base2int}
\alias{base2int,character-method}
\alias{base2int,data.frame-method}
\title{Replace bases with integers from Z4 and Z5}
\usage{
base2int(base, ...)

\S4method{base2int}{character}(
  base,
  group = c("Z4", "Z5", "Z64", "Z125", "Z4^3", "Z5^3", "Z2"),
  cube = c("ACGT", "AGCT", "TCGA", "TGCA", "CATG", "GTAC", "CTAG", "GATC", "ACTG",
    "ATCG", "GTCA", "GCTA", "CAGT", "TAGC", "TGAC", "CGAT", "AGTC", "ATGC", "CGTA",
    "CTGA", "GACT", "GCAT", "TACG", "TCAG"),
  phychem = list(A = NULL, T = NULL, C = NULL, G = NULL, N = NULL)
)

\S4method{base2int}{data.frame}(
  base,
  group = c("Z4", "Z5", "Z64", "Z125", "Z4^3", "Z5^3", "Z2"),
  cube = c("ACGT", "AGCT", "TCGA", "TGCA", "CATG", "GTAC", "CTAG", "GATC", "ACTG",
    "ATCG", "GTCA", "GCTA", "CAGT", "TAGC", "TGAC", "CGAT", "AGTC", "ATGC", "CGTA",
    "CTGA", "GACT", "GCAT", "TACG", "TCAG"),
  phychem = list(A = NULL, T = NULL, C = NULL, G = NULL, N = NULL)
)
}
\arguments{
\item{base}{A character vector, string , or a dataframe of letters from the
DNA/RNA alphabet.}

\item{...}{Not in use.}

\item{group}{A character string denoting the group representation for the
given base or codon as shown in reference (2-3).}

\item{cube}{A character string denoting one of the 24 Genetic-code cubes,
as given in references (2-3).}

\item{phychem}{Optional. Eventually, it could be useful to represent
DNA bases by numerical values of measured physicochemical properties. If
provided, then this argument must be a named numerical list. For example,
the \code{\link[base]{scale}} values of deoxyribonucleic acids proton
affinity (available at \url{https://www.wolframalpha.com/} and
in cell phone app: Wolfram Alpha):

\eqn{list('A' = 0.87, 'C' = 0.88, 'T' = 0.82, 'G' = 0.89, 'N' = NA)}

where symbol 'N' provide the value for any letter out of DNA base alphabet.
In this example, we could write NA or 0 (see example section).}
}
\value{
A numerical vector.
}
\description{
A simple function to represent DNA bases as elements from
the Abelian group of integers modulo 4 (Z4), 5 (Z5), or 2 (Z2).
}
\details{
For Z2 (binary representation of DNA bases), the cube bases are represented
in their order by: '00', '01', '10', and '11' (examples section).
}
\examples{
## A triplet with a letter not from DNA/RNA alphabet
## 'NA' is introduced by coercion!
base2int("UDG")

## The base replacement in cube "ACGT and group "Z4"
base2int("ACGT")

## The base replacement in cube "ACGT and group "Z5"
base2int("ACGT", group = "Z5")

## A vector of DNA base triplets
base2int(c("UTG", "GTA"))

##  A vector of DNA base triplets with different number of triplets.
##  Codon 'GTA' is recycled!
base2int(base = c("UTGGTA", "CGA"), group = "Z5")

## Data frames 

base2int(data.frame(x1 = c("UTG", "GTA"), x2 = c("UTG", "GTA")))


## Cube bases are represented n their order by: '00', '01', '10', and '11',
## For example for cube = "ACGT" we have mapping: A -> '00', C -> '01',
## G -> '11', and C -> '10'.

base2int("ACGT", group = "Z2", cube = "ACGT")

}
\references{
\enumerate{
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
}
\seealso{
\link{base_coord}, \link{codon_coord}, and \link{dna_phychem}.
}
\author{
Robersy Sanchez \url{https://genomaths.com}
}
