## Copyright (C) 2021 Robersy Sanchez <https://genomaths.com/>
## Author: Robersy Sanchez This file is part of the R package
## 'GenomAutomorphism'.  'GenomAutomorphism' is a free
## software: you can redistribute it and/or modify it under the
## terms of the GNU General Public License as published by the Free
## Software Foundation, either version 3 of the License, or (at
## your option) any later version.  This program is distributed in
## the hope that it will be useful, but WITHOUT ANY WARRANTY;
## without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for
## more details.  You should have received a copy of the GNU
## General Public License along with this program; if not, see
## <http://www.gnu.org/licenses/>.

#' Multiple Sequence Alignment (MSA) of Primate BRCA1 DNA repair genes.
#'
#' This is a \code{\link[Biostrings]{DNAMultipleAlignment}} carrying a MSA of
#' [BRCA1 DNA repair genes](https://bit.ly/3DimROD) to be used in the
#' examples provided for the package functions. The original file can be
#' downloaded from GitHub at: <https://bit.ly/3DimROD>
#' @usage 
#' data("brca1_aln", package = "GenomAutomorphism")
#'
#' @format \code{\link[Biostrings]{DNAMultipleAlignment}} class object.
#' @seealso [brca1_aln2], [brca1_autm], and [covid_aln].
#' @examples
#' data("brca1_aln", package = "GenomAutomorphism")
#' brca1_aln
#' 
"brca1_aln"
