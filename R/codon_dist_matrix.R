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

#' Compute Codon Distance Matrix
#' @rdname codon_dist_matrix
#' @aliases codon_dist_matrix
#' @description This function computes the codon distance matrix based on the 
#' weighted Manhattan distance between codons estimated with function
#' \code{\link{codon_dist}}.
#' @details By construction, a distance matrix is a symmetric matrix. Hence,
#' the knowledge of lower triangular matrix is enough for its application to
#' any dowstream analysis.
#' @param genetic_code A single string that uniquely identifies the genetic 
#' code to extract. Should be one of the values in the id or name2 columns of 
#' \code{\link[Biostrings]{GENETIC_CODE_TABLE}}.
#' @param group A character string denoting the group representation for the
#' given codon sequence as shown in reference (2-3).
#' @param cube A character string denoting one of the 24 Genetic-code cubes,
#' as given in references (2-3).
#' @param output Format of the returned lower triangular matrix: as a list of
#' 63 elements (labeled) or as a labeled vector using codons as labels.
#' @param num.cores 
#' @param verbose 
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster 
#' @import foreach
#' @export
#' @seealso \code{\link{codon_dist}}.
#' @return A lower triangular matrix excluding the diagonal.
#' @examples 
#' ## The distance matrix for codons for the standard genetic code,
#' ## cube "ACGT" with base-triplet represented on the group "Z4". 
#' ## Each coordinate of the returned numerical vector corresponds to the 
#' ## distance between codons given in the coordinate name. 
#' x <- codon_dist_matrix()
#' 
#' ## The two first elements of the list carrying the lower triangular 
#' ## distance matrix.
#' 
#' x[1:2]
#' 
#' ## The distance matrix for codons for the Invertebrate Mitochondrial.
#' x <- codon_dist_matrix(genetic_code = "5", cube = "TGCA", group = "Z5",
#'                     output = "vector")
#' 
#' head(x, 10)
#' 
codon_dist_matrix <- function(
    genetic_code = "1",
    group = c("Z4", "Z5"),
    cube = c("ACGT", "AGCT", "TCGA", "TGCA", "CATG", 
            "GTAC", "CTAG", "GATC", "ACTG", "ATCG", 
            "GTCA", "GCTA", "CAGT", "TAGC", "TGAC", 
            "CGAT", "AGTC", "ATGC", "CGTA", "CTGA", 
            "GACT", "GCAT", "TACG", "TCAG"),
    output = c("list", "vector"),
    num.cores = 1L,
    verbose = FALSE) {
    
    group <- match.arg(group)
    cube <- match.arg(cube)
    output <- match.arg(output)
    
    gc <- getGeneticCode(id_or_name2 = genetic_code)
    nms <- names(gc)
    
    
    cl <- makeCluster(num.cores, type = "FORK")
    registerDoParallel(cl)
    
    distm <- foreach(k = seq_len(63)) %dopar% {
        d <- as.vector(outer(nms[k], nms[seq((k + 1), 64, 1)], 
                            FUN = codon_dist, group = group, 
                            cube = cube))
        names(d) <- nms[seq((k + 1), 64, 1)]
        return(d)
    }
    stopCluster(cl)
    names(distm) <- nms[ seq(63)]
    if (output != "list")
        distm <- unlist(distm)
    return(distm)
}


