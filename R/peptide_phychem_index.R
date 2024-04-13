## Copyright (C) 2022-2024 Robersy Sanchez <https://genomaths.com/>
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


## =========================== Definition ==========================

#' @rdname peptide_phychem_index
#' @aliases peptide_phychem_index
#' @title Amino acid mutation matrix
#' @description 
#' This function applies the numerical indices representing various 
#' physicochemical and biochemical properties of amino acids and 
#' pairs of amino acids to whole DNA protein-coding or to aminoacid sequences.
#' As results, DNA protein-coding or the aminoacid sequences are represented
#' as numerical vectors which can be subject of further 
#' downstream statistical analysis and digital signal processing.
#' 
setGeneric("peptide_phychem_index",
    function(
        aa,  
        ...) standardGeneric("peptide_phychem_index"))



## =========================== Characters ==========================

#' @aliases peptide_phychem_index
#' @rdname peptide_phychem_index
#' @param acc Accession id for a specified mutation or contact potential 
#' matrix.
#' @param aaindex Database where the requested accession id is locate and from
#' where the aminoacid indices can be obtained. The possible values are:
#' "aaindex2" or "aaindex3".
#' @param userindex User provided aminoacid indices. This can be a numerical
#' vector or a matrix (20 x 20). If a numerical matrix is provided, then the 
#' aminoacid indices are computes as column averages. 
#' @param acc_list Logical. If TRUE, then the list of available matrices ids 
#' and index names is returned.
#' @return Depending on the user specifications, a mutation or contact 
#' potential matrix, a list of available matrices (indices) ids or index 
#' names can be returned. More specifically:
#' 
#' \itemize{
#'  \item{\strong{aa_mutmat}: }{Returns an aminoacid mutation matrix or
#'    a statistical protein contact potentials matrix.}
#'  \item{\strong{aa_index}: }{Returns the specified aminoacid physicochemical 
#'    indices.}
#' }
#' 
#' @param alphabet Whether the alphabet is from the 20 aminoacid (AA) or
#' four (DNA)/RNA base alphabet. This would prevent mistakes, i.e., 
#' the strings "ACG" would be a base-triplet on the DNA alphabet or simply
#' the amino acid sequence of alanine, cysteine, and glutamic acid.
#' 
#' @export
setMethod("peptide_phychem_index", signature(aa = "character"),
    function(
        aa, 
        acc = NULL,
        aaindex = NA,
        userindex = NULL,
        alphabet = c("AA", "DNA"),
        ...) {
        
        alphabet <- match.arg(alphabet)
        
        AA <- c("A","R","N","D","C","Q","E","G","H",
                "I","L","K","M","F","P","S","T","W","Y","V")
    
        if (alphabet == "DNA") {
            nc <- all(nchar(aa) != 3)
            
            if (length(aa) == 1 && nchar(aa) > 1 &&  nc) 
                aa <- base2codon(aa)
            
            if (all(nchar(aa) != 3))
                stop("*** The argument 'a' is not a DNA protein-coding",
                    " sequence, i.e., it is not a base-triplet sequence.")
            
            aa <- translation(aa)
            alphabet <- "AA"
        }

        if (length(aa) == 1 && nchar(aa) > 1)
            aa <- str2chr(aa)
        
        if (!is.null(userindex)) {
            if (is.matrix(userindex)) 
                phychem <- colMeans(userindex)
            
            phychem <- phychem[ match(aa, AA) ]
        }
        else {
            if (is.null(acc))
                stop("The accesion ID for the aaindex must be provided")
            phychem <- aa_phychem_index(acc = acc, aaindex = aaindex)
            if (is.matrix(phychem))
                phychem <- colMeans(phychem)
            
            phychem <- phychem[ match(aa, AA) ]
        }
        return(phychem)
    }
)


## ========================== DNAStringSet ========================


#' @aliases peptide_phychem_index
#' @rdname peptide_phychem_index
#' @param num.cores,tasks Parameters for parallel computation using package
#' \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
#' use, i.e. at most how many child processes will be run simultaneously
#' (see \code{\link[BiocParallel]{bplapply}} and the number of tasks per job
#' (only for Linux OS).
#' @param verbose If TRUE, prints the function log to stdout.
#' @param ... Not in use.
#' 
#' @export
setMethod("peptide_phychem_index", signature(aa = "DNAStringSet"),
    function(
        aa, 
        acc = NULL,
        aaindex = NA,
        userindex = NULL,
        num.cores = 1L,
        tasks = 0L,
        verbose = FALSE,
        ...) {
        
        aa <- translation(aa)
        l <- width(aa[1])
        seqs <- as.character(aa)
        
        phychem <- character()
        if (is.element(aaindex, c("aaindex1", "aaindex2", "aaindex3"))) {
            phychem <- aa_phychem_index(aaindex = aaindex, acc_list = TRUE)
            i <- grep(acc, phychem)
            if (length(i) > 0)
                phychem <- phychem[ i ]
            else
                stop("*** The accession provided '", acc, "' is not found in",
                    " the database '", aaindex, "'")
        }
        
        aa <- lapply(
                    seqs,
                    peptide_phychem_index, 
                    acc = acc, 
                    aaindex = aaindex,
                    userindex = userindex)
        
        
        nms <- names(seqs)
        aa <- matrix(unlist(aa, use.names = FALSE), ncol = l, byrow = TRUE)
        rownames(aa) <- paste0("S", seq(nrow(aa)))
        colnames(aa) <- paste0("A", seq(ncol(aa)))

        aa <- MatrixSeq(seqs = seqs,
                        matrix = aa, 
                        names = nms,
                        aaindex = if (is.na(aaindex)) "aaindex1" else aaindex,
                        phychem = phychem,
                        accession = acc)
        
        return(aa)
    }
)


