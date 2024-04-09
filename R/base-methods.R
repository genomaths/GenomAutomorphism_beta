## Copyright (C) 2021-2024 Robersy Sanchez <https://genomaths.com/>
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

## =========================== Docs ===========================

#' @rdname base-methods
#' @aliases seq2granges
#' @title DNA sequences to GRanges of bases and the reverse.
#' @description 
#' Functions 'seq2granges' and 'base_seq2string_set' are addressed to transform 
#' DNA sequences into \code{\link[GenomicRanges]{GRanges-class}} and the 
#' reverse transformation, respectively.
#' @param filepath A character vector containing the path to a file in
#' \emph{\strong{fasta}} format to be read. This argument must be given if
#' \emph{codon & base} arguments are not provided.
#' @param start,end,chr,strand Optional parameters required to build a
#' \code{\link[GenomicRanges]{GRanges-class}}. If not provided the default
#' values given for the function definition will be used.
#' @param ... Not in use.
#' @details
#' For the sake of brevety the metacolumns from the object returned by 
#' function 'seq2granges' are named as 'S1', 'S2', 'S3', and so on. The
#' original DNA sequence alias are stored in the slot named 'seq_alias'.
#' (see examples).
#' 
#' @returns 
#' 
#' ## seq2granges
#' 
#' This function returns a [BaseGroup] object carrying the DNA sequence(s), one
#' base per ranges. A [BaseGroup] class object inherits from  
#' \code{\link[GenomicRanges]{GRanges-class}}.
#' 
#' ## base_seq2string_set
#' 
#' This function returns a \code{\link[Biostrings]{DNAStringSet-class}}.
#' 
#' @seealso [Symmetric Group of the Genetic-Code Cubes.](
#' https://github.com/genomaths/GenomeAlgebra_SymmetricGroup)
#' @import S4Vectors
#' @import Biostrings
#' @importFrom methods new
#' @export
#' @author Robersy Sanchez <https://genomaths.com>
#' @seealso [base_coord] and [codon_coord].
#' @export
#' @examples
#' ## Load a multiple sequence alignment (MSA) of primate BRCA1 DNA repair 
#' ## genes 
#' data(brca1_aln2, package = "GenomAutomorphism")
#' brca1_aln2
#' 
#' ## Get BaseSeq-class object
#' gr <- seq2granges(brca1_aln2)
#' gr
#' 
#' ## Transform the BaseSeq-class object into a DNAStringSet-class object
#' str_set <- base_seq2string_set(gr)
#' str_set
#' 
#' ## Recovering the original MSA
#' require(Biostrings)
#' DNAMultipleAlignment(as.character(str_set))
#' 
setGeneric(
    "seq2granges",
    function(
        base = NULL,
        filepath = NULL,
        ...) {
        standardGeneric("seq2granges")
    }
)

## ===================== seq2granges ======================

#' @aliases seq2granges
#' @rdname base-methods
#' @import GenomicRanges
#' @importFrom methods new
#' @import Biostrings
#' @export
setMethod(
    "seq2granges", signature(base = "DNAStringSet_OR_NULL"),
    function(
        base = NULL,
        filepath = NULL,
        start = NA,
        end = NA,
        chr = 1L,
        strand = "+",
        seq_alias = NULL) {

    if (is.null(base) && is.null(filepath)) {
        stop(
            "*** Arguments 'base' & 'filepath' cannot be",
            " simultaneously NULL."
        )
    }
    
    if (!is.null(filepath) && is.character(filepath)) {
        base <- readDNAMultipleAlignment(filepath = filepath)
    }
    
    if (inherits(base, "DNAMultipleAlignment")) {
        base <- base@unmasked
        if (is.null(seq_alias))
            seq_alias <- names(base)
    }
    
    len <- min(width(base))
    
    if (!is.na(start) || !is.na(end)) {
        if (!is.na(start) && start > len) {
            stop(
                "*** The 'start' argument is greater than",
                " the 'base' length"
            )
        }
        if (!is.na(end) && end > len) {
            stop(
                "*** The 'end' argument is greater than",
                "   the 'base' length"
            )
        }
        base <- DNAStringSet(
            base,
            start = start,
            end = end
        )
    }
    
    if (length(base) > 1) {
        base <- t(as.matrix(base))
        colnames(base) <- NULL
        base <- data.frame(base)
    } else {
        base <- as.character(base)
        base <- strsplit(base, "")[[1]]
        base <- data.frame(base)
    }
    seq <- base
    if (is.na(start)) {
        start <- 1L
    }
    if (is.na(end)) {
        end <- len
    }
    
    pos <- seq(start, end, 1L)
    if (!is.null(dim(base))) 
        colnames(base) <- paste0("S", seq_len(ncol(base)))
    
    base <- data.frame(
        seqnames = chr,
        start = pos,
        end = pos,
        strand = strand,
        base
    )
    
    base <- makeGRangesFromDataFrame(base, keep.extra.columns = TRUE)
    
    base <- new(
        "BaseSeq",
        seqnames = seqnames(base),
        ranges = ranges(base),
        strand = strand(base),
        elementMetadata = base@elementMetadata,
        seqinfo = base@seqinfo,
        seq_alias = seq_alias
    )
    
    return(base)
    }
)

## ===================== base_seq2string_set ======================


#' @aliases base_seq2string_set
#' @rdname base-methods
#' @import GenomicRanges
#' @importFrom methods new
#' @import Biostrings
#' @param x A 'BaseSeq' class object. 
#' @export
setGeneric(
    "base_seq2string_set",
    function(
        x,
        ...) {
        standardGeneric("base_seq2string_set")
    }
)


#' @aliases base_seq2string_set
#' @rdname base-methods
#' @import GenomicRanges
#' @import Biostrings
#' @export
setMethod(
    "base_seq2string_set", signature(x = "BaseSeq"),
    function(x) {
        seq_alias <- x@seq_alias
        x <- mcols(x)
        x <- apply(x, 2, paste0, collapse = "")
        x <- DNAStringSet(x)
        names(x) <- seq_alias
        return(x)
    }
)


