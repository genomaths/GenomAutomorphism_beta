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

## ========================== BaseGroup =============================

#' @aliases BaseGroup
#' @rdname BaseGroup
#' @title A class definition to store codon automorphisms in given in the
#' Abelian group representation.
#' @importClassesFrom S4Vectors DataFrame
#' @seealso \code{\link{automorphisms}}
#' @keywords internal
#' @export
#' @return Given the slot values define a BaseGroup-class.
setClass("BaseGroup",
    slots = c(
        seqnames = "Rle",
        ranges = "IRanges_OR_IPos",
        strand = "Rle",
        elementMetadata = "DataFrame",
        seqinfo = "Seqinfo",
        colnames = "character",
        group = "character",
        cube = "character"
    ),
    contains = "GRanges"
)

# ====================  Validity BaseGroup ======================== #


#' @rdname BaseGroup
#' @title Valid BaseGroup mcols
#' @param x A 'BaseGroup' object
#' @keywords internal
valid.BaseGroup.elem <- function(x) {
    m1 <- paste0(
        "*** This is not a valid  BaseGroup-class object.",
        " Columns from the matacolumn have the wrong names"
    )
    m2 <- paste0(
        "*** This is not a valid  BaseGroup-class object.",
        " seq1 or seq2 columns is not a base sequence"
    )
    m3 <- paste0(
        "*** Argument 'x' is not a BaseGroup-class object.",
        " The slot 'group' is not present or wrong naming."
    )
    m4 <- paste0(
        "*** Argument 'x' is not a BaseGroup-class object.",
        " The slot 'cube' is not present or wrong naming."
    )
    r1 <- r2 <- r3 <- r4 <- FALSE
    if (length(x) > 0) {
        coln <- x@colnames
        if (any(!is.element(
            coln,
            c("seq1", "seq2", "coord1", "coord2")
        ))) {
            r1 <- TRUE
        }
        if (unique(nchar(x$seq1)) != 1 || unique(nchar(x$seq2)) != 1) {
            r2 <- TRUE
        }

        group <- try(x@group, silent = TRUE)
        elem <- is.element(group, c("Z4", "Z5", "Z4^3", "Z5^3"))
        if (inherits(group, "try-error") || !elem) {
            r3 <- TRUE
        }

        cube <- try(x@cube, silent = TRUE)
        celem <- is.element(cube, c(
            "ACGT", "AGCT", "TCGA", "TGCA", "CATG",
            "GTAC", "CTAG", "GATC", "ACTG", "ATCG",
            "GTCA", "GCTA", "CAGT", "TAGC", "TGAC",
            "CGAT", "AGTC", "ATGC", "CGTA", "CTGA",
            "GACT", "GCAT", "TACG", "TCAG"
        ))
        if (inherits(celem, "try-error") || !celem) {
            r4 <- TRUE
        }
    }
    if (any(c(r1, r2, r3, r4))) {
        res <- c(m1, m2, m3, m4)[c(r1, r2, r3, r4)]
        return(res[1])
    }
    NULL
}

#' @rdname BaseGroup
#' @title Valid 'BaseGroup' inheritance from 'GRanges' class
#' @param x A 'BaseGroup object'
#' @keywords internal
valid.GRanges <- function(x) {
    if (length(x) > 0) {
        if (!inherits(x, "GRanges")) {
            return("*** This is not a valid  Automorphism-class object.")
        }
    }
    NULL
}

#' @rdname valid.BaseGroup
#' @title Valid BaseGroup
#' @param x A 'BaseGroup object'
#' @keywords internal
valid.BaseGroup <- function(x) {
    c(valid.GRanges(x), valid.BaseGroup.elem(x))
}

S4Vectors:::setValidity2("BaseGroup", valid.BaseGroup)


## ========================== CodonGroup =============================

#' @aliases CodonGroup
#' @rdname CodonGroup
#' @title A class definition to store codon automorphisms in given in the
#' Abelian group representation.
#' @importClassesFrom S4Vectors DataFrame
#' @importClassesFrom GenomicRanges GRanges
#' @seealso \code{\link{automorphisms}}
#' @keywords internal
#' @export
#' @return Given the slot values define a CodonGroup-class.
setClass("CodonGroup",
    slots = c(
        seqnames = "Rle",
        ranges = "IRanges_OR_IPos",
        strand = "Rle",
        elementMetadata = "DataFrame",
        seqinfo = "Seqinfo",
        colnames = "character",
        group = "character",
        cube = "character"
    ),
    contains = "GRanges"
)

# ====================  Validity CodonGroup ======================== #


#' @rdname CodonGroup
#' @title Valid CodonGroup mcols
#' @param x A 'CodonGroup' object
#' @keywords internal
valid.CodonGroup.mcols <- function(x) {
    if (length(x) > 0) {
        coln <- x@colnames
        m1 <- paste0(
            "*** This is not a valid  CodonGroup-class object.",
            "Columns from the matacolumn have the wrong names"
        )
        m2 <- paste0(
            "*** This is not a valid  BaseGroup-class object.",
            "seq1 or seq2 columns is not a base-triplet sequence"
        )
        m3 <- paste0(
            "*** This is not a CodonGroup-class object.",
            " The slot 'cube' is not present or wrong naming."
        )
        m4 <- paste0(
            "*** Argument 'x' is not a CodonGroup-class object.",
            "The slot 'group' is not present or wrong naming."
        )
        r1 <- r2 <- r3 <- r4 <- FALSE

        if (any(!is.element(
            coln,
            c("seq1", "seq2", "coord1", "coord2")
        ))) {
            r1 <- TRUE
        }
        if (unique(nchar(x$seq1)) != 3 || unique(nchar(x$seq2)) != 3) {
            r2 <- TRUE
        }

        cube <- try(x@cube, silent = TRUE)
        celem <- is.element(cube, c(
            "ACGT", "AGCT", "TCGA", "TGCA", "CATG",
            "GTAC", "CTAG", "GATC", "ACTG", "ATCG",
            "GTCA", "GCTA", "CAGT", "TAGC", "TGAC",
            "CGAT", "AGTC", "ATGC", "CGTA", "CTGA",
            "GACT", "GCAT", "TACG", "TCAG"
        ))
        if (inherits(celem, "try-error") || !celem) {
            r3 <- TRUE
        }

        group <- try(x@group, silent = TRUE)
        elem <- is.element(group, c("Z4", "Z5", "Z4^3", "
                                    Z5^3", "Z64", "Z125"))
        if (inherits(group, "try-error") || !elem) {
            r4 <- TRUE
        }
    }
    if (any(c(r1, r2, r3, r4))) {
        res <- c(m1, m2, m3, m4)[c(r1, r2, r3, r4)]
        return(res[1])
    }
    NULL
}

#' @rdname valid.CodonGroup
#' @title Valid CodonGroup
#' @param x A 'CodonGroup object'
#' @keywords internal
valid.CodonGroup <- function(x) {
    c(valid.GRanges(x), valid.CodonGroup.mcols(x))
}

S4Vectors:::setValidity2("CodonGroup", valid.CodonGroup)


setClassUnion("BaseGroup_OR_CodonGroup", c("BaseGroup", "CodonGroup"))

## ========================== CodonSeq =============================

#' @rdname CodonSeq
#' @title A class definition to store codon coordinates given in the Abelian
#' group and the codon sequence.
#' @description An objects from 'CodonSeq' or 'MatrixList' class is returned by
#' function \code{\link{get_coord}}. This object will store the coordinate of
#' each sequence in a list of 3D-vectors or a list of vectors located in the
#' slot named 'CoordList'. The original codon sequence (if provided) will be
#' stored in the slot named 'SeqRanges'.
#' @importFrom S4Vectors setValidity2
#' @importFrom methods validObject
#' @keywords internal
#' @aliases CodonSeq
#' @export
#' @return Given the slot values define a CodonSeq-class.
setClass("CodonSeq",
    slots = c(
        CoordList = "list",
        SeqRanges = "GenomicRanges_OR_missing"
    )
)

#' @aliases coordList
#' @rdname CodonSeq
#' @title Method to extract 'CoordList' slot from a
#' \code{\link{CodonSeq-class}}
#' @param x An object from \code{\link{CodonSeq-class}}.
setGeneric(
    "coordList",
    function(x) standardGeneric("coordList")
)

#' @aliases coordList
#' @rdname CodonSeq
#' @export
setMethod(
    "coordList", signature(x = "CodonSeq"),
    function(x) x@CoordList
)

#' @aliases seqRanges
#' @rdname CodonSeq
#' @title Method to extract 'SeqRanges' slot from a
#' \code{\link{CodonSeq-class}}
#' @param x An object from \code{\link{CodonSeq-class}}.
setGeneric(
    "seqRanges",
    function(x) standardGeneric("seqRanges")
)

#' @aliases seqRanges
#' @rdname CodonSeq
#' @export
setMethod(
    "seqRanges", signature(x = "CodonSeq"),
    function(x) x@SeqRanges
)

## ========================== MatrixList class =============================

#' @aliases MatrixList
#' @rdname MatrixList
#' @title Definition of MatrixList-class
#' @description  A class denoting a list of matrices.
#' @keywords internal
#' @export
#' @return Given the slot values, it defines a MatrixList-class.
setClass("MatrixList",
    slots = c(
        matrices = "list",
        names = "character"
    )
)

## ======================== Validity MatrixList ======================= #
#' @rdname valid.MatrixList
#' @title Valid MatrixList
#' @param x A 'MatrixList object'
#' @keywords internal
valid.MatrixList <- function(x) {
    if (!all(slapply(x, function(y) inherits(y, "matrix")))) {
        return(
            "*** Not all the elements of the MatrixList object",
            " are from 'matrix' class."
        )
    }
    NULL
}


#' @importFrom S4Vectors setValidity2
#' @importClassesFrom Biostrings DNAMultipleAlignment DNAMultipleAlignment
#' @rdname valid.MatrixList
#' @keywords internal
#' @importClassesFrom Biostrings DNAMultipleAlignment DNAStringSet
#' @export
setClassUnion(
    "DNAStringSet_OR_NULL",
    c("DNAStringSet", "DNAMultipleAlignment", "NULL", "missing")
)

## ========================== Automorphism =============================

#' @rdname Automorphism
#' @aliases Automorphism
#' @title A class definition to store codon automorphisms in given as an
#' Abelian group representation.
#' @description Two classes are involved in to storing codon automorphisms:
#' \emph{\strong{Automorphism-class}} and
#' \emph{\strong{AutomorphismList-class}}.
#' @details \emph{\strong{Automorphism-class}} inherits from a
#' \code{\link[GenomicRanges]{GRanges}}.
#'
#' @section Methods:
#' ## as(from, "Automorphism")
#' Permits the transformation of a \code{\link[base]{data.frame}} or a
#' \code{\link[S4Vectors]{DataFrame-class}} object into
#' \emph{\strong{Automorphism-class}} object the proper columns are provided. 
#' An \emph{\strong{Automorphism-class}} object has six columns: "seq1", 
#' "seq2","coord1", "coord2", "autm", and "cube". See the examples for function
#' \code{\link{automorphisms}}. Observe that as the
#' \emph{\strong{Automorphism-class}} inherits from
#' \code{\link[GenomicRanges]{GRanges-class}} the transformation starting from a
#' \code{\link[GenomicRanges]{GRanges-class}} object into an
#' \emph{\strong{Automorphism-class}} is straightforward. However, the
#' transformation starting from a \code{\link[base]{data.frame}} or a
#' \code{\link[S4Vectors]{DataFrame-class}} object \eqn{"x"} requires for the
#' creation of an additional \code{\link[GenomicRanges]{GRanges-class}} object,
#' which by default will have the argument seqnames = "1", strand = "+"
#' start/end = 1:nrow(x), length = nrow(x). These details must be keep in mind
#' to prevent fundamental errors in the downstream analyses.
#'
#' ## \emph{\strong{AutomorphismList-class}} has the method
#' ### as.AutomorphismList(x)
#' \emph{\strong{as.AutomorphismList}} function transform a list of
#' \code{\link[GenomicRanges]{GRanges-class}}, a
#' \code{\link[GenomicRanges]{GRangesList-class}}, a list of
#' \code{\link[base]{data.frame}} or a \code{\link[S4Vectors]{DataFrame-class}}
#' objects into a \emph{\strong{AutomorphismList-class}} object.
#'
#' @seealso \code{\link{automorphisms}}
#' @keywords internal
#' @importClassesFrom S4Vectors DataFrame
#' @importClassesFrom GenomicRanges GRanges
#' @export
#' @return Given the slot values, it defines an Automorphism-class object.
setClass("Automorphism",
    slots = c(
        seqnames = "Rle",
        ranges = "IRanges_OR_IPos",
        strand = "Rle",
        elementMetadata = "DataFrame",
        seqinfo = "Seqinfo",
        colnames = "character"
    ),
    contains = "GRanges"
)

#' @rdname Automorphism
#' @importFrom S4Vectors DataFrame
#' @keywords internal
#' @export
setClassUnion(
    "DataFrame_OR_data.frame",
    c("DataFrame", "data.frame")
)

#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors mcols
#' @importFrom GenomeInfoDb Seqinfo
#' @importClassesFrom S4Vectors DataFrame
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
setAs(
    "DataFrame_OR_data.frame", "Automorphism",
    function(from) {
        nr <- nrow(from)
        pos <- seq(1, nr, 1)
        gr <- GRanges(
            seqnames = 1,
            ranges = IRanges(start = pos, end = pos),
            strand = "+"
        )
        mcols(gr) <- from

        new("Automorphism",
            seqnames = seqnames(gr),
            ranges = ranges(gr),
            strand = strand(gr),
            elementMetadata = from,
            seqinfo = seqinfo(gr),
            colnames = colnames(from)
        )
    }
)


# ======================== Validity Automorphism ======================= #
#' @rdname Automorphism
#' @title Valid Automorphism mcols
#' @param x A 'Automorphism object'
#' @keywords internal
valid.Automorphism.mcols <- function(x) {
    alf <- c("A", "C", "G", "T", "-")
    if (length(x) > 0) {
        m1 <- m2 <- FALSE
        if (inherits(x, "GRanges")) {
            coln <- colnames(x@elementMetadata)
        } else {
            coln <- x@colnames
        }
        if (any(!is.element(
            c(
                "seq1", "seq2", "coord1",
                "coord2", "autm", "cube"
            ),
            coln
        ))) {
            m1 <- TRUE
        }
        if (unique(nchar(x$seq1)) != 3 || unique(nchar(x$seq2)) != 3) {
            m2 <- TRUE
        }
        if (m2) {
            if (all(is.element(x$seq1, alf)) && all(is.element(x$seq1, alf))) {
                m2 <- FALSE
            }
        }

        if (m1 || m2) {
            return("*** This is not a valid Automorphism-class object.")
        }
    }
    NULL
}

#' @rdname Automorphism
#' @title Valid Automorphism
#' @param x A 'Automorphism object'
#' @keywords internal
valid.Automorphism <- function(x) {
    c(valid.GRanges(x), valid.Automorphism.mcols(x))
}

S4Vectors:::setValidity2("Automorphism", valid.Automorphism)


## ========================== AutomorphismList =============================

#' @rdname Automorphism
#' @title A class definition to store list of Automorphism class objects.
#' @description A class definition to store list of Automorphism class objects
#' derived from the pairwise automorphism estimation from pairwise
#' alignments.
#' @importFrom methods validObject setClass
#' @keywords internal
#' @export
#' @aliases AutomorphismList
setClass("AutomorphismList",
    slots = c(
        DataList = "list",
        SeqRanges = "GenomicRanges_OR_missing"
    )
)

#' @rdname Automorphism
#' @aliases as.AutomorphismList
#' @title AutomorphismList-class object constructor from a list.
#' @description The function build a AutomorphismList-class object from a
#' list of \code{\link[S4Vectors]{DataFrame}} or a \code{\link{automorphisms}}
#' class object.
#' @param x A \code{\link[S4Vectors]{DataFrame}} or a
#' \code{\link{automorphisms}} class object.
#' @param gr A \code{\link[GenomicRanges]{GRanges-class}} object.
#' @importFrom GenomicRanges GRanges GRangesList
#' @importFrom S4Vectors mcols DataFrame
#' @importFrom GenomicRanges GRanges GRangesList
#' @importFrom S4Vectors mcols
#' @importFrom methods setGeneric
#' @export
setGeneric(
    "as.AutomorphismList",
    function(x,
    grs = GRanges(),
    ...) {
        standardGeneric("as.AutomorphismList")
    }
)


setClassUnion("GRanges_OR_NULL", c("GRanges", "NULL", "missing"))


#' @rdname Automorphism
#' @aliases as.AutomorphismList
#' @export
setMethod(
    "as.AutomorphismList",
    signature(x = "GRangesList", grs = "GRanges_OR_NULL"),
    function(x,
    grs = GRanges(),
    ...) {
        if (length(grs) == 0) {
            grs <- x[[1]]
        }
        mcols(grs) <- NULL

        x <- lapply(x, function(y) {
            x <- as(x, "Automorphism")
            gr <- x
            mcols(gr) <- NULL
            if (gr != grs) {
                stop("*** The ranges from the GRanges-class objects
                    must equals.")
            }
            return(mcols(x))
        })

        new("AutomorphismList",
            DataList = x,
            SeqRanges = grs
        )
    }
)

#' @rdname Automorphism
#' @aliases as.AutomorphismList
#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors mcols
#' @export
setMethod(
    "as.AutomorphismList",
    signature(x = "list", grs = "GRanges_OR_NULL"),
    function(x,
    grs = GRanges(),
    ...) {
        if (length(grs) == 0) {
            if (inherits(x[[1]], "GRanges")) {
                grs <- x[[1]]
            } else {
                if (inherits(x[[1]], c("DataFrame", "data.frame"))) {
                    pos <- seq(1, nrow(x[[1]]), 1)
                    grs <- GRanges(
                        seqnames = 1,
                        ranges = IRanges(start = pos, end = pos),
                        strand = "+"
                    )
                } else {
                    stop(
                        "*** The argument of 'x' must be a list of ",
                        "objects from any of the classes: 'GRanges', ",
                        "'DataFrame', or 'data.frame'."
                    )
                }
            }
        }

        if (!is.null(mcols(grs))) {
            mcols(grs) <- NULL
        }

        if (all(slapply(x, function(y) inherits(y, "GRanges")))) {
            if (length(grs) == length(x)) {
                grs <- x
                mcols(grs) <- NULL
            }

            if (length(grs) == 0) {
                grs <- x[[1]]
                mcols(grs) <- NULL
            }

            x <- lapply(x, function(y) {
                y <- as(y, "Automorphism")
                return(mcols(y))
            })

            x <- new("AutomorphismList",
                DataList = x,
                SeqRanges = grs
            )
        }
        if (!inherits(x, "AutomorphismList")) {
            if (all(slapply(x, function(y) inherits(y, "DataFrame")))) {
                x <- new("AutomorphismList",
                    DataList = x,
                    SeqRanges = grs
                )
            }
        }
        return(x)
    }
)


#' @rdname Automorphism
#' @export
setMethod("names",
    signature = "AutomorphismList",
    function(x) names(x@DataList)
)

#' @rdname Automorphism
#' @export
setReplaceMethod(
    "names", "AutomorphismList",
    function(x, value) {
        names(x@DataList) <- value
        return(x)
    }
)

#' @rdname Automorphism
#' @export
setMethod("as.list",
    signature = "AutomorphismList",
    function(x) {
        x <- getAutomorphisms(x)
        return(x@DataList)
    }
)


setAs("AutomorphismList", "list", function(from) {
    from <- getAutomorphisms(from)
    return(from@DataList)
})


setAs("AutomorphismList", "GRangesList", function(from) {
    from <- getAutomorphisms(from)
    from <- as.list(from)
    return(as(from, "GRangesList"))
})


#' @importClassesFrom GenomicRanges GRangesList
setMethod("unlist",
    signature = "AutomorphismList",
    function(x) {
        x <- as(x, "GRangesList")
        return(unlist(x))
    }
)


## ====================== Validity AutomorphismList ================== #

#' @rdname Automorphism
#' @title Valid AutomorphismList mcols
#' @param x A 'AutomorphismList object'
#' @importFrom S4Vectors mcols
#' @keywords internal

valid.AutomorphismList <- function(x) {
    m1 <- FALSE
    if (!(inherits(x@DataList[[1]], "Automorphism") ||
        inherits(x@DataList[[1]], "DataFrame"))) {
        m1 <- TRUE
    }

    if (inherits(x@DataList[[1]], "Automorphism")) {
        if (any(!slapply(
            x@DataList,
            function(y) {
                return(inherits(y, "Automorphism") && validObject(y))
            }
        ))) {
            m1 <- TRUE
        }
    }

    if (inherits(x@DataList[[1]], "DataFrame")) {
        if (any(!slapply(
            x@DataList,
            function(y) {
                return(inherits(y, "DataFrame") && validObject(y))
            }
        ))) {
            m1 <- TRUE
        }
        if (!inherits(x@SeqRanges, "GRanges")) {
            m1 <- TRUE
        }
    }
    if (m1) {
        return("*** This is not a valid AutomorphismList class object.")
    }
    NULL
}

S4Vectors:::setValidity2("AutomorphismList", valid.AutomorphismList)

## ======================== Show AutomorphismList ==================== #

#' @rdname Automorphism
#' @aliases show-AutomorphismList
#' @title Show method for \code{\link{AutomorphismList-class}} object
#' @param object An object from \code{\link{AutomorphismList-class}}.
#' @importFrom methods show
#' @keywords internal
#' @export
setMethod(
    "show",
    signature = "AutomorphismList",
    definition = function(object) {
        nams <- names(object@DataList)
        l <- length(nams)
        if (l > 10) {
            nams <- nams[c(seq(4), l - 2, l - 1, l)]
            nams[4] <- "..."
        }
        cat(class(object), " object of length: ",
            length(object@DataList), "\n",
            sep = ""
        )
        cat(paste0("names(", l, "):"), nams, "\n")
        cat("------- \n")
        gr <- object@SeqRanges
        if (length(gr) > 0 && inherits(object@DataList[[1]], "DataFrame")) {
            mcols(gr) <- object@DataList[[1]]
        } else {
            gr <- object@DataList[[1]]
        }

        print(as(gr, "Automorphism"))
        cat("...\n")
        cat("<", l - 1, " more ",
            class(object@DataList[[1]])[1], " element(s)>\n",
            sep = ""
        )
        cat("Two slots: 'DataList' & 'SeqRanges'\n")
        cat("------- \n")
        invisible(object)
    }
)



## ========================== AutomorphismByCoef ===========================

#' @aliases AutomorphismByCoef
#' @rdname Automorphism
#' @title A class definition to store conserved gene/genomic regions found
#' in a MSA.
#' @keywords internal
#' @export
#' @return AutomorphismByCoef-class definition.
setClass("AutomorphismByCoef",
    contains = "GRanges"
)

# ======================== Validity AutomorphismByCoef ================== #
#' @rdname AutomorphismByCoef
#' @aliases valid.AutomorphismByCoef
#' @title Valid AutomorphismByCoef mcols
#' @param x A 'AutomorphismByCoef object'
#' @importFrom S4Vectors mcols
#' @keywords internal
valid.AutomorphismByCoef <- function(x) {
    coln <- colnames(mcols(x))
    if (!inherits(x, "GRanges") ||
        any(!is.element(c("autm", "cube"), coln))) {
        return("*** This is not a valid AutomorphismByCoef
                class object.")
    } else {
        NULL
    }
}

S4Vectors:::setValidity2("AutomorphismByCoef", valid.AutomorphismByCoef)

## ========================= AutomorphismByCoefList ======================

#' @aliases AutomorphismByCoefList
#' @rdname AutomorphismByCoef
#' @title A class definition for a list of AutomorphismByCoef class objects.
#' @keywords internal
#' @importClassesFrom S4Vectors DataFrame
#' @importFrom methods as
#' @details \strong{AutomorphismByCoefList-class} has the following methods:
#' ## as('from', "AutomorphismByCoefList")
#' Where 'from' is a list of \strong{AutomorphismByCoef-class}.
#'
#' ## unlist(x)
#' Where 'x' is a an \strong{AutomorphismByCoefList-class} object.
#' @export
#' @return AutomorphismByCoefList-class definition.
setClass(
    "AutomorphismByCoefList",
    slots = c(
        elementMetadata = "DataFrame",
        elementType = "character",
        metadata = "list",
        listData = "list"
    ),
    contains = "SimpleGRangesList"
)

as_list_of_AutomorphismByCoef <- function(from) {
    lapply(from, as, Class = "AutomorphismByCoef")
}

setAs("list", "AutomorphismByCoefList", function(from) {
    from <- as_list_of_AutomorphismByCoef(from)
    from <- new_SimpleList_from_list(
        Class = "SimpleGRangesList",
        x = from
    )
    new("AutomorphismByCoefList", from)
})

#' @importClassesFrom GenomicRanges GRangesList
setMethod("unlist",
    signature = "AutomorphismByCoefList",
    function(x) {
        x <- as(x, "GRangesList")
        return(unlist(x))
    }
)

# ===================== Validity AutomorphismByCoefList ================== #
#' @rdname AutomorphismByCoef
#' @aliases valid.AutomorphismByCoefList
#' @title Valid AutomorphismByCoefList mcols
#' @param x A 'AutomorphismByCoefList object'
#' @importFrom S4Vectors mcols
#' @keywords internal
valid.AutomorphismByCoefList <- function(x) {
    if (any(!slapply(x, validObject)) || any(slapply(x, function(y) {
        coln <- colnames(mcols(y))
        !is.element(c("autm", "cube"), coln)
    }))) {
        return("*** This is not a valid AutomorphismByCoefList
                class object.")
    } else {
        NULL
    }
}

S4Vectors:::setValidity2(
    "AutomorphismByCoefList",
    valid.AutomorphismByCoefList
)



## ========================== ConservedRegion ===========================

#' @aliases ConservedRegion
#' @rdname ConservedRegion
#' @title A class definition to store conserved gene/genomic regions found
#' in a MSA.
#' @keywords internal
#' @export
#' @return Definition of the \strong{ConservedRegion-class}.
setClass("ConservedRegion",
    contains = "GRanges"
)

# ======================== Validity ConservedRegion ================== #
#' @aliases valid.ConservedRegion
#' @rdname ConservedRegion
#' @title Valid ConservedRegion mcols
#' @param x A 'ConservedRegion object'
#' @importFrom S4Vectors mcols
#' @keywords internal

valid.ConservedRegion <- function(x) {
    coln <- colnames(mcols(x))
    if (!inherits(x, "GRanges") ||
        any(!is.element(coln, c("autm", "cube")))) {
        return("*** This is not a valid ConservedRegion
                class object.")
    } else {
        NULL
    }
}

S4Vectors:::setValidity2("ConservedRegion", valid.ConservedRegion)

## ========================= ConservedRegionList ======================

#' @aliases ConservedRegionList
#' @rdname ConservedRegion
#' @title A class definition for a list of ConservedRegion class objects.
#' @details \strong{ConservedRegionList-class} has the following method:
#' ## as('from', "ConservedRegionList")
#' Where 'from' is a list of \strong{ConservedRegion-class}.
#'
#' @keywords internal
#' @export
setClass(
    "ConservedRegionList",
    slots = c(
        elementMetadata = "DataFrame",
        elementType = "character",
        metadata = "list",
        listData = "list"
    ),
    contains = "SimpleGRangesList"
)

as_list_of_ConservedRegion <- function(from) {
    lapply(from, as, Class = "ConservedRegion")
}

setAs("list", "ConservedRegionList", function(from) {
    from <- as_list_of_ConservedRegion(from)
    from <- new_SimpleList_from_list(
        Class = "SimpleGRangesList",
        x = from
    )
    new("ConservedRegionList", from)
})

# ===================== Validity ConservedRegionList ================== #
#' @rdname ConservedRegion
#' @aliases valid.ConservedRegion
#' @title Valid ConservedRegionList mcols
#' @param x A 'ConservedRegionList object'
#' @importFrom S4Vectors mcols
#' @keywords internal

valid.ConservedRegionList <- function(x) {
    coln <- colnames(mcols(x))
    if (any(!slapply(x, validObject)) || any(slapply(x, function(y) {
        coln != "autm"
    }))) {
        return("*** This is not a valid ConservedRegionList
                class object.")
    } else {
        NULL
    }
}

S4Vectors:::setValidity2("ConservedRegionList", valid.ConservedRegionList)


## ======================= Show methods =================================


#' @rdname CodonSeq
#' @aliases show-CodonSeq
#' @title Show method for 'CodonSeq' class object
#' @param object An object from 'CodonSeq'.
#' @importFrom methods show
#' @keywords internal
#' @export
setMethod(
    "show",
    signature = "CodonSeq",
    definition = function(object) {
        nams <- names(object@CoordList)
        cat(class(object), " object of length: ",
            length(object@CoordList), "\n",
            sep = ""
        )
        cat(paste0("names(", length(object@CoordList), "):"), nams, "\n")
        cat("------- \n")
        print(.showMatrix(object@CoordList[[1]]))
        cat("...\n")
        cat("<",
            length(object@CoordList) - 1,
            " more ", class(object@CoordList[[1]])[1], " element(s)>\n",
            sep = ""
        )
        cat("Two slots: 'CoordList' & 'SeqRanges'\n")
        cat("------- \n")
        invisible(object)
    }
)

## ========================= Show MatrixList ============================= #

#' @rdname MatrixList
#' @aliases show-MatrixList
#' @title Show method for 'MatrixList' class object
#' @param object An object from 'MatrixList' class
#' @importFrom methods show
#' @keywords internal
#' @export
#' @return Print/show of a MatrixList-class object.
setMethod(
    "show",
    signature = "MatrixList",
    definition = function(object) {
        nams <- names(object@matrices)
        cat(class(object), " object of length: ",
            length(object@matrices), "\n",
            sep = ""
        )
        cat(paste0("names(", length(object@matrices), "):"), nams, "\n")
        cat("------- \n")
        print(.showMatrix(object@matrices[[1]]))
        cat("...\n")
        cat("<",
            length(object@matrices) - 1,
            " more ", class(object@matrices[[1]])[1], " element(s)>\n",
            sep = ""
        )
        cat("Two slots: 'matrices' & 'names'\n")
        cat("------- \n")
        invisible(object)
    }
)

.showMatrix <- function(x) {
    d <- dim(x)
    if (!is.null(d)) {
        cat("Matrix with", d[1], "rows and", d[2], "columns:\n")
        if (d[1] > 10) {
            r <- c()
            for (k in c(seq(5), seq(d[1] - 5, d[1]))) {
                r <- rbind(r, x[k, ])
            }
            r[6, ] <- "..."
        }
        r <- data.frame(r)
        rown <- paste0(c(seq(5), seq(d[1] - 5, d[1])), ":")
        rown[6] <- "..."
        rownames(r) <- rown
    } else {
        l <- length(x)
        cat("Vector of length:", l, "\n")
        if (l > 10) {
            r <- x[c(seq(5), seq(l - 5, l))]
            r[6] <- "..."
            r <- data.frame(matrix(r, ncol = length(r)))
            nms <- paste0(c(seq(5), seq(l - 5, l)), ":")
            nms[6] <- r[6]
            colnames(r) <- nms
        }
    }
    return(r)
}

## =========================== Auxiliary function ======================

make_zero_col_DataFrame <- function(nrow = 0L) {
    stopifnot(isSingleNumber(nrow))
    if (!is.integer(nrow)) {
        nrow <- as.integer(nrow)
    }
    stopifnot(nrow >= 0L)
    new2("DFrame", nrows = nrow, check = FALSE)
}

#' @importFrom methods extends is
new_SimpleList_from_list <- function(Class, x, type, ..., mcols) {
    if (!extends(Class, "SimpleList")) {
        stop("class ", Class, " must extend SimpleList")
    }
    if (!is.list(x)) {
        stop("'x' must be a list")
    }
    if (is.array(x)) {
        tmp_names <- names(x)
        dim(x) <- NULL
        names(x) <- tmp_names
    }
    class(x) <- "list"
    proto <- new(Class)
    if (missing(type)) {
        ans_elementType <- elementType(proto)
    } else {
        ans_elementType <- type
    }
    if (is(S4Vectors::mcols(proto, use.names = FALSE), "DataFrame")) {
        mcols <- make_zero_col_DataFrame(length(x))
    }
    if (!all(slapply(x, function(xi) extends(class(xi), ans_elementType)))) {
        stop(
            "all elements in 'x' must be ", ans_elementType,
            " objects"
        )
    }
    if (missing(mcols)) {
        return(new2(Class, listData = x, ..., check = FALSE))
    }
    new2(Class, listData = x, ..., elementMetadata = mcols, check = FALSE)
}
