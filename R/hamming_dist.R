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

#' @rdname hamming_dist
#' @title Hamming Distance
#' @param m1,m2 A \code{\link[IRanges]{AtomicList-class}} (NumericList), a
#' \code{\link{encNumSignalList-class}} or a data.frame where each row
#' contains a vector of 0s and 1s, and the row-names are the sequence names.
#' By default m2 is NULL.
#' @description This function compute de Hamming distance between binary
#' vectors.
#' @details A pairwise distance between vectors from m1 and m2 is computed when
#' both arguments, m1 and m2,  are provided. Otherwise, when m2 = NULL, the
#' distance matrix for vectors from m1 is computed.
#' @param return_matrix logical(1). Whether to return a lower distance matrix
#' when m2 = NULL. If TRUE, the diagonal and the upper matrix are fielded with
#' zeros and with lower matrix with the distance values. If FALSE, then a vector
#' with the distance values for the matrix cells with coordinates ((j + 1):l, j)
#' is returned (j = 1, ... l - 1, and l = nrow(m1)).
#' @param num.cores,tasks Parameters for parallel computation using package
#' \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
#' use, i.e. at most how many child processes will be run simultaneously
#' (see \code{\link[BiocParallel]{bplapply}} and the number of tasks per job
#' (only for Linux OS).
#' @param verbose if TRUE, prints the function log to stdout.
#' @return Depending on the provided data, it will return a lower distance
#' matrix, a vector with the distance values for the matrix cells with
#' coordinates ((j + 1):l, j) (j = 1, ... l - 1, and l = nrow(m1)) or a pairwise
#' distance vector between each vector k from m1 and the corresponding vector k
#' from m2.

#' @importFrom BiocParallel MulticoreParam bplapply SnowParam
#' @export
#' @author Robersy Sanchez (\url{https://genomaths.com}).
#' @examples
#' ### Generating a matrix carrying three binary sequence of 0s and 1s.
binary_seq <- matrix(sample(c(0,1), 750, replace = TRUE),
                      nrow = 3)
rownames(binary_seq) <- c("s1", "s2", "s3")
#'
#' ### Compute the distance matrix
#' hamming_dist(m1 = binary_seq, return_matrix = TRUE)
#'
#'
#' ### Compute the distance matrix and return the lower diagonal matrix
#' ### as a vector
#' hamming_dist(m1 = binary_seq, return_matrix = FALSE)

hamming_dist <- function(m1, m2 = NULL, return_matrix = FALSE, num.cores = 1L,
                         tasks = 0L, verbose = TRUE) {



if (!is.matrix(m1) && !inherits(m1, "SignalEncList")) {
    m1 <- try(as.matrix(m1), silent = TRUE)
    if (inherits(m1, "try-error"))
        stop("\n*** Argument for 'm1' cannot be coersed to a matrix")

    if (!all(is.element(na.omit(unique(as.vector(m1))), c(0,1))))
        stop("\n*** Not all vectors are binary vectors of 0s and 1s")
}

if (is.null(m2)) {
    l <- nrow(m1)
    if (return_matrix) {
        d <- matrix(0, nrow = l, ncol = l)
        for (j in seq_len(l - 1)) {
            for (k in (j + 1):l) {
              d[k,j] <- h_dist(m1[k,], m1[j,])
            }
        }
        rownames(d) <- nams1
        colnames(d) <- nams1
    } else {
        d <- numeric(length = l * (l - 1)/2)
        nams <- character(length = l * (l - 1)/2)
        i <- 1
        for (j in seq_len(l - 1)) {
            for (k in (j + 1):l) {
                d[i] <- h_dist(m1[k,], m1[j,])
                nams[i] <- paste(nams1[k], nams1[j], sep = "_")
                i <- i + 1
            }
        }
        names(d) <- nams
    }
} else {
    if (inherits(m2, c("NumericList", "encNumSignalList",
                      "data.frame", "matrix"))) {
        if (!all(dim(m1) == dim(m2)))
            stop("\n*** Datasets 'm1' and 'm2' must have the same dimensions")

        if (inherits(m2, "NumericList") || inherits(m2, "encNumSignalList")) {
            nams2 <- names(m2)
            valid_object <- TRUE
        }

        if (inherits(m2, "data.frame") || inherits(m2, "matrix")) {
            nams2 <- rownames(m2)
            valid_object <- TRUE
        }

        if (!valid_object)
            stop("\n*** Argument 'm2' must be an object from class:",
                " 'NumericList', 'encNumSignalList', 'data.frame' or 'matrix'")

        if (!is.matrix(m2) && !inherits(m2, "SignalEncList")) {
            m2 <- try(as.matrix(m2), silent = TRUE)
            if (inherits(m2, "try-error"))
                stop("\n*** Argument for 'm2' cannot be coersed to a matrix")

            if (!all(is.element(na.omit(unique(as.vector(m2))), c(0,1))))
                stop("\n*** Not all vectors are binary vectors of 0s and 1s")
        }

    }

    progressbar <- FALSE

    if (verbose) progressbar <- TRUE
    if (Sys.info()["sysname"] == "Linux") {
        bpparam <- MulticoreParam(workers = num.cores, tasks = tasks,
                                progressbar = progressbar)
    } else {
        bpparam <- SnowParam(workers = num.cores, type = "SOCK",
                            progressbar = progressbar)
    }

    if (inherits(m1, "SignalEncList") && inherits(m2, "SignalEncList")) {
        nams1 <- names(m1)
        nams2 <- names(m2)

        if (num.cores > 1) {
            d <- bplapply(seq_along(m1), function(k) {
                sg <- uniqueGRanges(list(m1[[k]], m2[[k]]),
                                    columns = 3L,
                                    missing = "000",
                                    ignore.strand = TRUE,
                                    type = "equal",
                                    verbose = FALSE)
                sg <- getBinary(sg)
                return(h_dist(sg[, 1], sg[, 2]))
            }, BPPARAM = bpparam)
            d <- unlist(d)
        } else {
            dists <- numeric(length(m1))
            for (k in seq_along(m1)) {
                cat("* Processing sample pair #", k, "\n")
                sg <- uniqueGRanges(list(m1[[k]], m2[[k]]),
                                    columns = 3L,
                                    missing = "000",
                                    ignore.strand = TRUE,
                                    type = "equal",
                                    verbose = verbose)
                sg <- getBinary(sg)
                dists[k] <- h_dist(sg[, 1], sg[, 2])
            }
        }

    }
    else {
        d <- bplapply(seq_len(nrow(m1)), function(k) h_dist(m1[k,], m2[k,]),
                    BPPARAM = bpparam)
    }

    if (any(nams1 != nams2)) names(d) <- paste(nams1, nams2, sep = "_")
    else names(d) <- nams1
}
    return(d)
}

### ========== Auxiliary functions ===============
h_dist <- function(x, y) sum((x + y) %% 2)


# =================== get encoding =========== #
getBinary <- function(sg) {
  sg1 <- sg$encoding
  sg2 <- sg$encoding.1

  sg1 <- as.integer(strsplit(paste0(sg1, collapse = ""), split = "")[[1]])
  sg2 <- as.integer(strsplit(paste0(sg2, collapse = ""), split = "")[[1]])
  return(cbind(sg1, sg2))
}
