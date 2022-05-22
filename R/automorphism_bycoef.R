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


#' @aliases automorphism_bycoef
#' @rdname automorphism_bycoef
#' @title Autmorphism Grouping by Coefficient
#' @description Automorphisms with the same automorphism's coefficients are
#' grouped.
#' @param x An \code{\link{Automorphism-class}} or an
#' \code{\link{AutomorphismList-class}} object returned by function
#' \code{\link{automorphisms}}.
#' @param mut.type Logical. Whether to include the mutation type as given by
#' function \code{\link{mut_type}}.
#' @param ... Not in use.
#' @return An \code{\link{AutomorphismByCoef}} class object. A coefficient with
#' 0 value is assigned to mutational events that are not automorphisms, e.g.,
#' indel mutations.
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @export
#' @examples
#' ## Load dataset
#' data(autm, package = "GenomAutomorphism")
#'
#' automorphism_bycoef(x = autm[1:2])
setGeneric(
    "automorphism_bycoef",
    function(x,
    ...) {
        standardGeneric("automorphism_bycoef")
    }
)

#' @aliases automorphism_bycoef
#' @rdname automorphism_bycoef
#' @param x An automorphism-class object returned by function
#' \code{\link{automorphisms}}.
#' @importFrom data.table data.table
#' @importFrom dplyr mutate lag %>%
#' @export
setMethod(
    "automorphism_bycoef", signature(x = "Automorphism"),
    function(x, mut.type = TRUE) {
        seq1 <- seq2 <- autm <- cube <- row_names <- NULL
        starts <- lagged <- NULL
        
        idx <- which(x$cube != "Trnl")
        x <- x[ idx ]
        nams <- names(x)
        x <- data.frame(x)
        
        if (mut.type) {
            x$mut_type <- slapply(
                seq_len(nrow(x)),
                function(k) mut_type(x$seq1[k], x$seq2[k])
            )
        }
        
        if (!is.null(nams)) {
            x$row_names <- nams
        }

        x$autm[which(is.na(x$autm))] <- 0
        x <- x %>% mutate(lagged = lag(autm)) %>% 
                mutate(starts = (autm != lagged))
        x$starts[1] <- TRUE
        x <- x %>% mutate(idx = cumsum(starts))
        x$idx <- as.factor(x$idx)

        x <- data.table(data.frame(x))
        if (mut.type)
            keys <- c("idx", "mut_type")
        else
            keys <- "idx"
        
        if (!is.null(nams)) {
            x <- x[, list(
                        seqnames = unique(seqnames),
                        start = min(start),
                        end = max(end),
                        strand = unique(strand),
                        seq1 = unique(seq1),
                        seq2 = unique(seq2),
                        autm = unique(autm),
                        cube = unique(cube),
                        row_names = unique(row_names)),
                    by = keys ]
            x <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
            names(x) <- x$row_names
        } else {
            x <- x[, list(
                        seqnames = unique(seqnames),
                        start = min(start),
                        end = max(end),
                        strand = unique(strand),
                        seq1 = unique(seq1),
                        seq2 = unique(seq2),
                        autm = unique(autm),
                        cube = unique(cube)),
                    by = keys
                ]
            x <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
        }
        x <- x[, c("seq1", "seq2", "autm", "mut_type", "cube")]
        x <- sortByChromAndStart(x)
        return(as(x, "AutomorphismByCoef"))
    }
)



#' @aliases automorphism_bycoef
#' @rdname automorphism_bycoef
#' @param min.len Minimum length of a range to be reported.
#' @param num.cores,tasks Integers. Argument \emph{num.cores} denotes the
#' number of cores to use, i.e. at most how many child processes will be run
#' simultaneously (see \code{\link[BiocParallel]{bplapply}} function from
#' BiocParallel package). Argument \emph{tasks} denotes the number of tasks per
#' job. value must be a scalar integer >= 0L. In this documentation a job is
#' defined as a single call to a function, such as
#' \code{\link[BiocParallel]{bplapply}}. A task is the division of the \eqn{X}
#' argument into chunks. When tasks == 0 (default), \eqn{X} is divided as evenly
#' as possible over the number of workers (see
#' \code{\link[BiocParallel]{MulticoreParam}} from BiocParallel package).
#' @param verbose logic(1). If TRUE, enable progress bar.
#' @importFrom GenomicRanges GRangesList
#' @importFrom parallel detectCores
#' @importFrom BiocParallel MulticoreParam bplapply SnowParam
#' @importFrom data.table data.table
#' @export
setMethod(
    "automorphism_bycoef", signature(x = "AutomorphismList"),
    function(x,
    min.len = 1L,
    mut.type = TRUE,
    num.cores = detectCores() - 1,
    tasks = 0L,
    verbose = TRUE) {
        gr <- try(x@SeqRanges, silent = TRUE)
        if (!inherits(gr, "try-error")) {
            x <- getAutomorphisms(x)
        }

        ## ---------------- Setting parallel distribution --------- ##

        progressbar <- FALSE
        if (verbose) progressbar <- TRUE
        if (Sys.info()["sysname"] == "Linux") {
            bpparam <- MulticoreParam(
                workers = num.cores, tasks = tasks,
                progressbar = progressbar
            )
        } else {
            bpparam <- SnowParam(
                workers = num.cores, type = "SOCK",
                progressbar = progressbar
            )
        }

        ## ------------------------------------------------------- ##

        if (length(gr) > 0) {
            nams <- names(x)
            x <- as.list(x)
            names(x) <- nams
            x0 <- try(bplapply(
                x,
                automorphism_bycoef,
                mut.type = mut.type,
                BPPARAM = bpparam
            ),
            silent = TRUE
            )

            if (inherits(x0, "try-error")) {
                x <- lapply(
                    x,
                    automorphism_bycoef,
                    mut.type = mut.type
                )
            } else {
                x <- x0
                rm(x0)
            }
        }

        idx <- which(slapply(x, function(x) length(x) > min.len))
        x <- x[idx]

        return(as(x, "AutomorphismByCoefList"))
        return(x)
    }
)
