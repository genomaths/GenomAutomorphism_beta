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

#' @rdname modlineq
#' @aliases modlineq
#' @title Modular System of Linear Equation Solver (MLE)
#' @param a An integer or a vector of integers. 
#' @param b An integer or a vector of integers.
#' @param n An integer or a vector of integers.
#' @description If \eqn{a, b}, and  \eqn{c} are integer vectors, this function 
#' try to find, at each coordinate, the solution of the MLE 
#' \eqn{a x = b}  mod \eqn{n}. If the MLE \eqn{a x = b mod n} has not 
#' solutions (see \code{\link[numbers]{modlin}}), the value reported for the 
#' coordinate will be 0 and the corresponding translation.
#' @details For \eqn{a, b}, and \eqn{c} integer scalars, it is just a 
#' wrapper function to call \code{\link[numbers]{modlin}}. 
#' @importFrom numbers modlin
#' @export
#' @return A numerical vector. If for some coordinate the equation has not 
#' solution in their definition domain it will return 0 for such coordinate.
#' @examples
#' ## The MLE 8x = 54 mod 64 is not solvable; hence the first vector 
#' ## coordinate returns 0.
#' modlineq(a = c(8,9), b = c(54, 34), n = c(64,64))
setGeneric(
    "modlineq",
    function(a, b, n) {
        a <- as.integer(a)
        b <- as.integer(b)
        n <- as.integer(n)
        
        if (length(a) < 2) 
            res <- modlin(a, b, n)
        else {
            if (length(n) == 1) {
                res <- mapply(function(x,y) modl(x, y, n),
                              a, b, USE.NAMES = FALSE)
            }
            else
                res <- mapply(function(x, y, z) {
                    modl(x, y, z)
                }, a, b, n, USE.NAMES = FALSE)
        }
        return(res)
    }
)

## ========= Auxiliary function ===================

modl <- function(x, y, m) {
    if (length(m) == 1 && length(x) > 1)
        m <- rep(m, length(x))
    if (x > 0 && y > 0) {
        res <- modlin(x, y, m)
        if (length(res) > 1) 
            res <- min(res)
    }
    else 
        res <- 0
    if (is.null(res)) 
        res <- 0
    
    return(res)
}



