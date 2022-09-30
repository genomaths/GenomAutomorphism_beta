#' Reexport useful functions to be available to users
## From S4Vectors ---------------------------------------

#' @importFrom S4Vectors mcols
#' @export
S4Vectors::mcols

#' @importFrom S4Vectors mcols<-
#' @export
S4Vectors::`mcols<-`

#' @importFrom S4Vectors setValidity2
#' @export
S4Vectors::setValidity2

## From Biostrings ---------------------------------------
#' @importFrom Biostrings DNAStringSet
#' @export
Biostrings::DNAStringSet

#' @importFrom Biostrings readDNAMultipleAlignment
#' @export
Biostrings::readDNAMultipleAlignment

#' @importFrom Biostrings getGeneticCode
#' @export
Biostrings::getGeneticCode

#' @importFrom Biostrings translate
#' @export
Biostrings::translate

#' @aliases translate
setMethod(translate, signature = "character", 
    function(
        x,
        genetic.code = getGeneticCode("1")) {
        
        if (length(x) == 1) {
            if (nchar(x) %% 3 != 0)
                stop("*** Argument 'x' must be a character vector or ",
                    "a character coercible to a character vector")
            x <- base2codon(x)
        }
        
        if (length(x) > 1) {
            if (all(nchar(x) == 1)) {
                x <- paste(x, collapse = "")
                x <- base2codon(x)
            }
        } 
        
        if (any((nchar(x) %% 3) != 0)) 
            stop("*** The number of characters in argument 'x' ,
                    must multiple of 3")
        
        x <- toupper(x)
        x <- gsub("U", "T", x)
        
        aa <- genetic.code[match(x, names(genetic.code))]
        aa[is.na(aa)] <- "-"
        return(aa)
    }
)

## From BiocGenerics ---------------------------------------
#' @importFrom BiocGenerics width
#' @export
BiocGenerics::width

#' @importFrom BiocGenerics start
#' @export
BiocGenerics::start

#' @importFrom BiocGenerics start<-
#' @export
BiocGenerics::`start<-`

#' @importFrom BiocGenerics end
#' @export
BiocGenerics::end

#' @importFrom BiocGenerics end<-
#' @export
BiocGenerics::`end<-`

#' @importFrom BiocGenerics strand
#' @export
BiocGenerics::strand

#' @importFrom BiocGenerics strand<-
#' @export
BiocGenerics::`strand<-`



