#' Reexport useful functions to be available to users
## From S4Vectors ---------------------------------------

#' @importFrom S4Vectors mcols
#' @export
S4Vectors::mcols

#' @importFrom S4Vectors mcols<-
#' @export
S4Vectors::`mcols<-`

## From Biostrings ---------------------------------------
#' @importFrom Biostrings DNAStringSet
#' @export
Biostrings::DNAStringSet

#' @importFrom Biostrings readDNAMultipleAlignment
#' @export
Biostrings::readDNAMultipleAlignment

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
