#' @export
setGeneric(
    "as.S3hclust",
    function(object)
    {
        standardGeneric("as.S3hclust")
    }
)

## coexList generics

#' @export
setGeneric("getPowerEstimate", function(x, ...) standardGeneric("getPowerEstimate"))

#' @export
setGeneric("isFiltered", function(x, ...) standardGeneric("isFiltered"))

#' @export
setGeneric("normCounts<-", function(x, ..., value) standardGeneric("normCounts<-"))
