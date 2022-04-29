
#' Convert a S4 hclust to S3 hclust
#'
#' @param object An S4 hclust object
#' @return An S3 hclust object
#'
#' @export
setGeneric(
    "as.S3hclust",
    function(object)
    {
        standardGeneric("as.S3hclust")
    }
)
