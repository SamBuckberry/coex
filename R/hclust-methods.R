
#' Convert a S4 hclust to S3 hclust
#'
#' @param object A S4 hclust object
#' @return An S3 hclust object
#'
#' @export
setMethod(
    "as.S3hclust",
    signature(object = "hclust"),
    function(object)
    {
        structure(list(merge = object@merge,
                       height = object@height,
                       order = object@order,
                       labels = object@labels,
                       method = object@method,
                       call = object@call,
                       dist.method = object@dist.method),
                  class = "hclust")
    }
)
