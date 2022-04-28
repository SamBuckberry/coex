
# Enable the conversion of S4 hclust to S3 hclust
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
