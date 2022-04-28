
#' An S4 class to represent the data used in a weighted-gene co-expression network analysis.
#'
#' @param counts A matrix of gene expression data.
#'
#' @return A coexList class object.
#'
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @examples
#' ngenes <- 1000
#' nsamples <- 16
#' edat <- matrix(rnorm(ngenes*nsamples,mean=5,sd=2),ngenes,nsamples)
#' rownames(edat) <- 1:ngenes
#' cl <- coexList(exprs = edat)
#' cl
coexList <- function(counts, ...){
    se <- SummarizedExperiment::SummarizedExperiment(list(counts=counts), ...)
    cl <- .coexList(se)

    # Set filtered to false for new data import
    cl@isFiltered <- FALSE
    return(cl)
}


## Defining the validity method
S4Vectors::setValidity2("coexList", function(object){
    msg <- NULL

    if (assayNames(object)[1] != "counts") {
        msg <- c(msg, "'counts' must be first assay")
    }

    if (is.null(msg)) {
        TRUE
    } else msg
})

## Creating Getter methods

#' @export
setMethod("getPowerEstimate", "coexList", function(x) {
    out <- x@powerEstimate
    out
})

#' @export
setMethod("isFiltered", "coexList", function(x) {
    out <- x@isFiltered
    out
})

## Creating a SHOW method

#' @export
#' @importMethodsFrom SummarizedExperiment show
setMethod("show", "coexList", function(object) {
    callNextMethod()
    cat(
        "exprs has ", ncol((object@exprs)), " columns\n",
        "Data filtered: ", object@isFiltered, "\n",
        "powerEstimate is ", object@powerEstimate, "\n",
        "networkType is", object@networkType, "\n",
        sep=""
    )
})


## Creating Setter methods

#' @export
setReplaceMethod("normCounts", "coexList", function(x, value) {
    x@normCounts <- value
    validObject(x)
    x
})
