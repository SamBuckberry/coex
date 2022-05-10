
#' An S4 class object to represent the data used in a weighted-gene co-expression network analysis.
#'
#' @param counts A matrix of gene expression data.
#' @param ... Parameters for SummarizedExperiment.
#'
#' @return A CoexList class object.
#'
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @examples
#' ngenes <- 1000
#' nsamples <- 16
#' edat <- matrix(rnorm(ngenes*nsamples,mean=5,sd=2),ngenes,nsamples)
#' rownames(edat) <- 1:ngenes
#' CoexList(counts = edat)
#'
CoexList <- function(counts, ...){
    se <- SummarizedExperiment::SummarizedExperiment(list(assays=counts), ...)
    cl <- .CoexList(se)

    # Set filtered to false for new data import
    cl@isFiltered <- FALSE
    return(cl)
}

## Defining the validity method
S4Vectors::setValidity2("CoexList", function(object){
    msg <- NULL

    if (assayNames(object)[1] != "assays") {
        msg <- c(msg, "'assays' must be first assay")
    }

    if (is.null(msg)) {
        TRUE
    } else msg
})

## Creating a SHOW method

#' Show method for when calling a CoexList object
#' @param object A CoexList object
#' @export
#' @importMethodsFrom SummarizedExperiment show
setMethod("show", "CoexList", function(object) {
    callNextMethod()
    cat(
        "Assay data has ", ncol(object), " columns\n",
        "Data filtered: ", object@isFiltered, "\n",
        "powerEstimate is ", object@powerEstimate, "\n",
        "networkType is", object@networkType, "\n",
        sep=""
    )
})
