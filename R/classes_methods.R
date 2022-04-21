#' An S4 class object for storing the S3 output of hclust, FlashClust or fastcluster
#'
#' @slot merge Matrix as represented in the hclust S3 object
#' @slot height Numeric  as represented in the hclust S3 object
#' @slot order Integer as represented in the hclust S3 object
#' @slot labels NULL as represented in the hclust S3 object
#' @slot method Character as represented in the hclust S3 object
#' @slot call Call represented in the hclust S3 object
#' @slot dist.method NULL as represented in the hclust S3 object
#'
setClass(Class = "hclust", representation(merge="matrix",
                                          height="numeric",
                                          order="integer",
                                          labels="NULL",
                                          method="character",
                                          call="call",
                                          dist.method = "NULL"))

setGeneric(
    "as.S3hclust",
    function(object)
    {
        standardGeneric("as.S3hclust")
    }
)

# Enable the conversion of S4 hclust to S3 hclust
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



#' An S4 class to represent the data used in a weighted-gene co-expression network analysis.
#'
#' @slot exprs A matrix of normalised counts.
#' @slot phenoData A data.frame with variable names (samples, libraries) as rows, description tags (e.g., unit of measurement) as columns.
#' @slot powerEstimate A numeric of length one. An estimate of an appropriate soft-thresholding power calculated usingWGCNA::pickSoftThreshold()
#' @slot fitIndices A data.frame containing the fit indices for scale free topology calculated using WGCNA::pickSoftThreshold().
#' @slot softPower A numeric of length one. The user-selected soft power to be used for calculating the co-expression adjacency matrix.
#' @slot networkType A character vector of length one. A network type. Allowed values are (unique abbreviations of) "unsigned", "signed", "signed hybrid", "distance".
#' @slot adjacencyMat Matrix. A WGCNA adjacency matrix generated with WGCNA::adjacency().
#' @slot dissTOM Matrix. A topological overlap matrix generated with WGCNA::TOMsimilarity().
#' @slot geneTree A S4 class hclust object as defined by tidycoex::hclust().
#'
setClass(Class = "coexList", slots=list(exprs = "matrix",
                                                    phenoData = "data.frame",
                                                    powerEstimate = "numeric",
                                                    fitIndices = "data.frame",
                                                    softPower = "numeric",
                                                    networkType = "character",
                                                    adjacencyMat = "matrix",
                                                    dissTOM = "matrix",
                                                    geneTree = "hclust"))

#' An S4 class to represent the data used in a weighted-gene co-expression network analysis.
#'
#' @param exprs A matrix of normalised gene expression data.
#'
#' @return A coexList class object.
#'
#' @export
#' @examples
#' ngenes <- 1000
#' nsamples <- 16
#' edat <- matrix(rnorm(ngenes*nsamples,mean=5,sd=2),ngenes,nsamples)
#' rownames(edat) <- 1:ngenes
#' cl <- coexList(exprs = edat)
#' cl
coexList <- function(exprs){

    methods::new(Class = "coexList",
                 exprs = exprs)

}

setMethod("show",
          "coexList",
          function(object) {
              message("A coexList object")
          }
)
