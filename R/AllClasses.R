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
.hclust <- setClass(Class = "hclust", representation(merge="matrix",
                                          height="numeric",
                                          order="integer",
                                          labels="NULL",
                                          method="character",
                                          call="call",
                                          dist.method = "NULL"))

#' An S4 class to represent the data used in a weighted-gene co-expression network analysis.
#'
#' @slot normCounts A matrix of normalised counts.
#' @slot normMethod Character. The normalisation method used.
#' @slot scaleMethod Character. The scaling method used.
#' @slot isFiltered Logical of length one. Indicates if data have been filtered on variance.
#' @slot powerEstimate A numeric of length one. An estimate of an appropriate
#' soft-thresholding power calculated usingWGCNA::pickSoftThreshold()
#' @slot fitIndices A data.frame containing the fit indices for scale free
#' topology calculated using WGCNA::pickSoftThreshold().
#' @slot softPower A numeric of length one. The user-selected soft power
#' to be used for calculating the co-expression adjacency matrix.
#' @slot networkType A character vector of length one. A network type.
#' Allowed values are (unique abbreviations of) "unsigned", "signed",
#' "signed hybrid", "distance".
#' @slot adjacencyMat Matrix. A WGCNA adjacency matrix generated
#' with WGCNA::adjacency().
#' @slot dissTOM Matrix. A topological overlap matrix generated
#' with WGCNA::TOMsimilarity().
#' @slot geneTree A S4 class hclust object as defined by tidycoex::hclust().
#' @slot moduleEigengenes A list of module eigengenes data as returned from
#' WGCNA::moduleEigengenes
#'
#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
.coexList <- setClass(Class = "coexList",
                      contains="SummarizedExperiment",
                      slots=representation(
                          normCounts = "matrix",
                          normMethod = "character",
                          scaleMethod = "character",
                          isFiltered = "logical",
                          powerEstimate = "numeric",
                          fitIndices = "data.frame",
                          softPower = "numeric",
                          networkType = "character",
                          adjacencyMat = "matrix",
                          dissTOM = "matrix",
                          geneTree = "hclust",
                          moduleEigengenes = "list"
                          ),
                      )
