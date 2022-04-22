library(devtools)
load_all()

# Packages used by tidycoex
usethis::use_package("WGCNA")
usethis::use_package("fastcluster")
usethis::use_package("matrixStats")
usethis::use_package("methods")
usethis::use_package("ggplot2")

# Files to ignore
usethis::use_build_ignore("setup.R")

devtools::document()

devtools::check()

devtools::install()

ngenes <- 1000
nsamples <- 16
edat <- matrix(rnorm(ngenes*nsamples,mean=5,sd=2),ngenes,nsamples)
rownames(edat) <- 1:ngenes
cl <- coexList(exprs = edat)
cl <- applyFilter(cl = )

cl <- calcSoftPower(cl)

#' An S4 class to represent the data used in a weighted-gene co-expression network analysis.
#'
#' @param exprs A matrix of normalised counts.
#' @param phenoData A data.frame with variable names (samples, libraries) as rows, description tags (e.g., unit of measurement) as columns.
#' @param powerEstimate A numeric of length one. An estimate of an appropriate soft-thresholding power calculated usingWGCNA::pickSoftThreshold()
#' @param fitIndices A data.frame containing the fit indices for scale free topology calculated using WGCNA::pickSoftThreshold().
#' @param softPower A numeric of length one. The user-selected soft power to be used for calculating the co-expression adjacency matrix.
#' @param networkType A character vector of length one. A network type. Allowed values are (unique abbreviations of) "unsigned", "signed", "signed hybrid", "distance".
#' @param adjacencyMat Matrix. A WGCNA adjacency matrix generated with WGCNA::adjacency().
#' @param dissTOM Matrix. A topological overlap matrix generated with WGCNA::TOMsimilarity().
#' @param geneTree A S4 class hclust object as defined by tidycoex::hclust().
#'
#' @return A coexList class object.
#'
#' @export
#'
coexList <- function(exprs, phenoData=NULL, powerEstimate=NULL, fitIndices=NULL,
                     softPower=NULL, networkType=NULL, adjacencyMat=NULL,
                     dissTOM=NULL, geneTree=NULL){

    methods::new(Class = "coexList",
                 exprs = exprs,
                 phenoData = phenoData,
                 powerEstimate = powerEstimate,
                 fitIndices = fitIndices,
                 softPower = softPower,
                 networkType = networkType,
                 adjacencyMat = adjacencyMat,
                 dissTOM = dissTOM,
                 geneTree = geneTree)

}
