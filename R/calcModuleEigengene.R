
#' Calculate an adjacency matrix
#' @param cl An object of class coexList.
#' @param deepSplit integer in the range 0 to 4.
#' Provides a rough control over sensitivity to cluster splitting.
#' The higher the value, the more and smaller clusters will be produced.
#' See ?dynamicTreeCut::cutreeHybrid for more details.
#' @param minClusterSize Numeric. The minimum cluster size.
#' See ?cutreeHybrid::dynamicTreeCut
#' @param pamStage logical. Default is TRUE. See ?cutreeHybrid::dynamicTreeCut
#' @param ... Parameters for cutreeHybrid::dynamicTreeCut
#' @return A coexList object.
#'
#' @export
#' @examples
#' ngenes <- 1000
#' nsamples <- 16
#' edat <- matrix(rpois(ngenes*nsamples, lambda=5), nrow=ngenes)
#' rownames(edat) <- paste0("gene_", 1:ngenes)
#' cl <- coexList(counts = edat)
#' cl <- normCounts(cl)
#' cl <- calcAdjacency(cl, softPower=6)
#' # To view adjacency matrix
#' cl@adjacencyMat[1:6, 1:6]

calcModuleEigengenes <- function(cl, deepSplit=2, minClusterSize=10,
                                pamStage=TRUE, ...){

    stopifnot("cl must be a coexList object" = class(cl)[1] == "coexList")
    stopifnot("deepSplit must be integer between 0 and 4" =
                  is.numeric(deepSplit) & deepSplit < 5 & deepSplit > 0)
    stopifnot("minClusterSize must be numeric > 1" =
                  is.numeric(minClusterSize) & minClusterSize > 1)
    stopifnot("pamStage must be logical of length 1" = is.logical(pamStage) &
                  length(pamStage) == 1)

    cat("Detecting clusters in dendrogram...")
    tree <- dynamicTreeCut::cutreeHybrid(dendro = cl@geneTree,
                                         pamStage = pamStage,
                                         minClusterSize = minClusterSize,
                                         deepSplit = deepSplit,
                                         distM = cl@dissTOM, ...)

    cat("Adding modules to rowData...")
    modules <- WGCNA::labels2colors(tree$labels)
    rowData(cl)$module <- modules

    cat("Calculating module eigengenes...")
    ME <- WGCNA::moduleEigengenes(t(cl@normCounts), colors=modules)

    cl@moduleEigengenes <- ME

    return(ME)
}
