
#' Calculate an adjacency matrix
#' @param cl An object of class CoexList.
#' @param deepSplit integer in the range 0 to 4.
#' Provides a rough control over sensitivity to cluster splitting.
#' The higher the value, the more and smaller clusters will be produced.
#' See ?dynamicTreeCut::cutreeHybrid for more details.
#' @param ... Parameters for cutreeHybrid::dynamicTreeCut
#' @return A CoexList object.
#'
#' @export
#' @examples
#' ngenes <- 1000
#' nsamples <- 16
#' edat <- matrix(rpois(ngenes*nsamples, lambda=5), nrow=ngenes)
#' rownames(edat) <- paste0("gene_", 1:ngenes)
#' cl <- CoexList(counts = edat)
#' cl <- normCounts(cl)
#' cl <- calcAdjacency(cl, softPower=6)
#' cl <- calcTOM(cl)
#' cl <- calcTree(cl)
#' cl <- calcModuleEigengenes(cl)
#' str(cl@moduleEigengenes)

calcModuleEigengenes <- function(cl, deepSplit=2, ...){

    stopifnot("cl must be a CoexList object" = class(cl)[1] == "CoexList")
    stopifnot("deepSplit must be integer between 0 and 4" =
                  is.numeric(deepSplit) & deepSplit < 5 & deepSplit > 0)
    #stopifnot("minClusterSize must be numeric > 1" =
    #              is.numeric(minClusterSize) & minClusterSize > 1)
    #stopifnot("pamStage must be logical of length 1" = is.logical(pamStage) &
    #              length(pamStage) == 1)

    cat("=== Running dynamicTreeCut::cutreeHybrid ===\n")
    tree <- dynamicTreeCut::cutreeHybrid(dendro = cl@geneTree,
                                         deepSplit = deepSplit,
                                         distM = cl@dissTOM, ...)

    modules <- WGCNA::labels2colors(tree$labels)
    #rowData(cl)$module <- modules

    cat("=== Running WGCNA::moduleEigengenes ===\n")
    cl@moduleEigengenes <- WGCNA::moduleEigengenes(t(cl@normCounts),
                                                   colors=modules)

    return(cl)
}
