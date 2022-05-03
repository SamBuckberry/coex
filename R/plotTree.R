
#' Plot a dendrogram with modules
#'
#' @param cl An object of class coexList.
#' @param ... Parameters for fastcluster::hclust()
#' @return A plot to the display device
#'
#' @importFrom stats as.dist
#' @export
#' @examples
#' ngenes <- 1000
#' nsamples <- 16
#' edat <- matrix(rpois(ngenes*nsamples, lambda=5), nrow=ngenes)
#' rownames(edat) <- paste0("gene_", 1:ngenes)
#' cl <- coexList(counts = edat)
#' cl <- normCounts(cl)
#' cl <- calcAdjacency(cl, softPower=6)
#' cl <- calcTOM(cl)
#' cl <- calcTree(cl)
#' plotTree(cl)
plotTree <- function(cl, ...){

    # Seperate list objects
    mColorh=NULL

    for (ds in 0:4){
        tree <- dynamicTreeCut::cutreeHybrid(dendro = cl@geneTree,
                                             deepSplit = ds,
                                             verbose = 0,
                                             distM = cl@dissTOM, ...)

        mColorh <- cbind(mColorh, WGCNA::labels2colors(tree$labels))
    }


    WGCNA::plotDendroAndColors(cl@geneTree, mColorh,
                               paste("dpSplt =",0:4), main = "",
                               dendroLabels=FALSE)
}
