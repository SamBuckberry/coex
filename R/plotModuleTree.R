
#' Plot a dendrogram with modules
#'
#' @param cl An object of class CoexList.
#' @param deepSplits integer in the range 0 to 4.
#' Provides a rough control over sensitivity to cluster splitting.
#' The higher the value, the more and smaller clusters will be produced.
#' See dynamicTreeCut::cutreeHybrid() for more information.
#' @param cores integer 1-5. How many cores to use?
#' Parameter passed to parallel::mclapply(mc.cores=cores)
#' @param ... Parameters for dynamicTreeCut::cutreeHybrid()
#' @return A plot to the display device
#'
#' @importFrom stats as.dist
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
#' plotModuleTree(cl)
plotModuleTree <- function(cl, deepSplits=c(0:4), cores=1, ...){


    stopifnot("cl must be a CoexList object" = class(cl)[1] == "CoexList")

    cat(paste0("=== Running dynamicTreeCut::cutreeHybrid for deepSplits ",
               min(deepSplits), "-", max(deepSplits), " ===\n"))
    cat("...This might take a while...\n")
    calc_start <- Sys.time()

    cut_tree <- function(x, ...){
        tree <- dynamicTreeCut::cutreeHybrid(dendro = cl@geneTree,
                                             deepSplit = x,
                                             verbose = 0,
                                             distM = cl@dissTOM, ...)

        mColor <- WGCNA::labels2colors(tree$labels)
        return(mColor)
    }

    #mColorh <- parallel::mclapply(X = deepSplits, FUN = cut_tree,
    #                              mc.cores = cores, ...)

    mColorh <- BiocParallel::bplapply(X = deepSplits, FUN = cut_tree, ...)

    mColorh <- do.call(cbind, mColorh)

    colnames(mColorh) <- paste("deepSpilt =",0:4)

    cat("...Plotting dendrogram with deepSplit modules...\n")
    WGCNA::plotDendroAndColors(dendro = cl@geneTree, colors = mColorh,
                               groupLabels = colnames(mColorh),
                               main = "cutreeHybrid",
                               dendroLabels=FALSE)

    calc_end <- Sys.time() - calc_start
    cat(paste0("Elapsed time: ",
               round(as.numeric(calc_end), digits = 2),
               " ", units(calc_end), "\n"))
}
