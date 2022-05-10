
#' Cluster topological overlap matrix using fastcluster::hclust
#'
#' @param cl An object of class CoexList.
#' @param ... Parameters for fastcluster::hclust()
#' @return A CoexList object.
#'
#' @description A fast wrapper function for running hclust on a CoexList object
#' on a CoexList object
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
#' # To view hclust object
#' cl@geneTree


calcTree <- function(cl, ...){

    cat("=== Clustering with fastcluster::hclust ===\n")

    calc_start <- Sys.time()

    dists <- as.dist(cl@dissTOM)
    geneTreeA1 <- fastcluster::hclust(dists, method="average", ...)

    cl@geneTree <- geneTreeA1

    # convert output to S4 object
    #cl@geneTree <- new(Class = "hclust",
    #                   merge=geneTreeA1$merge,
    #                   height=geneTreeA1$height,
    #                   order=geneTreeA1$order,
    #                   method=geneTreeA1$method,
    #                   call=geneTreeA1$call,
    #                   dist.method=geneTreeA1$dist.method)

    calc_end <- Sys.time() - calc_start
    cat(paste0("Elapsed time: ",
               round(as.numeric(calc_end), digits = 2),
               " ", units(calc_end), "\n"))

    return(cl)
}

