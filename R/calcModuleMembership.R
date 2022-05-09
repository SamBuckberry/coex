
#' Calculate co-expression module membership
#'
#' @param cl An object of class coexList.
#' @param ... Parameters for WGCNA::signedKME
#' @return A coexList object. This function adds data to the 'moduleMembership'
#' slot of a coexList object.
#'
#' @seealso \code{\link[WGCNA]{signedKME}}
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
#' cl <- calcModuleEigengenes(cl)
#' cl <- calcModuleMembership(cl)
#' mm <- getModuleMembership(cl)
#' mm[1:3, 1:3]

calcModuleMembership <- function(cl, ...){

    stopifnot("cl must be a coexList object" = class(cl)[1] == "coexList")

    cat("=== Running WGCNA::signedKME ===\n")
    kME <- WGCNA::signedKME(datExpr=t(cl@normCounts),
                                            datME=cl@moduleEigengenes$eigengenes,
                                            outputColumnName = "", ...)

    stopifnot(all(rownames(kME) == rownames(cl@normCounts)))

    cl@moduleMembership <- kME

    return(cl)
}



