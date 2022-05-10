
#' Calculate an adjacency matrix
#' @param cl An object of class CoexList.
#' @param ... Parameters for WGCNA::TOMsimilarity()
#' @return A CoexList object.
#'
#' @description A wrapper function for running the WGCNA::TOMsimilarity function
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
#' # To view TOM
#' cl@dissTOM[1:6, 1:6]

calcTOM <- function(cl, ...){

    cat("=== Running WGCNA::TOMsimilarity ===\n")
    cat("...This might take a while...\n")
    calc_start <- Sys.time()

    cl@dissTOM <- 1 - WGCNA::TOMsimilarity(cl@adjacencyMat, ...)

    calc_end <- Sys.time() - calc_start
    cat(paste0("Elapsed time: ",
               round(as.numeric(calc_end), digits = 2),
               " ", units(calc_end), "\n"))
    return(cl)
}


