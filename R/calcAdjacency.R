
#' Calculate an adjacency matrix
#' @param cl An object of class coexList.
#' @param method Character. The method to construct matrix.
#' Current options are 'wgcna', clr'. For 'wgcna', the function used is WGCNA::adjacency.
#' For 'clr', the functio used is
#' @param networkType character. The type of WGCNA network. Default is 'signed'.
#' See ?WGCNA::adjacency for options.
#' @param corFun character. The correlation method used.
#' Options are "pearson" or "spearman".
#' @param softPower Integer. The soft thresholding power used to construct the
#' WGCNA network. This can be determined using the determineSoftPowerWGCNA function.
#' @param ... Paramaters passed to adjacency function defined in method.
#' @return A coexList object.
#'
#' @description A wrapper function for constructing a co-expression matrix
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

calcAdjacency <- function(cl, method = "wgcna",
                          softPower=cl@powerEstimate,
                          networkType="signed", corFun="pearson", ...){

    # Check inputs
    stopifnot("method must be a character of length 1. Choose one of 'wgcna', 'clr'" =
                  length(method) == 1 & class(method) == "character")

    stopifnot("method not recognised. Choose one of 'wgcna', 'clr'" =
                  method %in% c("wgcna", "clr"))

    stopifnot("corFun must be one on 'pearson', 'spearman'" =
                  (corFun == "pearson") | (corFun == "spearman"))

    options(stringsAsFactors = FALSE)

    # Clear the adjcencyMat slot
    cl@adjacencyMat <- matrix(0,0,0)

    ### calculate adjacency matrix
    calc_start <- Sys.time()
    gc()
    if (method == "wgcna"){
        cat("=== Running WGCNA::adjacency ===\n")
        cat("...This may take a while...\n")
        cl@adjacencyMat <- WGCNA::adjacency(datExpr = t(cl@normCounts),
                                            corFnc = WGCNA::cor,
                                            type = networkType,
                                            power=softPower, ...)
        diag(cl@adjacencyMat) <- 0

    } else if (method == "clr"){
        cat("=== Running minet::clr ===\n")
        cat("...This may take a while...\n")
        cl@adjacencyMat <- minet::minet(dataset = t(cl@normCounts),
                                        estimator = "pearson",
                                        method = "clr", ...)
    }

    calc_end <- Sys.time() - calc_start
    cat(paste0("Elapsed time: ",
               round(as.numeric(calc_end), digits = 2),
               " ", units(calc_end), "\n"))

    return(cl)
}


