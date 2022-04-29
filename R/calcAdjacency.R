
#' Calculate an adjacency matrix
#' @param cl An object of class coexList.
#' @param method Character. The method to construct matrix.
#' Current options are 'wgcna', clr'.
#' @param softPower Integer. The soft thresholding power used to construct the
#' WGCNA network. This can be determined using the determineSoftPowerWGCNA function.
#' @param TOM Logical. Calculate Topological overlap matrix for clustering?
#' Default = TRUE
#' @param networkType A character vector of length one. A network type.
#' Allowed values are (unique abbreviations of) "unsigned", "signed",
#' "signed hybrid", "distance"
#' @param corFun Character. The correlation function to be used.
#' Can be one of pearson (default), or spearman.
#' @return A coexList object.
#'
#' @description A wrapper function for constructing a co-expression matrix
#' @export
#' @examples
#' ngenes <- 1000
#' nsamples <- 16
#' edat <- matrix(rnorm(ngenes*nsamples,mean=5,sd=2),ngenes,nsamples)
#' rownames(edat) <- 1:ngenes
#' cl <- coexList(counts = edat)
#' cl <- calcSoftPower(cl)
#' cl <- calcAdjacency(cl)
#' cl@adjacencyMat[1:6, 1:6]

calcAdjacency <- function(cl, method = "wgcna",
                          softPower=cl@powerEstimate, TOM=TRUE,
                          networkType="signed", corFun="pearson"){

    # Check inputs
    stopifnot("method must be a character of length 1. Choose one of 'wgcna', 'clr'" =
                  length(method) == 1 & class(method) == "character")

    stopifnot("method not recognised. Choose one of 'wgcna', 'clr'" =
                  method %in% c("wgcna", "clr"))

    stopifnot("corFun must be one on 'pearson', 'spearman'" =
                  (corFun == "pearson") | (corFun == "spearman"))

    options(stringsAsFactors = FALSE)

    msg <- NULL

    ### calculate adjacency matrix
    if (method == "wgcna"){
        cat("=== Running WGCNA::adjacency ===\n")
        cl@adjacencyMat <- WGCNA::adjacency(datExpr = t(SummarizedExperiment::assay(cl)),
                                            corOptions = list(use = 'p', method = corFun),
                                            power=softPower,
                                            type=networkType)
        diag(cl@adjacencyMat) <- 0

    } else if (method == "clr"){
        cat("=== Running minet::clr ===\n")
        cl@adjacencyMat <- minet::minet(dataset = t(SummarizedExperiment::assay(cl)),
                                        estimator = corFun, method = "clr")
    }

    #### Calculation of the topological overlap matrix
    if (TOM == TRUE){
        cat("=== Running WGCNA::TOMsimilarity ===\n")
        cl@dissTOM <- 1 - WGCNA::TOMsimilarity(cl@adjacencyMat,
                                               TOMType=networkType)
    }

    cat("=== Done! ===")
    return(cl)
}


