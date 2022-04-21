
#' Check to see how many genes or samples will be filtered for a given threshold
#'
#' @description A function to identify samples and genes to be filtered for network construction.
#'
#' @param cl An object of class coexList.
#' @param propGenes The proportion of genes to keep. Will rank and select genes with the highest variance.
#' @export
#'
#' @examples
#' \dontrun{
#' cl <- coexList(d1)
#' checkFilter(cl)
#' }
checkFilter <- function(cl, propGenes=1){

    good <- WGCNA::goodSamplesGenes(datExpr=t(cl@exprs), verbose = 0)

    # Filter genes based in variance
    propGenesCount <- round(propGenes * nrow(cl@exprs))
    vars <- matrixStats::rowVars(cl@exprs)
    varKeep <- rank(-vars) <= propGenesCount

    filterGenes <- (good$goodGenes == TRUE) & (varKeep == TRUE)
    filterSamples <- good$goodSamples

    message(paste0("Genes pre-filter: ", nrow(cl@exprs)))
    message(paste0("Genes post-filter: ", sum(filterGenes)))
    message(paste0("Samples pre-filter: ", ncol(cl@exprs)))
    message(paste0("Samples post-filter: ", sum(filterSamples)))

}

#' Filter genes and samples in coexList object for a given threshold
#'
#' @description A function to identify samples and genes to be filtered for network construction.
#'
#' @param cl An object of class coexList.
#' @param propGenes The proportion of genes to keep. Will rank and select genes with the highest variance.
#' @export
#'
#' @examples
#' \dontrun{
#' cl <- coexList(d1)
#' dim(cl@exprs)
#' cl <- applyFilter(cl, propGenes=0.5)
#' dim(cl@exprs)
#' }
applyFilter <- function(cl, propGenes=1){

    good <- WGCNA::goodSamplesGenes(datExpr=t(cl@exprs), verbose = 0)

    # Filter genes based in variance
    propGenesCount <- round(propGenes * nrow(cl@exprs))
    vars <- matrixStats::rowVars(cl@exprs)
    varKeep <- rank(-vars) <= propGenesCount

    filterGenes <- (good$goodGenes == TRUE) & (varKeep == TRUE)
    filterSamples <- good$goodSamples

    cl@exprs <- cl@exprs[filterGenes, filterSamples]

    return(cl)
}

