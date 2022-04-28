
#' Check to see how many genes or samples will be filtered for a given threshold
#'
#' @description A function to identify samples and genes to be filtered for
#' network construction.
#'
#' @param cl An object of class coexList.
#' @param propGenes The proportion of genes to keep. Will rank and select genes
#' with the highest variance.
#' @param ... Arguments passed to `WGCNA::goodSamplesGenes`
#' @export
#'
#' @seealso goodSamplesGenes
#'
#' @examples
#' \dontrun{
#' cl <- coexList(d1)
#' checkFilter(cl)
#' }
checkFilter <- function(cl, propGenes=1, ...){

    stopifnot("Data have already been filtered with applyFilter()" =
                  cl@isFiltered == FALSE)

    good <- WGCNA::goodSamplesGenes(datExpr=t(cl@normCounts),
                                    verbose = 1, ...)

    # Filter genes based in variance
    propGenesCount <- round(propGenes * nrow(cl@normCounts))
    vars <- matrixStats::rowVars(cl@normCounts)
    varKeep <- rank(-vars) <= propGenesCount

    filterGenes <- (good$goodGenes == TRUE) & (varKeep == TRUE)
    filterSamples <- good$goodSamples

    message(paste0("Genes pre-filter: ", nrow(cl@normCounts)))
    message(paste0("Genes post-filter: ", sum(filterGenes)))
    message(paste0("Samples pre-filter: ", ncol(cl@normCounts)))
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

    stopifnot("Data have already been filtered with applyFilter()" =
                  cl@isFiltered == FALSE)

    good <- WGCNA::goodSamplesGenes(datExpr=t(cl@normCounts), verbose = 0)

    # Filter genes based in variance
    propGenesCount <- round(propGenes * nrow(cl@normCounts))
    vars <- matrixStats::rowVars(cl@normCounts)
    varKeep <- rank(-vars) <= propGenesCount

    filterGenes <- (good$goodGenes == TRUE) & (varKeep == TRUE)
    filterSamples <- good$goodSamples

    cl@normCounts <- cl@normCounts[filterGenes, filterSamples]

    # Record the filtering in object
    cl@isFiltered <- TRUE

    return(cl)
}

