#' Normalise RNA-seq count data for co-expression analyses
#'
#' @param cl An object of class coexList.
#'
#' @export
#'
#' @examples
#' ngenes <- 1000
#' nsamples <- 16
#' edat <- matrix(rpois(ngenes*nsamples, lambda=5), nrow=ngenes)
#' rownames(edat) <- paste0("gene_", 1:ngenes)
#' cl <- coexList(counts = edat)
#' cl <- normCounts(cl)
#' cl@normCounts[1:6, 1:6]
#'
normCounts <- function(cl){

    # Normalise using the CTF method
    norm_factors <- edgeR::calcNormFactors(assay(cl),
                                    lib.size = colSums(assay(cl)),
                                    method = "TMM")

    dat_norm <- sweep(assay(cl), 2, norm_factors, "/")

    # Scale
    dat_norm <- asinh(dat_norm)

    # Add normalised data to object
    cl@normCounts <- dat_norm

    cl@normMethod <- "CTF"
    cl@scaleMethod <- "asinh"

    return(cl)

}
