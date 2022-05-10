#' Normalise RNA-seq count data for co-expression analyses
#'
#' @param cl An object of class CoexList.
#' @param normMethod character. The within-sample normalisation method.
#' One of "CTF", "CPM", "TPM"
#' @param scaleMethod character. One of "log2", "asinh",
#' @param geneLength A numeric vector of length nrow(cl)
#' that is the length of genes for length based normalisation such as TPM.
#' @export
#'
#' @examples
#' ngenes <- 1000
#' nsamples <- 16
#' edat <- matrix(rpois(ngenes*nsamples, lambda=5), nrow=ngenes)
#' rownames(edat) <- paste0("gene_", 1:ngenes)
#' cl <- CoexList(counts = edat)
#' cl <- normCounts(cl)
#' cl@normCounts[1:6, 1:6]
#' cl@normMethod
#' cl@scaleMethod
#'
normCounts <- function(cl, normMethod="CPM", scaleMethod="log2",
                       geneLength=NA){

    # Check CoexList object
    stopifnot("cl must be a CoexList object" = class(cl)[1] == "CoexList")

    # Check normalisation method is valid
    normMethod <- casefold(normMethod, upper = TRUE)
    stopifnot("normMethod must be one of 'CTF', 'CPM', 'TPM'" =
                  normMethod %in% c("CTF", "CPM", "TPM") &
                  length(normMethod) == 1 & class(normMethod) == "character")

    # Check scale method is valid
    scaleMethod <- casefold(scaleMethod, upper = FALSE)
    stopifnot("scaleMethod must be one of 'asinh', 'log2'" =
                  scaleMethod %in% c("asinh", "log2") &
                  length(scaleMethod) == 1 & class(scaleMethod) == "character")

    dat_norm <- NULL

    dat <- assay(cl)

    #################
    # Normalise data
    #################

    # Normalise using the CTF method
    if (normMethod == "CTF"){

        norm_factors <- edgeR::calcNormFactors(dat,
                                               lib.size = colSums(dat),
                                               method = "TMM")

        dat_norm <- sweep(dat, 2, norm_factors, "/")

    }

    # Normalise using the CPM method
    if (normMethod == "CPM"){

        dat_norm <- edgeR::cpm(assay(cl), log=TRUE)

    }

    # Normalise using the TPM method
    if (normMethod == "TPM"){

        rate <- log(dat) - log(geneLength)
        denom <- log(sum(exp(rate)))
        dat_norm <- exp(rate - denom + log(1e6))

    }

    ###############
    # Scale data
    ###############

    # Scale using asinh
    if (scaleMethod == "asinh"){

        dat_norm <- asinh(dat_norm)
    }

    # Scale using log2
    if (scaleMethod == "log2" & normMethod != "CPM"){

        dat_norm <- log2(dat_norm)

    }

    # Add normalised data to object
    cl@normCounts <- as.matrix(dat_norm)

    # Add methods record to CoexList object
    cl@normMethod <- normMethod
    cl@scaleMethod <- scaleMethod

    return(cl)

}
