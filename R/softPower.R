#' Analysis of scale free topology for soft-thresholding
#'
#' @description The aim of this function is to help choose an appropriate
#' soft-thresholding power for network construction. Use `plotSoftPower` to
#' inspect results
#'
#' @param cl A coexList class object.
#'
#' @return A coexList class object.
#'
#' @export
#'
#' @examples
#' ngenes <- 1000
#' nsamples <- 16
#' edat <- matrix(rnorm(ngenes*nsamples,mean=5,sd=2),ngenes,nsamples)
#' rownames(edat) <- 1:ngenes
#' cl <- coexList(exprs = edat)
#' cl <- calcSoftPower(cl)
#' cl@powerEstimate
#' head(cl@fitIndices)

calcSoftPower <- function(cl){

    # set options
    options(stringsAsFactors = FALSE)

    # Check inputs
    stopifnot("cl is not of class coexList" = class(cl) == "coexList")

    # Set powers
    powers <- c(c(1:10), seq(from = 12, to=20, by=2))

    # Call the network topology analysis function

    sft <- WGCNA::pickSoftThreshold(t(cl@exprs),
                                    powerVector=powers, verbose=5)

    cl@fitIndices <- sft$fitIndices
    cl@powerEstimate <- sft$powerEstimate
    message(paste0("Estimated soft power for network construction: ",
                   cl@powerEstimate))
    message("Run plotSoftPower() to inspect results.")

    return(cl)
}
