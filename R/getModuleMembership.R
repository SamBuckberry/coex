#' Get co-expression module membership data from CoexList object
#'
#' @param cl An object of class CoexList.
#' @return data.frame
#'
#' @seealso \code{\link[WGCNA]{signedKME}}
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
#' cl <- calcModuleEigengenes(cl)
#' cl <- calcModuleMembership(cl)
#' mm <- getModuleMembership(cl)
#' mm[1:3, 1:3]
#'
getModuleMembership <- function(cl){

    stopifnot("cl must be a CoexList object" = class(cl)[1] == "CoexList")

    return(cl@moduleMembership)

    }
