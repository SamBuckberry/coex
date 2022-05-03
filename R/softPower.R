#' Analysis of scale free topology for soft-thresholding
#'
#' @description The aim of this function is to help choose an appropriate
#' soft-thresholding power for network construction. Use `plotSoftPower` to
#' inspect results
#'
#' @param cl A coexList class object.
#' @param RsquaredCut Desired minimum scale free topology fitting index R^2.
#' @param ... Arguments passed to `WGCNA:::pickSoftThreshold()`
#' @return A coexList class object.
#'
#' @import SummarizedExperiment
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
#' cl <- calcSoftPower(cl)
#' cl@powerEstimate
#' @seealso WGCNA::pickSoftThreshold

calcSoftPower <- function(cl, RsquaredCut = 0.85, ...){

    stopifnot("RsquaredCut must be a numeric of length 1, with a value between 0-1." =
                  class(RsquaredCut) == "numeric" & length(RsquaredCut) == 1 &
                  RsquaredCut > 0 & RsquaredCut <=1)

    # set options
    options(stringsAsFactors = FALSE)

    # Check inputs
    stopifnot("cl is not of class coexList" = class(cl) == "coexList")

    # Set powers
    powers <- c(c(1:10), seq(from = 12, to=20, by=2))

    # Call the network topology analysis function
    sft <- WGCNA::pickSoftThreshold(t(cl@normCounts),
                                    RsquaredCut = RsquaredCut,
                                    powerVector=powers,
                                    verbose=1, ...)

    cl@fitIndices <- sft$fitIndices
    cl@powerEstimate <- sft$powerEstimate

    cat(paste0("Estimated soft power for network construction: ",
                   cl@powerEstimate))
    cat("Run plotSoftPower to inspect results.")

    return(cl)
}

#' Plot scale free topology and mean connectivity for soft power thresholds.
#'
#' @description Plot the results from `calcSoftPower`.
#' This function returns two plots: 1) Scale free topology model fit,
#' and 2) Mean connectivity against soft power thresholds.
#'
#' @param cl An object of class coexList.
#' @param RsquaredCut Desired minimum scale free topology fitting index R^2.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' ngenes <- 1000
#' nsamples <- 16
#' edat <- matrix(rpois(ngenes*nsamples, lambda=5), nrow=ngenes)
#' rownames(edat) <- paste0("gene_", 1:ngenes)
#' cl <- coexList(counts = edat)
#' cl <- normCounts(cl)
#' cl <- calcSoftPower(cl)
#' plotSoftPower(cl)
#'
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_point geom_line scale_x_continuous
#' geom_hline xlab ylab theme_bw
#' @importFrom cowplot plot_grid
#' @seealso calcSoftPower

plotSoftPower <- function(cl, RsquaredCut=0.85){

    stopifnot("RsquaredCut must be a numeric of length 1, with a value between 0-1." =
                  class(RsquaredCut) == "numeric" & length(RsquaredCut) == 1 &
                  RsquaredCut > 0 & RsquaredCut <=1)

    plot_dat <- cl@fitIndices
    plot_dat$Topology_Model <- -sign(plot_dat[ ,3])*plot_dat[ ,2]

    gg_si <- ggplot2::ggplot(data = plot_dat,
                    ggplot2::aes(x = .data$Power, y = .data$Topology_Model)) +
        ggplot2::geom_point() +
        ggplot2::geom_line(alpha=0.5) +
        ggplot2::scale_x_continuous(breaks=plot_dat$Power) +
        ggplot2::geom_hline(yintercept = RsquaredCut) +
        ggplot2::xlab(label = "Soft threshold (power)") +
        ggplot2::ylab(label = "Scale free topology model fit, signed R^2") +
        ggplot2::theme_bw()

    gg_mc <- ggplot2::ggplot(plot_dat,
                             ggplot2::aes(x = .data$Power, y = .data$mean.k.)) +
        ggplot2::geom_point() +
        ggplot2::geom_line(alpha=0.5) +
        ggplot2::scale_x_continuous(breaks=plot_dat$Power) +
        ggplot2::xlab(label = "Soft threshold (power)") +
        ggplot2::ylab(label = "Mean connectivity") +
        ggplot2::theme_bw()

    cp <- cowplot::plot_grid(gg_si, gg_mc, nrow = 1, ncol = 2)
    return(cp)
}

