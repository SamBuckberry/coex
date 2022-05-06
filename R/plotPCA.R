#' Plot PCA of samples from normalised data
#'
#' @param cl An object of class coexList.
#' @param colVar Character or numeric that refers to a
#' column name of colData(cl)
#' @param x_pc numeric of length 1. The principal component to plot on the x-axis
#' @param y_pc numeric of length 1. The principal component to plot on the y-axis
#' @param scale logical. Should the data be scaled before
#' calculating principal components? See ?prcomp for more details.
#' A column name of colData(cl) for point colours.
#' @return ggplot object
#' @export
#'
#' @examples
#' ngenes <- 1000
#' nsamples <- 16
#' edat <- matrix(rpois(ngenes*nsamples, lambda=5), nrow=ngenes)
#' rownames(edat) <- paste0("gene_", 1:ngenes)
#' cl <- coexList(counts = edat)
#' cl <- normCounts(cl)
#' plotPCA(cl)
#'

plotPCA <- function(cl, colVar="", x_pc=1, y_pc=2, scale=TRUE){

    stopifnot("cl must be a coexList object" = class(cl)[1] == "coexList")
    stopifnot("scale must be logical of length 1" = class(scale) == "logical")
    stopifnot("x_p1 must be numeric of length 1" = class(x_pc) == "numeric" &
                  length(x_pc) == 1)
    stopifnot("y_pc must be numeric of length 1" = class(y_pc) == "numeric" &
                  length(y_pc) == 1)

    # Remove incomplete cases
    mat <- cl@normCounts[stats::complete.cases(cl@normCounts), ]

    # Transpose
    mat <- t(mat)

    # Remove low variance features
    cat("Testing for near zero variance features...\n")
    lowVar <- caret::nearZeroVar(mat, saveMetrics = TRUE)
    cat(paste0("Low variance features removed = ", sum(lowVar$nzv)))
    mat <- mat[ ,!lowVar$nzv]

    # Calculate PC's
    pr <- stats::prcomp(x = mat, scale.=scale)
    pc1 <- round((summary(pr)$importance[2, x_pc] * 100),digits = 1)
    pc2 <- round((summary(pr)$importance[2, y_pc] * 100),digits = 1)

    pc1_dat <- pr$x[ ,x_pc]
    pc2_dat <- pr$x[ ,y_pc]
    #samples <- colnames(cl)

    pca_df <- data.frame(PC1=pc1_dat, PC2=pc2_dat)

    if (colVar != ""){
        pca_df <- cbind(pca_df, colData(cl)[ ,colVar])
        colnames(pca_df)[4] <- "colVar"
    }

    gg_pca <-  ggplot2::ggplot(data = pca_df,
                      mapping = ggplot2::aes(x = .data$PC1, y = .data$PC2,
                                    fill=colVar, colour=colVar)) +
        ggplot2::geom_point(alpha=0.8, size=1) +
        ggplot2::theme_linedraw() +
        ggplot2::theme(panel.grid =
                           ggplot2::element_line(colour = 'grey')) +
        ggplot2::xlab(paste0("PC", x_pc, " (", pc2, "%)")) +
        ggplot2::ylab(paste0("PC", y_pc, " (", pc1, "%)"))

    return(gg_pca)
}
