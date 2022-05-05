#' Plot dendrogram of samples from normalised data with colourbar
#'
#' @param cl An object of class coexList.
#' @param groupVar Character or numeric that refers to a
#' column name of colData(cl)
#' @param corMethod Correlation method. See ?cor
#' @param hclustMethod clustering method. See ?hclust
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
#' plotDendrogram(cl, colVar=1)
#'

plotDendrogram <- function(cl, colVar=NA, corMethod="spearman",
                           hclustMethod="average"){

    # Do the clustering
    cat("Calculating correlations...\n")
    dissimilarity <- 1 - cor(cl@normCounts, method = corMethod)
    cat("Clustering...\n")
    distance <- as.dist(dissimilarity)
    hc <- hclust(distance, method = hclustMethod)
    gg1 <- suppressMessages(
        ggdendro::ggdendrogram(hc, rotate = FALSE, labels = TRUE) +
            ggplot2::scale_x_discrete(expand=c(0,0)) +
            ggplot2::scale_y_continuous(expand = c(0,0)) +
            ggplot2::theme(axis.text.y = ggplot2::element_text(angle = 0),
                           plot.margin = ggplot2::unit(c(10,10,0,10), "points"))
    )

    # Get the colData for colour bar

    get_col_data <- function(x){

        df <- data.frame(id = rownames(colData(cl)),
                         variable = colnames(colData(cl))[x],
                         value = colData(cl)[ ,x],
                         order = hc$order)

        # Set order of dendrogram
        df$id <- factor(df$id, levels = hc$labels[hc$order])

        return(df)
    }

    col_dat <- get_col_data(x = colVar)

    # Plot the colour bar
    gg2 <- ggplot2::ggplot(col_dat, aes(x = .data$id,
                                        y = .data$variable,
                                        fill=.data$value)) +
        ggplot2::geom_tile() +
        ggplot2::scale_x_discrete(expand = c(0,0)) +
        ggplot2::scale_y_discrete(expand=c(0,0)) +
        ggplot2::scale_fill_viridis_d() +
        ggplot2::xlab("") +
        ggplot2::ylab("") +
        ggplot2::theme(plot.margin = unit(c(0,10,10,10), "points"),
                       axis.ticks = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_blank(),
                       axis.text.y = element_text(angle = 0),
                       legend.title = element_blank(),
                       legend.position="bottom")

    # Assemble the plot
    cp <- cowplot::plot_grid(gg1, gg2, axis = "l", align = "v",
                             nrow = 2,
                             rel_heights = c(0.4, 0.2))

    return(cp)
}
