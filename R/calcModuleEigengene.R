

calcModuleEigengenes <- function(cl, deepSplit=2, minClusterSize=10,
                                pamStage=TRUE, ...){

    stopifnot("cl must be a coexList object" = class(cl)[1] == "coexList")

    cat("Detecting clusters in dendrogram...")
    tree <- dynamicTreeCut::cutreeHybrid(dendro = cl@geneTree,
                                         pamStage=pamStage,
                                         minClusterSize = minClusterSize,
                                         deepSplit = deepSplit,
                                         distM = cl@dissTOM, ...)

    cat("Adding modules to rowData...")
    modules <- WGCNA::labels2colors(tree$labels)
    rowData(cl)$module <- modules

    cat("Calculating module eigengenes...")
    ME <- WGCNA::moduleEigengenes(t(cl@normCounts), colors=modules)

    cl@moduleEigengenes <- ME

    return(ME)
}
