########3
# To Do

#- Add stopifnot and unit tests for plotPCA
#- Check if plotPCA should be a generic like in DEseq2

#- Add spqn normalisation methods. See https://bioconductor.org/packages/release/bioc/vignettes/spqn/inst/doc/spqn.html

#- Add CPM/TPM normalisation methods
########

library(devtools)
load_all()

# Packages used by tidycoex
usethis::use_package("WGCNA")
usethis::use_package("edgeR")
usethis::use_package("recount")
usethis::use_package("magrittr")
usethis::use_package("stringr")
usethis::use_package("minet")
usethis::use_package("fastcluster")
usethis::use_package("matrixStats")
usethis::use_package("methods")
usethis::use_package("ggplot2")
usethis::use_package("dynamicTreeCut")
usethis::use_package("SummarizedExperiment")
usethis::use_package("caret")
usethis::use_testthat()
#usethis::use_vignette("my-vignette")

# Files to ignore
usethis::use_build_ignore("setup.R")

devtools::document()

devtools::check()

devtools::install()

#######


