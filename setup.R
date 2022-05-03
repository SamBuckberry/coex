library(devtools)
load_all()

# Packages used by tidycoex
usethis::use_package("WGCNA")
usethis::use_package("edgeR")
usethis::use_package("recount")
usethis::use_package("magrittr")
usethis::use_package("minet")
usethis::use_package("fastcluster")
usethis::use_package("matrixStats")
usethis::use_package("methods")
usethis::use_package("ggplot2")
usethis::use_package("SummarizedExperiment")
usethis::use_testthat()
#usethis::use_vignette("my-vignette")
usethis::use_data_raw()
# Files to ignore
usethis::use_build_ignore("setup.R")

devtools::document()

devtools::check()

devtools::install()

#######
