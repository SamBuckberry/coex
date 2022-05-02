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
usethis::use_vignette("my-vignette")
# Files to ignore
usethis::use_build_ignore("setup.R")

devtools::document()

devtools::check()

devtools::install()

#######

ngenes <- 1000
nsamples <- 16
edat <- matrix(rnorm(ngenes*nsamples,mean=5,sd=2),ngenes,nsamples)
rownames(edat) <- 1:ngenes
cl <- coexList(counts = edat)
cl <- calcSoftPower(cl)


