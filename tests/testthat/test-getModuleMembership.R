
## Setup the toy data for testing
ngenes <- 1000
nsamples <- 16
edat <- matrix(rpois(ngenes*nsamples, lambda=5), nrow=ngenes)
rownames(edat) <- paste0("gene_", 1:ngenes)
cl <- coexList(counts = edat)
cl <- normCounts(cl)
cl <- calcAdjacency(cl, softPower=6)
cl <- calcTOM(cl)
cl <- calcTree(cl)
cl <- calcModuleEigengenes(cl)
cl <- calcModuleMembership(cl)

testthat::test_that("calcModuleMembership returns data.frame", {
  testthat::expect_true(is.data.frame(getModuleMembership(cl)))
})


