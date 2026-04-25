# library(testthat); library(BiocNeighbors); source("setup.R"); source("test-AnnoyParam.R")

test_that("AnnoyParam construction behaves properly", {
    p <- AnnoyParam()
    expect_output(show(p), "AnnoyParam")
    expect_identical(bndistance(p), "Euclidean")

    p <- AnnoyParam(distance="Manhattan")
    expect_identical(bndistance(p), "Manhattan")
})

test_that("Cursory checks for AnnoyParam", {
    Y <- matrix(rnorm(10000), ncol=20)

    p <- AnnoyParam()
    out <- findKNN(Y, k=8, BNPARAM=p)
    expect_identical(ncol(out$distance), 8L)
    expect_identical(ncol(out$index), 8L)
})

test_that("AnnoyParam queries behave with cosine distance", {
    Y <- matrix(rnorm(10000), ncol=20)
    Z <- matrix(rnorm(2000), ncol=20)

    p <- AnnoyParam(distance="Cosine")

    out <- queryKNN(Y, Z, k=8, BNPARAM=p)
    Y1 <- Y/sqrt(rowSums(Y^2))
    Z1 <- Z/sqrt(rowSums(Z^2))
    ref <- queryKNN(Y1, Z1, k=8, BNPARAM=p)
    expect_equal(out, ref)
})

test_that("AnnoyParam accepts the shared threading interface", {
    Y <- matrix(rnorm(10000), ncol=20)
    Z <- matrix(rnorm(2000), ncol=20)

    idx <- buildIndex(Y, BNPARAM=AnnoyParam(), BPPARAM=BiocParallel::SnowParam(2))
    out <- queryKNN(idx, Z, k=5, BPPARAM=BiocParallel::SnowParam(2))

    expect_s4_class(idx, "AnnoyIndex")
    expect_identical(dim(out$index), c(nrow(Z), 5L))
    expect_identical(dim(out$distance), c(nrow(Z), 5L))
})
