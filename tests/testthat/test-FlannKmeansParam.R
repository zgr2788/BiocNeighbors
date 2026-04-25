# library(testthat); library(BiocNeighbors); source("setup.R"); source("test-FlannKmeansParam.R")

test_that("FlannKmeansParam construction behaves properly", {
    skip_if_not_installed("rflann")

    p <- FlannKmeansParam()
    expect_output(show(p), "FlannKmeansParam")
    expect_identical(bndistance(p), "Euclidean")

    p <- FlannKmeansParam(distance="Cosine")
    expect_identical(bndistance(p), "Cosine")
})

test_that("Cursory checks for FlannKmeansParam", {
    skip_if_not_installed("rflann")
    set.seed(1401)

    Y <- matrix(rnorm(10000), ncol=20)

    p <- FlannKmeansParam()
    out <- findKNN(Y, k=8, BNPARAM=p)
    expect_identical(ncol(out$distance), 8L)
    expect_identical(ncol(out$index), 8L)

    dist.only <- findDistance(Y, k=8, BNPARAM=p)
    expect_identical(length(dist.only), nrow(Y))

    transposed <- findKNN(Y, k=8, BNPARAM=p, get.index="transposed", get.distance=FALSE)
    expect_identical(dim(transposed$index), c(8L, nrow(Y)))
})

test_that("FlannKmeansParam queries behave with cosine distance", {
    skip_if_not_installed("rflann")
    set.seed(1402)

    Y <- matrix(rnorm(10000), ncol=20)
    Z <- matrix(rnorm(2000), ncol=20)

    p <- FlannKmeansParam(checks=64L, distance="Cosine")
    out <- queryKNN(Y, Z, k=8, BNPARAM=p)

    Y1 <- Y/sqrt(rowSums(Y^2))
    Z1 <- Z/sqrt(rowSums(Z^2))
    ref <- queryKNN(
        Y1,
        Z1,
        k=8,
        BNPARAM=FlannKmeansParam(checks=p[["checks"]], distance="Euclidean")
    )

    expect_equal(out, ref)
})

test_that("FlannKmeansParam can save and reload indices", {
    skip_if_not_installed("rflann")
    set.seed(1403)

    Y <- matrix(rnorm(10000), ncol=20)
    idx <- buildIndex(Y, BNPARAM=FlannKmeansParam())

    tmp <- tempfile()
    dir.create(tmp)
    saveIndex(idx, tmp)

    reloaded <- loadIndex(tmp)
    out <- findKNN(reloaded, k=5)
    expect_s4_class(reloaded, "FlannKmeansIndex")
    expect_identical(dim(out$index), c(nrow(Y), 5L))
    expect_identical(dim(out$distance), c(nrow(Y), 5L))
})

test_that("FlannKmeansParam works in findMutualNN", {
    skip_if_not_installed("rflann")
    set.seed(1404)

    A <- matrix(rnorm(10000), ncol=20)
    B <- matrix(rnorm(20000), ncol=20)
    p <- FlannKmeansParam()

    ref <- findMutualNN(A, B, k1=10, BNPARAM=p)
    expect_identical(length(ref$first), length(ref$second))
    expect_true(length(ref$first) > 0)

    B1 <- buildIndex(A, BNPARAM=p)
    B2 <- buildIndex(B, BNPARAM=p)
    out <- findMutualNN(A, B, k1=10, BNPARAM=p, BNINDEX1=B1, BNINDEX2=B2)
    expect_identical(length(out$first), length(out$second))
    expect_true(length(out$first) > 0)
})

test_that("FlannKmeansParam supports the shared threading interface", {
    skip_if_not_installed("rflann")
    set.seed(1405)

    Y <- matrix(rnorm(10000), ncol=20)
    Z <- matrix(rnorm(2000), ncol=20)
    p <- FlannKmeansParam()

    built.bp <- buildIndex(Y, BNPARAM=p, BPPARAM=BiocParallel::SnowParam(2))
    built.nt <- buildIndex(Y, BNPARAM=p, num.threads=2)
    out.bp <- queryKNN(built.bp, Z, k=5, BPPARAM=BiocParallel::SnowParam(2))
    out.nt <- queryKNN(built.bp, Z, k=5, num.threads=2)

    expect_s4_class(built.bp, "FlannKmeansIndex")
    expect_s4_class(built.nt, "FlannKmeansIndex")
    expect_identical(dim(out.bp$index), c(nrow(Z), 5L))
    expect_identical(dim(out.bp$distance), c(nrow(Z), 5L))
    expect_identical(dim(out.nt$index), c(nrow(Z), 5L))
    expect_identical(dim(out.nt$distance), c(nrow(Z), 5L))
    expect_true(all(out.bp$index >= 1L & out.bp$index <= nrow(Y)))
    expect_true(all(out.nt$index >= 1L & out.nt$index <= nrow(Y)))
    expect_true(all(is.finite(out.bp$distance)))
    expect_true(all(is.finite(out.nt$distance)))
})
