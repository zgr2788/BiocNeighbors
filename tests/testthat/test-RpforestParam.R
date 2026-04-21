# library(testthat); library(BiocNeighbors); source("setup.R"); source("test-RpforestParam.R")

test_that("RpforestParam construction behaves properly", {
    skip_if_not_installed("rnndescent")

    p <- RpforestParam()
    expect_output(show(p), "RpforestParam")
    expect_identical(bndistance(p), "Euclidean")

    p <- RpforestParam(distance="Manhattan")
    expect_identical(bndistance(p), "Manhattan")
})

test_that("Cursory checks for RpforestParam", {
    skip_if_not_installed("rnndescent")
    set.seed(1101)

    Y <- matrix(rnorm(10000), ncol=20)

    p <- RpforestParam()
    out <- findKNN(Y, k=8, BNPARAM=p)
    expect_identical(ncol(out$distance), 8L)
    expect_identical(ncol(out$index), 8L)

    dist.only <- findDistance(Y, k=8, BNPARAM=p)
    expect_identical(length(dist.only), nrow(Y))

    transposed <- findKNN(Y, k=8, BNPARAM=p, get.index="transposed", get.distance=FALSE)
    expect_identical(dim(transposed$index), c(8L, nrow(Y)))
})

test_that("RpforestParam queries behave with cosine distance", {
    skip_if_not_installed("rnndescent")
    set.seed(1102)

    Y <- matrix(rnorm(10000), ncol=20)
    Z <- matrix(rnorm(2000), ncol=20)

    set.seed(2102)
    p <- RpforestParam(distance="Cosine")
    out <- queryKNN(Y, Z, k=8, BNPARAM=p)

    Y1 <- Y/sqrt(rowSums(Y^2))
    Z1 <- Z/sqrt(rowSums(Z^2))
    set.seed(2102)
    ref <- queryKNN(
        Y1,
        Z1,
        k=8,
        BNPARAM=RpforestParam(
            nneighbors=p[["nneighbors"]],
            ntrees=p[["ntrees"]],
            leaf.size=p[["leaf.size"]],
            max.tree.depth=p[["max.tree.depth"]],
            distance="Euclidean"
        )
    )

    expect_equal(out, ref)
})

test_that("RpforestParam can save and reload indices", {
    skip_if_not_installed("rnndescent")
    set.seed(1103)

    Y <- matrix(rnorm(10000), ncol=20)
    idx <- buildIndex(Y, BNPARAM=RpforestParam())

    tmp <- tempfile()
    dir.create(tmp)
    saveIndex(idx, tmp)

    reloaded <- loadIndex(tmp)
    out <- findKNN(reloaded, k=5)
    expect_s4_class(reloaded, "RpforestIndex")
    expect_identical(dim(out$index), c(nrow(Y), 5L))
    expect_identical(dim(out$distance), c(nrow(Y), 5L))
})

test_that("RpforestParam works in findMutualNN", {
    skip_if_not_installed("rnndescent")
    set.seed(1104)

    A <- matrix(rnorm(10000), ncol=20)
    B <- matrix(rnorm(20000), ncol=20)
    p <- RpforestParam()

    ref <- findMutualNN(A, B, k1=10, BNPARAM=p)
    expect_identical(length(ref$first), length(ref$second))
    expect_true(length(ref$first) > 0)

    B1 <- buildIndex(A, BNPARAM=p)
    B2 <- buildIndex(B, BNPARAM=p)
    out <- findMutualNN(A, B, k1=10, BNPARAM=p, BNINDEX1=B1, BNINDEX2=B2)
    expect_identical(length(out$first), length(out$second))
    expect_true(length(out$first) > 0)
})
