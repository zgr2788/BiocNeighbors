# library(testthat); library(BiocNeighbors); source("setup.R"); source("test-FlannKdtreeParam.R")

test_that("FlannKdtreeParam construction behaves properly", {
    skip_if_not_installed("rflann")

    p <- FlannKdtreeParam()
    expect_output(show(p), "FlannKdtreeParam")
    expect_identical(bndistance(p), "Euclidean")

    p <- FlannKdtreeParam(distance="Cosine")
    expect_identical(bndistance(p), "Cosine")
})

test_that("Cursory checks for FlannKdtreeParam", {
    skip_if_not_installed("rflann")
    set.seed(1301)

    Y <- matrix(rnorm(10000), ncol=20)

    p <- FlannKdtreeParam()
    out <- findKNN(Y, k=8, BNPARAM=p)
    expect_identical(ncol(out$distance), 8L)
    expect_identical(ncol(out$index), 8L)

    dist.only <- findDistance(Y, k=8, BNPARAM=p)
    expect_identical(length(dist.only), nrow(Y))

    transposed <- findKNN(Y, k=8, BNPARAM=p, get.index="transposed", get.distance=FALSE)
    expect_identical(dim(transposed$index), c(8L, nrow(Y)))
})

test_that("FlannKdtreeParam queries behave with cosine distance", {
    skip_if_not_installed("rflann")
    set.seed(1302)

    Y <- matrix(rnorm(10000), ncol=20)
    Z <- matrix(rnorm(2000), ncol=20)

    p <- FlannKdtreeParam(distance="Cosine")
    out <- queryKNN(Y, Z, k=8, BNPARAM=p)

    Y1 <- Y/sqrt(rowSums(Y^2))
    Z1 <- Z/sqrt(rowSums(Z^2))
    ref <- queryKNN(
        Y1,
        Z1,
        k=8,
        BNPARAM=FlannKdtreeParam(checks=p[["checks"]], distance="Euclidean")
    )

    expect_equal(out, ref)
})

test_that("FlannKdtreeParam can save and reload indices", {
    skip_if_not_installed("rflann")
    set.seed(1303)

    Y <- matrix(rnorm(10000), ncol=20)
    idx <- buildIndex(Y, BNPARAM=FlannKdtreeParam())

    tmp <- tempfile()
    dir.create(tmp)
    saveIndex(idx, tmp)

    reloaded <- loadIndex(tmp)
    out <- findKNN(reloaded, k=5)
    expect_s4_class(reloaded, "FlannKdtreeIndex")
    expect_identical(dim(out$index), c(nrow(Y), 5L))
    expect_identical(dim(out$distance), c(nrow(Y), 5L))
})

test_that("FlannKdtreeParam works in findMutualNN", {
    skip_if_not_installed("rflann")
    set.seed(1304)

    A <- matrix(rnorm(10000), ncol=20)
    B <- matrix(rnorm(20000), ncol=20)
    p <- FlannKdtreeParam()

    ref <- findMutualNN(A, B, k1=10, BNPARAM=p)
    expect_identical(length(ref$first), length(ref$second))
    expect_true(length(ref$first) > 0)

    B1 <- buildIndex(A, BNPARAM=p)
    B2 <- buildIndex(B, BNPARAM=p)
    out <- findMutualNN(A, B, k1=10, BNPARAM=p, BNINDEX1=B1, BNINDEX2=B2)
    expect_identical(length(out$first), length(out$second))
    expect_true(length(out$first) > 0)
})

test_that("FlannKdtreeParam supports the shared threading interface", {
    skip_if_not_installed("rflann")
    set.seed(1305)

    Y <- matrix(rnorm(10000), ncol=20)
    Z <- matrix(rnorm(2000), ncol=20)
    p <- FlannKdtreeParam()

    built.bp <- buildIndex(Y, BNPARAM=p, BPPARAM=BiocParallel::SnowParam(2))
    built.nt <- buildIndex(Y, BNPARAM=p, num.threads=2)
    out.bp <- queryKNN(built.bp, Z, k=5, BPPARAM=BiocParallel::SnowParam(2))
    out.nt <- queryKNN(built.bp, Z, k=5, num.threads=2)

    expect_s4_class(built.bp, "FlannKdtreeIndex")
    expect_s4_class(built.nt, "FlannKdtreeIndex")
    expect_equal(out.bp, out.nt)
})
