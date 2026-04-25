# library(testthat); library(BiocNeighbors); source("setup.R"); source("test-MlpackLshParam.R")

.expect_valid_mlpack_lsh_knn <- function(out, nobs, k, nref) {
    expect_identical(dim(out$index), c(nobs, k))
    expect_identical(dim(out$distance), c(nobs, k))
    if (length(out$index)) {
        expect_true(all(out$index >= 1L & out$index <= nref))
        expect_true(all(is.finite(out$distance)))
    }
}

test_that("MlpackLshParam construction behaves properly", {
    skip_if_not_installed("mlpack")

    p <- MlpackLshParam()
    expect_output(show(p), "MlpackLshParam")
    expect_identical(bndistance(p), "Euclidean")
    expect_identical(p[["num.tables"]], 30L)
    expect_identical(p[["projections"]], 10L)
    expect_identical(p[["bucket.size"]], 500L)
    expect_identical(p[["second.hash.size"]], 99901L)
    expect_identical(p[["num.probes"]], 0L)
    expect_true(is.na(p[["seed"]]))

    p <- MlpackLshParam(distance="Cosine", search=2L, seed=123L)
    expect_identical(bndistance(p), "Cosine")
    expect_identical(p[["num.probes"]], 2L)
    expect_identical(p[["seed"]], 123L)

    expect_error(MlpackLshParam(num.tables=0L), "positive")
    expect_error(MlpackLshParam(projections=NA_integer_), "positive")
    expect_error(MlpackLshParam(bucket.size=-1L), "positive")
    expect_error(MlpackLshParam(second.hash.size=0L), "positive")
    expect_error(MlpackLshParam(num.probes=-1L), "non-negative")
    expect_error(MlpackLshParam(seed=-1L), "non-negative")
    expect_error(MlpackLshParam(num.probes=1L, search=2L), "different")
})

test_that("Cursory checks for MlpackLshParam", {
    skip_if_not_installed("mlpack")
    set.seed(1501)

    Y <- matrix(rnorm(10000), ncol=20)

    p <- MlpackLshParam(seed=1501L)
    out <- findKNN(Y, k=8, BNPARAM=p)
    expect_identical(ncol(out$distance), 8L)
    expect_identical(ncol(out$index), 8L)

    dist.only <- findDistance(Y, k=8, BNPARAM=p)
    expect_identical(length(dist.only), nrow(Y))

    transposed <- findKNN(Y, k=8, BNPARAM=p, get.index="transposed", get.distance=FALSE)
    expect_identical(dim(transposed$index), c(8L, nrow(Y)))
})

test_that("MlpackLshParam queries behave with cosine distance", {
    skip_if_not_installed("mlpack")
    set.seed(1502)

    Y <- matrix(rnorm(10000), ncol=20)
    Z <- matrix(rnorm(2000), ncol=20)

    p <- MlpackLshParam(distance="Cosine", seed=1502L)
    out <- queryKNN(Y, Z, k=8, BNPARAM=p)

    Y1 <- Y/sqrt(rowSums(Y^2))
    Z1 <- Z/sqrt(rowSums(Z^2))

    expected <- matrix(NA_real_, nrow=nrow(Z), ncol=8L)
    for (i in seq_len(nrow(Z))) {
        expected[i,] <- sqrt(rowSums((Y1[out$index[i,],,drop=FALSE] -
            matrix(Z1[i,], nrow=8L, ncol=ncol(Z1), byrow=TRUE))^2))
    }

    .expect_valid_mlpack_lsh_knn(out, nrow(Z), 8L, nrow(Y))
    expect_equal(out$distance, expected)
})

test_that("MlpackLshParam returns one-based indices and Euclidean distances", {
    skip_if_not_installed("mlpack")

    reference <- matrix(c(
        0, 0,
        3, 4,
        8, 9,
        10, 1,
        1, 11
    ), ncol=2, byrow=TRUE)
    query <- matrix(c(
        0, 0,
        8, 9
    ), ncol=2, byrow=TRUE)

    out <- queryKNN(reference, query, k=3, BNPARAM=MlpackLshParam(seed=1503L))
    .expect_valid_mlpack_lsh_knn(out, nrow(query), 3L, nrow(reference))

    expected <- matrix(NA_real_, nrow=nrow(query), ncol=3L)
    for (i in seq_len(nrow(query))) {
        expected[i,] <- sqrt(rowSums((reference[out$index[i,],,drop=FALSE] -
            matrix(query[i,], nrow=3L, ncol=ncol(query), byrow=TRUE))^2))
    }

    expect_equal(out$distance, expected)
})

test_that("MlpackLshParam can save and reload indices", {
    skip_if_not_installed("mlpack")
    set.seed(1504)

    Y <- matrix(rnorm(10000), ncol=20)
    idx <- buildIndex(Y, BNPARAM=MlpackLshParam(seed=1504L))

    tmp <- tempfile()
    dir.create(tmp)
    saveIndex(idx, tmp)

    reloaded <- loadIndex(tmp)
    out <- findKNN(reloaded, k=5)
    expect_s4_class(reloaded, "MlpackLshIndex")
    expect_identical(dim(out$index), c(nrow(Y), 5L))
    expect_identical(dim(out$distance), c(nrow(Y), 5L))
})

test_that("MlpackLshParam works in findMutualNN", {
    skip_if_not_installed("mlpack")
    set.seed(1505)

    A <- matrix(rnorm(10000), ncol=20)
    B <- A + matrix(rnorm(length(A), sd=0.01), nrow=nrow(A))
    p <- MlpackLshParam(seed=1505L)

    ref <- findMutualNN(A, B, k1=10, BNPARAM=p)
    expect_identical(length(ref$first), length(ref$second))
    expect_true(length(ref$first) > 0)

    B1 <- buildIndex(A, BNPARAM=p)
    B2 <- buildIndex(B, BNPARAM=p)
    out <- findMutualNN(A, B, k1=10, BNPARAM=p, BNINDEX1=B1, BNINDEX2=B2)
    expect_identical(length(out$first), length(out$second))
    expect_true(length(out$first) > 0)
})

test_that("MlpackLshParam supports the shared threading interface", {
    skip_if_not_installed("mlpack")
    set.seed(1506)

    Y <- matrix(rnorm(10000), ncol=20)
    Z <- matrix(rnorm(2000), ncol=20)
    p <- MlpackLshParam(seed=1506L)

    built.bp <- buildIndex(Y, BNPARAM=p, BPPARAM=BiocParallel::SnowParam(2))
    built.nt <- buildIndex(Y, BNPARAM=p, num.threads=2)
    out.bp <- queryKNN(built.bp, Z, k=5, BPPARAM=BiocParallel::SnowParam(2))
    out.nt <- queryKNN(built.bp, Z, k=5, num.threads=2)

    expect_s4_class(built.bp, "MlpackLshIndex")
    expect_s4_class(built.nt, "MlpackLshIndex")
    .expect_valid_mlpack_lsh_knn(out.bp, nrow(Z), 5L, nrow(Y))
    .expect_valid_mlpack_lsh_knn(out.nt, nrow(Z), 5L, nrow(Y))
})

test_that("MlpackLshParam does not support range searches", {
    skip_if_not_installed("mlpack")
    set.seed(1507)

    Y <- matrix(rnorm(1000), ncol=10)
    Z <- matrix(rnorm(200), ncol=10)
    p <- MlpackLshParam(seed=1507L)

    expect_error(findNeighbors(Y, threshold=1, BNPARAM=p), "range-based")
    expect_error(queryNeighbors(Y, Z, threshold=1, BNPARAM=p), "range-based")
})
