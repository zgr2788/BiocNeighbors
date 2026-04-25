# library(testthat); library(BiocNeighbors); source("setup.R"); source("test-approximate-backends.R")

.expect_scalar_knn <- function(out, nobs, k, nref, transposed=FALSE) {
    expected <- if (transposed) c(k, nobs) else c(nobs, k)

    if (!is.null(out$index)) {
        expect_identical(dim(out$index), expected)
        if (length(out$index)) {
            expect_true(all(out$index >= 1L & out$index <= nref))
        }
    }

    if (!is.null(out$distance)) {
        expect_identical(dim(out$distance), expected)
        if (length(out$distance)) {
            expect_true(all(is.finite(out$distance)))
        }
    }
}

.expect_variable_knn <- function(out, k, nref) {
    if (!is.null(out$index)) {
        expect_type(out$index, "list")
        expect_identical(length(out$index), length(k))
        expect_identical(lengths(out$index), as.integer(k))
        if (length(unlist(out$index, use.names=FALSE))) {
            expect_true(all(unlist(out$index, use.names=FALSE) >= 1L & unlist(out$index, use.names=FALSE) <= nref))
        }
    }

    if (!is.null(out$distance)) {
        expect_type(out$distance, "list")
        expect_identical(length(out$distance), length(k))
        expect_identical(lengths(out$distance), as.integer(k))
        if (length(unlist(out$distance, use.names=FALSE))) {
            expect_true(all(is.finite(unlist(out$distance, use.names=FALSE))))
        }
    }
}

.check_approximate_backend <- function(BNPARAM, package) {
    skip_if_not_installed(package)
    set.seed(9001)

    Y <- matrix(rnorm(120), ncol=4)
    rownames(Y) <- paste0("cell", seq_len(nrow(Y)))
    Z <- matrix(rnorm(40), ncol=4)
    rownames(Z) <- paste0("query", seq_len(nrow(Z)))

    built <- buildIndex(Y, BNPARAM=BNPARAM)
    scalar <- findKNN(built, k=3)
    .expect_scalar_knn(scalar, nrow(Y), 3L, nrow(Y))

    queried <- queryKNN(built, Z, k=3)
    .expect_scalar_knn(queried, nrow(Z), 3L, nrow(Y))

    fdist <- findDistance(built, k=3)
    expect_identical(length(fdist), nrow(Y))
    expect_true(all(is.finite(fdist)))

    qdist <- queryDistance(built, Z, k=3)
    expect_identical(length(qdist), nrow(Z))
    expect_true(all(is.finite(qdist)))

    by.integer <- findKNN(built, k=2, subset=c(1L, 5L, 7L))
    .expect_scalar_knn(by.integer, 3L, 2L, nrow(Y))

    by.logical <- findKNN(built, k=2, subset=seq_len(nrow(Y)) %% 3L == 0L)
    .expect_scalar_knn(by.logical, 10L, 2L, nrow(Y))

    by.character <- findKNN(built, k=2, subset=rownames(Y)[c(2L, 4L, 6L)])
    .expect_scalar_knn(by.character, 3L, 2L, nrow(Y))

    q.integer <- queryKNN(built, Z, k=2, subset=c(1L, 3L, 5L))
    .expect_scalar_knn(q.integer, 3L, 2L, nrow(Y))

    q.logical <- queryKNN(built, Z, k=2, subset=seq_len(nrow(Z)) %% 2L == 0L)
    .expect_scalar_knn(q.logical, 5L, 2L, nrow(Y))

    q.character <- queryKNN(built, Z, k=2, subset=rownames(Z)[c(2L, 4L, 6L)])
    .expect_scalar_knn(q.character, 3L, 2L, nrow(Y))

    forced <- findKNN(built, k=I(2L), subset=1L)
    .expect_variable_knn(forced, I(2L), nrow(Y))

    variable <- queryKNN(built, Z, k=I(c(1L, 3L, 2L)), subset=1:3)
    .expect_variable_knn(variable, I(c(1L, 3L, 2L)), nrow(Y))

    variable.unclassed <- findKNN(built, k=rep(c(1L, 2L), length.out=nrow(Y)))
    .expect_variable_knn(variable.unclassed, rep(c(1L, 2L), length.out=nrow(Y)), nrow(Y))

    no.index <- queryKNN(built, Z, k=3, get.index=FALSE)
    expect_null(no.index$index)
    .expect_scalar_knn(no.index, nrow(Z), 3L, nrow(Y))

    no.distance <- queryKNN(built, Z, k=3, get.distance=FALSE)
    expect_null(no.distance$distance)
    .expect_scalar_knn(no.distance, nrow(Z), 3L, nrow(Y))

    nothing <- queryKNN(built, Z, k=3, get.index=FALSE, get.distance=FALSE)
    expect_identical(nothing, list())

    normal <- queryKNN(built, Z, k=3, get.index="normal", get.distance="normal")
    .expect_scalar_knn(normal, nrow(Z), 3L, nrow(Y))

    transposed <- queryKNN(built, Z, k=3, get.index="transposed", get.distance="transposed")
    .expect_scalar_knn(transposed, nrow(Z), 3L, nrow(Y), transposed=TRUE)

    one <- suppressWarnings(buildIndex(Y[1,,drop=FALSE], BNPARAM=BNPARAM))
    expect_warning(one.find <- findKNN(one, k=3), "more neighbors")
    .expect_scalar_knn(one.find, 1L, 0L, 1L)
    expect_warning(one.query <- queryKNN(one, Z[1:2,,drop=FALSE], k=3), "more neighbors")
    .expect_scalar_knn(one.query, 2L, 1L, 1L)

    empty <- suppressWarnings(buildIndex(Y[0,,drop=FALSE], BNPARAM=BNPARAM))
    expect_warning(empty.find <- findKNN(empty, k=3), "more neighbors")
    .expect_scalar_knn(empty.find, 0L, 0L, 0L)
    expect_warning(empty.query <- queryKNN(empty, Z[1:2,,drop=FALSE], k=3), "more neighbors")
    .expect_scalar_knn(empty.query, 2L, 0L, 0L)

    expect_error(findKNN(built, k=NA_integer_), "NA")
    expect_error(queryKNN(built, Z, k=-1L), "non-negative")
    expect_error(findKNN(built, k=I(c(1L, 2L)), subset=1:3), "length of 'k'")
    expect_error(queryKNN(built, Z, k=I(c(1L, 2L)), subset=1:3), "length of 'k'")
    expect_error(findKNN(built, k=2, subset=c(1L, NA_integer_)), "finite")
    expect_error(queryKNN(built, Z, k=2, subset=c(TRUE, FALSE)), "logical 'subset'")
    expect_error(findKNN(built, k=2, subset="missing"), "failed to match")
}

test_that("NndescentParam supports the approximate backend contract", {
    .check_approximate_backend(NndescentParam(), "rnndescent")
})

test_that("RpforestParam supports the approximate backend contract", {
    .check_approximate_backend(RpforestParam(), "rnndescent")
})

test_that("FlannKdtreeParam supports the approximate backend contract", {
    .check_approximate_backend(FlannKdtreeParam(checks=128L), "rflann")
})

test_that("FlannKmeansParam supports the approximate backend contract", {
    .check_approximate_backend(FlannKmeansParam(checks=128L), "rflann")
})

test_that("MlpackLshParam supports the approximate backend contract", {
    .check_approximate_backend(MlpackLshParam(seed=9001L), "mlpack")
})

test_that("RP-forest fallback respects the requested distance", {
    reference <- matrix(c(0, 0, 2, 0, 0, 3, 10, 10), ncol=2, byrow=TRUE)
    query <- matrix(c(1, 2), ncol=2)
    ref <- queryKNN(reference, query, k=2, BNPARAM=KmknnParam(distance="Manhattan"))

    broken <- list(
        idx=matrix(NA_integer_, nrow=1L, ncol=2L),
        dist=matrix(NA_real_, nrow=1L, ncol=2L)
    )
    out <- BiocNeighbors:::.standardize_rpforest_hits(
        broken,
        reference=reference,
        query=query,
        k=2L,
        distance="Manhattan"
    )

    expect_identical(out$idx, ref$index)
    expect_identical(out$dist, ref$distance)
})

test_that("FLANN distances are returned on the Euclidean scale", {
    skip_if_not_installed("rflann")

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
    ref <- refQueryKNN(reference, query, k=3)

    out <- queryKNN(reference, query, k=3, BNPARAM=FlannKdtreeParam(checks=128L))
    expect_identical(out$index, ref$index)
    expect_equal(out$distance, ref$distance)

    out <- queryKNN(reference, query, k=3, BNPARAM=FlannKmeansParam(checks=128L))
    expect_identical(out$index, ref$index)
    expect_equal(out$distance, ref$distance)
})
