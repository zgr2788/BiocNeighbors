#' The FlannKmeansParam class
#'
#' A class to hold parameters for approximate nearest-neighbor searches with
#' FLANN hierarchical k-means trees.
#'
#' @param checks Integer scalar specifying the number of FLANN checks during
#' search. Larger values improve accuracy at the cost of slower searches.
#' The default is intentionally speed-oriented for approximate searches.
#' @inheritParams ExhaustiveParam
#'
#' @details
#' This backend uses the \pkg{rflann} package to query FLANN's hierarchical
#' k-means tree implementation. The current \pkg{rflann} wrapper does not
#' expose a reusable built-index object, so \code{\link{buildIndex}} stores the
#' processed reference data and reissues the FLANN search call for each query.
#'
#' Only Euclidean searches are directly supported by \pkg{rflann}. As with the
#' existing cosine implementations in \pkg{BiocNeighbors}, cosine distance is
#' handled by L2-normalizing each observation and running Euclidean search on
#' the normalized data.
#'
#' \pkg{rflann} reports squared Euclidean distances; these are square-rooted in
#' the wrapper so that the returned distances match the conventions used
#' elsewhere in \pkg{BiocNeighbors}.
#'
#' @return
#' The \code{FlannKmeansParam} constructor returns an instance of the
#' FlannKmeansParam class.
#'
#' The \code{\link{buildIndex}} method returns an instance of the
#' \code{FlannKmeansIndex} class.
#'
#' @author
#' Oz Beker
#'
#' @seealso
#' \linkS4class{BiocNeighborParam}, for the parent class and its available
#' methods.
#'
#' \url{https://github.com/YeeJeremy/rflann}, for details on the underlying
#' wrapper.
#'
#' @examples
#' (out <- FlannKmeansParam())
#' out[['checks']]
#'
#' out[['checks']] <- 128L
#' out
#'
#' @aliases
#' FlannKmeansParam-class
#' show,FlannKmeansParam-method
#' FlannKmeansIndex
#' FlannKmeansIndex-class
#'
#' @docType class
#'
#' @export
#' @importFrom methods new
FlannKmeansParam <- function(checks=8L, distance=c("Euclidean", "Cosine")) {
    new("FlannKmeansParam", checks=as.integer(checks), distance=match.arg(distance))
}

setValidity("FlannKmeansParam", function(object) {
    msg <- character(0)

    checks <- object[["checks"]]
    if (length(checks) != 1L || is.na(checks) || checks <= 0L) {
        msg <- c(msg, "'checks' should be a positive integer scalar")
    }

    if (length(msg)) return(msg)
    TRUE
})

#' @export
setMethod("show", "FlannKmeansParam", function(object) {
    callNextMethod()
    cat(sprintf("checks: %i\n", object[["checks"]]))
})

#' @export
FlannKmeansIndex <- function(data, names, param) {
    new("FlannKmeansIndex", data=data, names=names, param=param)
}

.flann_kmeans_algorithm_name <- "BiocNeighbors::FlannKmeans"

.run_flann_kmeans <- function(query, reference, k, checks, num.threads=1L, get.distance=TRUE) {
    nobs <- nrow(query)
    if (!nobs || !k) {
        return(list(
            index=matrix(integer(0), nobs, k),
            distance=if (get.distance) matrix(numeric(0), nobs, k) else NULL
        ))
    }

    if (nrow(reference) <= 1L) {
        out <- queryKNN(
            reference,
            query=query,
            k=k,
            BNPARAM=KmknnParam(distance="Euclidean"),
            get.distance=get.distance,
            num.threads=num.threads
        )
        return(list(index=out$index, distance=out$distance))
    }

    found <- rflann::Neighbour(
        query=query,
        ref=reference,
        k=k,
        build="kmeans",
        cores=max(1L, as.integer(num.threads)),
        checks=as.integer(checks)
    )

    idx <- .coerce_flann_matrix(found$indices, nobs, k)
    if (is.null(idx)) {
        stop("failed to coerce FLANN search indices into the expected matrix shape")
    }

    dist <- NULL
    if (get.distance) {
        dist <- .coerce_flann_matrix(found$distances, nobs, k)
        if (is.null(dist)) {
            stop("failed to coerce FLANN search distances into the expected matrix shape")
        }
        dist <- sqrt(dist)
    }

    list(index=idx, distance=dist)
}

#' @export
#' @rdname buildIndex
setMethod("buildIndex", "FlannKmeansParam", function(X, BNPARAM, transposed=FALSE, num.threads=1, BPPARAM=NULL, ..., .check.nonfinite=TRUE) {
    .require_rflann()
    .resolve_num_threads(num.threads, BPPARAM)

    X <- .coerce_flann_observation_rows(X, transposed=transposed)
    if (.check.nonfinite && any(!is.finite(as.matrix(X)))) {
        stop("cannot build an index from non-finite values")
    }

    processed <- .prepare_flann_reference(X, BNPARAM)
    FlannKmeansIndex(data=processed, names=rownames(processed), param=BNPARAM)
})

#' @export
#' @rdname findKNN
setMethod("findKnnFromIndex", "FlannKmeansIndex", function(BNINDEX, k, get.index=TRUE, get.distance=TRUE, num.threads=1, subset=NULL, ...) {
    .require_rflann()

    chosen <- .subset_flann_index(BNINDEX, subset)
    validated <- .validate_and_cap_k(k, length(chosen), max(nrow(BNINDEX@data) - 1L, 0L))
    k <- validated$k
    variable <- validated$variable
    max.k <- if (length(k)) max(k) else 0L

    report.distance <- !isFALSE(get.distance)

    if (!length(chosen)) {
        found <- list(
            index=matrix(integer(0), 0, max.k),
            distance=if (report.distance) matrix(numeric(0), 0, max.k) else NULL
        )
    } else if (!max.k) {
        found <- list(
            index=matrix(integer(0), length(chosen), 0),
            distance=if (report.distance) matrix(numeric(0), length(chosen), 0) else NULL
        )
    } else {
        found <- .run_flann_kmeans(
            query=BNINDEX@data[chosen,,drop=FALSE],
            reference=BNINDEX@data,
            k=min(max.k + 1L, nrow(BNINDEX@data)),
            checks=BNINDEX@param@checks,
            num.threads=num.threads,
            get.distance=report.distance
        )
        found <- .drop_flann_self_matches(found$index, found$distance, chosen, max.k)
    }

    if (variable) {
        sliced <- .extract_flann_variable_neighbors(found$index, found$distance, k)
        .format_flann_output(sliced$index, sliced$distance, get.index, get.distance, variable=TRUE)
    } else {
        found <- .orient_flann_output(found)
        .format_flann_output(found$index, found$distance, get.index, get.distance, variable=FALSE)
    }
})

#' @export
#' @rdname queryKNN
setMethod("queryKnnFromIndex", "FlannKmeansIndex", function(
    BNINDEX,
    query,
    k,
    get.index=TRUE,
    get.distance=TRUE,
    num.threads=1,
    subset=NULL,
    transposed=FALSE,
    ...,
    .check.nonfinite=TRUE
) {
    .require_rflann()

    query <- .coerce_flann_observation_rows(query, transposed=transposed, subset=subset)
    if (.check.nonfinite && any(!is.finite(as.matrix(query)))) {
        stop("cannot query an index with non-finite values")
    }

    query <- .prepare_flann_query(query, BNINDEX@param)
    validated <- .validate_and_cap_k(k, nrow(query), nrow(BNINDEX@data))
    k <- validated$k
    variable <- validated$variable
    max.k <- if (length(k)) max(k) else 0L

    report.distance <- !isFALSE(get.distance)

    found <- .run_flann_kmeans(
        query=query,
        reference=BNINDEX@data,
        k=max.k,
        checks=BNINDEX@param@checks,
        num.threads=num.threads,
        get.distance=report.distance
    )

    if (variable) {
        sliced <- .extract_flann_variable_neighbors(found$index, found$distance, k)
        .format_flann_output(sliced$index, sliced$distance, get.index, get.distance, variable=TRUE)
    } else {
        found <- .orient_flann_output(found)
        .format_flann_output(found$index, found$distance, get.index, get.distance, variable=FALSE)
    }
})

#' @export
#' @rdname findDistance
setMethod("findDistanceFromIndex", "FlannKmeansIndex", function(BNINDEX, k, num.threads=1, subset=NULL, ...) {
    found <- findKnnFromIndex(BNINDEX, k=k, get.index=FALSE, get.distance=TRUE, num.threads=num.threads, subset=subset)
    .last_distance_from_knn(found, k)
})

#' @export
#' @rdname queryDistance
setMethod("queryDistanceFromIndex", "FlannKmeansIndex", function(
    BNINDEX,
    query,
    k,
    num.threads=1,
    subset=NULL,
    transposed=FALSE,
    ...,
    .check.nonfinite=TRUE
) {
    found <- queryKnnFromIndex(
        BNINDEX,
        query=query,
        k=k,
        get.index=FALSE,
        get.distance=TRUE,
        num.threads=num.threads,
        subset=subset,
        transposed=transposed,
        ...,
        .check.nonfinite=.check.nonfinite
    )

    .last_distance_from_knn(found, k)
})

#' @export
#' @rdname findNeighbors
setMethod("findNeighborsFromIndex", "FlannKmeansIndex", function(BNINDEX, threshold, get.index=TRUE, get.distance=TRUE, num.threads=1, subset=NULL, ...) {
    stop("range-based searches are not supported for FlannKmeansParam()")
})

#' @export
#' @rdname queryNeighbors
setMethod("queryNeighborsFromIndex", "FlannKmeansIndex", function(BNINDEX, query, threshold, get.index=TRUE, get.distance=TRUE, num.threads=1, subset=NULL, transposed=FALSE, ..., .check.nonfinite=TRUE) {
    stop("range-based searches are not supported for FlannKmeansParam()")
})

#' @export
#' @rdname saveIndex
setMethod("saveIndex", "FlannKmeansIndex", function(BNINDEX, dir, ...) {
    writeLines(.flann_kmeans_algorithm_name, file.path(dir, "ALGORITHM"))
    saveRDS(
        list(
            data=BNINDEX@data,
            names=BNINDEX@names,
            param=BNINDEX@param
        ),
        file.path(dir, "index.rds")
    )
    invisible(NULL)
})

.load_flann_kmeans_index <- function(dir, ...) {
    .require_rflann()

    payload <- readRDS(file.path(dir, "index.rds"))
    data <- payload$data
    if (!is.null(payload$names)) {
        rownames(data) <- payload$names
    }

    out <- buildIndex(data, BNPARAM=payload$param, .check.nonfinite=FALSE)
    out@names <- payload$names
    out
}
