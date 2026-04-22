#' The FlannKdtreeParam class
#'
#' A class to hold parameters for approximate nearest-neighbor searches with
#' FLANN randomized kd-trees.
#'
#' @param checks Integer scalar specifying the number of FLANN checks during
#' search. Larger values improve accuracy at the cost of slower searches.
#' The default is intentionally speed-oriented for approximate searches.
#' @inheritParams ExhaustiveParam
#'
#' @details
#' This backend uses the \pkg{rflann} package to query FLANN's randomized
#' kd-tree implementation. The current \pkg{rflann} wrapper does not expose a
#' reusable built-index object, so \code{\link{buildIndex}} stores the
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
#' The \code{FlannKdtreeParam} constructor returns an instance of the
#' FlannKdtreeParam class.
#'
#' The \code{\link{buildIndex}} method returns an instance of the
#' \code{FlannKdtreeIndex} class.
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
#' (out <- FlannKdtreeParam())
#' out[['checks']]
#'
#' out[['checks']] <- 128L
#' out
#'
#' @aliases
#' FlannKdtreeParam-class
#' show,FlannKdtreeParam-method
#' FlannKdtreeIndex
#' FlannKdtreeIndex-class
#'
#' @docType class
#'
#' @export
#' @importFrom methods new
FlannKdtreeParam <- function(checks=8L, distance=c("Euclidean", "Cosine")) {
    new("FlannKdtreeParam", checks=as.integer(checks), distance=match.arg(distance))
}

setValidity("FlannKdtreeParam", function(object) {
    msg <- character(0)

    checks <- object[["checks"]]
    if (length(checks) != 1L || is.na(checks) || checks <= 0L) {
        msg <- c(msg, "'checks' should be a positive integer scalar")
    }

    if (length(msg)) return(msg)
    TRUE
})

#' @export
setMethod("show", "FlannKdtreeParam", function(object) {
    callNextMethod()
    cat(sprintf("checks: %i\n", object[["checks"]]))
})

#' @export
FlannKdtreeIndex <- function(data, names, param) {
    new("FlannKdtreeIndex", data=data, names=names, param=param)
}

.flann_kdtree_algorithm_name <- "BiocNeighbors::FlannKdtree"

.require_rflann <- function() {
    if (!requireNamespace("rflann", quietly=TRUE)) {
        stop(
            "the 'rflann' package must be installed to use FlannKdtreeParam(); ",
            "see https://github.com/YeeJeremy/rflann",
            call.=FALSE
        )
    }
}

.coerce_flann_observation_rows <- function(x, transposed=FALSE, subset=NULL) {
    if (transposed) {
        if (is.matrix(x)) {
            x <- t(x)
        } else {
            x <- t(as.matrix(x))
        }
    } else if (!is.matrix(x)) {
        x <- as.matrix(x)
    }

    if (!is.null(subset)) {
        x <- x[subset,,drop=FALSE]
    }

    x
}

.l2_normalize_rows_flann <- function(x) {
    x <- as.matrix(x)
    denom <- sqrt(rowSums(x^2))
    denom[denom == 0] <- 1
    x / denom
}

.prepare_flann_reference <- function(x, param) {
    if (bndistance(param) == "Cosine") {
        .l2_normalize_rows_flann(x)
    } else {
        as.matrix(x)
    }
}

.prepare_flann_query <- function(x, param) {
    if (bndistance(param) == "Cosine") {
        .l2_normalize_rows_flann(x)
    } else {
        as.matrix(x)
    }
}

.subset_flann_index <- function(index, subset) {
    if (is.null(subset)) {
        return(seq_len(nrow(index@data)))
    }

    if (is.character(subset)) {
        if (is.null(index@names)) {
            stop("cannot use character 'subset' when observation names are not available")
        }
        subset <- match(subset, index@names)
        if (anyNA(subset)) {
            stop("failed to match some entries in 'subset' to observation names")
        }
        return(as.integer(subset))
    }

    if (is.logical(subset)) {
        return(which(subset))
    }

    as.integer(subset)
}

.cap_flann_k <- function(k, limit) {
    capped <- pmin(as.integer(k), limit)
    if (any(as.integer(k) != capped)) {
        warning("more neighbors were requested than available and will be ignored")
    }
    if (is(k, "AsIs")) {
        I(capped)
    } else {
        capped
    }
}

.extract_flann_variable_neighbors <- function(index, distance, k) {
    output <- list()

    if (!is.null(index)) {
        output$index <- lapply(seq_len(nrow(index)), function(i) index[i,seq_len(k[i]),drop=TRUE])
    }
    if (!is.null(distance)) {
        output$distance <- lapply(seq_len(nrow(distance)), function(i) distance[i,seq_len(k[i]),drop=TRUE])
    }

    output
}

.orient_flann_output <- function(found) {
    output <- list()
    if (!is.null(found$index)) {
        output$index <- t(found$index)
    }
    if (!is.null(found$distance)) {
        output$distance <- t(found$distance)
    }
    output
}

.format_flann_output <- function(index, distance, get.index, get.distance, variable=FALSE) {
    output <- list()

    if (!isFALSE(get.index)) {
        output$index <- index
        if (!variable) {
            output <- .format_output(output, "index", get.index)
        }
    }

    if (!isFALSE(get.distance)) {
        output$distance <- distance
        if (!variable) {
            output <- .format_output(output, "distance", get.distance)
        }
    }

    output
}

.coerce_flann_matrix <- function(x, nobs, k) {
    if (!length(x)) {
        return(matrix(x, nrow=nobs, ncol=k))
    }

    if (is.null(dim(x))) {
        if (nobs == 1L) {
            return(matrix(x, nrow=1L))
        }
        return(NULL)
    }

    if (identical(dim(x), c(nobs, k))) {
        return(x)
    }

    if (identical(dim(x), c(k, nobs))) {
        return(t(x))
    }

    NULL
}

.run_flann_kdtree <- function(query, reference, k, checks, num.threads=1L, get.distance=TRUE) {
    nobs <- nrow(query)
    if (!k) {
        return(list(
            index=matrix(integer(0), nobs, 0),
            distance=if (get.distance) matrix(numeric(0), nobs, 0) else NULL
        ))
    }

    idx <- NULL
    dist <- NULL
    safe.threads <- max(1L, as.integer(num.threads))

    if (!get.distance && safe.threads == 1L) {
        idx <- .coerce_flann_matrix(rflann::FastKDNeighbour(query, reference, k=k), nobs, k)
        if (is.null(idx)) {
            stop("failed to coerce FLANN kd-tree indices into the expected matrix shape")
        }
    } else {
        found <- rflann::Neighbour(
            query=query,
            ref=reference,
            k=k,
            build="kdtree",
            cores=safe.threads,
            checks=as.integer(checks)
        )

        idx <- .coerce_flann_matrix(found$indices, nobs, k)
        if (is.null(idx)) {
            stop("failed to coerce FLANN search indices into the expected matrix shape")
        }

        if (get.distance) {
            dist <- .coerce_flann_matrix(found$distances, nobs, k)
            if (is.null(dist)) {
                stop("failed to coerce FLANN search distances into the expected matrix shape")
            }
            dist <- sqrt(dist)
        }
    }

    list(index=idx, distance=dist)
}

.drop_flann_self_matches <- function(index, distance, self.ids, keep) {
    nr <- nrow(index)
    out.index <- matrix(NA_integer_, nr, keep)
    out.distance <- if (is.null(distance)) NULL else matrix(NA_real_, nr, keep)

    for (i in seq_len(nr)) {
        chosen <- index[i,] != self.ids[i]
        current.index <- index[i,chosen]
        if (!is.null(distance)) {
            current.distance <- distance[i,chosen]
        }
        if (length(current.index) > keep) {
            current.index <- current.index[seq_len(keep)]
            if (!is.null(distance)) {
                current.distance <- current.distance[seq_len(keep)]
            }
        }

        if (keep && length(current.index)) {
            out.index[i,seq_along(current.index)] <- current.index
            if (!is.null(distance)) {
                out.distance[i,seq_along(current.distance)] <- current.distance
            }
        }
    }

    list(index=out.index, distance=out.distance)
}

#' @export
#' @rdname buildIndex
setMethod("buildIndex", "FlannKdtreeParam", function(X, BNPARAM, transposed=FALSE, ..., .check.nonfinite=TRUE) {
    .require_rflann()

    X <- .coerce_flann_observation_rows(X, transposed=transposed)
    if (.check.nonfinite && any(!is.finite(as.matrix(X)))) {
        stop("cannot build an index from non-finite values")
    }

    processed <- .prepare_flann_reference(X, BNPARAM)
    FlannKdtreeIndex(data=processed, names=rownames(processed), param=BNPARAM)
})

#' @export
#' @rdname findKNN
setMethod("findKnnFromIndex", "FlannKdtreeIndex", function(BNINDEX, k, get.index=TRUE, get.distance=TRUE, num.threads=1, subset=NULL, ...) {
    .require_rflann()

    chosen <- .subset_flann_index(BNINDEX, subset)
    k <- .cap_flann_k(k, max(nrow(BNINDEX@data) - 1L, 0L))
    variable <- is(k, "AsIs")
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
        found <- .run_flann_kdtree(
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
setMethod("queryKnnFromIndex", "FlannKdtreeIndex", function(
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
    k <- .cap_flann_k(k, nrow(BNINDEX@data))
    variable <- is(k, "AsIs")
    max.k <- if (length(k)) max(k) else 0L

    report.distance <- !isFALSE(get.distance)

    found <- .run_flann_kdtree(
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
setMethod("findDistanceFromIndex", "FlannKdtreeIndex", function(BNINDEX, k, num.threads=1, subset=NULL, ...) {
    found <- findKnnFromIndex(BNINDEX, k=k, get.index=FALSE, get.distance=TRUE, num.threads=num.threads, subset=subset)
    if (is(k, "AsIs")) {
        vapply(found$distance, function(x) if (length(x)) x[length(x)] else NA_real_, 0)
    } else if (ncol(found$distance)) {
        found$distance[,ncol(found$distance)]
    } else {
        rep(NA_real_, nrow(found$distance))
    }
})

#' @export
#' @rdname queryDistance
setMethod("queryDistanceFromIndex", "FlannKdtreeIndex", function(
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

    if (is(k, "AsIs")) {
        vapply(found$distance, function(x) if (length(x)) x[length(x)] else NA_real_, 0)
    } else if (ncol(found$distance)) {
        found$distance[,ncol(found$distance)]
    } else {
        rep(NA_real_, nrow(found$distance))
    }
})

#' @export
#' @rdname findNeighbors
setMethod("findNeighborsFromIndex", "FlannKdtreeIndex", function(BNINDEX, threshold, get.index=TRUE, get.distance=TRUE, num.threads=1, subset=NULL, ...) {
    stop("range-based searches are not supported for FlannKdtreeParam()")
})

#' @export
#' @rdname queryNeighbors
setMethod("queryNeighborsFromIndex", "FlannKdtreeIndex", function(BNINDEX, query, threshold, get.index=TRUE, get.distance=TRUE, num.threads=1, subset=NULL, transposed=FALSE, ..., .check.nonfinite=TRUE) {
    stop("range-based searches are not supported for FlannKdtreeParam()")
})

#' @export
#' @rdname saveIndex
setMethod("saveIndex", "FlannKdtreeIndex", function(BNINDEX, dir, ...) {
    writeLines(.flann_kdtree_algorithm_name, file.path(dir, "ALGORITHM"))
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

.load_flann_kdtree_index <- function(dir, ...) {
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
