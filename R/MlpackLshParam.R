#' The MlpackLshParam class
#'
#' A class to hold parameters for approximate nearest-neighbor searches with
#' locality-sensitive hashing from \pkg{mlpack}.
#'
#' @param num.tables Integer scalar specifying the number of hash tables.
#' @param projections Integer scalar specifying the number of hash functions
#' per table.
#' @param bucket.size Integer scalar specifying the size of the second-level
#' hash buckets.
#' @param second.hash.size Integer scalar specifying the size of the
#' second-level hash table.
#' @param num.probes Integer scalar specifying the number of additional
#' buckets to probe during multi-probe LSH. A value of zero uses standard LSH.
#' @param search Optional integer scalar alias for \code{num.probes}.
#' @param seed Integer scalar specifying the random seed to pass to
#' \pkg{mlpack}. If \code{NA}, the mlpack default is used.
#' @inheritParams ExhaustiveParam
#'
#' @details
#' This backend uses \code{\link[mlpack]{lsh}} to query mlpack's approximate
#' nearest-neighbor implementation based on locality-sensitive hashing.
#' \code{\link{buildIndex}} stores the processed reference data and reissues
#' the mlpack search call for each query. This avoids relying on mlpack's R
#' model pointer serialization or reuse behavior.
#'
#' Only Euclidean searches are directly supported by mlpack LSH. As with the
#' existing cosine implementations in \pkg{BiocNeighbors}, cosine distance is
#' handled by L2-normalizing each observation and running Euclidean search on
#' the normalized data.
#'
#' Very small references are searched exactly, as LSH is not useful in that
#' regime and can trigger backend-specific edge behavior.
#'
#' @return
#' The \code{MlpackLshParam} constructor returns an instance of the
#' MlpackLshParam class.
#'
#' The \code{\link{buildIndex}} method returns an instance of the
#' \code{MlpackLshIndex} class.
#'
#' @author
#' Oz Beker
#'
#' @seealso
#' \linkS4class{BiocNeighborParam}, for the parent class and its available
#' methods.
#'
#' \url{https://mlpack.org/}, for details on the underlying implementation.
#'
#' @examples
#' (out <- MlpackLshParam())
#' out[['num.tables']]
#'
#' out[['num.tables']] <- 50L
#' out
#'
#' @aliases
#' MlpackLshParam-class
#' show,MlpackLshParam-method
#' MlpackLshIndex
#' MlpackLshIndex-class
#'
#' @docType class
#'
#' @export
#' @importFrom methods new
MlpackLshParam <- function(
    num.tables=30L,
    projections=10L,
    bucket.size=500L,
    second.hash.size=99901L,
    num.probes=0L,
    search=NULL,
    seed=NA_integer_,
    distance=c("Euclidean", "Cosine")
) {
    if (!is.null(search)) {
        if (!missing(num.probes) && !identical(as.integer(num.probes), as.integer(search))) {
            stop("'search' and 'num.probes' should not specify different values")
        }
        num.probes <- search
    }

    new(
        "MlpackLshParam",
        num.tables=as.integer(num.tables),
        projections=as.integer(projections),
        bucket.size=as.integer(bucket.size),
        second.hash.size=as.integer(second.hash.size),
        num.probes=as.integer(num.probes),
        seed=as.integer(seed),
        distance=match.arg(distance)
    )
}

setValidity("MlpackLshParam", function(object) {
    msg <- character(0)

    distance <- bndistance(object)
    if (length(distance) != 1L || !distance %in% c("Euclidean", "Cosine")) {
        msg <- c(msg, "'distance' should be either 'Euclidean' or 'Cosine'")
    }

    num.tables <- object[["num.tables"]]
    if (length(num.tables) != 1L || is.na(num.tables) || num.tables <= 0L) {
        msg <- c(msg, "'num.tables' should be a positive integer scalar")
    }

    projections <- object[["projections"]]
    if (length(projections) != 1L || is.na(projections) || projections <= 0L) {
        msg <- c(msg, "'projections' should be a positive integer scalar")
    }

    bucket.size <- object[["bucket.size"]]
    if (length(bucket.size) != 1L || is.na(bucket.size) || bucket.size <= 0L) {
        msg <- c(msg, "'bucket.size' should be a positive integer scalar")
    }

    second.hash.size <- object[["second.hash.size"]]
    if (length(second.hash.size) != 1L || is.na(second.hash.size) || second.hash.size <= 0L) {
        msg <- c(msg, "'second.hash.size' should be a positive integer scalar")
    }

    num.probes <- object[["num.probes"]]
    if (length(num.probes) != 1L || is.na(num.probes) || num.probes < 0L) {
        msg <- c(msg, "'num.probes' should be a non-negative integer scalar")
    }

    seed <- object[["seed"]]
    if (length(seed) != 1L || (!is.na(seed) && seed < 0L)) {
        msg <- c(msg, "'seed' should be NA or a non-negative integer scalar")
    }

    if (length(msg)) return(msg)
    TRUE
})

#' @export
setMethod("show", "MlpackLshParam", function(object) {
    callNextMethod()
    cat(sprintf("num.tables: %i\n", object[["num.tables"]]))
    cat(sprintf("projections: %i\n", object[["projections"]]))
    cat(sprintf("bucket.size: %i\n", object[["bucket.size"]]))
    cat(sprintf("second.hash.size: %i\n", object[["second.hash.size"]]))
    cat(sprintf("num.probes: %i\n", object[["num.probes"]]))
    cat(sprintf("seed: %s\n", if (is.na(object[["seed"]])) "NA" else as.character(object[["seed"]])))
})

#' @export
MlpackLshIndex <- function(index, data, names, param) {
    new("MlpackLshIndex", index=index, data=data, names=names, param=param)
}

.mlpack_lsh_algorithm_name <- "BiocNeighbors::MlpackLsh"

.require_mlpack <- function() {
    if (!requireNamespace("mlpack", quietly=TRUE)) {
        stop(
            "the 'mlpack' package must be installed to use MlpackLshParam(); ",
            "see https://mlpack.org/",
            call.=FALSE
        )
    }
}

.as_mlpack_points <- function(x) {
    as.matrix(x)
}

.mlpack_lsh_build_args <- function(param) {
    args <- list(
        bucket_size=param@bucket.size,
        projections=param@projections,
        reference=NULL,
        second_hash_size=param@second.hash.size,
        tables=param@num.tables,
        num_probes=param@num.probes,
        k=0L
    )

    if (!is.na(param@seed)) {
        args$seed <- param@seed
    }

    args
}

.mlpack_lsh_query_args <- function(param, k) {
    args <- list(
        k=as.integer(k),
        num_probes=param@num.probes
    )

    if (!is.na(param@seed)) {
        args$seed <- param@seed
    }

    args
}

.safe_mlpack_threads <- function(num.threads) {
    num.threads <- suppressWarnings(as.integer(num.threads)[1L])
    if (is.na(num.threads) || num.threads < 1L) {
        num.threads <- 1L
    }
    num.threads
}

.with_mlpack_threads <- function(num.threads, expr) {
    safe.threads <- .safe_mlpack_threads(num.threads)
    old <- Sys.getenv("OMP_NUM_THREADS", unset=NA_character_)

    Sys.setenv(OMP_NUM_THREADS=as.character(safe.threads))
    on.exit({
        if (is.na(old)) {
            Sys.unsetenv("OMP_NUM_THREADS")
        } else {
            Sys.setenv(OMP_NUM_THREADS=old)
        }
    })

    force(expr)
}

.call_mlpack_lsh <- function(args, num.threads=1L) {
    formals <- names(formals(mlpack::lsh))
    if ("num_threads" %in% formals && is.null(args$num_threads)) {
        args$num_threads <- .safe_mlpack_threads(num.threads)
    }

    .with_mlpack_threads(num.threads, do.call(mlpack::lsh, args))
}

.normalize_mlpack_lsh_indices <- function(index, nref) {
    dims <- dim(index)
    index <- matrix(as.integer(index), nrow=dims[1L], ncol=dims[2L])

    if (!length(index)) {
        return(index)
    }

    if (!anyNA(index) && min(index) == 0L && max(index) <= nref - 1L) {
        index <- index + 1L
    }

    index
}

.compute_mlpack_lsh_distances <- function(query, reference, index) {
    output <- matrix(NA_real_, nrow=nrow(index), ncol=ncol(index))

    if (!length(index)) {
        return(output)
    }

    for (i in seq_len(nrow(index))) {
        current <- index[i,]
        ok <- !is.na(current) & current >= 1L & current <= nrow(reference)
        if (any(ok)) {
            ref <- reference[current[ok],,drop=FALSE]
            output[i,ok] <- sqrt(rowSums((ref - matrix(query[i,], nrow=nrow(ref), ncol=ncol(ref), byrow=TRUE))^2))
        }
    }

    output
}

.order_mlpack_lsh_hits <- function(index, distance) {
    if (is.null(distance) || !length(index)) {
        return(list(index=index, distance=distance))
    }

    for (i in seq_len(nrow(index))) {
        o <- order(distance[i,], na.last=TRUE)
        index[i,] <- index[i,o]
        distance[i,] <- distance[i,o]
    }

    list(index=index, distance=distance)
}

.repair_mlpack_lsh_hits <- function(index, distance, query, reference, k, num.threads=1L) {
    if (!length(index)) {
        return(list(index=index, distance=distance))
    }

    bad <- !is.finite(index) | index < 1L | index > nrow(reference)
    if (!is.null(distance)) {
        bad <- bad | !is.finite(distance)
    }
    bad.rows <- which(rowSums(bad) > 0L)

    if (length(bad.rows)) {
        exact <- queryKNN(
            reference,
            query=query[bad.rows,,drop=FALSE],
            k=k,
            BNPARAM=KmknnParam(distance="Euclidean"),
            get.index=TRUE,
            get.distance=!is.null(distance),
            num.threads=num.threads
        )

        index[bad.rows,] <- exact$index
        if (!is.null(distance)) {
            distance[bad.rows,] <- exact$distance
        }
    }

    list(index=index, distance=distance)
}

.run_mlpack_lsh <- function(query, reference, model=NULL, param, k, num.threads=1L, get.distance=TRUE) {
    nobs <- nrow(query)
    if (!nobs || !k) {
        return(list(
            index=matrix(integer(0), nobs, k),
            distance=if (get.distance) matrix(numeric(0), nobs, k) else NULL
        ))
    }

    if (nrow(reference) <= 100L) {
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

    args <- .mlpack_lsh_build_args(param)
    args$k <- as.integer(k)
    args$query <- .as_mlpack_points(query)
    args$reference <- .as_mlpack_points(reference)

    found <- .call_mlpack_lsh(args, num.threads=num.threads)

    idx <- .coerce_flann_matrix(found$neighbors, nobs, k)
    if (is.null(idx)) {
        stop("failed to coerce mlpack LSH indices into the expected matrix shape")
    }
    idx <- .normalize_mlpack_lsh_indices(idx, nrow(reference))

    dist <- NULL
    if (get.distance) {
        dist <- .compute_mlpack_lsh_distances(query, reference, idx)
    }

    out <- .repair_mlpack_lsh_hits(
        index=idx,
        distance=dist,
        query=query,
        reference=reference,
        k=k,
        num.threads=num.threads
    )
    .order_mlpack_lsh_hits(out$index, out$distance)
}

#' @export
#' @rdname buildIndex
setMethod("buildIndex", "MlpackLshParam", function(X, BNPARAM, transposed=FALSE, num.threads=1, BPPARAM=NULL, ..., .check.nonfinite=TRUE) {
    .require_mlpack()
    num.threads <- .resolve_num_threads(num.threads, BPPARAM)

    X <- .coerce_flann_observation_rows(X, transposed=transposed)
    if (.check.nonfinite && any(!is.finite(as.matrix(X)))) {
        stop("cannot build an index from non-finite values")
    }

    processed <- .prepare_flann_reference(X, BNPARAM)
    MlpackLshIndex(index=NULL, data=processed, names=rownames(processed), param=BNPARAM)
})

#' @export
#' @rdname findKNN
setMethod("findKnnFromIndex", "MlpackLshIndex", function(BNINDEX, k, get.index=TRUE, get.distance=TRUE, num.threads=1, subset=NULL, ...) {
    .require_mlpack()

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
        found <- .run_mlpack_lsh(
            query=BNINDEX@data[chosen,,drop=FALSE],
            reference=BNINDEX@data,
            model=BNINDEX@index,
            param=BNINDEX@param,
            k=min(max.k + 1L, nrow(BNINDEX@data)),
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
setMethod("queryKnnFromIndex", "MlpackLshIndex", function(
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
    .require_mlpack()

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

    found <- .run_mlpack_lsh(
        query=query,
        reference=BNINDEX@data,
        model=BNINDEX@index,
        param=BNINDEX@param,
        k=max.k,
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
setMethod("findDistanceFromIndex", "MlpackLshIndex", function(BNINDEX, k, num.threads=1, subset=NULL, ...) {
    found <- findKnnFromIndex(BNINDEX, k=k, get.index=FALSE, get.distance=TRUE, num.threads=num.threads, subset=subset)
    .last_distance_from_knn(found, k)
})

#' @export
#' @rdname queryDistance
setMethod("queryDistanceFromIndex", "MlpackLshIndex", function(
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
setMethod("findNeighborsFromIndex", "MlpackLshIndex", function(BNINDEX, threshold, get.index=TRUE, get.distance=TRUE, num.threads=1, subset=NULL, ...) {
    stop("range-based searches are not supported for MlpackLshParam()")
})

#' @export
#' @rdname queryNeighbors
setMethod("queryNeighborsFromIndex", "MlpackLshIndex", function(BNINDEX, query, threshold, get.index=TRUE, get.distance=TRUE, num.threads=1, subset=NULL, transposed=FALSE, ..., .check.nonfinite=TRUE) {
    stop("range-based searches are not supported for MlpackLshParam()")
})

#' @export
#' @rdname saveIndex
setMethod("saveIndex", "MlpackLshIndex", function(BNINDEX, dir, ...) {
    writeLines(.mlpack_lsh_algorithm_name, file.path(dir, "ALGORITHM"))
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

.load_mlpack_lsh_index <- function(dir, ...) {
    .require_mlpack()

    payload <- readRDS(file.path(dir, "index.rds"))
    data <- payload$data
    if (!is.null(payload$names)) {
        rownames(data) <- payload$names
    }

    out <- buildIndex(data, BNPARAM=payload$param, .check.nonfinite=FALSE)
    out@names <- payload$names
    out
}
