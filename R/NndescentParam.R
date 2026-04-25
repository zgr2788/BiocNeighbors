#' The NndescentParam class
#'
#' A class to hold parameters for approximate nearest-neighbor searches with the
#' nearest-neighbor descent algorithm.
#'
#' @param nneighbors Integer scalar specifying the number of neighbors used to
#' build the search graph.
#' @param ntrees Integer scalar specifying the number of random projection trees
#' used for initialization.
#' @param max.candidates Integer scalar specifying the maximum number of
#' candidate neighbors to consider for each point in each iteration.
#' @param niter Integer scalar specifying the number of refinement iterations.
#' @inheritParams ExhaustiveParam
#'
#' @details
#' This backend uses the \pkg{rnndescent} package to build a queryable search
#' graph. The \code{nneighbors} parameter controls the quality of the built
#' graph and acts as a lower bound on the neighborhood size that the index is
#' optimized for. Larger values improve recall at the cost of more memory and a
#' slower build.
#'
#' As with the existing cosine implementations in \pkg{BiocNeighbors}, cosine
#' distance is handled by L2-normalizing each observation and running the search
#' with Euclidean distances on the normalized data.
#'
#' @return
#' The \code{NndescentParam} constructor returns an instance of the
#' NndescentParam class.
#'
#' The \code{\link{buildIndex}} method returns an instance of the
#' \code{NndescentIndex} class.
#'
#' @author
#' Oz Beker
#'
#' @seealso
#' \linkS4class{BiocNeighborParam}, for the parent class and its available
#' methods.
#'
#' \url{https://search.r-project.org/CRAN/refmans/rnndescent/html/rnndescent-package.html},
#' for details on the underlying implementation.
#'
#' @examples
#' (out <- NndescentParam())
#' out[['nneighbors']]
#'
#' out[['nneighbors']] <- 40L
#' out
#'
#' @aliases
#' NndescentParam-class
#' show,NndescentParam-method
#' NndescentIndex
#' NndescentIndex-class
#'
#' @docType class
#'
#' @export
#' @importFrom methods new
NndescentParam <- function(nneighbors=20, ntrees=3, max.candidates=20, niter=3, distance=c("Euclidean", "Manhattan", "Cosine")) {
    new(
        "NndescentParam",
        nneighbors=as.integer(nneighbors),
        ntrees=as.integer(ntrees),
        max.candidates=as.integer(max.candidates),
        niter=as.integer(niter),
        distance=match.arg(distance)
    )
}

setValidity("NndescentParam", function(object) {
    msg <- character(0)

    nneighbors <- object[["nneighbors"]]
    if (length(nneighbors) != 1L || is.na(nneighbors) || nneighbors <= 0L) {
        msg <- c(msg, "'nneighbors' should be a positive integer scalar")
    }

    ntrees <- object[["ntrees"]]
    if (length(ntrees) != 1L || is.na(ntrees) || ntrees <= 0L) {
        msg <- c(msg, "'ntrees' should be a positive integer scalar")
    }

    max.candidates <- object[["max.candidates"]]
    if (length(max.candidates) != 1L || is.na(max.candidates) || max.candidates <= 0L) {
        msg <- c(msg, "'max.candidates' should be a positive integer scalar")
    }

    niter <- object[["niter"]]
    if (length(niter) != 1L || is.na(niter) || niter <= 0L) {
        msg <- c(msg, "'niter' should be a positive integer scalar")
    }

    if (length(msg)) return(msg)
    TRUE
})

#' @export
setMethod("show", "NndescentParam", function(object) {
    callNextMethod()
    cat(sprintf("nneighbors: %i\n", object[["nneighbors"]]))
    cat(sprintf("ntrees: %i\n", object[["ntrees"]]))
    cat(sprintf("max.candidates: %i\n", object[["max.candidates"]]))
    cat(sprintf("niter: %i\n", object[["niter"]]))
})

#' @export
NndescentIndex <- function(index, data, names, param) {
    new("NndescentIndex", index=index, data=data, names=names, param=param)
}

.nndescent_algorithm_name <- "BiocNeighbors::Nndescent"

.require_rnndescent <- function() {
    if (!requireNamespace("rnndescent", quietly=TRUE)) {
        stop("the 'rnndescent' package must be installed to use NndescentParam()", call.=FALSE)
    }
}

.coerce_observation_rows <- function(x, transposed=FALSE, subset=NULL) {
    if (transposed) {
        if (is.matrix(x)) {
            x <- t(x)
        } else {
            x <- t(as.matrix(x))
        }
    } else if (!is.matrix(x) && !inherits(x, "dgCMatrix")) {
        x <- as.matrix(x)
    }

    if (!is.null(subset)) {
        x <- x[.validate_subset_indices(subset, nrow(x), rownames(x)),,drop=FALSE]
    }

    x
}

.l2_normalize_rows <- function(x) {
    x <- as.matrix(x)
    denom <- sqrt(rowSums(x^2))
    denom[denom == 0] <- 1
    x / denom
}

.prepare_nndescent_reference <- function(x, param) {
    metric <- tolower(bndistance(param))
    if (metric == "cosine") {
        list(data=.l2_normalize_rows(x), metric="euclidean")
    } else {
        list(data=x, metric=metric)
    }
}

.prepare_nndescent_query <- function(x, param) {
    if (bndistance(param) == "Cosine") {
        .l2_normalize_rows(x)
    } else {
        x
    }
}

.subset_nndescent_index <- function(index, subset) {
    .validate_subset_indices(subset, nrow(index@data), index@names)
}

.extract_nndescent_graph <- function(BNINDEX) {
    graph <- BNINDEX@index$graph
    if (!is.list(graph) || is.null(graph$idx) || is.null(graph$dist)) {
        return(NULL)
    }
    graph
}

.can_reuse_nndescent_graph <- function(BNINDEX, k) {
    graph <- .extract_nndescent_graph(BNINDEX)
    if (is.null(graph)) {
        return(FALSE)
    }

    nstored <- ncol(graph$idx)
    if (.is_variable_k(k)) {
        if (!length(k)) {
            return(TRUE)
        }
        max(k) <= nstored
    } else {
        as.integer(k) <= nstored
    }
}

.cap_nndescent_k <- function(k, limit) {
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

.extract_variable_neighbors <- function(index, distance, k) {
    list(
        index=lapply(seq_len(nrow(index)), function(i) index[i,seq_len(k[i]),drop=TRUE]),
        distance=lapply(seq_len(nrow(distance)), function(i) distance[i,seq_len(k[i]),drop=TRUE])
    )
}

.orient_nndescent_output <- function(found) {
    list(
        index=t(found$index),
        distance=t(found$distance)
    )
}

.format_nndescent_output <- function(index, distance, get.index, get.distance, variable=FALSE) {
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

.drop_self_matches <- function(index, distance, self.ids, keep) {
    nr <- nrow(index)
    out.index <- matrix(NA_integer_, nr, keep)
    out.distance <- matrix(NA_real_, nr, keep)

    for (i in seq_len(nr)) {
        chosen <- index[i,] != self.ids[i]
        current.index <- index[i,chosen]
        current.distance <- distance[i,chosen]
        if (length(current.index) < keep) {
            current.index <- current.index[seq_len(length(current.index))]
            current.distance <- current.distance[seq_len(length(current.distance))]
        } else {
            current.index <- current.index[seq_len(keep)]
            current.distance <- current.distance[seq_len(keep)]
        }

        if (keep) {
            out.index[i,seq_along(current.index)] <- current.index
            out.distance[i,seq_along(current.distance)] <- current.distance
        }
    }

    list(index=out.index, distance=out.distance)
}

.fallback_nndescent_query <- function(reference, query, k, distance, num.threads=1L) {
    out <- queryKNN(
        reference,
        query=query,
        k=k,
        BNPARAM=KmknnParam(distance=distance),
        num.threads=num.threads
    )
    list(index=out$index, distance=out$distance)
}

.nndescent_graph_self <- function(BNINDEX, ids, k) {
    graph <- .extract_nndescent_graph(BNINDEX)
    if (!length(ids)) {
        return(list(index=matrix(integer(0), 0, k), distance=matrix(numeric(0), 0, k)))
    }

    if (!k) {
        return(list(index=matrix(integer(0), length(ids), 0), distance=matrix(numeric(0), length(ids), 0)))
    }

    list(
        index=graph$idx[ids, seq_len(k), drop=FALSE],
        distance=graph$dist[ids, seq_len(k), drop=FALSE]
    )
}

.nndescent_query_self <- function(BNINDEX, ids, k, num.threads) {
    if (!length(ids)) {
        return(list(index=matrix(integer(0), 0, k), distance=matrix(numeric(0), 0, k)))
    }

    if (!k) {
        return(list(index=matrix(integer(0), length(ids), 0), distance=matrix(numeric(0), length(ids), 0)))
    }

    query <- BNINDEX@data[ids,,drop=FALSE]
    found <- rnndescent::rnnd_query(
        index=BNINDEX@index,
        query=query,
        k=min(k + 1L, nrow(BNINDEX@data)),
        n_threads=as.integer(num.threads),
        verbose=FALSE
    )

    .drop_self_matches(found$idx, found$dist, ids, k)
}

#' @export
#' @rdname buildIndex
setMethod("buildIndex", "NndescentParam", function(X, BNPARAM, transposed=FALSE, num.threads=1, BPPARAM=NULL, ..., .check.nonfinite=TRUE) {
    .require_rnndescent()
    num.threads <- .resolve_num_threads(num.threads, BPPARAM)

    X <- .coerce_observation_rows(X, transposed=transposed)
    if (.check.nonfinite && any(!is.finite(as.matrix(X)))) {
        stop("cannot build an index from non-finite values")
    }

    processed <- .prepare_nndescent_reference(X, BNPARAM)
    safe.k <- .cap_nndescent_k(BNPARAM@nneighbors, max(nrow(processed$data) - 1L, 0L))

    if (nrow(processed$data) <= 1L) {
        idx <- list(
            graph=list(
                idx=matrix(integer(0), nrow(processed$data), 0),
                dist=matrix(numeric(0), nrow(processed$data), 0)
            ),
            empty=TRUE
        )
    } else {
        idx <- rnndescent::rnnd_build(
            data=processed$data,
            k=safe.k,
            metric=processed$metric,
            n_trees=BNPARAM@ntrees,
            n_iters=BNPARAM@niter,
            max_candidates=BNPARAM@max.candidates,
            n_threads=num.threads,
            verbose=FALSE,
            progress="bar"
        )
    }

    NndescentIndex(index=idx, data=processed$data, names=rownames(processed$data), param=BNPARAM)
})

#' @export
#' @rdname findKNN
setMethod("findKnnFromIndex", "NndescentIndex", function(BNINDEX, k, get.index=TRUE, get.distance=TRUE, num.threads=1, subset=NULL, ...) {
    .require_rnndescent()

    chosen <- .subset_nndescent_index(BNINDEX, subset)
    validated <- .validate_and_cap_k(k, length(chosen), max(nrow(BNINDEX@data) - 1L, 0L))
    k <- validated$k
    variable <- validated$variable

    if (variable) {
        max.k <- if (length(k)) max(k) else 0L
        if (.can_reuse_nndescent_graph(BNINDEX, k)) {
            found <- .nndescent_graph_self(BNINDEX, chosen, max.k)
        } else {
            found <- .nndescent_query_self(BNINDEX, chosen, max.k, num.threads)
        }
        sliced <- .extract_variable_neighbors(found$index, found$distance, k)
        .format_nndescent_output(sliced$index, sliced$distance, get.index, get.distance, variable=TRUE)
    } else {
        if (.can_reuse_nndescent_graph(BNINDEX, k)) {
            found <- .nndescent_graph_self(BNINDEX, chosen, as.integer(k))
        } else {
            found <- .nndescent_query_self(BNINDEX, chosen, as.integer(k), num.threads)
        }
        found <- .orient_nndescent_output(found)
        .format_nndescent_output(found$index, found$distance, get.index, get.distance, variable=FALSE)
    }
})

#' @export
#' @rdname queryKNN
setMethod("queryKnnFromIndex", "NndescentIndex", function(
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
    .require_rnndescent()

    query <- .coerce_observation_rows(query, transposed=transposed, subset=subset)
    if (.check.nonfinite && any(!is.finite(as.matrix(query)))) {
        stop("cannot query an index with non-finite values")
    }

    query <- .prepare_nndescent_query(query, BNINDEX@param)
    validated <- .validate_and_cap_k(k, nrow(query), nrow(BNINDEX@data))
    k <- validated$k
    variable <- validated$variable
    max.k <- if (length(k)) max(k) else 0L

    if (!nrow(query) || !max.k) {
        nr <- nrow(query)
        found <- list(index=matrix(integer(0), nr, max.k), distance=matrix(numeric(0), nr, max.k))
    } else if (isTRUE(BNINDEX@index$empty)) {
        found <- .fallback_nndescent_query(
            reference=BNINDEX@data,
            query=query,
            k=max.k,
            distance=bndistance(BNINDEX@param),
            num.threads=num.threads
        )
    } else {
        found <- rnndescent::rnnd_query(
            index=BNINDEX@index,
            query=query,
            k=max.k,
            n_threads=as.integer(num.threads),
            verbose=FALSE
        )
        found <- list(index=found$idx, distance=found$dist)
    }

    if (variable) {
        sliced <- .extract_variable_neighbors(found$index, found$distance, k)
        .format_nndescent_output(sliced$index, sliced$distance, get.index, get.distance, variable=TRUE)
    } else {
        found <- .orient_nndescent_output(found)
        .format_nndescent_output(found$index, found$distance, get.index, get.distance, variable=FALSE)
    }
})

#' @export
#' @rdname findDistance
setMethod("findDistanceFromIndex", "NndescentIndex", function(BNINDEX, k, num.threads=1, subset=NULL, ...) {
    found <- findKnnFromIndex(BNINDEX, k=k, get.index=FALSE, get.distance=TRUE, num.threads=num.threads, subset=subset)
    .last_distance_from_knn(found, k)
})

#' @export
#' @rdname queryDistance
setMethod("queryDistanceFromIndex", "NndescentIndex", function(
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
setMethod("findNeighborsFromIndex", "NndescentIndex", function(BNINDEX, threshold, get.index=TRUE, get.distance=TRUE, num.threads=1, subset=NULL, ...) {
    stop("range-based searches are not supported for NndescentParam()")
})

#' @export
#' @rdname queryNeighbors
setMethod("queryNeighborsFromIndex", "NndescentIndex", function(BNINDEX, query, threshold, get.index=TRUE, get.distance=TRUE, num.threads=1, subset=NULL, transposed=FALSE, ..., .check.nonfinite=TRUE) {
    stop("range-based searches are not supported for NndescentParam()")
})

#' @export
#' @rdname saveIndex
setMethod("saveIndex", "NndescentIndex", function(BNINDEX, dir, ...) {
    writeLines(.nndescent_algorithm_name, file.path(dir, "ALGORITHM"))
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

.load_nndescent_index <- function(dir, ...) {
    .require_rnndescent()

    payload <- readRDS(file.path(dir, "index.rds"))
    data <- payload$data
    if (!is.null(payload$names)) {
        rownames(data) <- payload$names
    }

    out <- buildIndex(data, BNPARAM=payload$param, .check.nonfinite=FALSE)
    out@names <- payload$names
    out
}
