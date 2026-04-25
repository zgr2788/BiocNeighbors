#' The RpforestParam class
#'
#' A class to hold parameters for approximate nearest-neighbor searches with a
#' random projection forest.
#'
#' @param nneighbors Integer scalar specifying the number of neighbors to cache
#' for self-searches on the built index.
#' @param ntrees Integer scalar specifying the number of random projection trees
#' used in the forest.
#' @param leaf.size Integer scalar specifying the target leaf size during tree
#' construction.
#' @param max.tree.depth Integer scalar specifying the maximum tree depth.
#' @inheritParams ExhaustiveParam
#'
#' @details
#' This backend uses the \pkg{rnndescent} package to build a random projection
#' forest that can be queried for approximate nearest neighbors.
#'
#' The \code{nneighbors} parameter controls the size of the cached neighbor
#' graph used for self-searches on the built index. Larger values improve the
#' quality of cached \code{\link{findKNN}} results at the cost of more work
#' during index construction.
#'
#' As with the existing cosine implementations in \pkg{BiocNeighbors}, cosine
#' distance is handled by L2-normalizing each observation and running the search
#' with Euclidean distances on the normalized data.
#'
#' @return
#' The \code{RpforestParam} constructor returns an instance of the
#' RpforestParam class.
#'
#' The \code{\link{buildIndex}} method returns an instance of the
#' \code{RpforestIndex} class.
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
#' (out <- RpforestParam())
#' out[['ntrees']]
#'
#' out[['ntrees']] <- 20L
#' out
#'
#' @aliases
#' RpforestParam-class
#' show,RpforestParam-method
#' RpforestIndex
#' RpforestIndex-class
#'
#' @docType class
#'
#' @export
#' @importFrom methods new
RpforestParam <- function(nneighbors=20, ntrees=10, leaf.size=20, max.tree.depth=200, distance=c("Euclidean", "Manhattan", "Cosine")) {
    new(
        "RpforestParam",
        nneighbors=as.integer(nneighbors),
        ntrees=as.integer(ntrees),
        leaf.size=as.integer(leaf.size),
        max.tree.depth=as.integer(max.tree.depth),
        distance=match.arg(distance)
    )
}

setValidity("RpforestParam", function(object) {
    msg <- character(0)

    nneighbors <- object[["nneighbors"]]
    if (length(nneighbors) != 1L || is.na(nneighbors) || nneighbors <= 0L) {
        msg <- c(msg, "'nneighbors' should be a positive integer scalar")
    }

    ntrees <- object[["ntrees"]]
    if (length(ntrees) != 1L || is.na(ntrees) || ntrees <= 0L) {
        msg <- c(msg, "'ntrees' should be a positive integer scalar")
    }

    leaf.size <- object[["leaf.size"]]
    if (length(leaf.size) != 1L || is.na(leaf.size) || leaf.size <= 0L) {
        msg <- c(msg, "'leaf.size' should be a positive integer scalar")
    }

    max.tree.depth <- object[["max.tree.depth"]]
    if (length(max.tree.depth) != 1L || is.na(max.tree.depth) || max.tree.depth <= 0L) {
        msg <- c(msg, "'max.tree.depth' should be a positive integer scalar")
    }

    if (length(msg)) return(msg)
    TRUE
})

#' @export
setMethod("show", "RpforestParam", function(object) {
    callNextMethod()
    cat(sprintf("nneighbors: %i\n", object[["nneighbors"]]))
    cat(sprintf("ntrees: %i\n", object[["ntrees"]]))
    cat(sprintf("leaf.size: %i\n", object[["leaf.size"]]))
    cat(sprintf("max.tree.depth: %i\n", object[["max.tree.depth"]]))
})

#' @export
RpforestIndex <- function(index, data, names, param) {
    new("RpforestIndex", index=index, data=data, names=names, param=param)
}

.rpforest_algorithm_name <- "BiocNeighbors::Rpforest"

.fallback_rpforest_query <- function(reference, query, k, distance, num.threads=1L) {
    out <- queryKNN(
        reference,
        query=query,
        k=k,
        BNPARAM=KmknnParam(distance=distance),
        num.threads=num.threads
    )
    list(idx=out$index, dist=out$distance)
}

.coerce_rpforest_matrix <- function(x, nobs, k) {
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

.standardize_rpforest_hits <- function(found, reference, query, k, distance, num.threads=1L) {
    nobs <- nrow(query)
    if (!k) {
        return(list(
            idx=matrix(integer(0), nobs, 0),
            dist=matrix(numeric(0), nobs, 0)
        ))
    }

    idx <- .coerce_rpforest_matrix(found$idx, nobs, k)
    dist <- .coerce_rpforest_matrix(found$dist, nobs, k)

    if (is.null(idx) || is.null(dist)) {
        return(.fallback_rpforest_query(reference, query, k, distance=distance, num.threads=num.threads))
    }

    bad.index <- anyNA(idx) || any(idx < 1L) || any(idx > nrow(reference))
    bad.distance <- anyNA(dist) || any(!is.finite(dist))
    if (bad.index || bad.distance) {
        return(.fallback_rpforest_query(reference, query, k, distance=distance, num.threads=num.threads))
    }

    list(idx=idx, dist=dist)
}

.extract_rpforest_graph <- function(BNINDEX) {
    graph <- BNINDEX@index$graph
    if (!is.list(graph) || is.null(graph$idx) || is.null(graph$dist)) {
        return(NULL)
    }
    graph
}

.can_reuse_rpforest_graph <- function(BNINDEX, k) {
    graph <- .extract_rpforest_graph(BNINDEX)
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

.rpforest_graph_self <- function(BNINDEX, ids, k) {
    graph <- .extract_rpforest_graph(BNINDEX)
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

.rpforest_query_self <- function(BNINDEX, ids, k, num.threads) {
    if (!length(ids)) {
        return(list(index=matrix(integer(0), 0, k), distance=matrix(numeric(0), 0, k)))
    }

    if (!k) {
        return(list(index=matrix(integer(0), length(ids), 0), distance=matrix(numeric(0), length(ids), 0)))
    }

    query <- BNINDEX@data[ids,,drop=FALSE]
    found <- rnndescent::rpf_knn_query(
        query=query,
        reference=BNINDEX@data,
        forest=BNINDEX@index$forest,
        k=min(k + 1L, nrow(BNINDEX@data)),
        n_threads=as.integer(num.threads),
        verbose=FALSE
    )
    found <- .standardize_rpforest_hits(
        found,
        reference=BNINDEX@data,
        query=query,
        k=min(k + 1L, nrow(BNINDEX@data)),
        distance=bndistance(BNINDEX@param),
        num.threads=num.threads
    )

    .drop_self_matches(found$idx, found$dist, ids, k)
}

#' @export
#' @rdname buildIndex
setMethod("buildIndex", "RpforestParam", function(X, BNPARAM, transposed=FALSE, num.threads=1, BPPARAM=NULL, ..., .check.nonfinite=TRUE) {
    .require_rnndescent()
    num.threads <- .resolve_num_threads(num.threads, BPPARAM)

    X <- .coerce_observation_rows(X, transposed=transposed)
    if (.check.nonfinite && any(!is.finite(as.matrix(X)))) {
        stop("cannot build an index from non-finite values")
    }

    processed <- .prepare_nndescent_reference(X, BNPARAM)
    safe.k <- .cap_nndescent_k(BNPARAM@nneighbors, max(nrow(processed$data) - 1L, 0L))

    if (safe.k > 0L) {
        built <- rnndescent::rpf_knn(
            processed$data,
            k=safe.k,
            metric=processed$metric,
            n_trees=BNPARAM@ntrees,
            leaf_size=BNPARAM@leaf.size,
            max_tree_depth=BNPARAM@max.tree.depth,
            include_self=FALSE,
            ret_forest=TRUE,
            n_threads=num.threads,
            verbose=FALSE
        )
        forest <- built$forest
        graph <- .standardize_rpforest_hits(
            built,
            reference=processed$data,
            query=processed$data,
            k=safe.k,
            distance=bndistance(BNPARAM),
            num.threads=num.threads
        )
        idx <- list(
            forest=forest,
            graph=list(idx=graph$idx, dist=graph$dist),
            metric=processed$metric
        )
    } else if (nrow(processed$data) > 1L) {
        idx <- list(
            forest=rnndescent::rpf_build(
                processed$data,
                metric=processed$metric,
                n_trees=BNPARAM@ntrees,
                leaf_size=BNPARAM@leaf.size,
                max_tree_depth=BNPARAM@max.tree.depth,
                n_threads=num.threads,
                verbose=FALSE
            ),
            graph=list(
                idx=matrix(integer(0), nrow(processed$data), 0),
                dist=matrix(numeric(0), nrow(processed$data), 0)
            ),
            metric=processed$metric
        )
    } else {
        idx <- list(
            forest=NULL,
            graph=list(
                idx=matrix(integer(0), nrow(processed$data), 0),
                dist=matrix(numeric(0), nrow(processed$data), 0)
            ),
            metric=processed$metric,
            empty=TRUE
        )
    }

    RpforestIndex(index=idx, data=processed$data, names=rownames(processed$data), param=BNPARAM)
})

#' @export
#' @rdname findKNN
setMethod("findKnnFromIndex", "RpforestIndex", function(BNINDEX, k, get.index=TRUE, get.distance=TRUE, num.threads=1, subset=NULL, ...) {
    .require_rnndescent()

    chosen <- .subset_nndescent_index(BNINDEX, subset)
    validated <- .validate_and_cap_k(k, length(chosen), max(nrow(BNINDEX@data) - 1L, 0L))
    k <- validated$k
    variable <- validated$variable

    if (variable) {
        max.k <- if (length(k)) max(k) else 0L
        if (.can_reuse_rpforest_graph(BNINDEX, k)) {
            found <- .rpforest_graph_self(BNINDEX, chosen, max.k)
        } else {
            found <- .rpforest_query_self(BNINDEX, chosen, max.k, num.threads)
        }
        sliced <- .extract_variable_neighbors(found$index, found$distance, k)
        .format_nndescent_output(sliced$index, sliced$distance, get.index, get.distance, variable=TRUE)
    } else {
        if (.can_reuse_rpforest_graph(BNINDEX, k)) {
            found <- .rpforest_graph_self(BNINDEX, chosen, as.integer(k))
        } else {
            found <- .rpforest_query_self(BNINDEX, chosen, as.integer(k), num.threads)
        }
        found <- .orient_nndescent_output(found)
        .format_nndescent_output(found$index, found$distance, get.index, get.distance, variable=FALSE)
    }
})

#' @export
#' @rdname queryKNN
setMethod("queryKnnFromIndex", "RpforestIndex", function(
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
    } else if (is.null(BNINDEX@index$forest)) {
        found <- .fallback_rpforest_query(
            reference=BNINDEX@data,
            query=query,
            k=max.k,
            distance=bndistance(BNINDEX@param),
            num.threads=num.threads
        )
        found <- list(index=found$idx, distance=found$dist)
    } else {
        found <- rnndescent::rpf_knn_query(
            query=query,
            reference=BNINDEX@data,
            forest=BNINDEX@index$forest,
            k=max.k,
            n_threads=as.integer(num.threads),
            verbose=FALSE
        )
        found <- .standardize_rpforest_hits(
            found,
            reference=BNINDEX@data,
            query=query,
            k=max.k,
            distance=bndistance(BNINDEX@param),
            num.threads=num.threads
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
setMethod("findDistanceFromIndex", "RpforestIndex", function(BNINDEX, k, num.threads=1, subset=NULL, ...) {
    found <- findKnnFromIndex(BNINDEX, k=k, get.index=FALSE, get.distance=TRUE, num.threads=num.threads, subset=subset)
    .last_distance_from_knn(found, k)
})

#' @export
#' @rdname queryDistance
setMethod("queryDistanceFromIndex", "RpforestIndex", function(
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
setMethod("findNeighborsFromIndex", "RpforestIndex", function(BNINDEX, threshold, get.index=TRUE, get.distance=TRUE, num.threads=1, subset=NULL, ...) {
    stop("range-based searches are not supported for RpforestParam()")
})

#' @export
#' @rdname queryNeighbors
setMethod("queryNeighborsFromIndex", "RpforestIndex", function(BNINDEX, query, threshold, get.index=TRUE, get.distance=TRUE, num.threads=1, subset=NULL, transposed=FALSE, ..., .check.nonfinite=TRUE) {
    stop("range-based searches are not supported for RpforestParam()")
})

#' @export
#' @rdname saveIndex
setMethod("saveIndex", "RpforestIndex", function(BNINDEX, dir, ...) {
    writeLines(.rpforest_algorithm_name, file.path(dir, "ALGORITHM"))
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

.load_rpforest_index <- function(dir, ...) {
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
