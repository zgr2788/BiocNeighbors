.integerize_subset <- function(index, subset) {
    if (is.null(subset)) {
        subset
    } else if (!is.numeric(subset)) {
        dummy <- seq_len(generic_num_obs(index@ptr))
        names(dummy) <- index@names
        dummy[subset]
    } else if (!is.integer(subset)) {
        as.integer(subset)
    } else {
        subset
    }
}

.resolve_num_threads <- function(num.threads, BPPARAM=NULL) {
    if (!is.null(BPPARAM)) {
        num.threads <- BiocParallel::bpnworkers(BPPARAM)
    }

    as.integer(num.threads)
}

.is_variable_k <- function(k) {
    is(k, "AsIs") || length(k) != 1L
}

.validate_and_cap_k <- function(k, nobs, limit) {
    variable <- .is_variable_k(k)
    as.is <- is(k, "AsIs")
    converted <- tryCatch(suppressWarnings(as.integer(k)), error=function(e) NULL)
    if (is.null(converted)) {
        stop("'k' should be an integer scalar or vector")
    }

    if (!length(converted)) {
        if (variable && nobs == 0L) {
            return(list(k=if (as.is) I(converted) else converted, variable=variable))
        }
        stop("length of 'k' must be equal to the number of observations")
    }

    if (anyNA(converted)) {
        stop("'k' should not contain NA values")
    }
    if (any(converted < 0L)) {
        stop("'k' should be non-negative")
    }

    if (variable && length(converted) != nobs) {
        stop("length of 'k' must be equal to the number of observations")
    }

    limit <- as.integer(limit)
    capped <- pmin(converted, limit)
    if (any(converted != capped)) {
        warning("more neighbors were requested than available and will be ignored")
    }

    list(k=if (as.is) I(capped) else capped, variable=variable)
}

.validate_subset_indices <- function(subset, nobs, names=NULL) {
    if (is.null(subset)) {
        return(seq_len(nobs))
    }

    if (is.character(subset)) {
        if (is.null(names)) {
            stop("cannot use character 'subset' when observation names are not available")
        }
        chosen <- match(subset, names)
        if (anyNA(chosen)) {
            stop("failed to match some entries in 'subset' to observation names")
        }
        return(as.integer(chosen))
    }

    if (is.logical(subset)) {
        if (length(subset) != nobs) {
            stop("logical 'subset' should have length equal to the number of observations")
        }
        if (anyNA(subset)) {
            stop("'subset' should not contain NA values")
        }
        return(which(subset))
    }

    if (!is.numeric(subset)) {
        stop("'subset' should be an integer, logical, or character vector")
    }

    if (anyNA(subset) || any(!is.finite(subset))) {
        stop("'subset' should contain finite values")
    }

    chosen <- as.integer(subset)
    if (any(chosen < 1L | chosen > nobs)) {
        stop("'subset' contains out-of-range indices")
    }

    chosen
}

.last_distance_from_knn <- function(found, k) {
    if (.is_variable_k(k)) {
        vapply(found$distance, function(x) if (length(x)) x[length(x)] else NA_real_, 0)
    } else if (ncol(found$distance)) {
        found$distance[,ncol(found$distance)]
    } else {
        rep(NA_real_, nrow(found$distance))
    }
}

.format_output <- function(output, name, to.get) {
    if (isFALSE(to.get)) {
        output[[name]] <- NULL
    } else if (identical(to.get, "transposed")) {
        ; # no-op
    } else if (isTRUE(to.get) || identical(to.get, "normal")) {
        if (is.matrix(output[[name]])) {
            output[[name]] <- t(output[[name]])
        }
    } else {
        stop("unsupported option '", to.get, "'")
    }
    output
}

.transpose_and_subset <- function(x, transposed, subset) {
    if (is.matrix(x)) {
        # For back-compatibility, and to avoid loading beachmat and Matrix,
        # we just apply these operations directly if a matrix is provided.
        if (!transposed) {
            x <- t(x)
        }
    } else {
        x <- DelayedArray::DelayedArray(x)
        if (!transposed) {
            x <- Matrix::t(x)
        }
    }

    if (!is.null(subset)) {
        x <- x[,subset,drop=FALSE]
    }

    x
}
