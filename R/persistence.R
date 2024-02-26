#' @title Efficient Format for Persistence Data (Diagrams)
#' @description Transform persistence data to the format used by **plt**.
#'
#' @details
#'
#' This class and function will soon be spun off to a lower-level package that
#' can be imported by others as needed. Refer to
#' <http://adv-r.had.co.nz/S3.html> for guidance on methods to make available.
#'
#' @name as_persistence
#' @aliases persistence diagram persistence_diagram
#' @param x Persistence data to be transformed; either a \eqn{\geq 3}-column
#'   matrix (or object coercible to one) with dimension/degree, start/birth, and
#'   end/death columns, a list whose first element is such an object, an object
#'   of class 'PHom' as returned by `[ripserr::vietoris_rips()]`, or (a list as
#'   returned by a `*Diag()` function in **TDA** (e.g. `[TDA::ripsDiag()]`)
#'   whose first element is) an object of class 'diagram'.
#' @param ... Parameters passed to methods.
#' @param degree Non-negative integer; the homology degree for which to recover
#'   a matrix of persistence pairs.
#' @param modulus,max_dim,threshold Possibly missing parameters of the
#'   calculation of the persistence data `x` to be included in the 'persistence'
#'   object.
#' @inheritParams base::as.data.frame

#' @return Persistence data in the form of a list with named entries:
#'   \item{`pairs`}{
#'   A list of 2-column matrices containing birth-death pairs
#'   in each dimension/degree, starting with zero.
#'   }
#'   \item{`modulus`}{
#'   The (integer) modulus of the prime field in which
#'   homology coefficientts were calculated.
#'   }
#'   \item{`max_dim`}{
#'   The (integer) maximum dimension of calculated features.
#'   }
#'   \item{`threshold`}{
#'   The (real) maximum length of an edge in the simplicial filtration.
#'   }

#' @example inst/examples/ex-persistence.r
#' @export
as_persistence <- function(x, ...) UseMethod("as_persistence")

#' @rdname as_persistence
#' @export
as_persistence.default <- function(
    x, modulus = NULL, max_dim = NULL, threshold = NULL, ...
) {
  x <- as.matrix(x)
  
  # initialize list
  pd <- list()
  pd$pairs <- list()
  
  # concatenate a matrix for each dimension
  if (is.null(max_dim)) max_dim <- if (nrow(x) == 0L) {
    -Inf
  } else {
    max(x[, 1L, drop = TRUE])
  }
  if (max_dim >= 0) for (i in seq(0L, max_dim)) {
    pd$pairs[[i + 1L]] <- x[x[, 1L] == i, c(2L, 3L), drop = FALSE]
    dimnames(pd$pairs[[i + 1L]]) <- NULL
  }
  
  # assign parameters
  pd$modulus <- modulus %||% NA_integer_
  pd$max_dim <- max_dim %||% as.integer(max(x[, 1L]))
  pd$threshold <- threshold %||% NA_real_
  
  class(pd) <- "persistence"
  pd
}

#' @rdname as_persistence
#' @export
as_persistence.persistence <- function(x, ...) x

#' @rdname as_persistence
#' @export
as_persistence.list <- function(x, ...) {
  # if list contains one object, reroute to default method
  # (handles `TDA::*Diag()` output with 'diagram' object as list singleton)
  # TODO: check all possible outputs of `TDA::*Diag()`
  if (length(x) == 1L) {
    
    return(as_persistence(x[[1L]], ...))
    
  } else if (all(vapply(x, is.matrix, NA)) && all(vapply(x, ncol, 0L) == 2L)) {
    
    # collect recognized parameters
    dots <- list(...)
    params <- dots[intersect(c("modulus", "max_dim", "threshold"), names(dots))]
    
    # interpret as pairs
    pd <- c(list(pairs = x), params)
    
    class(pd) <- "persistence"
    pd
  }
}

#' @rdname as_persistence
#' @export
as_persistence.diagram <- function(x, ...) {
  
  # get 3-column matrix
  m <- matrix(as.vector(x), nrow = dim(x)[[1L]], ncol = dim(x)[[2L]])
  
  # reconcile dots
  dots <- list(...)
  params <- list(
    x = m,
    max_dim = attr(x, "call")$maxdimension %||% attr(x, "maxdimension"),
    threshold = attr(x, "call")$maxscale %||% {
      # discern maximum scale
      x_scale <- attr(x, "scale")
      if (! is.null(x_scale) && length(x_scale) == 2L) x_scale <- x_scale[[2L]]
      x_scale
    }
  )
  
  # reroute to default method with extracted parameters
  do.call(as_persistence.default, utils::modifyList(params, dots))
}

#' @rdname as_persistence
#' @export
as_persistence.PHom <- function(x, ...) {
  # coerce to matrix and reroute to default method
  # (will need to change if 'PHom' class changes)
  as_persistence.default(as.matrix(as.data.frame(x)), ...)
}

#' @rdname as_persistence
#' @export
print.persistence <- function(x, ...) {
  cat(format(x, ...), "\n")
  invisible(x)
}

#' @rdname as_persistence
#' @export
format.persistence <- function(x, ...) {
  # empty persistence data
  if (pd$max_dim == -Inf && length(pd$pairs) == 0L)
    return("empty 'persistence' data")
  
  # parameters
  if (is.na(x$max_dim)) x$max_dim <- length(x$pairs) - 1L
  
  # header
  fmt1 <- sprintf(
    "'persistence' data computed up to degree %i:",
    x$max_dim %||% "<unknown>"
  )
  
  # features
  num_deg <- length(x$pairs) - 1L
  num_feat <- vapply(x$pairs, nrow, 0L)[seq(min(6L, num_deg + 1L))]
  wid_feat <- floor(log(max(num_feat), 10)) + 1L
  lines2 <- c(
    paste(
      "* ",
      format(seq_along(num_feat) - 1L, width = 1L),
      "-degree features: ",
      format(num_feat, width = wid_feat),
      sep = ""
    ),
    if (num_deg > 6L) "..."
  )
  fmt2 <- paste(lines2, collapse = "\n")
  
  # other parameters, if any
  fmt3 <- do.call(
    paste,
    c(
      if (! is.na(x$modulus)) list(sprintf("modulus = %i", x$modulus)),
      if (! is.na(x$threshold)) list(sprintf("threshold = %f", x$threshold)),
      list(sep = "; ")
    )
  )
  
  if (length(fmt3) == 0L) {
    paste(fmt1, fmt2, sep = "\n\n")
  } else {
    paste(fmt1, fmt2, fmt3, sep = "\n\n")
  }
}

# `get_pairs()`

#' @rdname as_persistence
#' @export
get_pairs <- function(x, degree, ...) {
  stopifnot(inherits(x, "persistence"))
  if (length(x$pairs) > degree) {
    x$pairs[[degree + 1L]]
  } else {
    matrix(NA_real_, nrow = 0L, ncol = 2L)
  }
}

# `as.data.frame()`

#' @rdname as_persistence
#' @export
as.data.frame.persistence <- function(
    x, row.names = NULL, optional = TRUE, ...
) {
  features <- vapply(x$pairs, nrow, 0L)
  degree <- rep(seq_along(x$pairs) - 1L, features)
  pairs <- Reduce(rbind, x$pairs)
  df <- data.frame(degree, pairs)
  names(df) <- c("degree", "birth", "death")
  if (! is.null(row.names)) {
    rownames(df) <- row.names
  } else if (! optional) {
    id <- unlist(lapply(features, seq))
    rownames(df) <- paste(degree, id, sep = ".")
  }
  df
}
