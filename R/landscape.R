#' @title Persistence Landscapes
#' @description Compute persistence landscapes from persistence data.
#'
#' @details `landscape()` is a wrapper around the S4 class constructor
#'   `[methods:new()]`. The `pl_*()` helper functions query a persistence
#'   landscape as returned by `landscape()` for specific information or
#'   manipulate its internal representation.
#'
#'   Use `pl_is_exact()` and `pl_type()` to get a landscape's internal
#'   representation, `pl_num_levels()` its number of levels, `pl_limits()` the
#'   endpoints of its internal representation (excluding infinities), and
#'   `pl_support()` the infimum and supremum of its support (the points at which
#'   its value is nonzero).
#'   
#'   Use `pl_delimit()` to change the limits of a PL and `pl_discretize()` to
#'   convert an exact landscape to a discrete one (using its internally-stored
#'   range and resolution).
#'
#' @name landscape
#' @include plt-package.R
#' @include PersistenceLandscape.R
#' @param pd Persistence data (or diagram), stored as a 2-column matrix, as a
#'   '[persistence]' object, or in a format coercible to 'persistence'.
#' @param degree Non-negative integer; if input is a persistence diagram object,
#'   then the dimension for which to compute a landscape. (For degree \eqn{d},
#'   the \eqn{(d+1)}th matrix in the list will be selected.)
#' @param exact Set to `TRUE` for exact representation, `FALSE` (default) for
#'   discrete.
#' @param xmin,xmax Domain thresholds for discrete PL; if not specified, then
#'   taken to be the support of the PL constructed from the data or the internal
#'   values of the 'Rcpp_PersistenceLandscape' object.
#' @param xby Domain grid diameter for discrete PL; if not specified, then set
#'   to the power of 10 that yields between 100 and 1000 intervals.
#' @param pl A persistence landscape as returned by `landscape()`.
#' @return `landscape()` returns a persistence landscape (an object of S4 class
#'   'Rcpp_PersistenceLandscape'). Other functions return summary information
#'   about such an object.
#' @seealso [Rcpp_PersistenceLandscape-class] for the exported C++ class.
#' @example inst/examples/ex-landscape-exact.R
#' @example inst/examples/ex-landscape-discrete.R
#' @export
landscape <- function(
    pd, degree = NULL,
    exact = FALSE,
    xmin = NULL, xmax = NULL, xby = NULL
) {
  
  # birth-death pairs matrix `diagram`
  if (is.atomic(pd)) {
    stopifnot(ncol(pd) >= 2L, is.numeric(pd))
  } else {
    pd <- try(as_persistence(pd))
    if (inherits(pd, "try-error"))
      stop("There is no `as_persistence()` method for object `pd`.")
    if (is.null(degree))
      stop("`landscape()` requires a homological degree (`degree = <int>`).")
    pd <- if (length(pd$pairs) < degree + 1L) {
      # TODO: fix this on the C++ side
      # matrix(NA_real_, nrow = 0L, ncol = 2L)
      matrix(0, nrow = 1L, ncol = 2L)
    } else {
      pd$pairs[[degree + 1L]]
    }
  }
  
  # remove infinities
  pd <- pd[! apply(pd, 1L, \(f) any(is.infinite(f))), , drop = FALSE]
  
  # content check
  if (is.null(pd) || (nrow(pd) > 0L && all(is.na(pd)))) {
    stop("`landscape()` requires non-missing finite persistence data.")
  }
  
  # infer any missing parameters
  xmin <- xmin %||% 0
  xmax <- xmax %||% if (nrow(pd) == 0L) 1 else max(pd)
  # grid of between 100 and 1000 intervals of length a power of 10
  # TODO: Make this a non-default global option.
  xby <- xby %||% (10 ^ (floor(log(xmax - xmin, 10)) - 2L))
  
  # construct persistence landscape
  new(PersistenceLandscape, pd, exact, xmin, xmax, xby)
}

#' @rdname landscape
#' @export
pl_is_exact <- function(pl) {
  stopifnot(inherits(pl, "Rcpp_PersistenceLandscape"))
  pl$isExact()
}

#' @rdname landscape
#' @export
pl_type <- function(pl) {
  stopifnot(inherits(pl, "Rcpp_PersistenceLandscape"))
  # if (is.atomic(pl$getInternal())) "discrete" else "exact"
  if (pl$isExact()) "exact" else "discrete"
}

#' @rdname landscape
#' @export
pl_num_levels <- function(pl) {
  stopifnot(inherits(pl, "Rcpp_PersistenceLandscape"))
  switch (
    pl_type(pl),
    exact = return(length(pl$getInternal())),
    discrete = return(dim(pl$getInternal())[[1L]])
  )
}

#' @rdname landscape
#' @export
pl_limits <- function(pl) {
  stopifnot(inherits(pl, "Rcpp_PersistenceLandscape"))
  switch (
    pl_type(pl),
    exact = {
      xvec <- Reduce(rbind, pl$getInternal())[, 1L, drop = TRUE]
      return(range(xvec[! is.infinite(xvec)]))
    },
    discrete = {
      return(range(pl$getInternal()[1L, , 1L]))
    }
  )
}

#' @rdname landscape
#' @export
pl_support <- function(pl) {
  stopifnot(inherits(pl, "Rcpp_PersistenceLandscape"))
  switch (
    pl_type(pl),
    exact = {
      xvec <- Reduce(rbind, pl$getInternal())[, 1L, drop = TRUE]
      return(range(xvec[! is.infinite(xvec)]))
    },
    discrete = {
      nz <- which(apply(pl$getInternal()[, , 2L, drop = FALSE], 2L, max) > 0)
      supp <- intersect(
        Reduce(union, list(nz - 1L, nz, nz + 1L)),
        seq(dim(pl$getInternal())[[2L]])
      )
      return(range(pl$getInternal()[1L, supp, 1L]))
    }
  )
}

#' @rdname landscape
#' @export
pl_delimit <- function(pl, xmin = NULL, xmax = NULL, xby = NULL) {
  if (is.null(xmin)) xmin <- pl$xMin()
  if (is.null(xmax)) xmax <- pl$xMax()
  if (is.null(xby)) xby <- pl$xBy()
  pl$delimit(xmin, xmax, xby)
}

#' @rdname landscape
#' @export
pl_discretize <- function(pl) pl$discretize()
