#' @title Persistence Landscapes
#' @description Compute persistence landscapes from persistence data.
#'
#' @details `pl_new()` is a wrapper around the S4 class constructor
#'   `[methods:new()]`. The `pl_*()` helper functions query a persistence
#'   landscape as returned by `pl_new()` for specific information or manipulate
#'   its internal representation.
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
#' @name pl_new
#' @include plt-package.R
#' @include PersistenceLandscape.R
#' @param x Persistence data (or diagram), stored as a 2-column matrix, as a
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
#' @param pl A persistence landscape as returned by `pl_new()`.
#' @return `pl_new()` returns a persistence landscape (an object of S4 class
#'   'Rcpp_PersistenceLandscape'). Other functions return summary information
#'   about such an object.
#' @seealso [Rcpp_PersistenceLandscape-class] for the exported C++ class.
#' @example inst/examples/ex-landscape-exact.R
#' @example inst/examples/ex-landscape-discrete.R
#' @export
pl_new <- function(
    x, degree = NULL,
    exact = FALSE, xmin = NULL, xmax = NULL, xby = NULL
) UseMethod("pl_new")

#' @export
pl_new.matrix <- function(
    x, degree = NULL,
    exact = FALSE, xmin = NULL, xmax = NULL, xby = NULL
) {
  
  # check format
  if (! is.numeric(x) || ncol(x) < 2L || ncol(x) > 3L)
    stop("`x` must be a 2- or 3-column numeric matrix.")
  if (ncol(x) == 2L && ! is.null(degree)) {
    warning("`x` has only 2 columns; `degree` will be ignored.")
  }
  if (ncol(x) == 3L && is.null(degree))
    stop("Please specify a homological degree.")
  
  # remove infinities
  x <- x[! apply(x, 1L, \(f) any(is.infinite(f))), , drop = FALSE]
  
  # check content
  if (is.null(x) || (nrow(x) > 0L && all(is.na(x))))
    stop("`pl_new()` requires non-missing finite persistence data.")
  
  # infer any missing parameters
  xmin <- xmin %||% 0
  xmax <- xmax %||% if (nrow(x) == 0L) 1 else max(x)
  # grid of between 100 and 1000 intervals of length a power of 10
  # TODO: Make this a non-default global option.
  xby <- xby %||% (10 ^ (floor(log(xmax - xmin, 10)) - 2L))
  
  # construct persistence landscape
  new(PersistenceLandscape, x, exact, xmin, xmax, xby)
}

#' @export
pl_new.persistence <- function(x, degree = NULL, ...) {
  
  if (is.null(degree))
    stop("Please specify a homological degree.")
  
  x <- if (length(x$pairs) < degree + 1L) {
    # TODO: fix this on the C++ side
    # matrix(NA_real_, nrow = 0L, ncol = 2L)
    matrix(0, nrow = 1L, ncol = 2L)
  } else {
    x$pairs[[degree + 1L]]
  }
  
  pl_new.matrix(x, degree = NULL, ...)
}

#' @export
pl_new.default <- function(x, ...) {
  
  x <- try(as_persistence(x))
  if (inherits(x, "try-error")) {
    stop("No `as_persistence()` method for class '", class(x)[[1L]], "'.")
  }
  
  pl_new.persistence(x, ...)
}

#' @rdname pl_new
#' @export
pl_is_exact <- function(pl) {
  stopifnot(inherits(pl, "Rcpp_PersistenceLandscape"))
  pl$isExact()
}

#' @rdname pl_new
#' @export
pl_type <- function(pl) {
  stopifnot(inherits(pl, "Rcpp_PersistenceLandscape"))
  # if (is.atomic(pl$getInternal())) "discrete" else "exact"
  if (pl$isExact()) "exact" else "discrete"
}

#' @rdname pl_new
#' @export
pl_num_levels <- function(pl) {
  stopifnot(inherits(pl, "Rcpp_PersistenceLandscape"))
  switch (
    pl_type(pl),
    exact = return(length(pl$getInternal())),
    discrete = return(dim(pl$getInternal())[[1L]])
  )
}

#' @rdname pl_new
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

#' @rdname pl_new
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

#' @rdname pl_new
#' @export
pl_delimit <- function(pl, xmin = NULL, xmax = NULL, xby = NULL) {
  if (is.null(xmin)) xmin <- pl$xMin()
  if (is.null(xmax)) xmax <- pl$xMax()
  if (is.null(xby)) xby <- pl$xBy()
  pl$delimit(xmin, xmax, xby)
}

#' @rdname pl_new
#' @export
pl_discretize <- function(pl) pl$discretize()
