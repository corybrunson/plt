#' @title Linear Algebra with Persistence Landscapes
#' @description Convert between persistence landscapes and their vectorizations,
#'   and perform linear algebra on vectorized persistence landscapes.
#'   
#' @details `pl_t()` extracts the common grid of \eqn{t}-values, discretizing
#' exact landscapes on the fly.
#'
#' `pl_to_vector()` pad/truncates levels and flattens a persistence landscape to a 
#' numeric vector; the grid is stored in the `"t_vals"` attribute.
#'
#' `pl_to_matrix()` – standardizes persistence landscapes by coordinating 
#' resolution, limits, and level count across multiple landscapes and then bind 
#' their vectorizations row-wise into a matrix (grid in `"t_vals"` attribute).
#'
#' `pl_from_vector()` rebuilds one landscape from a flattened vector + grid, 
#' optionally dropping empty levels.
#'
#' `pl_from_matrix()` applies `pl_from_vector()` row-wise to turn a matrix 
#' back into a list of landscapes.
#' 
#' @name algebra
#' @include arithmetic.R
#' @inheritParams pl_new
#' @inheritParams arithmetic
#' @param ... Persistence landscapes or lists of persistence landscapes.
#' @param x A numeric vector, optionally with an attribute `"t_vals"` whose
#'   length divides that of `x`, but with no other attributes.
#' @param t_vals A numeric vector of equally-spaced values assumed to be the
#'   grid support of a discrete persistence landscape.
#' @param num_levels Number of levels to vectorize; may me more or less than
#'   `pl_num_levels(pl)`.
#' @param drop_levels Logical; whether to omit levels that are empty in all PLs.
#' @returns `pl_t()` returns a numeric vector containing the grid of a 
#' discrete persistence landscape. `pl_to_vector()` returns a numeric 
#' vector with a `"t_vals"` attribute to reconstruct the discrete 
#' landscape. `pl_to_matrix()` returns a matrix of such vectorizations.  
#' `pl_from_vector()` returns a persistence landscape with class 
#' "Rcpp_PersistenceLandscape", and `pl_from_matrix()` returns a list of 
#' such persistence landscapes.
#' @seealso [arithmetic], [analysis], [inference] for other landscape functions.
#' @example inst/examples/ex-algebra.R
NULL

#' @rdname algebra
#' @export
pl_t <- function(pl) {
  pl_internal <- if (pl$isExact()) {
    pl$discretize()$getInternal()
  } else {
    pl$getInternal()
  }
  pl_internal[1L, , 1L]
}

#' @rdname algebra
#' @export
pl_to_vector <- function(pl, num_levels = pl_num_levels(pl)) {
  if (pl$isExact()) pl <- pl_discretize(pl)
  # get matrix of levels from discrete representation
  # x <- pl$getInternal()[, , 2L, drop = 3L]
  x <- pl$getInternal()[, , 2L, drop = FALSE]
  x <- array(x, dim = dim(x)[seq(2L)])
  # remove any superfluous levels
  if (num_levels < nrow(x)) x <- x[seq(num_levels), , drop = FALSE]
  # augment any additional empty levels
  if (num_levels > nrow(x))
    x <- rbind(x, matrix(0, nrow = num_levels - nrow(x), ncol = ncol(x)))
  # concatenate levels
  pl_vec <- as.vector(t(x), mode = "any")
  # assign grid as attribute
  attr(pl_vec, "t_vals") <- pl_t(pl)
  
  return(pl_vec)
}

#' @rdname algebra
#' @export
pl_to_matrix <- function(...) {
  
  # aggregate arguments into a single list of PLs
  args <- list(...)
  args <- unlist(args, recursive = TRUE)
  if (! all(vapply(args, inherits, FALSE, what = "Rcpp_PersistenceLandscape")))
    stop("`pl_to_matrix()` accepts only PLs and lists of PLs.", call. = FALSE)
  
  # unique resolution of discrete PLs or else first resolution of exact PLs
  xbys <- vapply(args, \(x) x$xBy(), 0)
  exacts <- vapply(args, \(x) x$isExact(), FALSE)
  if (! almostUnique(xbys[! exacts]))
    stop("All discrete PLs must have (almost) equal resolution `dx`.")
  xby1 <- if (any(! exacts)) args[! exacts][[1L]]$xBy() else args[[1L]]$xBy()
  args[exacts] <- lapply(args[exacts], pl_delimit, xby = xby1)
  
  # union of limits of all PLs
  xmin <- min(vapply(args, \(x) x$xMin(), 0))
  xmax <- max(vapply(args, \(x) x$xMax(), 0))
  args <- lapply(args, pl_delimit, xmin = xmin, xmax = xmax)
  
  # use greatest number of levels
  max_levels <- max(vapply(args, pl_num_levels, 0L))
  # concatenate vectorizations (`pl_to_vector()` discretizes any exacts)
  pl_mat <- sapply(args, pl_to_vector, num_levels = max_levels, simplify = TRUE)
  # encode vectorized landscapes as rows
  pl_mat <- t(pl_mat)
  # assign grid from first PL as attribute
  attr(pl_mat, "t_vals") <- pl_t(args[[1L]])
  
  return(pl_mat)
}

#' @rdname algebra
#' @export
pl_from_vector <- function(x, t_vals = NULL, drop_levels = FALSE) {
  
  # run checks
  if (is.null(t_vals) && is.null(attr(x, "t_vals")))
    stop("'t_vals' must be an attribute of `x` or provided separately.")
  if (! is.null(t_vals) && ! is.null(attr(x, "t_vals")))
    warning("`t_vals` was provided, so the `x` attribute will be ignored.")
  t_vals <- t_vals %||% attr(x, "t_vals")
  attr(x, "t_vals") <- NULL
  
  # check regularity of `t_vals` up to `epsi`
  stopifnot(almostUnique(diff(t_vals)))
  
  # reformat `x` if necessary; wrap levels to rows
  x <- as.vector(x)
  x <- matrix(x, ncol = length(t_vals), byrow = TRUE)
  
  if (drop_levels) {
    # remove empty levels
    empty_levs <- apply(x, 1L, \(r) all(r == 0))
    x <- x[seq(max(which(! empty_levs))), , drop = FALSE]
  }
  
  # construct persistence landscape
  new(PersistenceLandscape, t_vals, x)
}

#' @rdname algebra
#' @export
pl_from_matrix <- function(x, t_vals = NULL, drop_levels = FALSE) {
  
  # run checks
  if (! is.matrix(x)) x <- t(as.matrix(x))
  if (is.null(t_vals) && is.null(attr(x, "t_vals")))
    stop("'t_vals' must be an attribute of `x` or provided separately.")
  if (! is.null(t_vals) && ! is.null(attr(x, "t_vals")))
    warning("`t_vals` was provided, so the `x` attribute will be ignored.")
  t_vals <- t_vals %||% attr(x, "t_vals")
  attr(x, "t_vals") <- NULL
  
  # apply `pl_from_vector()` rowwise
  vs <- lapply(seq(nrow(x)), \(i) x[i, , drop = TRUE])
  res <- lapply(vs, pl_from_vector, t_vals = t_vals, drop_levels = drop_levels)
  if (nrow(x) == 1L) res <- res[[1L]]
  res
}
