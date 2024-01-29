#' @title Linear Algebra with Persistence Landscapes
#' @description Convert between persistence landscapes and their vectorizations,
#'   and perform linear algebra on vectorized persistence landscapes.
#'
#' @name algebra
#' @include arithmetic.R
#' @inheritParams landscape
#' @inheritParams arithmetic
#' @param ... Persistence landscapes or lists of persistence landscapes.
#' @param x A numeric vector, optionally with an attribute `"t_vals"` whose
#'   length divides that of `x`, but with no other attributes.
#' @param t_vals A numeric vector of equally-spaced values assumed to be the
#'   grid support of a discrete persistence landscape.
#' @example inst/examples/ex-algebra.R
NULL

#' @rdname algebra
#' @export
pl_t <- function(pl) {
  pl_internal <- if (pl$isExact()) {
    pl$toDiscrete()$getInternal()
  } else {
    pl$getInternal()
  }
  pl_internal[1L, , 1L]
}

#' @rdname algebra
#' @export
pl_to_vector <- function(pl, num_levels = pl_num_levels(pl)) {
  # get matrix of levels from discrete representation
  x <- pl_discretize(pl)$getInternal()[, , 2L, drop = 3L]
  # remove any superfluous levels
  if (num_levels < nrow(x)) x <- x[seq(num_levels), , drop = FALSE]
  # augment any additional empty levels
  if (num_levels > nrow(x))
    x <- rbind(x, matrix(0, nrow = num_levels - nrow(x), ncol = ncol(x)))
  # concatenate levels
  as.vector(t(x), mode = "any")
}

#' @rdname algebra
#' @export
pl_vectorize <- function(...) {
  
  # aggregate arguments into a single list of PLs
  args <- list(...)
  args <- unlist(args, recursive = TRUE)
  if (! all(vapply(args, inherits, FALSE, what = "Rcpp_PersistenceLandscape")))
    stop("`pl_vectorize()` accepts only PLs and lists of PLs.", call. = FALSE)
  
  # ensure that all PLs have the same `t_vals`
  xmins <- vapply(args, \(x) x$xMin(), 0)
  xmaxs <- vapply(args, \(x) x$xMax(), 0)
  xbys <- vapply(args, \(x) x$xBy(), 0)
  if (! almostUnique(xmins) || ! almostUnique(xmaxs) || ! almostUnique(xbys))
    stop("PL fields `xmin`, `xmax`, `dx` must be equal.")
  
  # ensure that all PLs are discretized
  args_exact <- vapply(args, \(pl) pl$isExact(), FALSE)
  if (any(args_exact))
    args[args_exact] <- lapply(args[args_exact], \(pl) pl$toDiscrete())
  
  # use greatest number of levels
  max_levels <- max(vapply(args, pl_num_levels, 0L))
  pl_mat <- sapply(args, pl_to_vector, num_levels = max_levels, simplify = TRUE)
  # encode vectorized landscapes as rows
  pl_mat <- t(pl_mat)
  
  # construct a matrix of vectorized PLs
  attr(pl_mat, "t_vals") <- pl_t(args[[1L]])
  return(pl_mat)
}

#' @rdname algebra
#' @export
pl_devectorize <- function(x, t_vals = NULL, drop_levels = FALSE) {
  
  # run checks
  if (is.null(t_vals) && is.null(attr(x, "t_vals")))
    stop("'t_vals' must be an attribute of `x` or provided separately.")
  if (! is.null(t_vals) && ! is.null(attr(x, "t_vals")))
    warning("`t_vals` was provided, so the `x` attribute will be ignored.")
  t_vals <- t_vals %||% attr(x, "t_vals")
  attr(x, "t_vals") <- NULL
  
  # recurse a matrix to vectors
  if (is.matrix(x) && nrow(x) > 1L) {
    vs <- lapply(seq(nrow(x)), \(i) x[i, , drop = TRUE])
    return(lapply(vs, pl_devectorize,
                  t_vals = t_vals, drop_levels = drop_levels))
  }
  
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
