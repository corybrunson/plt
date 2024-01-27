#' @title Linear Algebra with Persistence Landscapes
#' @description Convert between persistence landscapes and their vectorizations,
#'   and perform linear algebra on vectorized persistence landscapes.
#'
#' @name algebra
#' @include arithmetic.R
#' @inheritParams landscape
#' @inheritParams arithmetic
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
pl_vectorize <- function(pl) {
  if (pl$isExact()) pl <- pl$toDiscrete()
  pl_vec <- as.vector(pl)
  attr(pl_vec, "t_vals") <- pl_t(pl)
  pl_vec
}

#' @rdname algebra
#' @export
pl_devectorize <- function(x, t_vals = NULL) {
  
  # run checks
  if (is.null(t_vals) && is.null(attr(x, "t_vals")))
    stop("'t_vals' must be an attribute of `x` or provided separately.")
  if (! is.null(t_vals) && ! is.null(attr(x, "t_vals")))
    warning("`t_vals` was provided, so the `x` attribute will be ignored.")
  t_vals <- t_vals %||% attr(x, "t_vals")
  attr(x, "t_vals") <- NULL
  
  # check regularity of `t_vals` up to `epsi`
  stopifnot(almostUnique(diff(t_vals)))
  
  # reformat `x` if necessary
  stopifnot(is.vector(x) || is.matrix(x))
  if (! is.matrix(x)) x <- matrix(x, ncol = length(t_vals), byrow = TRUE)
  
  # construct persistence landscape
  new(PersistenceLandscape, t_vals, x)
}
