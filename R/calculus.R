#' @title Norms, Distances, and Integrals of Persistence Landscapes
#' @description Calculate distances and perform infinitesimal calculus on
#'   persistence landscapes. See Section 3 of Bubenik (2015).
#'
#' @name calculus
#' @include arithmetic.R
#' @inheritParams landscape
#' @inheritParams arithmetic
#' @param supports List of support intervals for landscape levels.
#' @param r Non-negative number; the power of the coefficient \eqn{1/k} in the
#'   indicator linear form.
#' @param p Positive integer or infinity; the power used to compute an integral.
#' @return A persistence landscape (an object of S4 class
#'   'Rcpp_PersistenceLandscape').
#' @seealso PersistenceLandscape-methods
#' @example inst/examples/ex-calculus.R
NULL

#' @rdname calculus
#' @export
pl_integrate <- function(pl, p = 1) {
  p <- ensure_p(p)
  pl$integrate(p)
}

#' @rdname calculus
#' @export
pl_distance <- function(pl1, pl2, p = 2) {
  p <- ensure_p(p)
  if (p == Inf) p <- 0
  pl1$distance(pl2, p)
}

#' @rdname calculus
#' @export
pl_dist <- function(pl_list, p = 2) {
  p <- ensure_p(p)
  PLdist(pl_list, p)
}

#' @rdname calculus
#' @export
pl_norm <- function(pl, p = 2) {
  p <- ensure_p(p)
  if (p == Inf) p <- 0
  pl$norm(p)
}

#' @rdname calculus
#' @export
pl_indicator <- function(pl, supports, r = 0) {
  pl$indicator(supports, r = r)
}

#' @rdname calculus
#' @export
pl_indicator_form <- function(pl, supports, r = 0, p = 1) {
  p <- ensure_p(p)
  pl$indicator_form(supports, r = r, p = p)
}
