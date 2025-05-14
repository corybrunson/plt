#' @title Integrals, Distances, and Norms of Persistence Landscapes
#' @description Calculate distances and perform infinitesimal calculus on
#'   persistence landscapes. See Section 3 of Bubenik (2015).
#'
#' @details `pl_integrate()` computes the \eqn{L^p} integral of a persistence
#'   landscape. The `p` parameter denotes the integral's power;
#'   `p = 1` computes the area under the curves, while `p = 2` corresponds to the
#'    \eqn{L^2} norm. Letting `p = Inf` returns the maximum landscape height.
#'
#'   `pl_distance()` calculates the \eqn{L^p} distance between two persistence
#'   landscapes. 
#'
#'   `pl_dist()` utilizes `pl_distance()` to compute all pairwise distances
#'   between a list of landscapes. It returns a symmetric matrix of distances.
#'
#'   `pl_norm()` returns the \eqn{L^p} norm of a single landscape, used as a
#'   scalar summary of its magnitude. Similar to `pl_integrate()`, it supports
#'   all finite powers `p >= 1`, as well as `p = Inf` for the supremum norm.
#'
#'   `pl_indicator()` computes the linear form defined by a sum of indicator
#'   functions over user-specified support intervals, weighted by \eqn{1/k^r}
#'   where \eqn{k} indexes the landscape levels.
#'
#'   `pl_indicator_form()` evaluates the integral of the indicator function form
#'   described above with respect to the \eqn{L^p} norm. It returns a real number
#'   representing the landscape's integral onto the chosen intervals.
#'
#' @name analysis
#' @include arithmetic.R
#' @inheritParams pl_new
#' @inheritParams arithmetic
#' @param supports List of support intervals for landscape levels.
#' @param r Non-negative number; the power of the coefficient \eqn{1/k} in the
#'   indicator linear form.
#' @param p Positive integer or infinity; the power used to compute an integral.
#' @return A persistence landscape (an object of S4 class
#'   'Rcpp_PersistenceLandscape'), a real number, or a vector of real numbers.
#' @seealso [arithmetic], [algebra], [inference] for other landscape functions.
#' @example inst/examples/ex-analysis.R
NULL

#' @rdname analysis
#' @export
pl_integrate <- function(pl, p = 1) {
  p <- ensure_p(p)
  pl$integrate(p)
}

#' @rdname analysis
#' @export
pl_distance <- function(pl1, pl2, p = 2) {
  p <- ensure_p(p)
  if (p == Inf) p <- 0
  pl1$distance(pl2, p)
}

#' @rdname analysis
#' @export
pl_dist <- function(pl_list, p = 2) {
  if (! is.list(pl_list)) return(pl_dist(list(pl_list)))
  p <- ensure_p(p)
  PLdist(pl_list, p)
}

#' @rdname analysis
#' @export
pl_norm <- function(pl, p = 2) {
  p <- ensure_p(p)
  if (p == Inf) p <- 0
  pl$norm(p)
}

#' @rdname analysis
#' @export
pl_indicator <- function(pl, supports, r = 0) {
  pl$indicator(supports, r = r)
}

#' @rdname analysis
#' @export
pl_indicator_form <- function(pl, supports, r = 0, p = 1) {
  p <- ensure_p(p)
  pl$indicator_form(supports, r = r, p = p)
}
