#' @title Arithmetic and Statistical Operations on Persistence Landscapes
#' @description Calculate sums, scalar multiples, absolute values, means, inner
#'   products, extrema, and moments of persistent landscapes. These operations
#'   arise from the Hilbert space structure on persistence landscapes (Bubenik,
#'   2015).
#'
#' @details These functions are prefixed `pl_*()` to help users access them via
#'   tab-completion. Some take their names from the underlying S4 class methods
#'   and are only provided to enable composition via pipes: `add`, `scale`,
#'   `abs`, `inner`, `min` (`minimum`), `max` (`maximum`), and `moment`; `range`
#'   combines `min` and `max`. Others mimic classic R functions to handle lists
#'   of persistence landscapes: `sum`, `diff`, `mean`, `var`, and `sd`. Finally,
#'   some are vectorizations of the preceding: `vmin`, `vmax`, `vrange`, and
#'   `vmoment`.
#'
#' @name arithmetic
#' @include PersistenceLandscape.R
#' @param pl Persistence landscapes.
#' @param pl1,pl2 Persistence landscapes.
#' @param pl_list A list of persistence landscapes.
#' @param mult Double; a real-valued scale factor.
#' @param level Positive integer; the level of the persistence landscape (up to)
#'   whose moment to calculate.
#' @param p Positive integer or infinity; the power used to compute a norm or
#'   moment.
#' @param center Double; where to center the moment.
#' @return A persistence landscape (an object of S4 class
#'   'Rcpp_PersistenceLandscape'), a list of persistence landscapes, a real
#'   number, or a vector of real numbers.
#' @seealso [Rcpp_PersistenceLandscape-methods] for prefix and infix syntax.
#' @seealso [algebra], [analysis], [inference] for other landscape functions.
#' @example inst/examples/ex-arithmetic.R
NULL

#' @rdname arithmetic
#' @export
pl_add <- function(pl1, pl2) {
  pl1$add(pl2)
}

#' @rdname arithmetic
#' @export
pl_sum <- function(pl_list) {
  PLsum(pl_list)
}

#' @rdname arithmetic
#' @export
pl_scale <- function(pl, mult) {
  pl$scale(mult)
}

#' @rdname arithmetic
#' @export
pl_abs <- function(pl) {
  pl$abs()
}

#' @rdname arithmetic
#' @export
pl_diff <- function(pl_list) {
  PLdiff(pl_list)
}

#' @rdname arithmetic
#' @export
pl_mean <- function(pl_list) {
  PLmean(pl_list)
}

#' @rdname arithmetic
#' @export
pl_var <- function(pl_list, p = 2) {
  if (length(pl_list) <= 1L) return(NA_real_)
  p <- ensure_p(p)
  if (p == Inf) p <- 0
  PLvar(pl_list, p)
}

#' @rdname arithmetic
#' @export
pl_sd <- function(pl_list, p = 2) {
  if (length(pl_list) <= 1L) return(NA_real_)
  p <- ensure_p(p)
  if (p == Inf) p <- 0
  PLsd(pl_list, p)
}

#' @rdname arithmetic
#' @export
pl_inner <- function(pl1, pl2) {
  pl1$inner(pl2)
}

#' @rdname arithmetic
#' @export
pl_min <- function(pl, level = 1L) {
  pl$minimum(level)
}

#' @rdname arithmetic
#' @export
pl_max <- function(pl, level = 1L) {
  pl$maximum(level)
}

#' @rdname arithmetic
#' @export
pl_range <- function(pl, level = 1L) {
  c(pl$minimum(level), pl$maximum(level))
}

#' @rdname arithmetic
#' @export
pl_vmin <- function(pl, level = pl_num_levels(pl)) {
  if (level < 1L) return(numeric(0L))
  vapply(seq(level), function(l) pl$minimum(l), 0.)
}

#' @rdname arithmetic
#' @export
pl_vmax <- function(pl, level = pl_num_levels(pl)) {
  if (level < 1L) return(numeric(0L))
  vapply(seq(level), function(l) pl$maximum(l), 0.)
}

#' @rdname arithmetic
#' @export
pl_vrange <- function(pl, level = pl_num_levels(pl)) {
  if (level < 1L) return(cbind(numeric(0L), numeric(0L)))
  cbind(
    vapply(seq(level), function(l) pl$minimum(l), 0.),
    vapply(seq(level), function(l) pl$maximum(l), 0.)
  )
}

#' @rdname arithmetic
#' @export
pl_moment <- function(pl, p = 1L, center = 0, level = 1L) {
  p <- ensure_p(p)
  pl$moment(p, center, level)
}

#' @rdname arithmetic
#' @export
pl_vmoment <- function(pl, p = 1L, center = 0, level = pl_num_levels(pl)) {
  p <- ensure_p(p)
  vapply(
    seq(level),
    function(l) pl$moment(p, center, l),
    0.
  )
}
