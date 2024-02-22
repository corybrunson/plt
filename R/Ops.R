#' @title Prefix and Infix Operators for Persistence Landscapes
#' @description Perform arithmetic on persistence landscapes.
#'
#' @details Most of these operators extend selected members of the [S4 Arith
#'   group generic][base::Arithmetic] to the 'Rcpp_PersistenceLandscape' class:
#'   unary and binary `+` and `-` for persistence landscapes, `*`, and `/` for
#'   one persistence landscape (either factor or the numerator) and one numeric
#'   (either factor or the denominator). The exception is the binary `%*%`,
#'   which extends [matrix multiplication][base::matmult] to the inner product
#'   on persistence landscapes.
#'
#' @name Rcpp_PersistenceLandscape-methods
#' @docType methods
#' @include PersistenceLandscape.R

#' @param e1,e2,x,y Arguments of unary and binary operators.
#' @return A persistence landscape (an object of S4 class
#'   'Rcpp_PersistenceLandscape').
#' @seealso [base::Arithmetic] and [base::matmult] for base methods;
#'   [arithmetic] for function syntax and extensions.
#' @example inst/examples/ex-Ops.R
NULL

#' @rdname Rcpp_PersistenceLandscape-methods
#' @aliases Rcpp_PersistenceLandscape,missing-method
#' @export
setMethod(
  `+`,
  signature(e1 = "Rcpp_PersistenceLandscape", e2 = "missing"),
  # include ', e2' to avoid the following NOTE:
  # Error: no comma in argument list following \S4method
  function(e1, e2) e1
)

#' @rdname Rcpp_PersistenceLandscape-methods
#' @aliases Rcpp_PersistenceLandscape,Rcpp_PersistenceLandscape-method
#' @export
setMethod(
  `+`,
  signature(e1 = "Rcpp_PersistenceLandscape", e2 = "Rcpp_PersistenceLandscape"),
  function(e1, e2) e1$add(e2)
)

#' @rdname Rcpp_PersistenceLandscape-methods
#' @aliases Rcpp_PersistenceLandscape,missing-method
#' @export
setMethod(
  `-`,
  signature(e1 = "Rcpp_PersistenceLandscape", e2 = "missing"),
  # include ', e2' to avoid the following NOTE:
  # Error: no comma in argument list following \S4method
  function(e1, e2) e1$scale(-1)
)

#' @rdname Rcpp_PersistenceLandscape-methods
#' @aliases Rcpp_PersistenceLandscape,Rcpp_PersistenceLandscape-method
#' @export
setMethod(
  `-`,
  signature(e1 = "Rcpp_PersistenceLandscape", e2 = "Rcpp_PersistenceLandscape"),
  function(e1, e2) e1$add(e2$scale(-1))
)

#' @rdname Rcpp_PersistenceLandscape-methods
#' @aliases numeric,Rcpp_PersistenceLandscape-method
#' @export
setMethod(
  `*`,
  signature(e1 = "numeric", e2 = "Rcpp_PersistenceLandscape"),
  function(e1, e2) e2$scale(e1)
)

#' @rdname Rcpp_PersistenceLandscape-methods
#' @aliases Rcpp_PersistenceLandscape,numeric-method
#' @export
setMethod(
  `*`,
  signature(e1 = "Rcpp_PersistenceLandscape", e2 = "numeric"),
  function(e1, e2) e1$scale(e2)
)

#' @rdname Rcpp_PersistenceLandscape-methods
#' @aliases Rcpp_PersistenceLandscape,numeric-method
#' @export
setMethod(
  `/`,
  signature(e1 = "Rcpp_PersistenceLandscape", e2 = "numeric"),
  function(e1, e2) e1$scale(1/e2)
)

#' @rdname Rcpp_PersistenceLandscape-methods
#' @aliases Rcpp_PersistenceLandscape,Rcpp_PersistenceLandscape-method
#' @export
setMethod(
  `%*%`,
  signature(x = "Rcpp_PersistenceLandscape", y = "Rcpp_PersistenceLandscape"),
  function(x, y) x$inner(y)
)
