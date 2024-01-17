#' @title Prefix and Infix Operators for Persistence Landscapes
#' @description Perform arithmetic on persistence landscapes.
#' 
#' @docType methods
#' @name persistence_landscape
#' @rdname Rcpp_PersistenceLandscape-methods
#' @include PersistenceLandscape.R
#' @param e1,e2,x,y Arguments of unary and binary operators.
#' @return A persistence landscape (an object of S4 class
#'   'Rcpp_PersistenceLandscape').
#' @seealso arithmetic-operations
#' @example inst/examples/ex-Ops.R
NULL

#' @rdname Rcpp_PersistenceLandscape-methods
#' @aliases Rcpp_PersistenceLandscape,missing-method
#' @export
setMethod(
  `+`,
  c(e1 = "Rcpp_PersistenceLandscape", e2 = "missing"),
  # include ', e2' to avoid the following NOTE:
  # Error: no comma in argument list following \S4method
  function(e1, e2) e1
)

#' @rdname Rcpp_PersistenceLandscape-methods
#' @aliases Rcpp_PersistenceLandscape,Rcpp_PersistenceLandscape-method
#' @export
setMethod(
  `+`,
  c(e1 = "Rcpp_PersistenceLandscape", e2 = "Rcpp_PersistenceLandscape"),
  function(e1, e2) e1$add(e2)
)

#' @rdname Rcpp_PersistenceLandscape-methods
#' @aliases Rcpp_PersistenceLandscape,missing-method
#' @export
setMethod(
  `-`,
  c(e1 = "Rcpp_PersistenceLandscape", e2 = "missing"),
  # include ', e2' to avoid the following NOTE:
  # Error: no comma in argument list following \S4method
  function(e1, e2) e1$scale(-1)
)

#' @rdname Rcpp_PersistenceLandscape-methods
#' @aliases Rcpp_PersistenceLandscape,Rcpp_PersistenceLandscape-method
#' @export
setMethod(
  `-`,
  c(e1 = "Rcpp_PersistenceLandscape", e2 = "Rcpp_PersistenceLandscape"),
  function(e1, e2) e1$add(e2$scale(-1))
)

#' @rdname Rcpp_PersistenceLandscape-methods
#' @aliases numeric,Rcpp_PersistenceLandscape-method
#' @export
setMethod(
  `*`,
  c(e1 = "numeric", e2 = "Rcpp_PersistenceLandscape"),
  function(e1, e2) e2$scale(e1)
)

#' @rdname Rcpp_PersistenceLandscape-methods
#' @aliases Rcpp_PersistenceLandscape,numeric-method
#' @export
setMethod(
  `*`,
  c(e1 = "Rcpp_PersistenceLandscape", e2 = "numeric"),
  function(e1, e2) e1$scale(e2)
)

#' @rdname Rcpp_PersistenceLandscape-methods
#' @aliases Rcpp_PersistenceLandscape,numeric-method
#' @export
setMethod(
  `/`,
  c(e1 = "Rcpp_PersistenceLandscape", e2 = "numeric"),
  function(e1, e2) e1$scale(1/e2)
)

#' @rdname Rcpp_PersistenceLandscape-methods
#' @aliases Rcpp_PersistenceLandscape,Rcpp_PersistenceLandscape-method
#' @export
setMethod(
  `%*%`,
  c(x = "Rcpp_PersistenceLandscape", y = "Rcpp_PersistenceLandscape"),
  function(x, y) x$inner(y)
)
