#' @title 'plt' package
#' @description An Rcpp interface to Dłotko's Persistence Landscapes Toolbox
#'
#' @docType package
#' @name plt-package
#' @aliases plt
#' @import methods Rcpp
#' @useDynLib plt, .registration = TRUE
#' @details This package provides an R interface to the Persistence Landscapes
#'   Toolbox by Paweł Dłotko and Peter Bubenik. It is adapted from the
#'   'tdatools' package prepared by Jose Bouza.
"_PACKAGE"

# load modules (can be done anywhere, will activate at load time)
Rcpp::loadModule("persistence_landscape_module", TRUE)

# promise `PersistenceLandscape` to namespace (lazily) before load time
plt_ns <- topenv(environment())
if (! exists("PersistenceLandscape", envir = plt_ns)) {
  delayedAssign(
    "PersistenceLandscape",
    Rcpp::Module("persistence_landscape_module")$PersistenceLandscape,
    assign.env = plt_ns
  )
}
