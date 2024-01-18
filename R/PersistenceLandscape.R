Rcpp_PersistenceLandscape <- setClass("Rcpp_PersistenceLandscape")

#' @rdname PersistenceLandscape
#' @inheritParams base::summary
#' @export
summary.Rcpp_PersistenceLandscape <- function(object, ...) {
  res <- list(
    str = pl_type(object),
    n.levels = pl_num_levels(object),
    limits = pl_limits(object),
    resolution = if (object$isExact()) NA_real_ else object$xBy(),
    support = pl_support(object)
    # range = c(pl_min(object), pl_max(object)),
    # magnitude = object %*% object,
    # integral = pl_integrate(object)
  )
  class(res) <- "summary.Rcpp_PersistenceLandscape"
  res
}

#' @rdname PersistenceLandscape
#' @inheritParams base::print
#' @export
print.summary.Rcpp_PersistenceLandscape <- function(
    x, digits = max(1L, getOption("digits") - 3L), ...
) {
  fmt <- function(s) {
    if (is.na(s) || is.infinite(s)) return(format(s))
    format(round(s, max(0, digits - log10(abs(s)))))
  }
  cat("Internal representation: ", x$str, "\n")
  cat("Number of levels: ", x$n.levels, "\n")
  cat(
    "Representation limits: (",
    fmt(x$limits[[1L]]), ",", fmt(x$limits[[2L]]), ")",
    if (x$str == "discrete")
      paste0(" at resolution ", format(round(x$resolution))),
    "\n"
  )
  # cat("Landscape range: (",
  #     fmt(x$range[[1L]]), ",", fmt(x$range[[2L]]), ")", "\n")
  # cat("Magnitude: ", fmt(x$magnitude), "\n")
  # cat("Integral:  ", fmt(x$integral), "\n")
}

#' @rdname PersistenceLandscape
#' @inheritParams base::as.vector
#' @export
as.vector.Rcpp_PersistenceLandscape <- function(x, mode = "any") {
  # get discrete representation for consistent abscissa values
  internal <- pl_discretize(x)$getInternal()
  
  # export concatenation of levels
  as.vector(t(internal[, , 2L]), mode = mode)
}

#' @rdname PersistenceLandscape
#' @inheritParams base::as.data.frame
#' @export
as.data.frame.Rcpp_PersistenceLandscape <- function(
    x, row.names = NULL, optional = FALSE, ...
) {
  internal <- x$getInternal()
  
  if (x$isExact()) {
    
    levs <- vapply(internal, nrow, 0L)
    lev <- rep(seq_along(internal), levs)
    yfy <- Reduce(rbind, internal)
    df <- data.frame(lev, yfy)
    names(df) <- c("level", "x", "fx")
    if (! is.null(row.names)) {
      rownames(df) <- row.names
    } else if (! optional) {
      id <- unlist(lapply(levs, seq))
      rownames(df) <- paste(lev, id, sep = ".")
    }
    
  } else {
    
    lev <- rep(seq(dim(internal)[[2L]]), each = dim(internal)[[1L]])
    y <- as.vector(t(internal[, , 1L]))
    fy <- as.vector(t(internal[, , 2L]))
    df <- data.frame(level = lev, x = y, fx = fy)
    if (! is.null(row.names)) {
      # adapted from `as.data.frame.matrix()`
      .rowNamesDF(df, make.names = TRUE) <- row.names
    } else if (! optional) {
      id <- rep(seq(dim(internal)[[1L]]), times = dim(internal)[[2L]])
      attr(df, "row.names") <- paste(lev, id, sep = ".")
    }
  }
  
  df
}
