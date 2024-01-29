#' @title Plot Method for Persistence Landscapes.
#' @description A [base::plot()] method for persistence landscape objects.
#'
#' @name plot
#' @aliases plot,Rcpp_PersistenceLandscape-method
#' @include PersistenceLandscape.R
#' @param x A persistence landscape object of class 'Rcpp_PersistenceLandscape'.
#' @param replace_inf When using an exact representation of a persistence
#'   landscape, infinite values can appear. If not `NULL`, this value will
#'   replace `Inf` (and its negative will replace `-Inf`) in the plot.
#' @param n_levels Integer; number of levels to plot. If `NULL` (the default),
#'   determined to be the number of levels in `pl` or `x`.
#' @param palette Character; either a color palette from
#'   [grDevices::hcl.pals()], or the name of a palette function documented there
#'   (e.g. `'rainbow()'`), or a vector of colors for [grDevices::colorRamp()] to
#'   interpolate.
#' @param alpha,rev Parameters passed to [grDevices::hcl.colors()],
#'   [grDevices::colorRampPalette()].
#' @param ... Additional parameters passed to [base::plot()]. Values passed to
#'   `type` or `col` will be ignored with a message.
#' @param silent Logical; whether to silence messages.
#' @example inst/examples/ex-plot.R
#' @export
setMethod(
  "plot",
  c(x = "Rcpp_PersistenceLandscape"),
  function(
    x, replace_inf = NULL, n_levels = NULL,
    palette = "viridis", alpha = NULL, rev = FALSE, ...,
    silent = TRUE
  ) {
    # pre-process parameters
    n_env <- pl_num_levels(x)
    if (is.null(n_levels)) {
      n_levels <- n_env
    # } else if (n_levels > n_env) {
    } else if (n_levels < n_env && ! silent) {
      warning("`", deparse(substitute(x)), "` has more than `n_levels = ",
              n_levels,
              "` levels; lower levels will not be drawn.")
    }
    
    # pre-process internal representation of landscape
    x_exact <- pl_is_exact(x)
    internal <- x$getInternal()
    if (! is.null(replace_inf)) {
      if (! x_exact) {
        internal[internal == Inf] <- replace_inf
        internal[internal == -Inf] <- -replace_inf
      } else {
        for (level in seq(n_levels)) {
          internal[[level]][internal[[level]] == Inf] <- replace_inf
          internal[[level]][internal[[level]] == -Inf] <- replace_inf
        }
      }
    }
    
    # handle dots
    dots <- list(...)
    
    # create color palette
    cols <- colorLevels(
      n = n_levels,
      palette = palette, alpha = alpha, rev = rev
    )
    # reset `n_levels` if necessary
    if (n_levels > n_env) n_levels <- n_env
    
    # reconcile defaults with dots
    win_dots <- intersect(
      c("main", "sub", "xlab", "ylab", "asp", "xlim", "ylim"),
      names(dots)
    )
    rem_dots <- intersect(names(dots), c("type", "col"))
    aes_dots <- setdiff(names(dots), c(win_dots, rem_dots))
    # alert unused dots
    if (length(rem_dots) > 0L) {
      message(
        "Unrecognized or immutable parameters were ignored: ",
        gsub("(^list\\(|\\)$)", "", utils::capture.output(dput(dots[rem_dots])))
      )
    }
    
    # plot last level
    xyran <- apply(
      accessLevel(internal, 1L),
      2L,
      function(v) range(v[! is.infinite(v)])
    )
    def_dots <- list(
      type = "l",
      xlim = c(x$xMin(), x$xMax()), ylim = xyran[, 2L],
      xlab = "x", ylab = NA, lwd = 2
    )
    level_n <- accessLevel(internal, n_levels)
    plot_dots <- c(
      list(x = level_n[, 1L], y = level_n[, 2L], col = cols[[n_levels]]),
      utils::modifyList(def_dots, dots[c(win_dots, aes_dots)])
    )
    do.call(plot, plot_dots)
    # plot (any) remaining levels
    if (n_levels > 1) {
      for (level in seq(n_levels - 1L, 1L)) {
        level_d <- accessLevel(internal, level)
        line_dots <- c(
          list(x = level_d[, 1L], y = level_d[, 2L], col = cols[[level]]),
          utils::modifyList(list(lwd = 2), dots[aes_dots])
        )
        do.call(graphics::lines, line_dots)
      }
    }
  }
)

accessLevel <- function(internal, level) {
  # Trick to check if internal is vector or list.
  # From: https://stackoverflow.com/a/19501276/4556798
  if (is.atomic(internal)) {
    return(internal[level, , ])
  }
  else {
    return(internal[[level]])
  }
}

colorLevels <- function(n, palette, alpha, rev) {
  
  if (length(palette) == 1L) {
    
    # palette name
    pal <- try(grDevices::hcl.colors(
      n, palette = palette,
      alpha = alpha, rev = rev
    ), silent = TRUE)
    if (! inherits(pal, "try-error")) return(pal)
    
    # palette function
    pal_fun <- try(utils::getFromNamespace(palette, "grDevices"))
    if (! inherits(pal_fun, "try-error"))
      return(pal_fun(n = n, alpha = alpha, rev = rev))
    
    # single color
    pal <- try(grDevices::col2rgb(palette))
    if (! inherits(pal, "try-error")) return(rep(palette, n))
    
  } else {
    
    # multiple colors
    # TODO: unify error messages regardless of `length(palette)`
    pal <- grDevices::colorRampPalette(colors = palette)(n = n)
    if (! is.null(alpha)) pal <- grDevices::adjustcolor(pal, alpha.f = alpha)
    if (rev) pal <- rev(pal)
    return(pal)
    
  }
  
  # oops
  stop("`palette = '", palette, "'` is not a valid palette or color.")
  
}
