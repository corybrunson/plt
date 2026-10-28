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
#' @param palette Character; either the name of a qualitative palette from
#'   [grDevices::palette.pals()] (passed to [grDevices::palette.colors()]) or
#'   a sequential/diverging palette from [grDevices::hcl.pals()] (passed to
#'   [grDevices::hcl.colors()]), the name of a palette function in
#'   [grDevices] (e.g. `'rainbow'` or `'terrain.colors'`), a single color, or a
#'   vector of colors for [grDevices::colorRampPalette()] to interpolate. The
#'   default `palette = "default"` uses [grDevices::palette.colors()] with
#'   `palette = "Okabe-Ito"` (the colorblind-friendly default introduced in R
#'   4.0.0), falling back to the current [grDevices::palette()].
#' @param alpha,rev Parameters passed to [grDevices::palette.colors()],
#'   [grDevices::hcl.colors()], or [grDevices::colorRampPalette()]. For
#'   qualitative palettes returned by [grDevices::palette.colors()], `rev`
#'   reverses the order. Note that `recycle = TRUE` is set internally when
#'   calling [grDevices::palette.colors()].
#' @param ... Additional parameters passed to [base::plot()]. Values passed to
#'   `type` or `col` will be ignored with a message.
#' @param silent Logical; whether to silence messages.
#' @returns None (invisible NULL).
#' @seealso [Rcpp_PersistenceLandscape-class] for the exported C++ class and
#'   [pl_new()] for the R wrapper.
#' @example inst/examples/ex-plot.R
#' @export
setMethod(
  "plot",
  c(x = "Rcpp_PersistenceLandscape"),
  function(
    x, replace_inf = NULL, n_levels = NULL,
    palette = "default", alpha = NULL, rev = FALSE, ...,
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
    
    # range across all levels
    xyran <- apply(
      # accessLevel(internal, 1L),
      do.call(
        rbind,
        lapply(seq(n_env), accessLevel, internal = internal)
      ),
      2L,
      function(v) range(v[! is.infinite(v)])
    )
    # plot last level
    def_dots <- list(
      type = "l",
      xlim = xyran[, 1L],
      # xlim = c(x$xMin(), x$xMax()),
      ylim = xyran[, 2L],
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
  } else {
    return(internal[[level]])
  }
}

colorLevels <- function(n, palette, alpha, rev) {
  
  if (length(palette) > 1L) {
    
    # multiple colors: interpolate with a color ramp
    pal <- grDevices::colorRampPalette(colors = palette)(n = n)
    if (! is.null(alpha)) pal <- grDevices::adjustcolor(pal, alpha.f = alpha)
    if (isTRUE(rev)) pal <- rev(pal)
    return(pal)
    
  } else if (length(palette) == 1L) {
    
    # default palette (depends on R version)
    if (identical(palette, "default")) {
      pal <- grDevices::palette.colors(
        n = n,
        alpha = alpha,
        recycle = TRUE
      )
      if (isTRUE(rev)) {
        pal <- rev(pal)
      }
      return(pal)
    }
    
    # try qualitative palettes from palette.colors() (R 4.0+)
    pal <- try(
      grDevices::palette.colors(
        n = n,
        palette = palette,
        alpha = alpha,
        recycle = TRUE
      ),
      silent = TRUE
    )
    if (! inherits(pal, "try-error")) {
      if (isTRUE(rev)) {
        pal <- rev(pal)
      }
      return(pal)
    }
    
    # try HCL palettes from hcl.colors() (R 4.0+)
    pal <- try(
      grDevices::hcl.colors(
        n = n,
        palette = palette,
        alpha = alpha,
        rev = rev
      ),
      silent = TRUE
    )
    if (! inherits(pal, "try-error")) {
      return(pal)
    }
    
    # try palette function from grDevices (e.g. rainbow, topo.colors)
    pal_fun <- try(utils::getFromNamespace(palette, "grDevices"), silent = TRUE)
    if (! inherits(pal_fun, "try-error")) {
      pal_args <- list(n = n)
      if (! is.null(alpha)) {
        pal_args$alpha <- alpha
      }
      if (isTRUE(rev)) {
        pal_args$rev <- rev
      }
      pal <- try(do.call(pal_fun, pal_args), silent = TRUE)
      if (inherits(pal, "try-error")) {
        pal <- try(pal_fun(n = n), silent = TRUE)
      }
      if (! inherits(pal, "try-error")) {
        return(pal)
      }
    }
    
    # try as a single color
    pal <- try(grDevices::col2rgb(palette), silent = TRUE)
    if (! inherits(pal, "try-error")) {
      pal <- grDevices::rgb(t(pal / 255))
      if (! is.null(alpha)) {
        alpha_hex <- format(
          as.hexmode(round(alpha * 255 + 1e-04)),
          width = 2L,
          upper.case = TRUE
        )
        pal <- rep_len(paste0(substr(pal, 1L, 7L), alpha_hex), n)
      } else {
        pal <- substr(pal, 1L, 7L)
        pal <- rep_len(pal, n)
      }
      if (isTRUE(rev)) {
        pal <- rev(pal)
      }
      return(pal)
    }
    
  }
  
  # oops
  stop("`palette = '", palette, "'` is not recognized.")
  
}
