#' @title Statistical Inference with Persistence Landscapes
#' @description Compute test statistics and conduct null hypothesis tests for
#'   persistence data using persistence landscapes. See Section 3 of Bubenik
#'   (2015).
#'
#' @name inference
#' @include PersistenceLandscape.R
#' @inheritParams landscape
#' @inheritParams analysis
#' @param x,y Lists of persistence landscapes.
#' @param complete Logical; whether to compute averages between all combinations
#'   from the two lists of landscapes.
#' @param max_iter Positive integer; the maximum number of combinations using
#'   which to estimate the null distance between mean landscapes.
#' @inheritParams stats::t.test
#' @return A list with class `"htest"` containing the following components:
#'   \item{`statistic`}{
#'   (z-test only) the value of the test statistic.
#'   }
#'   \item{`parameter`}{
#'   (z-test only) the degrees of freedom of the test,
#'   \eqn{\lvert x \rvert + \lvert y \rvert - 2}.
#'   }
#'   \item{`p.value`}{
#'   the p-value for the test.
#'   }
#'   \item{`estimate`}{
#'   the estimated difference difference in means.
#'   }
#'   \item{`null.value`}{
#'   the difference in means under the null hypothesis, always \eqn{0}.
#'   }
#'   \item{`alternative`}{
#'   a character string describing the alternative hypothesis.
#'   }
#'   \item{`method`}{
#'   a character string indicating the test performed.
#'   }
#'   \item{`conf.int`}{
#'   (z-test only) a confidence interval for the estimated difference in means.
#'   Depends on the choice of `conf.level`.
#'   }
#' @seealso [arithmetic], [algebra], [analysis] for other landscape functions.
#' @example inst/examples/ex-inference.R
NULL

#' @rdname inference
#' @export
pl_z_test <- function(
    x, y,
    alternative = c("two.sided", "less", "greater"), conf.level = 0.95,
    supports = NULL, r = 0, p = 1
) {
  # calculate summary statistics
  if (is.null(supports)) {
    xstat <- vapply(x, pl_integrate, 0, p = p)
    ystat <- vapply(y, pl_integrate, 0, p = p)
  } else {
    xstat <- lapply(x, pl_indicator_form, supports = supports, r = r, p = p)
    ystat <- lapply(y, pl_indicator_form, supports = supports, r = r, p = p)
  }
  # calculate z-statistic
  xbar <- mean(xstat)
  ybar <- mean(ystat)
  xvar <- stats::var(xstat)
  yvar <- stats::var(ystat)
  z <- (xbar - ybar) / sqrt(xvar + yvar)
  # REVIEW: Ensure consistency across test parameters.
  # compute p-value
  alternative <- match.arg(alternative, c("two.sided", "less", "greater"))
  if (alternative == "two-sided") {
    pval <- stats::pnorm(abs(z), lower.tail = FALSE) * 2
  } else {
    pval <- stats::pnorm(z, lower.tail = alternative == "less")
  }
  # compute confidence interval
  conf <- xbar - ybar + c(-1, 1) * stats::qnorm(conf.level) * sqrt(xvar + yvar)
  attr(conf, "conf.level") <- conf.level
  # format output
  res <- list(
    statistic = c(z = z),
    parameter = c(df = length(x) + length(y) - 2L),
    p.value = pval,
    estimate = c(`mean of x` = xbar, `mean of y` = ybar),
    null.value = c(`difference in means` = 0),
    # stderr = c(),
    alternative = alternative,
    method = "z-test",
    conf.int = conf
  )
  class(res) <- "htest"
  return(res)
}

#' @rdname inference
#' @export
pd_z_test <- function(
    x, y,
    degree = NULL, exact = FALSE,
    xmin = NULL, xmax = NULL, by = NULL,
    alternative = c("two.sided", "less", "greater"), conf.level = 0.95,
    supports = NULL, r = 0, p = 1
) {
  pl_x <- lapply(x, landscape, degree, exact, xmin, xmax, by)
  pl_y <- lapply(x, landscape, degree, exact, xmin, xmax, by)
  pl_z_test(pl_x, pl_y, alternative, conf.level, supports, r, p)
}

#' @rdname inference
#' @export
pl_perm_test <- function(x, y, p = 1, complete = FALSE, max_iter = 1000L) {
  p <- ensure_p(p)
  # calculate sample means
  xbar <- pl_mean(x)
  ybar <- pl_mean(y)
  # calculate distance between means
  xydist <- pl_distance(xbar, ybar)
  
  # prepare permutations
  if (complete && choose(length(x) + length(y), length(x)) < max_iter) {
    xypairs <- utils::combn(length(x) + length(y), length(x))
  } else {
    if (complete)
      warning("More than `max_iter` combinations; sampling instead.")
    xypairs <- replicate(
      max_iter,
      matrix(sample(seq(length(x) + length(y)), length(x)))
    )[, 1L, ]
  }
  # calculate distances
  xydists0 <- rep(NA_real_, ncol(xypairs))
  for (i in seq(ncol(xypairs))) {
    xi <- xypairs[, i][xypairs[, i] <= length(x)]
    yi <- xypairs[, i][xypairs[, i] > length(x)] - length(x)
    xbar0 <- pl_mean(c(x[xi], y[yi]))
    ybar0 <- pl_mean(c(x[setdiff(seq(length(x)), xi)],
                       y[setdiff(seq(length(y)), yi)]))
    xydists0[[i]] <- pl_distance(xbar0, ybar0)
  }
  # calculate p-value
  pval <- length(which(xydists0 >= xydist)) / length(xydists0)
  
  # format output
  res <- list(
    # statistic = c(),
    # parameter = c(df = df),
    p.value = pval,
    # conf.int = conf,
    estimate = c(`distance between mean landscapes` = xydist),
    null.value = c(`distance between mean landscapes` = 0),
    # stderr = c(),
    alternative = "greater",
    method = "permutation test"
  )
  class(res) <- "htest"
  return(res)
}
