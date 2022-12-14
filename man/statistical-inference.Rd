% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inference.r
\name{statistical-inference}
\alias{statistical-inference}
\alias{pl_integral}
\alias{pl_indicator}
\alias{pl_indicator_form}
\alias{pl_z_test}
\alias{pd_z_test}
\alias{pl_perm_test}
\title{Statistical Inference with Persistence Landscapes}
\usage{
pl_integral(pl, p = 1)

pl_indicator(pl, supports, r = 0)

pl_indicator_form(pl, supports, r = 0, p = 1)

pl_z_test(
  x,
  y,
  alternative = c("two.sided", "less", "greater"),
  conf.level = 0.95,
  supports = NULL,
  r = 0,
  p = 1
)

pd_z_test(
  x,
  y,
  degree = NULL,
  exact = FALSE,
  min_x = NULL,
  max_x = NULL,
  by = NULL,
  max_y = NULL,
  alternative = c("two.sided", "less", "greater"),
  conf.level = 0.95,
  supports = NULL,
  r = 0,
  p = 1
)

pl_perm_test(x, y, p = 1, complete = FALSE, max_iter = 1000L)
}
\arguments{
\item{pl}{A persistent landscape.}

\item{p}{Positive integer or infinity; the power used to compute an integral.}

\item{supports}{List of support intervals for landscape levels.}

\item{r}{Non-negative number; the power of the coefficient \eqn{1/k} in the
indicator linear form.}

\item{x, y}{Lists of persistence landscapes.}

\item{alternative}{a character string specifying the alternative
    hypothesis, must be one of \code{"two.sided"} (default),
    \code{"greater"} or \code{"less"}.  You can specify just the initial
    letter.}

\item{conf.level}{confidence level of the interval.}

\item{degree}{Non-negative integer; if input is a persistence diagram object,
then the dimension for which to compute a landscape. (For degree $d$, the
$(d+1)$th matrix in the list will be selected.)}

\item{exact}{Set to \code{TRUE} for exact computation, \code{FALSE} (default) for
discrete.}

\item{min_x, max_x}{Domain thresholds for discrete PL; if not specified, then
taken to be the support of the PL constructed from the data or the internal
values of the 'Rcpp_PersistenceLandscape' object.}

\item{by}{Domain grid diameter for discrete PL; if not specified, then set to
the power of 10 that yields between 100 and 1000 intervals.}

\item{max_y}{Numeric; the threshold used to compute the persistence diagram
(could be infinite).}

\item{complete}{Logical; whether to compute averages between all combinations
from the two lists of landscapes.}

\item{max_iter}{Positive integer; the maximum number of combinations using
which to estimate the null distance between mean landscapes.}
}
\value{
A persistence landscape (an object of S4 class
'Rcpp_PersistenceLandscape').
}
\description{
Compute test statistics and conduct null hypothesis tests for
persistence data using persistence landscapes. See Section 3 of Bubenik
(2015).
}
\examples{
# multiply by an indicator function
pd <- as_persistence(data.frame(dim = 1, birth = c(0, 1), death = c(2, 4)))
f <- list(c(.5, 3), c(1.5, 2), c(1.75, 2))
# exact landscape
pl <- landscape(pd, degree = 1, exact = TRUE)
par(mfcol = c(2L, 2L), mar = c(0, 2, 0, 2))
plot(pl, xlim = c(0, 5))
# plf <- pl$indicator(f, 0)
plf <- pl_indicator(pl, f)
plot(plf, xlim = c(0, 5))
# discrete landscape
pl <- landscape(pd, degree = 1, min_x = 0, max_x = 5, by = .2)
plot(pl)
# plf <- pl$indicator(f, 0)
plf <- pl_indicator(pl, f)
plot(plf, xlim = c(0, 5))
par(mfrow = c(1L, 1L), mar = c(5.1, 4.1, 4.1, 2.1))

# apply z-test to two sets of landscapes
set.seed(711018L)
circlescapes <- replicate(
  6,
  tdaunif::sample_circle(n = rpois(n = 1, lambda = 24)) |>
    ripserr::vietoris_rips(max_dim = 2L, threshold = 2) |>
    as_persistence() |>
    landscape(degree = 1, exact = TRUE)
)
toruscapes <- replicate(
  6,
  tdaunif::sample_torus_tube(n = rpois(n = 1, lambda = 24)) |>
    ripserr::vietoris_rips(max_dim = 2L, threshold = 2) |>
    as_persistence() |>
    landscape(degree = 1, exact = TRUE)
)
pl_z_test(circlescapes, toruscapes)
pl_perm_test(circlescapes, toruscapes)

\dontrun{
# benchmark one- and two-step computation of indicator-based linear form
bench::mark(
  pl$indicator_form(f, 0, p = 1),
  pl$indicator(f, 0)$integral(p = 1)
)
bench::mark(
  pl$indicator_form(f, r = 1, p = 1),
  pl$indicator(f, r = 1)$integral(p = 1)
)
}
}
\seealso{
PersistenceLandscape-methods
}
