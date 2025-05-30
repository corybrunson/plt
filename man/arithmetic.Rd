% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/arithmetic.R
\name{arithmetic}
\alias{arithmetic}
\alias{pl_add}
\alias{pl_sum}
\alias{pl_scale}
\alias{pl_abs}
\alias{pl_diff}
\alias{pl_mean}
\alias{pl_var}
\alias{pl_sd}
\alias{pl_inner}
\alias{pl_min}
\alias{pl_max}
\alias{pl_range}
\alias{pl_vmin}
\alias{pl_vmax}
\alias{pl_vrange}
\alias{pl_moment}
\alias{pl_vmoment}
\title{Arithmetic and Statistical Operations on Persistence Landscapes}
\usage{
pl_add(pl1, pl2)

pl_sum(pl_list)

pl_scale(pl, mult)

pl_abs(pl)

pl_diff(pl_list)

pl_mean(pl_list)

pl_var(pl_list, p = 2)

pl_sd(pl_list, p = 2)

pl_inner(pl1, pl2)

pl_min(pl, level = 1L)

pl_max(pl, level = 1L)

pl_range(pl, level = 1L)

pl_vmin(pl, level = pl_num_levels(pl))

pl_vmax(pl, level = pl_num_levels(pl))

pl_vrange(pl, level = pl_num_levels(pl))

pl_moment(pl, p = 1L, center = 0, level = 1L)

pl_vmoment(pl, p = 1L, center = 0, level = pl_num_levels(pl))
}
\arguments{
\item{pl1, pl2}{Persistence landscapes.}

\item{pl_list}{A list of persistence landscapes. A non-list object will first
be wrapped in a list.}

\item{pl}{Persistence landscapes.}

\item{mult}{Double; a real-valued scale factor.}

\item{p}{Positive integer or infinity; the power used to compute a norm or
moment.}

\item{level}{Positive integer; the level of the persistence landscape (up to)
whose moment to calculate. This value is internally decreased by \code{1L} to
prevent off-by-one errors when passing to C++ code.}

\item{center}{Double; where to center the moment.}
}
\value{
A persistence landscape (an object of S4 class
'Rcpp_PersistenceLandscape'), a list of persistence landscapes, a real
number, or a vector of real numbers.
}
\description{
Calculate sums, scalar multiples, absolute values, means, inner
products, extrema, and moments of persistent landscapes. These operations
arise from the Hilbert space structure on persistence landscapes (Bubenik,
2015).
}
\details{
These functions are prefixed \verb{pl_*()} to help users access them via
tab-completion. Some take their names from the underlying S4 class methods
and are only provided to enable composition via pipes: \code{add}, \code{scale},
\code{abs}, \code{inner}, \code{min} (\code{minimum}), \code{max} (\code{maximum}), and \code{moment}; \code{range}
combines \code{min} and \code{max}. Others mimic classic R functions to handle lists
of persistence landscapes: \code{sum}, \code{diff}, \code{mean}, \code{var}, and \code{sd}. Finally,
some are vectorizations of the preceding: \code{vmin}, \code{vmax}, \code{vrange}, and
\code{vmoment}.
}
\examples{

# scale a landscape
x <- tdaunif::sample_torus_tube(40, 2.5)
pd <- ripserr::vietoris_rips(x, dim = 1L, threshold = 2)
pl <- pl_new(pd, degree = 1, exact = FALSE, xmax = 2.5, xby = 0.1)
print(pl$getInternal()[2, , ])
print(pl$scale(0.5)$getInternal()[2, , ])

# create two landscapes from the same sampling distribution
x <- tdaunif::sample_torus_tube(40, 2.5)
y <- tdaunif::sample_torus_tube(40, 2.5)
x_pd <- ripserr::vietoris_rips(x, dim = 1, threshold = 2)
y_pd <- ripserr::vietoris_rips(y, dim = 1, threshold = 2)
x_pl <- pl_new(x_pd, degree = 1, exact = FALSE, xmax = 2.5, xby = 0.1)
y_pl <- pl_new(y_pd, degree = 1, exact = FALSE, xmax = 2.5, xby = 0.1)

# compare landscapes and calculate mean landscape
print(x_pl$getInternal())
print(y_pl$getInternal())
print(pl_mean(list(x_pl, y_pl))$getInternal())

\dontrun{

set.seed(492869L)

# compute landscape for a large sample
pt <- tdaunif::sample_torus_tube(1000, 5)
pd <- ripserr::vietoris_rips(pt, dim = 2, threshold = 2)
pl <- pl_new(pd, degree = 1, exact = FALSE, xby = 0.1, xmin = 0, xmax = 2)

# compute landscapes for a large sample of small samples
pl_list <- c()
for (i in seq(100)) {
  pti <- tdaunif::sample_torus_tube(100, 5)
  pdi <- ripserr::vietoris_rips(pti, dim = 2, threshold = 2)
  pli <- pl_new(pdi, degree = 1, exact = FALSE, xby = 0.1, xmin = 0, xmax = 2)
  pl_list <- c(pl_list, pli)
}

# compute the mean landscape
pl_avg <- pl_mean(pl_list)

# compute the distance between the landscapes
pl_diff <- pl$add(pl_avg$scale(-1))
print(pl_inner(pl_diff, pl_diff))

}
}
\seealso{
\link{Rcpp_PersistenceLandscape-methods} for prefix and infix syntax.

\link{algebra}, \link{analysis}, \link{inference} for other landscape functions.
}
