
# TODO: When `epsi` in 'PersistenceLandscape.h' is made into a user-settable
# option, set this equal to that option.
epsi <- 0.000005

# insert critical points (y = 0)
insert_critical <- function(x) {
  wz <- which(sign(x[-1L, 2L]) != 0 & sign(x[-nrow(x), 2L]) != 0 &
                sign(x[-1L, 2L]) != sign(x[-nrow(x), 2L]))
  if (length(wz) == 0L) return(x)
  xz <- x[wz, 1L] +
    x[wz, 2L] / (x[wz, 2L] - x[wz + 1L, 2L]) * (x[wz + 1L, 1L] - x[wz, 1L])
  x1 <- matrix(0, nrow = nrow(x) + length(wz), ncol = 2L)
  w1 <- wz + seq(length(wz))
  x1[-w1, ] <- x[, ]
  x1[w1, 1L] <- xz; x1[w1, 2L] <- 0
  x1
}

# calculate the p-norm, via integration, of a piecewise linear function
integrate_power <- function(x, p) {
  # ignore flatlines (assume y = 0 from `-Inf` and to `Inf`)
  if (any(diff(x[, 1L]) <= 0)) warning("Zero-length interval encountered.")
  if (x[1L, 1L] == -Inf) x <- x[-1L, , drop = FALSE]
  if (x[nrow(x), 1L] == Inf) x <- x[-nrow(x), , drop = FALSE]
  # insert critical points if necessary
  x <- insert_critical(x)
  # calculate integral of y ^ p analytically
  integrand <- ifelse(
    # consider intervals with near-constant height as rectangles
    abs( ( x[-1L, 2L] - x[-nrow(x), 2L] ) / (p+1) ) < epsi,
    x[-1L, 2L] ^ p,
    ( x[-1L, 2L] ^ (p+1) - x[-nrow(x), 2L] ^ (p+1) ) /
      ( x[-1L, 2L] - x[-nrow(x), 2L] ) / (p+1)
  )
  # print(abs(integrand) * ( x[-1L, 1L] - x[-nrow(x), 1L] ))
  sum(abs(integrand) * ( x[-1L, 1L] - x[-nrow(x), 1L] ))
}

# norm for exact landscapes
pl_norm_exact <- function(pl, p = 2) {
  p <- plt:::ensure_p(p)
  # calculate norm
  if (p == Inf) {
    max(vapply(pl$getInternal(), function(x) max(abs(x[, 2L])), 0))
  } else {
    z <- vapply(pl$getInternal(), function(x) {
      integrate_power(x, p)
    }, 0)
    sum(z) ^ (1/p)
  }
}

# norm for discrete landscapes
pl_norm_discrete <- function(pl, p = 2) {
  p <- plt:::ensure_p(p)
  # calculate norm
  if (p == Inf) {
    max(abs(pl$getInternal()[, , 2L]))
  } else {
    sum(apply(pl$getInternal(), 1L, integrate_power, p = p)) ^ (1/p)
  }
}

# norm for landscapes
pl_norm <- function(pl, p = 2) {
  stopifnot(is.numeric(p), p >= 1)
  switch(
    pl_type(pl),
    exact = pl_norm_exact(pl = pl, p = p),
    discrete = pl_norm_discrete(pl = pl, p = p)
  )
}

# L^p distance between two comparable landscapes
pl_dist <- function(pl1, pl2, p = 2) {
  pl <- pl2 - pl1
  pl_norm(pl = pl, p = p)
}

# persistence data for two point clouds
x1 <- tdaunif::sample_torus_tube(120, 2.5)
x2 <- tdaunif::sample_torus_tube(60, 2.5)
pd1 <- as_persistence(ripserr::vietoris_rips(x1, dim = 1, threshold = 2))
pd2 <- as_persistence(ripserr::vietoris_rips(x2, dim = 1, threshold = 2))

# two exact persistence landscapes
ple1 <- landscape(pd1, degree = 1, exact = TRUE)
ple2 <- landscape(pd2, degree = 1, exact = TRUE)

# inf norm
pl_norm_exact(ple1, p = Inf)
# 2 norm
pl_norm_exact(ple1, p = 2)
# 1 norm
pl_norm_exact(ple1, p = 1)

# two discrete persistence landscapes
pl1 <- landscape(pd1, degree = 1, xmin = 0, xmax = 1.5, by = 0.01)
pl2 <- landscape(pd2, degree = 1, xmin = 0, xmax = 1.5, by = 0.01)

# inf norm
pl_norm_discrete(pl1, p = Inf)
# 2 norm
pl_norm_discrete(pl1, p = 2)
# 1 norm
pl_norm_discrete(pl1, p = 1)

# arithmetic differences
ple12 <- ple2 - ple1
plot(ple12)
pl12 <- pl2 - pl1
plot(pl12)
dim(pl12$getInternal())

# L^inf distance
pl_dist(pl1, pl2, p = Inf)
# L^2 distance
pl_dist(pl1, pl2, p = 2)
# L^1 distance
pl_dist(pl1, pl2, p = 1)

# benchmark comparison of distance functions
ps <- c(Inf, 1 + 2 ^ seq(8L, 0L, -2L), 1)
# exact
bench::mark(
  Cpp = vapply(ps, function(p) pl_distance(ple1, ple2, p), 0),
  R = vapply(ps, function(p) pl_dist(ple1, ple2, p), 0),
  check = \(x, y) all(abs(y - x) < epsi)
)
# discrete
bench::mark(
  Cpp = vapply(ps, function(p) pl_distance(pl1, pl2, p), 0),
  R = vapply(ps, function(p) pl_dist(pl1, pl2, p), 0),
  check = \(x, y) all(abs(y - x) < epsi)
)

# toy example from Bubenik & Dłotko (2017)

n <- 5L

pls_n <- list()
for (i in seq(11L)) {
  A <- do.call(rbind, lapply(seq(0L, n), function(m) {
    tdaunif::sample_circle(n = 50L, sd = .15/2) +
      matrix(rep(c(2 * m, 0), each = 50L), ncol = 2L)
  }))
  pd <- as_persistence(ripserr::vietoris_rips(A, dim = 1, threshold = n * 2))
  pl <- landscape(pd, degree = 1L, exact = TRUE)
  pls_n <- c(pls_n, list(pl))
}

pls_n_avg <- pl_mean(pls_n)
plot(pls_n_avg)

pl_norm(pls_n_avg, 1)
pl_norm(pls_n_avg, 2)
pl_norm(pls_n_avg, Inf)
