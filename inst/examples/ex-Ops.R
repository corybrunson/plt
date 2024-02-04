
# infix arithmetic with two landscapes
pt1 <- tdaunif::sample_torus_tube(100, 2.5)
pt2 <- tdaunif::sample_torus_tube(100, 2.5)
pd1 <- ripserr::vietoris_rips(pt1, dim = 2L, threshold = 3)
pd2 <- ripserr::vietoris_rips(pt2, dim = 2L, threshold = 3)
pl1 <- landscape(pd1, degree = 1L, exact = TRUE)
pl2 <- landscape(pd2, degree = 1L, exact = TRUE)
# plot
par(mfrow = c(3L, 2L), mar = c(0, 2, 0, 2))
plot(pl1)
plot(pl2)
# negation
plot(-pl1)
# scalar multiplication
plot(pl2 * 3)
# addition
plot(pl1 + 2 * pl2)
# subtraction
plot(pl1 - pl2 / 2)
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

# add an exact landscape to a discrete one
pl2d <- landscape(pd2, degree = 1, exact = FALSE,
                  xmin = 0, xmax = 1.2, by = 0.001)
# plot the summand landscapes and their sum with consistent parameters;
# the exact landscape is automatically converted to a discrete one
n_lev <- pl_num_levels(pl1)
par(mfrow = c(3, 1), mar = c(0, 2, 0, 2))
plot(pl1, xlim = c(0, 1.2), n_levels = n_lev)
plot(pl2d, xlim = c(0, 1.2), n_levels = n_lev)
plot(pl1 + pl2d, xlim = c(0, 1.2), n_levels = n_lev)
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

# inner products of similar landscapes
pl2 %*% pl2
pl2d %*% pl2d
# FIXME: This should not be much smaller than the above two.
pl2 %*% pl2d

\dontrun{

set.seed(031537L)

# compute exact landscape for a large sample
pt <- tdaunif::sample_torus_tube(600, 2.5)
pd <- ripserr::vietoris_rips(pt, dim = 2, threshold = 3)
pl <- landscape(pd, degree = 1, exact = TRUE)

# compute exact landscapes for a large sample of small samples
pl_list <- c()
for (i in seq(100)) {
  pti <- tdaunif::sample_torus_tube(100, 2.5)
  pdi <- ripserr::vietoris_rips(pti, dim = 2, threshold = 3)
  pli <- landscape(pdi, degree = 1L, exact = TRUE)
  pl_list <- c(pl_list, pli)
}

# compute the mean exact landscape
pl_avg <- Reduce(`+`, pl_list) / length(pl_list)

# compute the distance between the exact landscapes
pl_diff <- pl - pl_avg
print(pl_diff %*% pl_diff)

}
